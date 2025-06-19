"""Tests for progress tracking and checkpointing."""

import pytest
import json
import time
from pathlib import Path
from unittest.mock import patch, MagicMock

from chimeric_detector.utils.progress import ProgressTracker, ProgressReporter, ProgressState


class TestProgressState:
    """Test ProgressState data class."""
    
    def test_progress_state_creation(self):
        """Test ProgressState creation and basic functionality."""
        state = ProgressState(
            total_contigs=10,
            processed_contigs=3,
            current_contig="contig3",
            start_time=time.time(),
            last_update=time.time(),
            completed_contigs=["contig1", "contig2"],
            failed_contigs=["contig0"],
            partial_results={"contig1": "result1"}
        )
        
        assert state.total_contigs == 10
        assert state.processed_contigs == 3
        assert state.current_contig == "contig3"
        assert len(state.completed_contigs) == 2
        assert len(state.failed_contigs) == 1
        assert "contig1" in state.partial_results
    
    def test_progress_state_serialization(self):
        """Test ProgressState to_dict and from_dict methods."""
        original_state = ProgressState(
            total_contigs=5,
            processed_contigs=2,
            current_contig="contig2",
            start_time=1000.0,
            last_update=1100.0,
            completed_contigs=["contig1"],
            failed_contigs=[],
            partial_results={"contig1": ["result"]}
        )
        
        # Convert to dict and back
        state_dict = original_state.to_dict()
        restored_state = ProgressState.from_dict(state_dict)
        
        assert restored_state.total_contigs == original_state.total_contigs
        assert restored_state.processed_contigs == original_state.processed_contigs
        assert restored_state.current_contig == original_state.current_contig
        assert restored_state.completed_contigs == original_state.completed_contigs
        assert restored_state.failed_contigs == original_state.failed_contigs
        assert restored_state.partial_results == original_state.partial_results


class TestProgressTracker:
    """Test ProgressTracker functionality."""
    
    def test_progress_tracker_initialization(self, temp_dir):
        """Test ProgressTracker initialization."""
        checkpoint_file = temp_dir / "checkpoint.json"
        
        tracker = ProgressTracker(total_contigs=5, checkpoint_file=checkpoint_file)
        
        assert tracker.state.total_contigs == 5
        assert tracker.state.processed_contigs == 0
        assert tracker.state.current_contig is None
        assert len(tracker.state.completed_contigs) == 0
        assert len(tracker.state.failed_contigs) == 0
        assert tracker.checkpoint_file == checkpoint_file
    
    def test_progress_tracker_no_checkpoint(self):
        """Test ProgressTracker without checkpoint file."""
        tracker = ProgressTracker(total_contigs=3)
        
        assert tracker.state.total_contigs == 3
        assert tracker.checkpoint_file is None
    
    def test_contig_processing_lifecycle(self, temp_dir):
        """Test complete contig processing lifecycle."""
        tracker = ProgressTracker(total_contigs=3)
        
        # Start first contig
        tracker.start_contig("contig1")
        assert tracker.state.current_contig == "contig1"
        
        # Complete first contig
        tracker.complete_contig("contig1", "result1")
        assert tracker.state.processed_contigs == 1
        assert "contig1" in tracker.state.completed_contigs
        assert tracker.state.partial_results["contig1"] == "result1"
        assert tracker.state.current_contig is None
        
        # Start and fail second contig
        tracker.start_contig("contig2")
        tracker.fail_contig("contig2", "Error message")
        assert tracker.state.processed_contigs == 2
        assert "contig2" in tracker.state.failed_contigs
        assert "contig2" not in tracker.state.partial_results
        
        # Complete third contig
        tracker.start_contig("contig3")
        tracker.complete_contig("contig3", "result3")
        assert tracker.state.processed_contigs == 3
        assert tracker.is_complete()
    
    def test_contig_mismatch_warning(self, temp_dir):
        """Test warning when contig names don't match."""
        tracker = ProgressTracker(total_contigs=2)
        
        tracker.start_contig("contig1")
        
        with patch('chimeric_detector.utils.progress.logging.getLogger') as mock_logger:
            # Complete different contig than started
            tracker.complete_contig("contig2", "result")
            
            # Should log warning
            mock_logger.return_value.warning.assert_called()
    
    def test_get_remaining_contigs(self):
        """Test getting remaining contigs."""
        tracker = ProgressTracker(total_contigs=5)
        
        all_contigs = ["contig1", "contig2", "contig3", "contig4", "contig5"]
        
        # Initially all contigs remaining
        remaining = tracker.get_remaining_contigs(all_contigs)
        assert len(remaining) == 5
        
        # Process some contigs
        tracker.start_contig("contig1")
        tracker.complete_contig("contig1", "result1")
        tracker.start_contig("contig3")
        tracker.fail_contig("contig3", "error")
        
        remaining = tracker.get_remaining_contigs(all_contigs)
        assert len(remaining) == 3
        assert "contig1" not in remaining
        assert "contig3" not in remaining
        assert "contig2" in remaining
        assert "contig4" in remaining
        assert "contig5" in remaining
    
    def test_get_results(self):
        """Test getting accumulated results."""
        tracker = ProgressTracker(total_contigs=3)
        
        # Initially empty
        results = tracker.get_results()
        assert len(results) == 0
        
        # Add some results
        tracker.start_contig("contig1")
        tracker.complete_contig("contig1", ["bp1", "bp2"])
        tracker.start_contig("contig2")
        tracker.complete_contig("contig2", ["bp3"])
        
        results = tracker.get_results()
        assert len(results) == 2
        assert results["contig1"] == ["bp1", "bp2"]
        assert results["contig2"] == ["bp3"]
    
    def test_get_summary(self):
        """Test getting progress summary."""
        tracker = ProgressTracker(total_contigs=5)
        
        # Process some contigs
        tracker.start_contig("contig1")
        tracker.complete_contig("contig1", "result1")
        tracker.start_contig("contig2")
        tracker.fail_contig("contig2", "error")
        tracker.start_contig("contig3")
        
        summary = tracker.get_summary()
        
        assert summary["total_contigs"] == 5
        assert summary["completed"] == 1
        assert summary["failed"] == 1
        assert summary["remaining"] == 3
        assert summary["success_rate"] == 0.5  # 1 success out of 2 processed
        assert summary["current_contig"] == "contig3"
        assert "elapsed_time" in summary
        assert "estimated_remaining" in summary
    
    def test_time_estimation(self):
        """Test time estimation functionality."""
        tracker = ProgressTracker(total_contigs=10)
        
        # Mock time to control elapsed time
        with patch('time.time') as mock_time:
            # Start at time 0
            mock_time.return_value = 0
            tracker.state.start_time = 0
            
            # Process 2 contigs in 10 seconds (5 seconds each)
            mock_time.return_value = 10
            tracker.state.processed_contigs = 2
            
            summary = tracker.get_summary()
            
            # Should estimate 5 seconds per remaining contig
            # 8 remaining * 5 seconds = 40 seconds
            assert summary["estimated_remaining"] == 40.0
    
    def test_format_time(self):
        """Test time formatting utility."""
        tracker = ProgressTracker(total_contigs=1)
        
        # Test seconds
        assert tracker._format_time(30) == "30s"
        
        # Test minutes
        assert tracker._format_time(90) == "1.5m"
        
        # Test hours
        assert tracker._format_time(7200) == "2.0h"


class TestCheckpointing:
    """Test checkpointing functionality."""
    
    def test_save_checkpoint(self, temp_dir):
        """Test saving checkpoint to file."""
        checkpoint_file = temp_dir / "checkpoint.json"
        tracker = ProgressTracker(total_contigs=3, checkpoint_file=checkpoint_file)
        
        # Process some contigs
        tracker.start_contig("contig1")
        tracker.complete_contig("contig1", "result1")
        
        # Save checkpoint
        tracker._save_checkpoint()
        
        # Verify file exists and contains expected data
        assert checkpoint_file.exists()
        
        with open(checkpoint_file) as f:
            data = json.load(f)
        
        assert "timestamp" in data
        assert "state" in data
        assert data["state"]["total_contigs"] == 3
        assert data["state"]["processed_contigs"] == 1
        assert "contig1" in data["state"]["completed_contigs"]
    
    def test_load_checkpoint(self, temp_dir):
        """Test loading checkpoint from file."""
        checkpoint_file = temp_dir / "checkpoint.json"
        
        # Create checkpoint data
        checkpoint_data = {
            "timestamp": "2024-01-01T00:00:00",
            "state": {
                "total_contigs": 5,
                "processed_contigs": 2,
                "current_contig": None,
                "start_time": 1000.0,
                "last_update": 1100.0,
                "completed_contigs": ["contig1", "contig2"],
                "failed_contigs": [],
                "partial_results": {"contig1": "result1", "contig2": "result2"}
            }
        }
        
        with open(checkpoint_file, 'w') as f:
            json.dump(checkpoint_data, f)
        
        # Create tracker with existing checkpoint
        tracker = ProgressTracker(total_contigs=5, checkpoint_file=checkpoint_file)
        
        # Should load existing state
        assert tracker.state.processed_contigs == 2
        assert len(tracker.state.completed_contigs) == 2
        assert "contig1" in tracker.state.completed_contigs
        assert "contig2" in tracker.state.completed_contigs
        assert tracker.state.partial_results["contig1"] == "result1"
    
    def test_load_invalid_checkpoint(self, temp_dir):
        """Test loading invalid checkpoint file."""
        checkpoint_file = temp_dir / "invalid_checkpoint.json"
        
        # Create invalid JSON
        with open(checkpoint_file, 'w') as f:
            f.write("invalid json content")
        
        with patch('chimeric_detector.utils.progress.logging.getLogger') as mock_logger:
            # Should handle invalid checkpoint gracefully
            tracker = ProgressTracker(total_contigs=3, checkpoint_file=checkpoint_file)
            
            # Should start fresh
            assert tracker.state.processed_contigs == 0
            # Should log warning
            mock_logger.return_value.warning.assert_called()
    
    def test_checkpoint_save_error(self, temp_dir):
        """Test handling of checkpoint save errors."""
        checkpoint_file = temp_dir / "readonly" / "checkpoint.json"
        # Don't create the parent directory to cause an error
        
        tracker = ProgressTracker(total_contigs=3, checkpoint_file=checkpoint_file)
        
        with patch('chimeric_detector.utils.progress.logging.getLogger') as mock_logger:
            # Should handle save error gracefully
            tracker._save_checkpoint()
            
            # Should log warning
            mock_logger.return_value.warning.assert_called()
    
    def test_automatic_checkpointing(self, temp_dir):
        """Test automatic checkpointing during processing."""
        checkpoint_file = temp_dir / "auto_checkpoint.json"
        tracker = ProgressTracker(total_contigs=15, checkpoint_file=checkpoint_file)
        
        # Process 10 contigs (should trigger checkpoint)
        for i in range(10):
            tracker.start_contig(f"contig{i}")
            tracker.complete_contig(f"contig{i}", f"result{i}")
        
        # Should have saved checkpoint
        assert checkpoint_file.exists()
        
        with open(checkpoint_file) as f:
            data = json.load(f)
        
        assert data["state"]["processed_contigs"] == 10
    
    def test_checkpoint_cleanup(self, temp_dir):
        """Test checkpoint cleanup after completion."""
        checkpoint_file = temp_dir / "cleanup_checkpoint.json"
        tracker = ProgressTracker(total_contigs=2, checkpoint_file=checkpoint_file)
        
        # Create checkpoint
        tracker.start_contig("contig1")
        tracker.complete_contig("contig1", "result1")
        tracker._save_checkpoint()
        assert checkpoint_file.exists()
        
        # Complete all processing
        tracker.start_contig("contig2")
        tracker.complete_contig("contig2", "result2")
        
        # Cleanup checkpoint
        tracker.cleanup_checkpoint()
        assert not checkpoint_file.exists()
    
    def test_checkpoint_cleanup_error(self, temp_dir):
        """Test handling of checkpoint cleanup errors."""
        checkpoint_file = temp_dir / "cleanup_error_checkpoint.json"
        tracker = ProgressTracker(total_contigs=1, checkpoint_file=checkpoint_file)
        
        # Create checkpoint
        checkpoint_file.write_text("{}")
        
        with patch.object(Path, 'unlink', side_effect=OSError("Permission denied")):
            with patch('chimeric_detector.utils.progress.logging.getLogger') as mock_logger:
                tracker.cleanup_checkpoint()
                
                # Should log warning but not raise
                mock_logger.return_value.warning.assert_called()


class TestProgressReporter:
    """Test ProgressReporter functionality."""
    
    def test_progress_reporter_creation(self):
        """Test ProgressReporter creation."""
        reporter = ProgressReporter(total=100, description="Test Progress")
        
        assert reporter.total == 100
        assert reporter.current == 0
        assert reporter.description == "Test Progress"
    
    def test_progress_reporter_updates(self):
        """Test ProgressReporter update functionality."""
        reporter = ProgressReporter(total=10, description="Test")
        
        with patch('chimeric_detector.utils.progress.logging.getLogger') as mock_logger:
            # Update progress
            reporter.update(1)
            assert reporter.current == 1
            
            # Should report after significant progress
            reporter.update(5)  # Now at 6/10
            assert reporter.current == 6
            
            # Should have logged progress
            mock_logger.return_value.info.assert_called()
    
    def test_progress_reporter_completion(self):
        """Test ProgressReporter completion."""
        reporter = ProgressReporter(total=5, description="Test")
        
        with patch('chimeric_detector.utils.progress.logging.getLogger') as mock_logger:
            # Complete all at once
            reporter.update(5)
            assert reporter.current == 5
            
            # Should report completion
            mock_logger.return_value.info.assert_called()
    
    def test_progress_reporter_large_total(self):
        """Test ProgressReporter with large total."""
        reporter = ProgressReporter(total=1000, description="Large Test")
        
        with patch('chimeric_detector.utils.progress.logging.getLogger') as mock_logger:
            # Update by small amounts
            for i in range(50):
                reporter.update(1)
            
            # Should have reported progress (every 10 items)
            assert mock_logger.return_value.info.call_count >= 4


class TestProgressIntegration:
    """Test progress tracking integration with main workflow."""
    
    def test_progress_tracking_resumption(self, temp_dir):
        """Test resuming progress from checkpoint."""
        checkpoint_file = temp_dir / "resume_checkpoint.json"
        
        # Create initial tracker and process some contigs
        tracker1 = ProgressTracker(total_contigs=5, checkpoint_file=checkpoint_file)
        tracker1.start_contig("contig1")
        tracker1.complete_contig("contig1", ["bp1"])
        tracker1.start_contig("contig2")
        tracker1.complete_contig("contig2", ["bp2"])
        tracker1._save_checkpoint()
        
        # Simulate restart with new tracker
        tracker2 = ProgressTracker(total_contigs=5, checkpoint_file=checkpoint_file)
        
        # Should have loaded previous state
        assert tracker2.state.processed_contigs == 2
        
        all_contigs = ["contig1", "contig2", "contig3", "contig4", "contig5"]
        remaining = tracker2.get_remaining_contigs(all_contigs)
        
        assert len(remaining) == 3
        assert "contig1" not in remaining
        assert "contig2" not in remaining
        
        # Get existing results
        results = tracker2.get_results()
        assert "contig1" in results
        assert "contig2" in results
    
    def test_progress_tracking_with_failures(self, temp_dir):
        """Test progress tracking with mixed success/failure."""
        checkpoint_file = temp_dir / "failure_checkpoint.json"
        tracker = ProgressTracker(total_contigs=5, checkpoint_file=checkpoint_file)
        
        # Mix of successes and failures
        tracker.start_contig("contig1")
        tracker.complete_contig("contig1", ["bp1"])
        
        tracker.start_contig("contig2")
        tracker.fail_contig("contig2", "Processing error")
        
        tracker.start_contig("contig3")
        tracker.complete_contig("contig3", ["bp3"])
        
        tracker.start_contig("contig4")
        tracker.fail_contig("contig4", "Another error")
        
        tracker.start_contig("contig5")
        tracker.complete_contig("contig5", ["bp5"])
        
        summary = tracker.get_summary()
        assert summary["completed"] == 3
        assert summary["failed"] == 2
        assert summary["success_rate"] == 0.6
        assert tracker.is_complete()
        
        # Results should only include successful contigs
        results = tracker.get_results()
        assert len(results) == 3
        assert "contig1" in results
        assert "contig3" in results
        assert "contig5" in results
        assert "contig2" not in results
        assert "contig4" not in results