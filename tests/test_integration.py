"""Integration tests for the main workflow."""

import pytest
import json
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
import sys

import detect_chimeras
from chimeric_detector.config.config import Config
from chimeric_detector.core.bam_parser import BAMParser
from chimeric_detector.core.window_analysis import WindowAnalyzer, Breakpoint
from chimeric_detector.utils.output import OutputFormatter
from chimeric_detector.utils.progress import ProgressTracker


class TestMainWorkflow:
    """Test the main chimera detection workflow."""
    
    @patch('detect_chimeras.BAMParser')
    @patch('detect_chimeras.WindowAnalyzer')
    @patch('detect_chimeras.OutputFormatter')
    def test_detect_chimeras_in_contig_success(self, mock_formatter, mock_analyzer, mock_bam_parser, sample_config):
        """Test successful chimera detection in a single contig."""
        # Setup mocks
        mock_bam_instance = MagicMock()
        mock_bam_parser.return_value.__enter__.return_value = mock_bam_instance
        mock_bam_instance.estimate_insert_size_distribution.return_value = {
            "median": 300, "mad": 50, "mean": 310, "std": 75
        }
        mock_bam_instance.iterate_windows.return_value = [
            (0, 1000, []),  # Empty pairs for simplicity
            (500, 1500, []),
            (1000, 2000, [])
        ]
        mock_bam_instance.calculate_coverage_in_window.return_value = 10.0
        
        mock_analyzer_instance = MagicMock()
        mock_analyzer.return_value = mock_analyzer_instance
        mock_analyzer_instance.calculate_window_metrics.return_value = MagicMock()
        mock_analyzer_instance.detect_breakpoints.return_value = [
            Breakpoint(
                contig="test_contig", position=750, confidence=0.85,
                evidence={"proper_pair_drop": 0.4}, 
                left_metrics=None, right_metrics=None
            )
        ]
        
        # Test the function
        breakpoints = detect_chimeras.detect_chimeras_in_contig(
            mock_bam_instance, mock_analyzer_instance, "test_contig", sample_config
        )
        
        assert len(breakpoints) == 1
        assert breakpoints[0].contig == "test_contig"
        assert breakpoints[0].position == 750
        assert breakpoints[0].confidence == 0.85
    
    @patch('detect_chimeras.BAMParser')
    def test_detect_chimeras_in_contig_insert_size_error(self, mock_bam_parser, sample_config):
        """Test chimera detection when insert size estimation fails."""
        # Setup mocks
        mock_bam_instance = MagicMock()
        mock_bam_parser.return_value.__enter__.return_value = mock_bam_instance
        mock_bam_instance.estimate_insert_size_distribution.side_effect = Exception("Insert size error")
        mock_bam_instance.iterate_windows.return_value = []
        
        mock_analyzer = MagicMock()
        
        with patch('detect_chimeras.logging.getLogger') as mock_logger:
            breakpoints = detect_chimeras.detect_chimeras_in_contig(
                mock_bam_instance, mock_analyzer, "test_contig", sample_config
            )
            
            # Should return empty list and log error
            assert len(breakpoints) == 0
            mock_logger.return_value.error.assert_called()
    
    @patch('detect_chimeras.BAMParser')
    def test_detect_chimeras_in_contig_window_iteration_error(self, mock_bam_parser, sample_config):
        """Test chimera detection when window iteration fails."""
        # Setup mocks
        mock_bam_instance = MagicMock()
        mock_bam_parser.return_value.__enter__.return_value = mock_bam_instance
        mock_bam_instance.estimate_insert_size_distribution.return_value = {"median": 300, "mad": 50}
        mock_bam_instance.iterate_windows.side_effect = Exception("Window iteration error")
        
        mock_analyzer = MagicMock()
        
        with patch('detect_chimeras.logging.getLogger') as mock_logger:
            breakpoints = detect_chimeras.detect_chimeras_in_contig(
                mock_bam_instance, mock_analyzer, "test_contig", sample_config
            )
            
            # Should return empty list and log error
            assert len(breakpoints) == 0
            mock_logger.return_value.error.assert_called()
    
    @patch('detect_chimeras.BAMParser')
    def test_detect_chimeras_in_contig_partial_window_failures(self, mock_bam_parser, sample_config):
        """Test chimera detection with partial window processing failures."""
        # Setup mocks
        mock_bam_instance = MagicMock()
        mock_bam_parser.return_value.__enter__.return_value = mock_bam_instance
        mock_bam_instance.estimate_insert_size_distribution.return_value = {"median": 300, "mad": 50}
        mock_bam_instance.iterate_windows.return_value = [
            (0, 1000, []),
            (500, 1500, []),
            (1000, 2000, [])
        ]
        
        # Make coverage calculation fail for some windows
        coverage_results = [10.0, Exception("Coverage error"), 8.0]
        mock_bam_instance.calculate_coverage_in_window.side_effect = coverage_results
        
        mock_analyzer = MagicMock()
        mock_analyzer.calculate_window_metrics.return_value = MagicMock()
        mock_analyzer.detect_breakpoints.return_value = []
        
        with patch('detect_chimeras.logging.getLogger') as mock_logger:
            breakpoints = detect_chimeras.detect_chimeras_in_contig(
                mock_bam_instance, mock_analyzer, "test_contig", sample_config
            )
            
            # Should process successfully and log warnings for failed windows
            assert isinstance(breakpoints, list)
            mock_logger.return_value.warning.assert_called()


class TestMainFunction:
    """Test the main function and command line interface."""
    
    @patch('detect_chimeras.validate_bam_file')
    @patch('detect_chimeras.validate_output_path')
    @patch('detect_chimeras.BAMParser')
    @patch('detect_chimeras.OutputFormatter')
    @patch('detect_chimeras.Logger.setup_logging')
    @patch('sys.argv', ['detect_chimeras.py', 'test.bam', 'output.json'])
    def test_main_success(self, mock_logging, mock_formatter, mock_bam_parser, 
                         mock_validate_output, mock_validate_bam):
        """Test successful main function execution."""
        # Setup path validation
        mock_validate_bam.return_value = Path("test.bam")
        mock_validate_output.return_value = Path("output.json")
        
        # Setup BAM parser
        mock_bam_instance = MagicMock()
        mock_bam_parser.return_value.__enter__.return_value = mock_bam_instance
        mock_bam_instance.get_contigs.return_value = ["contig1", "contig2"]
        
        # Setup formatter
        mock_formatter_instance = MagicMock()
        mock_formatter.return_value = mock_formatter_instance
        
        # Mock the contig processing function to avoid complex setup
        with patch('detect_chimeras.detect_chimeras_in_contig') as mock_detect:
            mock_detect.return_value = []
            
            # Should not raise any exceptions
            detect_chimeras.main()
            
            # Verify key operations were called
            mock_logging.assert_called_once()
            mock_formatter_instance.write_results.assert_called_once()
    
    @patch('detect_chimeras.validate_bam_file')
    @patch('sys.argv', ['detect_chimeras.py', '../../../etc/passwd', 'output.json'])
    def test_main_security_error(self, mock_validate_bam):
        """Test main function with security error."""
        from chimeric_detector.utils.security import SecurityError
        
        mock_validate_bam.side_effect = SecurityError("Path traversal detected")
        
        with pytest.raises(SystemExit):
            detect_chimeras.main()
    
    @patch('detect_chimeras.validate_bam_file')
    @patch('detect_chimeras.validate_output_path') 
    @patch('detect_chimeras.Config.from_file')
    @patch('sys.argv', ['detect_chimeras.py', 'test.bam', 'output.json', '-c', 'config.yaml'])
    def test_main_with_config_file(self, mock_config_from_file, mock_validate_output, mock_validate_bam):
        """Test main function with configuration file."""
        from chimeric_detector.utils.security import validate_config_file
        
        # Setup mocks
        mock_validate_bam.return_value = Path("test.bam")
        mock_validate_output.return_value = Path("output.json")
        
        mock_config = MagicMock()
        mock_config_from_file.return_value = mock_config
        
        with patch('detect_chimeras.validate_config_file') as mock_validate_config:
            mock_validate_config.return_value = Path("config.yaml")
            
            with patch('detect_chimeras.BAMParser'), \
                 patch('detect_chimeras.OutputFormatter'), \
                 patch('detect_chimeras.Logger.setup_logging'), \
                 patch('detect_chimeras.detect_chimeras_in_contig') as mock_detect:
                
                mock_detect.return_value = []
                
                # Should load config from file
                detect_chimeras.main()
                mock_config_from_file.assert_called_once()
    
    @patch('sys.argv', ['detect_chimeras.py', 'test.bam', 'output.json', '--window-size', '50'])
    def test_main_invalid_numeric_parameter(self):
        """Test main function with invalid numeric parameter."""
        with pytest.raises(SystemExit):
            detect_chimeras.main()
    
    @patch('sys.argv', ['detect_chimeras.py', 'test.bam', 'output.json', '--contigs', 'con/tig'])
    def test_main_invalid_contig_name(self):
        """Test main function with invalid contig name."""
        with pytest.raises(SystemExit):
            detect_chimeras.main()


class TestProgressTracking:
    """Test progress tracking integration."""
    
    @patch('detect_chimeras.validate_bam_file')
    @patch('detect_chimeras.validate_output_path')
    @patch('detect_chimeras.BAMParser')
    @patch('detect_chimeras.OutputFormatter')
    @patch('detect_chimeras.Logger.setup_logging')
    @patch('sys.argv', ['detect_chimeras.py', 'test.bam', 'output.json', '--checkpoint', 'checkpoint.json'])
    def test_main_with_checkpoint(self, mock_logging, mock_formatter, mock_bam_parser,
                                 mock_validate_output, mock_validate_bam, temp_dir):
        """Test main function with checkpoint functionality."""
        # Setup path validation
        mock_validate_bam.return_value = Path("test.bam")
        mock_validate_output.return_value = Path("output.json")
        
        # Setup BAM parser
        mock_bam_instance = MagicMock()
        mock_bam_parser.return_value.__enter__.return_value = mock_bam_instance
        mock_bam_instance.get_contigs.return_value = ["contig1", "contig2", "contig3"]
        
        # Setup formatter
        mock_formatter_instance = MagicMock()
        mock_formatter.return_value = mock_formatter_instance
        
        checkpoint_file = temp_dir / "checkpoint.json"
        
        with patch('detect_chimeras.detect_chimeras_in_contig') as mock_detect:
            mock_detect.return_value = [
                Breakpoint(
                    contig="test", position=100, confidence=0.8,
                    evidence={}, left_metrics=None, right_metrics=None
                )
            ]
            
            with patch('detect_chimeras.Path') as mock_path:
                mock_path.return_value = checkpoint_file
                
                # Should not raise any exceptions
                detect_chimeras.main()
                
                # Progress tracking should be used
                assert mock_detect.call_count == 3  # Three contigs
    
    def test_progress_tracker_functionality(self, temp_dir):
        """Test ProgressTracker functionality in isolation."""
        checkpoint_file = temp_dir / "test_checkpoint.json"
        
        # Test creating and using progress tracker
        tracker = ProgressTracker(total_contigs=3, checkpoint_file=checkpoint_file)
        
        # Start processing contigs
        tracker.start_contig("contig1")
        tracker.complete_contig("contig1", ["result1"])
        
        tracker.start_contig("contig2")
        tracker.fail_contig("contig2", "Test error")
        
        tracker.start_contig("contig3")
        tracker.complete_contig("contig3", ["result3"])
        
        # Check final state
        assert tracker.is_complete()
        
        summary = tracker.get_summary()
        assert summary["completed"] == 2
        assert summary["failed"] == 1
        assert summary["success_rate"] == 2/3
        
        results = tracker.get_results()
        assert "contig1" in results
        assert "contig3" in results
        assert "contig2" not in results


class TestValidationPipeline:
    """Test validation pipeline integration."""
    
    @patch('detect_chimeras.validate_bam_file')
    @patch('detect_chimeras.validate_output_path')
    @patch('detect_chimeras.validate_assembly_file')
    @patch('detect_chimeras.BAMParser')
    @patch('detect_chimeras.ValidationPipeline')
    @patch('detect_chimeras.OutputFormatter')
    @patch('detect_chimeras.Logger.setup_logging')
    @patch('sys.argv', ['detect_chimeras.py', 'test.bam', 'output.json', 
                        '--assembly', 'assembly.fasta', '--validate-taxonomy'])
    def test_main_with_validation(self, mock_logging, mock_formatter, mock_validation_pipeline,
                                 mock_bam_parser, mock_validate_assembly, mock_validate_output,
                                 mock_validate_bam):
        """Test main function with validation enabled."""
        # Setup path validation
        mock_validate_bam.return_value = Path("test.bam")
        mock_validate_output.return_value = Path("output.json")
        mock_validate_assembly.return_value = Path("assembly.fasta")
        
        # Setup BAM parser
        mock_bam_instance = MagicMock()
        mock_bam_parser.return_value.__enter__.return_value = mock_bam_instance
        mock_bam_instance.get_contigs.return_value = ["contig1"]
        
        # Setup validation pipeline
        mock_pipeline_instance = MagicMock()
        mock_validation_pipeline.return_value = mock_pipeline_instance
        mock_pipeline_instance.run_validation.return_value = {"validation": "results"}
        
        # Setup formatter
        mock_formatter_instance = MagicMock()
        mock_formatter.return_value = mock_formatter_instance
        
        with patch('detect_chimeras.detect_chimeras_in_contig') as mock_detect:
            mock_detect.return_value = []
            
            # Mock config to enable validation
            with patch('detect_chimeras.Config') as mock_config_class:
                mock_config = MagicMock()
                mock_config.validation.enabled = True
                mock_config.validation.taxonomy_db = "dummy"
                mock_config_class.return_value = mock_config
                
                # Should not raise any exceptions
                detect_chimeras.main()
                
                # Validation should be called
                mock_pipeline_instance.run_validation.assert_called_once()


class TestOutputGeneration:
    """Test output generation integration."""
    
    def test_output_formatter_integration(self, temp_dir, sample_breakpoints):
        """Test OutputFormatter integration with real files."""
        from chimeric_detector.config.config import OutputConfig
        
        # Test JSON output
        json_config = OutputConfig(format="json", include_read_evidence=True)
        formatter = OutputFormatter(json_config)
        
        json_output = temp_dir / "test_output.json"
        metadata = {"test": "metadata"}
        
        formatter.write_results(sample_breakpoints, json_output, metadata)
        
        # Verify file was created and contains expected content
        assert json_output.exists()
        
        with open(json_output) as f:
            data = json.load(f)
        
        assert "metadata" in data
        assert "breakpoints" in data
        assert len(data["breakpoints"]) == len(sample_breakpoints)
        assert data["metadata"]["test"] == "metadata"
        
        # Test TSV output
        tsv_config = OutputConfig(format="tsv", include_read_evidence=False)
        formatter = OutputFormatter(tsv_config)
        
        tsv_output = temp_dir / "test_output.tsv"
        formatter.write_results(sample_breakpoints, tsv_output)
        
        assert tsv_output.exists()
        
        # Check TSV format
        lines = tsv_output.read_text().strip().split('\n')
        assert len(lines) == len(sample_breakpoints) + 1  # +1 for header
        assert 'contig' in lines[0]  # Header
        assert 'test_contig' in lines[1]  # Data
    
    def test_atomic_write_integration(self, temp_dir, sample_breakpoints):
        """Test atomic write functionality."""
        from chimeric_detector.config.config import OutputConfig
        
        formatter = OutputFormatter(OutputConfig(format="json"))
        output_file = temp_dir / "atomic_test.json"
        
        # Simulate interruption during write by patching shutil.move
        with patch('chimeric_detector.utils.output.shutil.move') as mock_move:
            mock_move.side_effect = Exception("Simulated interruption")
            
            with pytest.raises(RuntimeError):
                formatter.write_results(sample_breakpoints, output_file)
            
            # Original file should not exist (atomic operation failed)
            assert not output_file.exists()
        
        # Normal operation should work
        formatter.write_results(sample_breakpoints, output_file)
        assert output_file.exists()


class TestErrorPropagation:
    """Test error propagation through the system."""
    
    @patch('detect_chimeras.validate_bam_file')
    @patch('detect_chimeras.validate_output_path')
    @patch('detect_chimeras.BAMParser')
    @patch('detect_chimeras.Logger.setup_logging')
    @patch('sys.argv', ['detect_chimeras.py', 'test.bam', 'output.json', '--debug'])
    def test_debug_mode_error_propagation(self, mock_logging, mock_bam_parser,
                                         mock_validate_output, mock_validate_bam):
        """Test that errors are propagated in debug mode."""
        # Setup mocks
        mock_validate_bam.return_value = Path("test.bam")
        mock_validate_output.return_value = Path("output.json")
        
        mock_bam_instance = MagicMock()
        mock_bam_parser.return_value.__enter__.return_value = mock_bam_instance
        mock_bam_instance.get_contigs.return_value = ["contig1"]
        
        with patch('detect_chimeras.detect_chimeras_in_contig') as mock_detect:
            mock_detect.side_effect = Exception("Test error")
            
            # Mock config to enable debug
            with patch('detect_chimeras.Config') as mock_config_class:
                mock_config = MagicMock()
                mock_config.output.debug = True
                mock_config_class.return_value = mock_config
                
                # Should raise the exception in debug mode
                with pytest.raises(Exception, match="Test error"):
                    detect_chimeras.main()
    
    @patch('detect_chimeras.validate_bam_file')
    @patch('detect_chimeras.validate_output_path')
    @patch('detect_chimeras.BAMParser')
    @patch('detect_chimeras.Logger.setup_logging')
    @patch('sys.argv', ['detect_chimeras.py', 'test.bam', 'output.json'])
    def test_normal_mode_error_handling(self, mock_logging, mock_bam_parser,
                                       mock_validate_output, mock_validate_bam):
        """Test that errors are handled gracefully in normal mode."""
        # Setup mocks
        mock_validate_bam.return_value = Path("test.bam")
        mock_validate_output.return_value = Path("output.json")
        
        mock_bam_instance = MagicMock()
        mock_bam_parser.return_value.__enter__.return_value = mock_bam_instance
        mock_bam_instance.get_contigs.return_value = ["contig1"]
        
        mock_formatter = MagicMock()
        
        with patch('detect_chimeras.detect_chimeras_in_contig') as mock_detect, \
             patch('detect_chimeras.OutputFormatter') as mock_formatter_class:
            
            mock_detect.side_effect = Exception("Test error")
            mock_formatter_class.return_value = mock_formatter
            
            # Mock config to disable debug
            with patch('detect_chimeras.Config') as mock_config_class:
                mock_config = MagicMock()
                mock_config.output.debug = False
                mock_config_class.return_value = mock_config
                
                # Should handle error gracefully and continue
                detect_chimeras.main()
                
                # Should still write output (empty results)
                mock_formatter.write_results.assert_called_once()


class TestCommandLineInterface:
    """Test command line interface parsing and validation."""
    
    def test_parse_arguments_minimal(self):
        """Test argument parsing with minimal arguments."""
        with patch('sys.argv', ['detect_chimeras.py', 'input.bam', 'output.json']):
            args = detect_chimeras.parse_arguments()
            
            assert args.bam_file == 'input.bam'
            assert args.output == 'output.json'
            assert args.config is None
            assert args.contigs is None
            assert not args.verbose
            assert not args.debug
    
    def test_parse_arguments_full(self):
        """Test argument parsing with all arguments."""
        argv = [
            'detect_chimeras.py', 'input.bam', 'output.json',
            '--config', 'config.yaml',
            '--contigs', 'contig1', 'contig2',
            '--window-size', '2000',
            '--window-step', '200',
            '--min-coverage', '5.0',
            '--min-proper-pair-drop', '0.4',
            '--insert-zscore-threshold', '3.0',
            '--min-confidence', '0.7',
            '--output-format', 'tsv',
            '--include-evidence',
            '--assembly', 'assembly.fasta',
            '--validate-taxonomy',
            '--evaluate-splits',
            '--verbose',
            '--debug',
            '--log-file', 'debug.log',
            '--checkpoint', 'checkpoint.json'
        ]
        
        with patch('sys.argv', argv):
            args = detect_chimeras.parse_arguments()
            
            assert args.bam_file == 'input.bam'
            assert args.output == 'output.json'
            assert args.config == 'config.yaml'
            assert args.contigs == ['contig1', 'contig2']
            assert args.window_size == 2000
            assert args.window_step == 200
            assert args.min_coverage == 5.0
            assert args.min_proper_pair_drop == 0.4
            assert args.insert_zscore_threshold == 3.0
            assert args.min_confidence == 0.7
            assert args.output_format == 'tsv'
            assert args.include_evidence
            assert args.assembly == 'assembly.fasta'
            assert args.validate_taxonomy
            assert args.evaluate_splits
            assert args.verbose
            assert args.debug
            assert args.log_file == 'debug.log'
            assert args.checkpoint == 'checkpoint.json'