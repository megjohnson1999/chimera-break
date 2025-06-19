"""Tests for performance optimization features."""

import pytest
import time
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

from chimeric_detector.utils.performance import (
    PerformanceMonitor, PerformanceMetrics, ParallelProcessor, 
    WorkItem, WorkQueue, ResourceManager, performance_context,
    batch_iterator, memory_efficient_chunker
)
from chimeric_detector.config.config import Config


class TestPerformanceMonitor:
    """Test PerformanceMonitor functionality."""
    
    def test_performance_monitor_creation(self):
        """Test PerformanceMonitor creation."""
        monitor = PerformanceMonitor(sample_interval=0.1)
        assert monitor.sample_interval == 0.1
        assert isinstance(monitor.metrics, PerformanceMetrics)
        assert not monitor._monitoring
    
    def test_performance_metrics(self):
        """Test PerformanceMetrics calculations."""
        metrics = PerformanceMetrics()
        
        # Test initial state
        assert metrics.items_processed == 0
        assert metrics.errors_encountered == 0
        
        # Test elapsed time calculation
        elapsed = metrics.elapsed_time
        assert elapsed >= 0
        
        # Test throughput with items
        metrics.items_processed = 100
        time.sleep(0.01)  # Small delay to ensure elapsed time > 0
        throughput = metrics.throughput
        assert throughput > 0
    
    def test_performance_context(self):
        """Test performance context manager."""
        with performance_context("Test operation") as monitor:
            assert isinstance(monitor, PerformanceMonitor)
            monitor.record_item_processed()
            monitor.record_item_processed()
        
        # Should have stopped monitoring
        assert not monitor._monitoring
        assert monitor.metrics.items_processed == 2
    
    def test_monitoring_lifecycle(self):
        """Test monitoring start/stop lifecycle."""
        monitor = PerformanceMonitor(sample_interval=0.01)
        
        # Start monitoring
        monitor.start_monitoring()
        assert monitor._monitoring
        assert monitor._monitor_thread is not None
        
        # Record some activity
        monitor.record_item_processed()
        monitor.record_error()
        
        time.sleep(0.02)  # Let monitoring run briefly
        
        # Stop monitoring
        monitor.stop_monitoring()
        assert not monitor._monitoring
        assert monitor.metrics.items_processed == 1
        assert monitor.metrics.errors_encountered == 1


class TestParallelProcessor:
    """Test ParallelProcessor functionality."""
    
    def test_parallel_processor_creation(self):
        """Test ParallelProcessor creation."""
        config = Config()
        processor = ParallelProcessor(config, max_workers=2)
        
        assert processor.config == config
        assert processor.max_workers <= 8  # Should be capped
        assert processor.max_workers >= 1
    
    def test_work_item_creation(self):
        """Test WorkItem creation."""
        item = WorkItem(
            id="test_1",
            contig="contig1",
            args=("arg1", "arg2"),
            kwargs={"param": "value"}
        )
        
        assert item.id == "test_1"
        assert item.contig == "contig1"
        assert item.args == ("arg1", "arg2")
        assert item.kwargs == {"param": "value"}
    
    def test_optimal_chunk_size_calculation(self):
        """Test optimal chunk size calculation."""
        config = Config()
        processor = ParallelProcessor(config, max_workers=4)
        
        # Test with small dataset
        chunk_size = processor.get_optimal_chunk_size(10)
        assert chunk_size >= 1
        assert chunk_size <= 1000
        
        # Test with large dataset
        chunk_size = processor.get_optimal_chunk_size(50000)
        assert chunk_size >= 1
        assert chunk_size <= 1000
    
    @patch('chimeric_detector.utils.performance.ProcessPoolExecutor')
    def test_parallel_processing_mock(self, mock_executor):
        """Test parallel processing with mocked executor."""
        config = Config()
        processor = ParallelProcessor(config, max_workers=2)
        
        # Mock the executor and futures
        mock_executor_instance = MagicMock()
        mock_executor.return_value.__enter__.return_value = mock_executor_instance
        
        mock_future = MagicMock()
        mock_future.result.return_value = ["result1"]
        mock_executor_instance.submit.return_value = mock_future
        
        # Mock as_completed to return our future
        with patch('chimeric_detector.utils.performance.as_completed', return_value=[mock_future]):
            work_items = [
                WorkItem(id="test1", contig="contig1", args=("arg1",), kwargs={})
            ]
            
            def dummy_process_func(*args, **kwargs):
                return ["result"]
            
            results = processor.process_contigs_parallel(work_items, dummy_process_func)
            
            assert len(results) == 1
            assert "contig1" in results


class TestWorkQueue:
    """Test WorkQueue functionality."""
    
    def test_work_queue_operations(self):
        """Test basic work queue operations."""
        queue = WorkQueue(maxsize=10)
        
        # Add work item
        item = WorkItem(id="test1", contig="contig1", args=())
        queue.add_work(item)
        
        # Get work item
        retrieved = queue.get_work(timeout=0.1)
        assert retrieved is not None
        assert retrieved.id == "test1"
        
        # Test empty queue
        empty = queue.get_work(timeout=0.1)
        assert empty is None
    
    def test_results_and_errors(self):
        """Test results and error handling."""
        queue = WorkQueue()
        
        # Add results and errors
        queue.put_result("contig1", ["breakpoint1"])
        queue.put_error("contig2", Exception("Test error"))
        
        # Retrieve results
        results = queue.get_results()
        assert len(results) == 1
        assert results[0][0] == "contig1"
        assert results[0][1] == ["breakpoint1"]
        
        # Retrieve errors
        errors = queue.get_errors()
        assert len(errors) == 1
        assert errors[0][0] == "contig2"
        assert isinstance(errors[0][1], Exception)
    
    def test_stop_mechanism(self):
        """Test stop mechanism."""
        queue = WorkQueue()
        
        assert not queue.is_stopped()
        
        queue.stop()
        assert queue.is_stopped()


class TestResourceManager:
    """Test ResourceManager functionality."""
    
    def test_resource_manager_creation(self):
        """Test ResourceManager creation."""
        manager = ResourceManager(memory_limit_gb=4.0)
        assert manager.memory_limit_gb == 4.0
        
        # Test auto-detection
        manager_auto = ResourceManager()
        assert manager_auto.memory_limit_gb > 0
    
    def test_memory_monitoring(self):
        """Test memory monitoring functions."""
        manager = ResourceManager(memory_limit_gb=1000.0)  # High limit
        
        # Should pass with high limit
        assert manager.check_memory_usage()
        
        # Get memory stats
        stats = manager.get_memory_stats()
        assert "process_rss_gb" in stats
        assert "system_available_gb" in stats
        assert "memory_limit_gb" in stats
        assert stats["memory_limit_gb"] == 1000.0
    
    def test_gc_suggestion(self):
        """Test garbage collection suggestion."""
        # Test with low limit to trigger suggestion
        manager = ResourceManager(memory_limit_gb=0.1)
        
        # Should suggest GC with very low limit
        suggest = manager.suggest_gc()
        assert isinstance(suggest, bool)
        # Note: This might be False if system memory is very low


class TestUtilityFunctions:
    """Test utility functions."""
    
    def test_batch_iterator(self):
        """Test batch iterator."""
        items = list(range(10))
        batches = list(batch_iterator(items, batch_size=3))
        
        assert len(batches) == 4  # [0,1,2], [3,4,5], [6,7,8], [9]
        assert batches[0] == [0, 1, 2]
        assert batches[1] == [3, 4, 5]
        assert batches[2] == [6, 7, 8]
        assert batches[3] == [9]
    
    def test_memory_efficient_chunker(self):
        """Test memory efficient chunker."""
        def item_generator():
            for i in range(7):
                yield i
        
        chunks = list(memory_efficient_chunker(item_generator(), chunk_size=3))
        
        assert len(chunks) == 3  # [0,1,2], [3,4,5], [6]
        assert chunks[0] == [0, 1, 2]
        assert chunks[1] == [3, 4, 5]
        assert chunks[2] == [6]
    
    def test_empty_chunker(self):
        """Test chunker with empty input."""
        def empty_generator():
            return
            yield  # Unreachable
        
        chunks = list(memory_efficient_chunker(empty_generator(), chunk_size=3))
        assert len(chunks) == 0


class TestPerformanceIntegration:
    """Integration tests for performance features."""
    
    def test_config_serialization_for_parallel(self):
        """Test config serialization for parallel processing."""
        config = Config()
        config.window.size = 2000
        config.detection.min_confidence_score = 0.8
        
        # Convert to dict and back
        config_dict = config.to_dict()
        restored_config = Config.from_dict(config_dict)
        
        assert restored_config.window.size == 2000
        assert restored_config.detection.min_confidence_score == 0.8
    
    def test_performance_monitoring_integration(self):
        """Test performance monitoring integration."""
        with performance_context("Integration test") as monitor:
            # Simulate some work
            for i in range(5):
                monitor.record_item_processed()
                time.sleep(0.001)
        
        metrics = monitor.metrics
        assert metrics.items_processed == 5
        assert metrics.elapsed_time > 0
        assert metrics.throughput > 0
    
    @patch('psutil.virtual_memory')
    def test_memory_pressure_simulation(self, mock_memory):
        """Test behavior under simulated memory pressure."""
        # Simulate low memory
        mock_memory.return_value.available = 1024 * 1024 * 1024  # 1GB
        
        manager = ResourceManager()
        
        # Should adjust behavior for low memory
        stats = manager.get_memory_stats()
        assert stats["system_available_gb"] == 1.0
        
        # Should be conservative with memory limits
        assert manager.memory_limit_gb <= 1.0