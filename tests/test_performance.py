"""Performance benchmarks and stress tests."""

import pytest
import time
import numpy as np
from unittest.mock import MagicMock, patch

from chimeric_detector.core.window_analysis import WindowAnalyzer, RobustStatistics
from chimeric_detector.core.bam_parser import ReadPairInfo
from chimeric_detector.utils.progress import ProgressTracker, ProgressReporter
from chimeric_detector.utils.output import OutputFormatter
from chimeric_detector.config.config import Config


@pytest.mark.performance
class TestPerformanceBenchmarks:
    """Performance benchmarks for key operations."""
    
    def test_window_metrics_calculation_performance(self, sample_config, performance_timer, large_read_dataset):
        """Benchmark window metrics calculation performance."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        with performance_timer() as timer:
            for i in range(100):  # 100 windows
                chunk = large_read_dataset[i*100:(i+1)*100]  # 100 reads per window
                metrics = analyzer.calculate_window_metrics(chunk, coverage=15.0)
        
        elapsed = timer.elapsed
        assert elapsed is not None
        
        # Should process 100 windows in reasonable time
        assert elapsed < 5.0  # Less than 5 seconds
        
        # Calculate throughput
        windows_per_second = 100 / elapsed
        assert windows_per_second > 20  # At least 20 windows per second
    
    def test_statistical_method_performance(self, performance_timer):
        """Benchmark statistical method performance."""
        method = RobustStatistics()
        
        # Generate large datasets
        values = np.random.normal(300, 50, 1000)
        reference_values = np.random.normal(280, 45, 5000)
        
        with performance_timer() as timer:
            for _ in range(1000):
                score = method.calculate_anomaly_score(values, reference_values)
        
        elapsed = timer.elapsed
        assert elapsed is not None
        
        # Should handle 1000 calculations quickly
        assert elapsed < 2.0  # Less than 2 seconds
        
        calculations_per_second = 1000 / elapsed
        assert calculations_per_second > 500  # At least 500 calculations per second
    
    def test_large_breakpoint_detection_performance(self, sample_config, performance_timer):
        """Benchmark breakpoint detection with large datasets."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Create many window metrics
        windows = []
        for i in range(1000):
            windows.append(MagicMock(
                start=i*100, end=(i+1)*100,
                proper_pair_rate=0.8 + np.random.normal(0, 0.1),
                insert_size_median=300 + np.random.normal(0, 20),
                insert_size_mad=25,
                insert_size_zscore=0
            ))
        
        insert_size_stats = {"median": 300, "mad": 25}
        
        with performance_timer() as timer:
            breakpoints = analyzer.detect_breakpoints(windows, "large_contig", insert_size_stats)
        
        elapsed = timer.elapsed
        assert elapsed is not None
        
        # Should process 1000 windows quickly
        assert elapsed < 1.0  # Less than 1 second
        
        windows_per_second = 1000 / elapsed
        assert windows_per_second > 1000  # At least 1000 windows per second
    
    def test_output_generation_performance(self, temp_dir, sample_config, performance_timer):
        """Benchmark output generation performance."""
        from chimeric_detector.core.window_analysis import Breakpoint
        
        formatter = OutputFormatter(sample_config.output)
        
        # Create many breakpoints
        breakpoints = []
        for i in range(10000):
            breakpoints.append(Breakpoint(
                contig=f"contig_{i % 100}",
                position=i * 100,
                confidence=0.5 + np.random.random() * 0.5,
                evidence={"test": np.random.random()},
                left_metrics=None,
                right_metrics=None
            ))
        
        output_file = temp_dir / "performance_output.json"
        
        with performance_timer() as timer:
            formatter.write_results(breakpoints, output_file)
        
        elapsed = timer.elapsed
        assert elapsed is not None
        
        # Should write 10k breakpoints quickly
        assert elapsed < 5.0  # Less than 5 seconds
        
        breakpoints_per_second = 10000 / elapsed
        assert breakpoints_per_second > 2000  # At least 2000 breakpoints per second
        
        # Verify file was created correctly
        assert output_file.exists()
        assert output_file.stat().st_size > 1000  # Should be substantial
    
    def test_progress_tracking_performance(self, temp_dir, performance_timer):
        """Benchmark progress tracking overhead."""
        checkpoint_file = temp_dir / "perf_checkpoint.json"
        tracker = ProgressTracker(total_contigs=1000, checkpoint_file=checkpoint_file)
        
        with performance_timer() as timer:
            for i in range(1000):
                tracker.start_contig(f"contig_{i}")
                tracker.complete_contig(f"contig_{i}", f"result_{i}")
        
        elapsed = timer.elapsed
        assert elapsed is not None
        
        # Progress tracking should have minimal overhead
        assert elapsed < 2.0  # Less than 2 seconds
        
        operations_per_second = 2000 / elapsed  # 2 operations per contig
        assert operations_per_second > 1000  # At least 1000 operations per second


@pytest.mark.stress
class TestStressTests:
    """Stress tests for edge cases and resource limits."""
    
    def test_memory_usage_large_dataset(self, sample_config):
        """Test memory usage with large datasets."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Generate very large read pair dataset
        large_pairs = []
        for i in range(50000):  # 50k read pairs
            large_pairs.append(ReadPairInfo(
                position=i,
                mate_position=i + 300,
                insert_size=300 + int(np.random.normal(0, 50)),
                is_proper_pair=np.random.random() > 0.2,
                is_discordant=np.random.random() < 0.1,
                is_near_end=i < 1000 or i > 48000,
                contig="large_contig",
                mate_contig="large_contig"
            ))
        
        # Should handle large dataset without excessive memory usage
        metrics = analyzer.calculate_window_metrics(large_pairs, coverage=25.0)
        assert metrics is not None
        
        # Clean up
        del large_pairs
    
    def test_extreme_values_handling(self, sample_config):
        """Test handling of extreme values."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Create pairs with extreme values
        extreme_pairs = [
            ReadPairInfo(
                position=0, mate_position=1000000, insert_size=1000000,
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="test", mate_contig="test"
            ),
            ReadPairInfo(
                position=1, mate_position=2, insert_size=1,
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="test", mate_contig="test"
            ),
            ReadPairInfo(
                position=2, mate_position=3, insert_size=0,
                is_proper_pair=False, is_discordant=True, is_near_end=False,
                contig="test", mate_contig="other"
            )
        ]
        
        # Should handle extreme values gracefully
        metrics = analyzer.calculate_window_metrics(extreme_pairs, coverage=5.0)
        assert metrics is not None
        assert np.isfinite(metrics.insert_size_median)
        assert np.isfinite(metrics.insert_size_mad)
    
    def test_many_small_contigs_performance(self, sample_config, performance_timer):
        """Test performance with many small contigs."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        with performance_timer() as timer:
            total_breakpoints = 0
            
            # Process 1000 small contigs
            for contig_id in range(1000):
                # Small windows for each contig
                windows = [
                    MagicMock(
                        start=0, end=1000,
                        proper_pair_rate=0.8,
                        insert_size_median=300,
                        insert_size_mad=25,
                        insert_size_zscore=0.5
                    ),
                    MagicMock(
                        start=500, end=1500,
                        proper_pair_rate=0.4,  # Potential breakpoint
                        insert_size_median=400,
                        insert_size_mad=25,
                        insert_size_zscore=3.5
                    )
                ]
                
                breakpoints = analyzer.detect_breakpoints(
                    windows, f"contig_{contig_id}", {"median": 300, "mad": 25}
                )
                total_breakpoints += len(breakpoints)
        
        elapsed = timer.elapsed
        assert elapsed is not None
        
        # Should process many small contigs efficiently
        assert elapsed < 5.0  # Less than 5 seconds
        
        contigs_per_second = 1000 / elapsed
        assert contigs_per_second > 200  # At least 200 contigs per second
    
    def test_checkpoint_file_size_scaling(self, temp_dir):
        """Test checkpoint file size with large datasets."""
        checkpoint_file = temp_dir / "large_checkpoint.json"
        tracker = ProgressTracker(total_contigs=10000, checkpoint_file=checkpoint_file)
        
        # Process many contigs with results
        for i in range(5000):
            tracker.start_contig(f"contig_{i}")
            # Simulate some breakpoints
            results = [f"breakpoint_{j}" for j in range(10)]
            tracker.complete_contig(f"contig_{i}", results)
        
        # Force checkpoint save
        tracker._save_checkpoint()
        
        # Checkpoint file should exist but not be excessively large
        assert checkpoint_file.exists()
        file_size = checkpoint_file.stat().st_size
        
        # Should be reasonable size (less than 100MB for 5k contigs)
        assert file_size < 100 * 1024 * 1024
        
        # Should be able to load checkpoint quickly
        start_time = time.time()
        tracker2 = ProgressTracker(total_contigs=10000, checkpoint_file=checkpoint_file)
        load_time = time.time() - start_time
        
        assert load_time < 1.0  # Less than 1 second to load
        assert tracker2.state.processed_contigs == 5000


@pytest.mark.scaling
class TestScalingTests:
    """Test scaling behavior with different dataset sizes."""
    
    @pytest.mark.parametrize("dataset_size", [100, 1000, 10000])
    def test_window_processing_scaling(self, sample_config, dataset_size, performance_timer):
        """Test scaling of window processing with dataset size."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Generate dataset of specified size
        pairs = []
        for i in range(dataset_size):
            pairs.append(ReadPairInfo(
                position=i * 10,
                mate_position=i * 10 + 300,
                insert_size=300,
                is_proper_pair=True,
                is_discordant=False,
                is_near_end=False,
                contig="test",
                mate_contig="test"
            ))
        
        with performance_timer() as timer:
            metrics = analyzer.calculate_window_metrics(pairs, coverage=10.0)
        
        elapsed = timer.elapsed
        assert elapsed is not None
        assert metrics is not None
        
        # Calculate throughput
        reads_per_second = dataset_size / elapsed
        
        # Should scale reasonably (allow some degradation for larger datasets)
        if dataset_size == 100:
            assert reads_per_second > 10000
        elif dataset_size == 1000:
            assert reads_per_second > 5000
        else:  # 10000
            assert reads_per_second > 1000
    
    @pytest.mark.parametrize("num_breakpoints", [10, 100, 1000, 10000])
    def test_output_scaling(self, temp_dir, sample_config, num_breakpoints, performance_timer):
        """Test output generation scaling with number of breakpoints."""
        from chimeric_detector.core.window_analysis import Breakpoint
        
        formatter = OutputFormatter(sample_config.output)
        
        # Generate specified number of breakpoints
        breakpoints = []
        for i in range(num_breakpoints):
            breakpoints.append(Breakpoint(
                contig=f"contig_{i % 10}",
                position=i * 100,
                confidence=0.8,
                evidence={"test": 0.5},
                left_metrics=None,
                right_metrics=None
            ))
        
        output_file = temp_dir / f"scaling_{num_breakpoints}.json"
        
        with performance_timer() as timer:
            formatter.write_results(breakpoints, output_file)
        
        elapsed = timer.elapsed
        assert elapsed is not None
        
        # Calculate throughput
        breakpoints_per_second = num_breakpoints / elapsed
        
        # Should maintain reasonable throughput across scales
        if num_breakpoints <= 100:
            assert breakpoints_per_second > 5000
        elif num_breakpoints <= 1000:
            assert breakpoints_per_second > 2000
        else:  # 10000
            assert breakpoints_per_second > 1000
        
        # Verify output
        assert output_file.exists()


@pytest.mark.regression
class TestRegressionTests:
    """Regression tests to ensure performance doesn't degrade."""
    
    def test_baseline_performance_metrics(self, sample_config, performance_timer):
        """Establish baseline performance metrics."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Standard test dataset
        pairs = []
        for i in range(1000):
            pairs.append(ReadPairInfo(
                position=i * 10,
                mate_position=i * 10 + 300,
                insert_size=300 + int(np.random.normal(0, 50)),
                is_proper_pair=np.random.random() > 0.2,
                is_discordant=np.random.random() < 0.1,
                is_near_end=False,
                contig="baseline_test",
                mate_contig="baseline_test"
            ))
        
        # Measure window metrics calculation
        with performance_timer() as timer:
            metrics = analyzer.calculate_window_metrics(pairs, coverage=15.0)
        
        window_metrics_time = timer.elapsed
        assert window_metrics_time is not None
        assert metrics is not None
        
        # Baseline: should process 1000 reads in < 0.1 seconds
        assert window_metrics_time < 0.1
        
        # Measure breakpoint detection
        windows = [MagicMock(
            start=i*1000, end=(i+1)*1000,
            proper_pair_rate=0.8 + np.random.normal(0, 0.1),
            insert_size_median=300,
            insert_size_mad=25,
            insert_size_zscore=np.random.normal(0, 1)
        ) for i in range(100)]
        
        with performance_timer() as timer:
            breakpoints = analyzer.detect_breakpoints(
                windows, "baseline_test", {"median": 300, "mad": 25}
            )
        
        breakpoint_detection_time = timer.elapsed
        assert breakpoint_detection_time is not None
        
        # Baseline: should process 100 windows in < 0.05 seconds
        assert breakpoint_detection_time < 0.05
        
        # Record baseline metrics for future comparison
        baseline_metrics = {
            "window_metrics_time": window_metrics_time,
            "breakpoint_detection_time": breakpoint_detection_time,
            "total_time": window_metrics_time + breakpoint_detection_time
        }
        
        # These can be used for performance regression testing
        assert baseline_metrics["total_time"] < 0.15


class TestMemoryUsage:
    """Test memory usage patterns."""
    
    def test_memory_cleanup_after_processing(self, sample_config):
        """Test that memory is properly cleaned up after processing."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Process large dataset
        large_pairs = []
        for i in range(10000):
            large_pairs.append(ReadPairInfo(
                position=i, mate_position=i + 300, insert_size=300,
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="memory_test", mate_contig="memory_test"
            ))
        
        # Process and then clean up
        metrics = analyzer.calculate_window_metrics(large_pairs, coverage=10.0)
        assert metrics is not None
        
        # Explicit cleanup
        del large_pairs, metrics
        
        # Should be able to process another large dataset
        # (implicit test that memory was freed)
        large_pairs2 = []
        for i in range(10000):
            large_pairs2.append(ReadPairInfo(
                position=i, mate_position=i + 300, insert_size=300,
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="memory_test2", mate_contig="memory_test2"
            ))
        
        metrics2 = analyzer.calculate_window_metrics(large_pairs2, coverage=10.0)
        assert metrics2 is not None
    
    def test_streaming_vs_batch_processing(self, sample_config, performance_timer):
        """Compare streaming vs batch processing approaches."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Generate test data
        all_pairs = []
        for i in range(10000):
            all_pairs.append(ReadPairInfo(
                position=i, mate_position=i + 300, insert_size=300,
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="stream_test", mate_contig="stream_test"
            ))
        
        # Test batch processing
        with performance_timer() as timer:
            batch_metrics = analyzer.calculate_window_metrics(all_pairs, coverage=10.0)
        batch_time = timer.elapsed
        
        # Test streaming (simulated by processing in chunks)
        with performance_timer() as timer:
            chunk_results = []
            chunk_size = 1000
            for i in range(0, len(all_pairs), chunk_size):
                chunk = all_pairs[i:i+chunk_size]
                if len(chunk) >= analyzer.window_config.min_reads:
                    chunk_metrics = analyzer.calculate_window_metrics(chunk, coverage=10.0)
                    if chunk_metrics:
                        chunk_results.append(chunk_metrics)
        streaming_time = timer.elapsed
        
        # Both should complete successfully
        assert batch_metrics is not None
        assert len(chunk_results) > 0
        
        # Performance comparison (streaming might be slower due to overhead)
        assert batch_time is not None
        assert streaming_time is not None