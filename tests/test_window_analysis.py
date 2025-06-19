"""Tests for window analysis with edge cases."""

import pytest
import numpy as np
from unittest.mock import Mock, patch, MagicMock

from chimeric_detector.core.window_analysis import (
    WindowAnalyzer, WindowMetrics, Breakpoint, StatisticalMethod,
    RobustStatistics, ParametricStatistics, NonParametricStatistics
)
from chimeric_detector.config.config import WindowConfig, DetectionConfig, QualityConfig
from chimeric_detector.core.bam_parser import ReadPairInfo


class TestStatisticalMethods:
    """Test statistical method implementations."""
    
    def test_robust_statistics_normal_case(self):
        """Test robust statistics with normal data."""
        method = RobustStatistics()
        
        values = np.array([1, 2, 3, 4, 5])
        reference_values = np.array([2, 3, 4, 5, 6])
        
        score = method.calculate_anomaly_score(values, reference_values)
        assert isinstance(score, float)
        assert score >= 0
    
    def test_robust_statistics_empty_arrays(self):
        """Test robust statistics with empty arrays."""
        method = RobustStatistics()
        
        # Empty values
        score = method.calculate_anomaly_score(np.array([]), np.array([1, 2, 3]))
        assert score == 0.0
        
        # Empty reference
        score = method.calculate_anomaly_score(np.array([1, 2, 3]), np.array([]))
        assert score == 0.0
        
        # Both empty
        score = method.calculate_anomaly_score(np.array([]), np.array([]))
        assert score == 0.0
    
    def test_robust_statistics_identical_values(self):
        """Test robust statistics with identical values."""
        method = RobustStatistics()
        
        # All values identical in reference
        values = np.array([1, 2, 3])
        reference_values = np.array([5, 5, 5, 5])  # All identical
        
        score = method.calculate_anomaly_score(values, reference_values)
        # Should handle zero MAD case
        assert score == float('inf')  # Expected behavior for varying current vs constant reference
        
        # All values identical in both
        values = np.array([5, 5, 5])
        reference_values = np.array([5, 5, 5, 5])
        
        score = method.calculate_anomaly_score(values, reference_values)
        assert score == 0.0  # Both constant, no difference
    
    def test_parametric_statistics_normal_case(self):
        """Test parametric statistics with normal data."""
        method = ParametricStatistics()
        
        values = np.array([1, 2, 3, 4, 5])
        reference_values = np.array([2, 3, 4, 5, 6])
        
        score = method.calculate_anomaly_score(values, reference_values)
        assert isinstance(score, float)
        assert score >= 0
    
    def test_parametric_statistics_insufficient_data(self):
        """Test parametric statistics with insufficient data."""
        method = ParametricStatistics()
        
        # Too few values
        values = np.array([1, 2])
        reference_values = np.array([2, 3, 4, 5, 6])
        
        score = method.calculate_anomaly_score(values, reference_values)
        assert score == 0.0
        
        # Too few reference values
        values = np.array([1, 2, 3, 4, 5])
        reference_values = np.array([2, 3])
        
        score = method.calculate_anomaly_score(values, reference_values)
        assert score == 0.0
    
    def test_nonparametric_statistics_normal_case(self):
        """Test non-parametric statistics with normal data."""
        method = NonParametricStatistics()
        
        values = np.array([1, 2, 3, 4, 5])
        reference_values = np.array([6, 7, 8, 9, 10])
        
        score = method.calculate_anomaly_score(values, reference_values)
        assert isinstance(score, float)
        assert score >= 0
    
    def test_nonparametric_statistics_insufficient_data(self):
        """Test non-parametric statistics with insufficient data."""
        method = NonParametricStatistics()
        
        # Too few values
        values = np.array([1, 2])
        reference_values = np.array([2, 3, 4, 5, 6])
        
        score = method.calculate_anomaly_score(values, reference_values)
        assert score == 0.0


class TestWindowAnalyzerInitialization:
    """Test WindowAnalyzer initialization."""
    
    def test_analyzer_creation(self, sample_config):
        """Test WindowAnalyzer creation."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        assert analyzer.window_config == sample_config.window
        assert analyzer.detection_config == sample_config.detection
        assert analyzer.quality_config == sample_config.quality
        assert analyzer.stat_method is not None
    
    def test_statistical_method_selection(self, sample_config):
        """Test statistical method selection."""
        # Test robust method (default)
        config = sample_config.detection
        config.statistical_method = "robust"
        analyzer = WindowAnalyzer(sample_config.window, config)
        assert isinstance(analyzer.stat_method, RobustStatistics)
        
        # Test parametric method
        config.statistical_method = "parametric"
        analyzer = WindowAnalyzer(sample_config.window, config)
        assert isinstance(analyzer.stat_method, ParametricStatistics)
        
        # Test non-parametric method
        config.statistical_method = "nonparametric"
        analyzer = WindowAnalyzer(sample_config.window, config)
        assert isinstance(analyzer.stat_method, NonParametricStatistics)
        
        # Test unknown method (should default to robust)
        config.statistical_method = "unknown"
        analyzer = WindowAnalyzer(sample_config.window, config)
        assert isinstance(analyzer.stat_method, RobustStatistics)


class TestWindowMetricsCalculation:
    """Test window metrics calculation."""
    
    def test_calculate_window_metrics_normal_case(self, sample_config, sample_read_pairs):
        """Test normal window metrics calculation."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        metrics = analyzer.calculate_window_metrics(sample_read_pairs, coverage=10.0)
        
        assert metrics is not None
        assert isinstance(metrics, WindowMetrics)
        assert metrics.coverage == 10.0
        assert metrics.read_count == len(sample_read_pairs)
        assert 0 <= metrics.proper_pair_rate <= 1
        assert 0 <= metrics.discordant_pair_rate <= 1
        assert metrics.insert_size_median > 0
    
    def test_calculate_window_metrics_empty_pairs(self, sample_config):
        """Test window metrics calculation with empty pairs."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        metrics = analyzer.calculate_window_metrics([], coverage=10.0)
        assert metrics is None
    
    def test_calculate_window_metrics_insufficient_reads(self, sample_config):
        """Test window metrics calculation with insufficient reads."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Create minimal read pairs (less than min_reads)
        minimal_pairs = [
            ReadPairInfo(
                position=100, mate_position=200, insert_size=150,
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="test", mate_contig="test"
            )
        ]
        
        # Assuming min_reads > 1 in config
        if sample_config.window.min_reads > 1:
            metrics = analyzer.calculate_window_metrics(minimal_pairs, coverage=10.0)
            assert metrics is None
    
    def test_calculate_window_metrics_low_coverage(self, sample_config, sample_read_pairs):
        """Test window metrics calculation with low coverage."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Use coverage below min_coverage threshold
        low_coverage = sample_config.window.min_coverage - 1
        metrics = analyzer.calculate_window_metrics(sample_read_pairs, coverage=low_coverage)
        assert metrics is None
    
    def test_calculate_window_metrics_all_end_reads(self, sample_config):
        """Test window metrics calculation with all end reads."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Create pairs that are all near ends
        end_pairs = [
            ReadPairInfo(
                position=i * 10, mate_position=i * 10 + 200, insert_size=300,
                is_proper_pair=True, is_discordant=False, is_near_end=True,
                contig="test", mate_contig="test"
            ) for i in range(20)
        ]
        
        metrics = analyzer.calculate_window_metrics(end_pairs, coverage=10.0)
        
        # Should handle end reads differently
        assert metrics is not None
        assert metrics.read_count == len(end_pairs)
    
    def test_calculate_window_metrics_no_reliable_insert_sizes(self, sample_config):
        """Test window metrics calculation with no reliable insert sizes."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Create pairs with unreliable insert sizes (all end reads with large inserts)
        unreliable_pairs = [
            ReadPairInfo(
                position=i * 10, mate_position=i * 10 + 5000, insert_size=5000,
                is_proper_pair=True, is_discordant=False, is_near_end=True,
                contig="test", mate_contig="test"
            ) for i in range(20)
        ]
        
        metrics = analyzer.calculate_window_metrics(unreliable_pairs, coverage=10.0)
        
        assert metrics is not None
        # Should handle case with no reliable insert sizes
        assert metrics.insert_size_median == 0
        assert metrics.insert_size_mad == 0


class TestBreakpointDetection:
    """Test breakpoint detection algorithms."""
    
    def test_detect_breakpoints_normal_case(self, sample_config, sample_window_metrics):
        """Test normal breakpoint detection."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        insert_size_stats = {"median": 180, "mad": 15}
        
        breakpoints = analyzer.detect_breakpoints(
            sample_window_metrics, "test_contig", insert_size_stats
        )
        
        assert isinstance(breakpoints, list)
        # Should detect breakpoint in second window (low proper pair rate, high z-score)
        assert len(breakpoints) >= 1
        
        for bp in breakpoints:
            assert isinstance(bp, Breakpoint)
            assert bp.contig == "test_contig"
            assert 0 <= bp.confidence <= 1
            assert isinstance(bp.evidence, dict)
    
    def test_detect_breakpoints_empty_windows(self, sample_config):
        """Test breakpoint detection with empty windows."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        insert_size_stats = {"median": 180, "mad": 15}
        
        breakpoints = analyzer.detect_breakpoints([], "test_contig", insert_size_stats)
        assert len(breakpoints) == 0
    
    def test_detect_breakpoints_insufficient_windows(self, sample_config):
        """Test breakpoint detection with insufficient windows."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Only one window
        single_window = [
            WindowMetrics(
                start=0, end=1000, coverage=10.0, read_count=20,
                proper_pair_rate=0.8, discordant_pair_rate=0.1,
                insert_size_median=180, insert_size_mad=15, insert_size_zscore=0.5
            )
        ]
        
        insert_size_stats = {"median": 180, "mad": 15}
        
        breakpoints = analyzer.detect_breakpoints(
            single_window, "test_contig", insert_size_stats
        )
        assert len(breakpoints) == 0
    
    def test_detect_breakpoints_zero_mad(self, sample_config, sample_window_metrics):
        """Test breakpoint detection with zero MAD in insert size stats."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Zero MAD means all insert sizes are identical
        insert_size_stats = {"median": 180, "mad": 0}
        
        breakpoints = analyzer.detect_breakpoints(
            sample_window_metrics, "test_contig", insert_size_stats
        )
        
        # Should handle zero MAD case
        assert isinstance(breakpoints, list)
        # Some windows should have infinite z-scores
        for window in sample_window_metrics:
            if window.insert_size_median != 180:
                assert window.insert_size_zscore == float('inf')
    
    def test_detect_breakpoints_zero_proper_pair_rate(self, sample_config):
        """Test breakpoint detection with zero proper pair rate."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Create windows with zero proper pair rate
        windows = [
            WindowMetrics(
                start=0, end=1000, coverage=10.0, read_count=20,
                proper_pair_rate=0.0, discordant_pair_rate=0.9,
                insert_size_median=180, insert_size_mad=15, insert_size_zscore=0.5
            ),
            WindowMetrics(
                start=100, end=1100, coverage=10.0, read_count=20,
                proper_pair_rate=0.0, discordant_pair_rate=0.9,
                insert_size_median=180, insert_size_mad=15, insert_size_zscore=0.5
            ),
            WindowMetrics(
                start=200, end=1200, coverage=10.0, read_count=20,
                proper_pair_rate=0.0, discordant_pair_rate=0.9,
                insert_size_median=180, insert_size_mad=15, insert_size_zscore=0.5
            )
        ]
        
        insert_size_stats = {"median": 180, "mad": 15}
        
        breakpoints = analyzer.detect_breakpoints(windows, "test_contig", insert_size_stats)
        
        # Should handle zero proper pair rate gracefully
        assert isinstance(breakpoints, list)


class TestBreakpointMerging:
    """Test breakpoint merging functionality."""
    
    def test_merge_nearby_breakpoints_normal_case(self, sample_config):
        """Test normal breakpoint merging."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Create nearby breakpoints
        breakpoints = [
            Breakpoint(
                contig="test", position=1000, confidence=0.7,
                evidence={}, left_metrics=None, right_metrics=None
            ),
            Breakpoint(
                contig="test", position=1100, confidence=0.8,  # Higher confidence
                evidence={}, left_metrics=None, right_metrics=None
            ),
            Breakpoint(
                contig="test", position=2000, confidence=0.9,  # Far away
                evidence={}, left_metrics=None, right_metrics=None
            )
        ]
        
        merged = analyzer._merge_nearby_breakpoints(breakpoints, merge_distance=500)
        
        # Should merge first two, keep third separate
        assert len(merged) == 2
        # Should keep the one with higher confidence
        assert any(bp.position == 1100 and bp.confidence == 0.8 for bp in merged)
        assert any(bp.position == 2000 and bp.confidence == 0.9 for bp in merged)
    
    def test_merge_nearby_breakpoints_empty_list(self, sample_config):
        """Test breakpoint merging with empty list."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        merged = analyzer._merge_nearby_breakpoints([])
        assert len(merged) == 0
    
    def test_merge_nearby_breakpoints_single_breakpoint(self, sample_config):
        """Test breakpoint merging with single breakpoint."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        breakpoints = [
            Breakpoint(
                contig="test", position=1000, confidence=0.8,
                evidence={}, left_metrics=None, right_metrics=None
            )
        ]
        
        merged = analyzer._merge_nearby_breakpoints(breakpoints)
        assert len(merged) == 1
        assert merged[0].position == 1000
    
    def test_merge_nearby_breakpoints_all_close(self, sample_config):
        """Test breakpoint merging when all breakpoints are close."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # All breakpoints within merge distance
        breakpoints = [
            Breakpoint(
                contig="test", position=1000, confidence=0.7,
                evidence={}, left_metrics=None, right_metrics=None
            ),
            Breakpoint(
                contig="test", position=1200, confidence=0.8,
                evidence={}, left_metrics=None, right_metrics=None
            ),
            Breakpoint(
                contig="test", position=1400, confidence=0.9,  # Highest confidence
                evidence={}, left_metrics=None, right_metrics=None
            )
        ]
        
        merged = analyzer._merge_nearby_breakpoints(breakpoints, merge_distance=500)
        
        # Should merge all into one
        assert len(merged) == 1
        # Should keep the one with highest confidence
        assert merged[0].position == 1400
        assert merged[0].confidence == 0.9


class TestConfidenceScoring:
    """Test confidence scoring calculations."""
    
    def test_calculate_confidence_normal_case(self, sample_config):
        """Test normal confidence calculation."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        evidence = {
            "proper_pair_drop": 0.5,
            "insert_size_zscore": 3.0,
            "discordant_rate": 0.3
        }
        
        confidence = analyzer._calculate_confidence(evidence)
        
        assert 0 <= confidence <= 1
        # Should be relatively high given the evidence
        assert confidence > 0.5
    
    def test_calculate_confidence_no_evidence(self, sample_config):
        """Test confidence calculation with no evidence."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        evidence = {}
        confidence = analyzer._calculate_confidence(evidence)
        assert confidence == 0.0
    
    def test_calculate_confidence_weak_evidence(self, sample_config):
        """Test confidence calculation with weak evidence."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        evidence = {
            "proper_pair_drop": 0.1,  # Low drop
            "insert_size_zscore": 1.0,  # Low z-score
            "discordant_rate": 0.05  # Low discordant rate
        }
        
        confidence = analyzer._calculate_confidence(evidence)
        
        assert 0 <= confidence <= 1
        # Should be relatively low given weak evidence
        assert confidence < 0.5
    
    def test_calculate_confidence_strong_evidence(self, sample_config):
        """Test confidence calculation with strong evidence."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        evidence = {
            "proper_pair_drop": 0.8,  # High drop
            "insert_size_zscore": 5.0,  # High z-score
            "discordant_rate": 0.5  # High discordant rate
        }
        
        confidence = analyzer._calculate_confidence(evidence)
        
        assert 0 <= confidence <= 1
        # Should be high given strong evidence
        assert confidence > 0.8
    
    def test_calculate_confidence_capped_at_one(self, sample_config):
        """Test confidence calculation is capped at 1.0."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        evidence = {
            "proper_pair_drop": 1.0,  # Maximum drop
            "insert_size_zscore": 10.0,  # Very high z-score
            "discordant_rate": 1.0,  # Maximum discordant rate
            "extra_evidence": 2.0  # Additional evidence
        }
        
        confidence = analyzer._calculate_confidence(evidence)
        
        assert confidence <= 1.0
        assert confidence == 1.0  # Should be capped


class TestEdgeCases:
    """Test various edge cases in window analysis."""
    
    def test_window_metrics_with_nan_values(self, sample_config):
        """Test window metrics handling NaN values."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        # Create pairs with some invalid insert sizes
        pairs = [
            ReadPairInfo(
                position=100, mate_position=200, insert_size=float('nan'),
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="test", mate_contig="test"
            ),
            ReadPairInfo(
                position=150, mate_position=250, insert_size=180,
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="test", mate_contig="test"
            )
        ]
        
        metrics = analyzer.calculate_window_metrics(pairs, coverage=10.0)
        
        # Should handle NaN values gracefully
        assert metrics is not None
        assert not np.isnan(metrics.insert_size_median)
        assert not np.isnan(metrics.insert_size_mad)
    
    def test_window_metrics_with_infinite_values(self, sample_config):
        """Test window metrics handling infinite values."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        pairs = [
            ReadPairInfo(
                position=100, mate_position=200, insert_size=float('inf'),
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="test", mate_contig="test"
            ),
            ReadPairInfo(
                position=150, mate_position=250, insert_size=180,
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="test", mate_contig="test"
            )
        ]
        
        metrics = analyzer.calculate_window_metrics(pairs, coverage=10.0)
        
        # Should handle infinite values gracefully
        assert metrics is not None
        assert np.isfinite(metrics.insert_size_median)
        assert np.isfinite(metrics.insert_size_mad)
    
    def test_very_large_insert_sizes(self, sample_config):
        """Test handling of very large insert sizes."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        pairs = [
            ReadPairInfo(
                position=100, mate_position=200, insert_size=1000000,  # Very large
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="test", mate_contig="test"
            ),
            ReadPairInfo(
                position=150, mate_position=250, insert_size=180,
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="test", mate_contig="test"
            )
        ]
        
        metrics = analyzer.calculate_window_metrics(pairs, coverage=10.0)
        
        # Should handle large values
        assert metrics is not None
        assert metrics.insert_size_median > 0
    
    def test_negative_insert_sizes(self, sample_config):
        """Test handling of negative insert sizes."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        pairs = [
            ReadPairInfo(
                position=100, mate_position=200, insert_size=-180,  # Negative
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="test", mate_contig="test"
            ),
            ReadPairInfo(
                position=150, mate_position=250, insert_size=180,
                is_proper_pair=True, is_discordant=False, is_near_end=False,
                contig="test", mate_contig="test"
            )
        ]
        
        metrics = analyzer.calculate_window_metrics(pairs, coverage=10.0)
        
        # Should handle negative values
        assert metrics is not None
        assert metrics.insert_size_median > 0  # Should use absolute values
    
    def test_all_discordant_pairs(self, sample_config):
        """Test window with all discordant pairs."""
        analyzer = WindowAnalyzer(
            sample_config.window,
            sample_config.detection,
            sample_config.quality
        )
        
        pairs = [
            ReadPairInfo(
                position=i * 10, mate_position=i * 10 + 200, insert_size=300,
                is_proper_pair=False, is_discordant=True, is_near_end=False,
                contig="test", mate_contig="other"
            ) for i in range(20)
        ]
        
        metrics = analyzer.calculate_window_metrics(pairs, coverage=10.0)
        
        assert metrics is not None
        assert metrics.proper_pair_rate == 0.0
        assert metrics.discordant_pair_rate == 1.0