"""Sliding window analysis for breakpoint detection."""

import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Callable
from scipy import stats
import logging

from ..config.config import WindowConfig, DetectionConfig
from .bam_parser import ReadPairInfo


@dataclass
class WindowMetrics:
    """Metrics calculated for each window."""
    start: int
    end: int
    coverage: float
    read_count: int
    proper_pair_rate: float
    discordant_pair_rate: float
    insert_size_median: float
    insert_size_mad: float
    insert_size_zscore: float = 0.0
    

@dataclass
class Breakpoint:
    """Detected breakpoint with confidence score."""
    contig: str
    position: int
    confidence: float
    evidence: Dict[str, float]
    left_metrics: WindowMetrics
    right_metrics: WindowMetrics
    

class StatisticalMethod:
    """Base class for pluggable statistical methods."""
    
    def calculate_anomaly_score(self, values: np.ndarray, reference_values: np.ndarray) -> float:
        """Calculate anomaly score for a set of values."""
        raise NotImplementedError
        

class RobustStatistics(StatisticalMethod):
    """Robust statistical method using median and MAD."""
    
    def calculate_anomaly_score(self, values: np.ndarray, reference_values: np.ndarray) -> float:
        if len(values) == 0 or len(reference_values) == 0:
            return 0.0
            
        median_ref = np.median(reference_values)
        mad_ref = np.median(np.abs(reference_values - median_ref))
        
        if mad_ref == 0:
            mad_ref = 1.0
            
        current_median = np.median(values)
        z_score = (current_median - median_ref) / (1.4826 * mad_ref)
        
        return abs(z_score)
        

class ParametricStatistics(StatisticalMethod):
    """Traditional parametric statistics (assumes normality)."""
    
    def calculate_anomaly_score(self, values: np.ndarray, reference_values: np.ndarray) -> float:
        if len(values) < 3 or len(reference_values) < 3:
            return 0.0
            
        _, p_value = stats.ttest_ind(values, reference_values)
        return -np.log10(p_value + 1e-10)
        

class NonParametricStatistics(StatisticalMethod):
    """Non-parametric statistical method."""
    
    def calculate_anomaly_score(self, values: np.ndarray, reference_values: np.ndarray) -> float:
        if len(values) < 3 or len(reference_values) < 3:
            return 0.0
            
        _, p_value = stats.mannwhitneyu(values, reference_values, alternative='two-sided')
        return -np.log10(p_value + 1e-10)


class WindowAnalyzer:
    """Performs sliding window analysis for breakpoint detection."""
    
    def __init__(self, window_config: WindowConfig, detection_config: DetectionConfig, quality_config = None):
        self.window_config = window_config
        self.detection_config = detection_config
        self.quality_config = quality_config
        self.logger = logging.getLogger(__name__)
        
        # Select statistical method
        self.stat_method = self._get_statistical_method(detection_config.statistical_method)
        
    def _get_statistical_method(self, method_name: str) -> StatisticalMethod:
        """Get the appropriate statistical method."""
        methods = {
            "robust": RobustStatistics(),
            "parametric": ParametricStatistics(),
            "nonparametric": NonParametricStatistics()
        }
        return methods.get(method_name, RobustStatistics())
    
    def calculate_window_metrics(self, pairs: List[ReadPairInfo], coverage: float) -> Optional[WindowMetrics]:
        """Calculate metrics for a single window."""
        if len(pairs) < self.window_config.min_reads:
            return None
            
        if coverage < self.window_config.min_coverage:
            return None
            
        # Separate end reads from internal reads for different treatment
        end_reads = [p for p in pairs if p.is_near_end]
        internal_reads = [p for p in pairs if not p.is_near_end]
        
        # Calculate metrics, being more lenient with end reads
        proper_pairs = 0
        discordant_pairs = 0
        
        # For internal reads, use standard classification
        for p in internal_reads:
            if p.is_proper_pair:
                proper_pairs += 1
            elif p.is_discordant:
                discordant_pairs += 1
                
        # For end reads, be more lenient about proper pair classification if configured
        for p in end_reads:
            if self.quality_config and self.quality_config.lenient_end_reads:
                if p.is_proper_pair or (not p.is_discordant and not p.mate_unmapped):
                    # Consider non-discordant end reads as "acceptable" pairs
                    proper_pairs += 1
                elif p.is_discordant or p.mate_unmapped:
                    discordant_pairs += 1
            else:
                # Use standard classification for end reads too
                if p.is_proper_pair:
                    proper_pairs += 1
                elif p.is_discordant:
                    discordant_pairs += 1
        
        proper_pair_rate = proper_pairs / len(pairs) if len(pairs) > 0 else 0
        discordant_pair_rate = discordant_pairs / len(pairs) if len(pairs) > 0 else 0
        
        # Calculate insert size statistics, excluding end reads with unreliable sizes
        reliable_insert_sizes = []
        for p in pairs:
            if p.is_proper_pair and p.insert_size > 0:
                # For end reads, only include if insert size seems reasonable
                if p.is_near_end:
                    if p.insert_size < 2000:  # Conservative threshold for end reads
                        reliable_insert_sizes.append(p.insert_size)
                else:
                    reliable_insert_sizes.append(p.insert_size)
        
        if len(reliable_insert_sizes) > 0:
            insert_size_median = np.median(reliable_insert_sizes)
            insert_size_mad = np.median(np.abs(reliable_insert_sizes - insert_size_median))
        else:
            insert_size_median = 0
            insert_size_mad = 0
            
        return WindowMetrics(
            start=pairs[0].position if pairs else 0,
            end=pairs[-1].position if pairs else 0,
            coverage=coverage,
            read_count=len(pairs),
            proper_pair_rate=proper_pair_rate,
            discordant_pair_rate=discordant_pair_rate,
            insert_size_median=insert_size_median,
            insert_size_mad=insert_size_mad
        )
    
    def detect_breakpoints(self, windows: List[WindowMetrics], contig: str, 
                          insert_size_stats: Dict[str, float]) -> List[Breakpoint]:
        """Detect breakpoints from window metrics."""
        breakpoints = []
        
        if len(windows) < 3:
            return breakpoints
            
        # Calculate insert size z-scores
        ref_median = insert_size_stats.get("median", 500)
        ref_mad = insert_size_stats.get("mad", 100)
        
        for window in windows:
            if ref_mad > 0:
                window.insert_size_zscore = abs(window.insert_size_median - ref_median) / (1.4826 * ref_mad)
        
        # Scan for breakpoints
        for i in range(1, len(windows) - 1):
            prev_window = windows[i-1]
            curr_window = windows[i]
            next_window = windows[i+1]
            
            # Check for proper pair rate drop
            if prev_window.proper_pair_rate > 0:
                pair_rate_drop = (prev_window.proper_pair_rate - curr_window.proper_pair_rate) / prev_window.proper_pair_rate
            else:
                pair_rate_drop = 0
                
            # Check for insert size anomaly
            insert_size_anomaly = curr_window.insert_size_zscore > self.detection_config.insert_size_zscore_threshold
            
            # Check for high discordant pairs
            high_discordant = curr_window.discordant_pair_rate > self.detection_config.discordant_pair_threshold
            
            # Calculate confidence score
            evidence = {
                "proper_pair_drop": pair_rate_drop,
                "insert_size_zscore": curr_window.insert_size_zscore,
                "discordant_rate": curr_window.discordant_pair_rate
            }
            
            confidence = self._calculate_confidence(evidence)
            
            if confidence >= self.detection_config.min_confidence_score:
                breakpoint = Breakpoint(
                    contig=contig,
                    position=(curr_window.start + curr_window.end) // 2,
                    confidence=confidence,
                    evidence=evidence,
                    left_metrics=prev_window,
                    right_metrics=next_window
                )
                breakpoints.append(breakpoint)
                
        return self._merge_nearby_breakpoints(breakpoints)
    
    def _calculate_confidence(self, evidence: Dict[str, float]) -> float:
        """Calculate confidence score from multiple evidence sources."""
        scores = []
        
        # Proper pair drop score
        if evidence["proper_pair_drop"] >= self.detection_config.min_proper_pair_drop:
            scores.append(min(evidence["proper_pair_drop"] / 0.5, 1.0))
            
        # Insert size anomaly score
        if evidence["insert_size_zscore"] > self.detection_config.insert_size_zscore_threshold:
            scores.append(min(evidence["insert_size_zscore"] / 5.0, 1.0))
            
        # Discordant pair score
        if evidence["discordant_rate"] > self.detection_config.discordant_pair_threshold:
            scores.append(min(evidence["discordant_rate"] / 0.4, 1.0))
            
        if not scores:
            return 0.0
            
        # Weighted average with bonus for multiple evidence types
        base_score = np.mean(scores)
        evidence_bonus = 0.1 * (len(scores) - 1)
        
        return min(base_score + evidence_bonus, 1.0)
    
    def _merge_nearby_breakpoints(self, breakpoints: List[Breakpoint], 
                                 merge_distance: int = 500) -> List[Breakpoint]:
        """Merge breakpoints that are close together."""
        if not breakpoints:
            return breakpoints
            
        # Sort by position
        breakpoints.sort(key=lambda x: x.position)
        
        merged = []
        current = breakpoints[0]
        
        for bp in breakpoints[1:]:
            if bp.position - current.position <= merge_distance:
                # Merge: keep the one with higher confidence
                if bp.confidence > current.confidence:
                    current = bp
            else:
                merged.append(current)
                current = bp
                
        merged.append(current)
        return merged