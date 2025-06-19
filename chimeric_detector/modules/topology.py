"""Topology detection module for identifying circular genomes and other structural features."""

import numpy as np
import logging
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
from collections import defaultdict

from ..core.window_analysis import WindowMetrics, Breakpoint
from ..core.bam_parser import BAMParser
from ..config.config import Config


@dataclass
class TopologyAnalysis:
    """Results of topology analysis for a contig."""
    contig: str
    likely_circular: bool
    confidence: float
    evidence: Dict[str, Any]
    breakpoint_adjustments: Dict[int, float]  # position -> confidence_penalty


class CircularGenomeDetector:
    """Detects circular genomes that create systematic insert size artifacts."""
    
    def __init__(self, config: Config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
    def analyze_contig_topology(self, contig: str, window_metrics: List[WindowMetrics], 
                               breakpoints: List[Breakpoint], bam_parser: BAMParser) -> TopologyAnalysis:
        """Analyze a contig for circular topology indicators."""
        
        if len(window_metrics) < 5:  # Need minimum windows for analysis
            return TopologyAnalysis(
                contig=contig,
                likely_circular=False,
                confidence=0.0,
                evidence={"insufficient_data": True},
                breakpoint_adjustments={}
            )
        
        evidence = {}
        circular_indicators = []
        
        # 1. Check for systematic insert size anomalies
        insert_size_pattern = self._analyze_insert_size_pattern(window_metrics)
        evidence["insert_size_analysis"] = insert_size_pattern
        circular_indicators.append(insert_size_pattern["circular_score"])
        
        # 2. Check breakpoint distribution
        breakpoint_pattern = self._analyze_breakpoint_distribution(breakpoints, window_metrics)
        evidence["breakpoint_distribution"] = breakpoint_pattern
        circular_indicators.append(breakpoint_pattern["circular_score"])
        
        # 3. Check for terminal read patterns (DTR-like)
        terminal_pattern = self._analyze_terminal_patterns(contig, bam_parser)
        evidence["terminal_analysis"] = terminal_pattern
        circular_indicators.append(terminal_pattern["circular_score"])
        
        # 4. Overall assessment
        overall_score = np.mean(circular_indicators)
        likely_circular = overall_score > 0.6  # Threshold for circular detection
        
        # 5. Calculate breakpoint confidence adjustments
        adjustments = {}
        if likely_circular:
            # Penalize all breakpoints if genome appears circular
            penalty = min(0.8, overall_score)  # Cap penalty at 0.8
            for bp in breakpoints:
                adjustments[bp.position] = -penalty
                
        evidence["overall_score"] = overall_score
        evidence["circular_indicators"] = circular_indicators
        
        return TopologyAnalysis(
            contig=contig,
            likely_circular=likely_circular,
            confidence=overall_score,
            evidence=evidence,
            breakpoint_adjustments=adjustments
        )
    
    def _analyze_insert_size_pattern(self, window_metrics: List[WindowMetrics]) -> Dict[str, Any]:
        """Analyze insert size patterns for circular indicators."""
        
        # Extract insert size z-scores from windows
        z_scores = []
        positions = []
        
        for wm in window_metrics:
            if hasattr(wm, 'insert_size_zscore') and wm.insert_size_zscore is not None:
                z_scores.append(abs(wm.insert_size_zscore))
                positions.append((wm.start + wm.end) / 2)
        
        if len(z_scores) < 3:
            return {"circular_score": 0.0, "pattern": "insufficient_data"}
        
        z_scores = np.array(z_scores)
        positions = np.array(positions)
        
        # Check for systematic elevation of z-scores
        high_zscore_fraction = np.mean(z_scores > 2.0)
        
        # Check for even distribution (not clustered)
        # Calculate coefficient of variation of positions with high z-scores
        high_zscore_positions = positions[z_scores > 2.0]
        
        if len(high_zscore_positions) < 3:
            distribution_score = 0.0
        else:
            # Check if high z-scores are evenly distributed across contig
            pos_range = np.max(positions) - np.min(positions)
            high_pos_range = np.max(high_zscore_positions) - np.min(high_zscore_positions)
            distribution_score = high_pos_range / pos_range if pos_range > 0 else 0.0
        
        # Circular score based on systematic and distributed anomalies
        circular_score = min(1.0, (high_zscore_fraction * 2.0) * distribution_score)
        
        return {
            "circular_score": circular_score,
            "high_zscore_fraction": high_zscore_fraction,
            "distribution_score": distribution_score,
            "mean_zscore": np.mean(z_scores),
            "pattern": "systematic_distributed" if circular_score > 0.5 else "clustered_or_sparse"
        }
    
    def _analyze_breakpoint_distribution(self, breakpoints: List[Breakpoint], 
                                       window_metrics: List[WindowMetrics]) -> Dict[str, Any]:
        """Analyze breakpoint distribution patterns."""
        
        if len(breakpoints) < 5:  # Need multiple breakpoints for pattern analysis
            return {"circular_score": 0.0, "pattern": "too_few_breakpoints"}
        
        positions = [bp.position for bp in breakpoints]
        confidences = [bp.confidence for bp in breakpoints]
        
        # Get contig span from window metrics
        if window_metrics:
            contig_start = min(wm.start for wm in window_metrics)
            contig_end = max(wm.end for wm in window_metrics)
            contig_length = contig_end - contig_start
        else:
            return {"circular_score": 0.0, "pattern": "no_window_data"}
        
        # Check for even distribution across contig
        position_gaps = np.diff(sorted(positions))
        gap_cv = np.std(position_gaps) / np.mean(position_gaps) if len(position_gaps) > 0 else 0
        
        # Even distribution = low coefficient of variation
        distribution_evenness = max(0, 1.0 - gap_cv)
        
        # Check for high density of breakpoints
        breakpoint_density = len(breakpoints) / (contig_length / 1000)  # per kb
        density_score = min(1.0, breakpoint_density / 10.0)  # Normalize (10+ per kb = max score)
        
        # Circular genomes often have many evenly distributed "breakpoints"
        circular_score = min(1.0, distribution_evenness * density_score)
        
        return {
            "circular_score": circular_score,
            "distribution_evenness": distribution_evenness,
            "density_score": density_score,
            "breakpoint_density": breakpoint_density,
            "gap_coefficient_variation": gap_cv,
            "pattern": "evenly_distributed_dense" if circular_score > 0.6 else "clustered_or_sparse"
        }
    
    def _analyze_terminal_patterns(self, contig: str, bam_parser: BAMParser) -> Dict[str, Any]:
        """Analyze read patterns at contig termini for DTR-like signatures."""
        
        try:
            contig_length = bam_parser._bam.get_reference_length(contig)
            terminal_region = 1000  # Check first/last 1kb
            
            # Get reads from terminal regions
            start_reads = list(bam_parser.stream_read_pairs_in_window(
                contig, 0, min(terminal_region, contig_length)
            ))
            
            end_reads = list(bam_parser.stream_read_pairs_in_window(
                contig, max(0, contig_length - terminal_region), contig_length
            ))
            
            # Look for reads with mates mapping to opposite terminus
            # (suggesting circular structure)
            cross_terminal_pairs = 0
            total_terminal_reads = len(start_reads) + len(end_reads)
            
            for read_pair in start_reads:
                if (read_pair.mate_position is not None and 
                    read_pair.mate_position > contig_length - terminal_region):
                    cross_terminal_pairs += 1
                    
            for read_pair in end_reads:
                if (read_pair.mate_position is not None and 
                    read_pair.mate_position < terminal_region):
                    cross_terminal_pairs += 1
            
            # Calculate score based on cross-terminal pairs
            if total_terminal_reads > 0:
                cross_terminal_fraction = cross_terminal_pairs / total_terminal_reads
                circular_score = min(1.0, cross_terminal_fraction * 10.0)  # Scale up
            else:
                circular_score = 0.0
                cross_terminal_fraction = 0.0
            
            return {
                "circular_score": circular_score,
                "cross_terminal_pairs": cross_terminal_pairs,
                "cross_terminal_fraction": cross_terminal_fraction,
                "total_terminal_reads": total_terminal_reads,
                "contig_length": contig_length,
                "pattern": "DTR_like" if circular_score > 0.3 else "linear_like"
            }
            
        except Exception as e:
            self.logger.warning(f"Error analyzing terminal patterns for {contig}: {e}")
            return {
                "circular_score": 0.0,
                "error": str(e),
                "pattern": "analysis_failed"
            }


class TopologyValidator:
    """Validation module that adjusts breakpoint confidence based on topology."""
    
    def __init__(self, config: Config):
        self.config = config
        self.detector = CircularGenomeDetector(config)
        self.logger = logging.getLogger(__name__)
    
    def validate_with_topology(self, breakpoints: List[Breakpoint], 
                              window_metrics_by_contig: Dict[str, List[WindowMetrics]],
                              bam_parser: BAMParser) -> Dict[str, Any]:
        """Validate breakpoints considering contig topology."""
        
        results = {}
        topology_analyses = {}
        adjusted_breakpoints = []
        
        # Group breakpoints by contig
        breakpoints_by_contig = defaultdict(list)
        for bp in breakpoints:
            breakpoints_by_contig[bp.contig].append(bp)
        
        # Analyze each contig
        for contig, contig_breakpoints in breakpoints_by_contig.items():
            window_metrics = window_metrics_by_contig.get(contig, [])
            
            # Run topology analysis
            topology = self.detector.analyze_contig_topology(
                contig, window_metrics, contig_breakpoints, bam_parser
            )
            
            topology_analyses[contig] = topology
            
            # Adjust breakpoint confidences
            for bp in contig_breakpoints:
                if bp.position in topology.breakpoint_adjustments:
                    adjustment = topology.breakpoint_adjustments[bp.position]
                    new_confidence = max(0.0, bp.confidence + adjustment)
                    
                    # Create adjusted breakpoint
                    adjusted_bp = Breakpoint(
                        contig=bp.contig,
                        position=bp.position,
                        confidence=new_confidence,
                        evidence={
                            **bp.evidence,
                            "topology_adjustment": adjustment,
                            "original_confidence": bp.confidence,
                            "topology_analysis": topology.likely_circular
                        },
                        left_metrics=bp.left_metrics,
                        right_metrics=bp.right_metrics
                    )
                    adjusted_breakpoints.append(adjusted_bp)
                else:
                    adjusted_breakpoints.append(bp)
        
        results = {
            "topology_analyses": {
                contig: {
                    "likely_circular": bool(analysis.likely_circular),
                    "confidence": float(analysis.confidence),
                    "evidence": self._make_json_serializable(analysis.evidence)
                }
                for contig, analysis in topology_analyses.items()
            },
            "adjusted_breakpoints": adjusted_breakpoints,  # Will be serialized by ValidationPipeline
            "original_breakpoint_count": len(breakpoints),
            "adjusted_breakpoint_count": len([bp for bp in adjusted_breakpoints if bp.confidence >= self.config.detection.min_confidence_score])
        }
        
        self.logger.info(f"Topology analysis: {len(topology_analyses)} contigs analyzed")
        circular_contigs = [c for c, a in topology_analyses.items() if a.likely_circular]
        if circular_contigs:
            self.logger.info(f"Detected likely circular contigs: {circular_contigs}")
        
        return results
    
    def _make_json_serializable(self, obj):
        """Convert objects to JSON-serializable format."""
        if isinstance(obj, dict):
            return {k: self._make_json_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._make_json_serializable(item) for item in obj]
        elif isinstance(obj, (bool, int, float, str, type(None))):
            return obj
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.integer, np.floating)):
            return obj.item()
        elif hasattr(obj, '__dict__'):
            return self._make_json_serializable(obj.__dict__)
        else:
            return str(obj)