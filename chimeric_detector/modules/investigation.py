"""Breakpoint investigation module for detailed analysis of detected chimeric junctions."""

import numpy as np
import logging
import json
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, asdict
from pathlib import Path
import subprocess
from collections import defaultdict

from ..core.window_analysis import Breakpoint
from ..core.bam_parser import BAMParser, ReadPairInfo
from ..config.config import Config


@dataclass
class BreakpointEvidence:
    """Detailed evidence for a specific breakpoint."""
    breakpoint: Breakpoint
    flanking_sequences: Dict[str, str]
    coverage_profile: Dict[str, List[Tuple[int, float]]]
    read_alignments: Dict[str, List[Dict[str, Any]]]
    split_reads: List[Dict[str, Any]]
    discordant_pairs: List[Dict[str, Any]]
    local_assembly: Optional[str]
    junction_confidence: float


class BreakpointInvestigator:
    """Detailed investigation of specific breakpoints."""
    
    def __init__(self, config: Config, bam_path: str, assembly_path: Optional[str] = None):
        self.config = config
        self.bam_path = bam_path
        self.assembly_path = assembly_path
        self.logger = logging.getLogger(__name__)
        
        # Load assembly if provided
        self.assembly = None
        if assembly_path:
            try:
                import pysam
                self.assembly = pysam.FastaFile(assembly_path)
                self.logger.info(f"Loaded assembly from {assembly_path}")
            except Exception as e:
                self.logger.warning(f"Could not load assembly: {e}")
    
    def investigate_breakpoint(self, breakpoint: Breakpoint, 
                             flank_size: int = 2000) -> BreakpointEvidence:
        """Perform detailed investigation of a single breakpoint."""
        
        self.logger.info(f"Investigating breakpoint {breakpoint.contig}:{breakpoint.position}")
        
        evidence = BreakpointEvidence(
            breakpoint=breakpoint,
            flanking_sequences={},
            coverage_profile={},
            read_alignments={},
            split_reads=[],
            discordant_pairs=[],
            local_assembly=None,
            junction_confidence=0.0
        )
        
        try:
            with BAMParser(self.bam_path, self.config.quality, use_streaming=True) as bam_parser:
                
                # 1. Extract flanking sequences
                if self.assembly:
                    evidence.flanking_sequences = self._extract_flanking_sequences(
                        breakpoint, flank_size
                    )
                
                # 2. Analyze coverage around breakpoint
                evidence.coverage_profile = self._analyze_coverage_profile(
                    breakpoint, bam_parser, flank_size
                )
                
                # 3. Extract and analyze read alignments
                evidence.read_alignments = self._extract_read_alignments(
                    breakpoint, bam_parser, flank_size
                )
                
                # 4. Find split reads (if any)
                evidence.split_reads = self._find_split_reads(
                    breakpoint, bam_parser, flank_size
                )
                
                # 5. Analyze discordant read pairs
                evidence.discordant_pairs = self._analyze_discordant_pairs(
                    breakpoint, bam_parser, flank_size
                )
                
                # 6. Calculate junction confidence
                evidence.junction_confidence = self._calculate_junction_confidence(evidence)
                
                # 7. Attempt local assembly (optional)
                if len(evidence.split_reads) > 5:  # Only if we have sufficient evidence
                    evidence.local_assembly = self._attempt_local_assembly(
                        breakpoint, evidence.split_reads
                    )
                
        except Exception as e:
            self.logger.error(f"Error investigating breakpoint {breakpoint.contig}:{breakpoint.position}: {e}")
        
        return evidence
    
    def _extract_flanking_sequences(self, breakpoint: Breakpoint, 
                                  flank_size: int) -> Dict[str, str]:
        """Extract DNA sequences flanking the breakpoint."""
        
        if not self.assembly:
            return {}
        
        try:
            # Extract left and right flanking sequences
            left_start = max(0, breakpoint.position - flank_size)
            left_end = breakpoint.position
            right_start = breakpoint.position
            right_end = breakpoint.position + flank_size
            
            left_seq = self.assembly.fetch(breakpoint.contig, left_start, left_end)
            right_seq = self.assembly.fetch(breakpoint.contig, right_start, right_end)
            
            return {
                "left_flank": left_seq.upper(),
                "right_flank": right_seq.upper(),
                "left_coords": f"{breakpoint.contig}:{left_start}-{left_end}",
                "right_coords": f"{breakpoint.contig}:{right_start}-{right_end}",
                "junction_position": breakpoint.position
            }
            
        except Exception as e:
            self.logger.warning(f"Could not extract flanking sequences: {e}")
            return {}
    
    def _analyze_coverage_profile(self, breakpoint: Breakpoint, bam_parser: BAMParser,
                                flank_size: int) -> Dict[str, List[Tuple[int, float]]]:
        """Analyze read coverage around the breakpoint."""
        
        try:
            # Define analysis windows
            start = max(0, breakpoint.position - flank_size)
            end = breakpoint.position + flank_size
            window_size = 100  # 100bp windows for coverage analysis
            
            coverage_data = []
            
            # Calculate coverage in sliding windows
            for pos in range(start, end, window_size):
                window_end = min(pos + window_size, end)
                coverage = bam_parser.calculate_coverage_in_window(
                    breakpoint.contig, pos, window_end
                )
                coverage_data.append((pos + window_size // 2, coverage))
            
            return {
                "profile": coverage_data,
                "breakpoint_position": breakpoint.position,
                "window_size": window_size,
                "analysis_region": f"{start}-{end}"
            }
            
        except Exception as e:
            self.logger.warning(f"Could not analyze coverage profile: {e}")
            return {}
    
    def _extract_read_alignments(self, breakpoint: Breakpoint, bam_parser: BAMParser,
                               flank_size: int) -> Dict[str, List[Dict[str, Any]]]:
        """Extract detailed read alignment information around breakpoint."""
        
        try:
            # Get reads in expanded region
            start = max(0, breakpoint.position - flank_size)
            end = breakpoint.position + flank_size
            
            read_pairs = list(bam_parser.stream_read_pairs_in_window(
                breakpoint.contig, start, end
            ))
            
            # Categorize reads
            spanning_reads = []
            left_reads = []
            right_reads = []
            anomalous_reads = []
            
            for rp in read_pairs:
                read_info = {
                    "position": rp.position,
                    "mate_position": rp.mate_position,
                    "insert_size": rp.insert_size,
                    "is_proper_pair": rp.is_proper_pair,
                    "is_discordant": rp.is_discordant,
                    "mapping_quality": rp.mapping_quality,
                    "distance_to_breakpoint": abs(rp.position - breakpoint.position)
                }
                
                # Categorize based on position relative to breakpoint
                if rp.position < breakpoint.position and rp.mate_position and rp.mate_position > breakpoint.position:
                    spanning_reads.append(read_info)
                elif rp.position < breakpoint.position:
                    left_reads.append(read_info)
                elif rp.position > breakpoint.position:
                    right_reads.append(read_info)
                
                # Flag anomalous reads
                if rp.is_discordant or abs(rp.insert_size) > self.config.quality.max_insert_size:
                    anomalous_reads.append(read_info)
            
            return {
                "spanning_reads": spanning_reads,
                "left_reads": left_reads,
                "right_reads": right_reads,
                "anomalous_reads": anomalous_reads,
                "total_reads": len(read_pairs),
                "analysis_region": f"{start}-{end}"
            }
            
        except Exception as e:
            self.logger.warning(f"Could not extract read alignments: {e}")
            return {}
    
    def _find_split_reads(self, breakpoint: Breakpoint, bam_parser: BAMParser,
                         flank_size: int) -> List[Dict[str, Any]]:
        """Find reads that might be split across the breakpoint."""
        
        # This is a simplified implementation - in practice you'd look for
        # supplementary alignments, soft-clipped reads, etc.
        
        try:
            start = max(0, breakpoint.position - 200)  # Narrow window for split reads
            end = breakpoint.position + 200
            
            split_candidates = []
            
            # In a full implementation, you would examine the BAM file directly
            # looking for reads with supplementary alignments (SA tag) or
            # significant soft clipping near the breakpoint
            
            # For now, return placeholder structure
            return split_candidates
            
        except Exception as e:
            self.logger.warning(f"Could not find split reads: {e}")
            return []
    
    def _analyze_discordant_pairs(self, breakpoint: Breakpoint, bam_parser: BAMParser,
                                flank_size: int) -> List[Dict[str, Any]]:
        """Analyze discordant read pairs around the breakpoint."""
        
        try:
            start = max(0, breakpoint.position - flank_size)
            end = breakpoint.position + flank_size
            
            read_pairs = list(bam_parser.stream_read_pairs_in_window(
                breakpoint.contig, start, end
            ))
            
            discordant_analysis = []
            
            for rp in read_pairs:
                if rp.is_discordant:
                    discordant_info = {
                        "position": rp.position,
                        "mate_position": rp.mate_position,
                        "mate_contig": rp.mate_contig,
                        "insert_size": rp.insert_size,
                        "mapping_quality": rp.mapping_quality,
                        "distance_to_breakpoint": abs(rp.position - breakpoint.position),
                        "discordant_type": self._classify_discordant_type(rp)
                    }
                    discordant_analysis.append(discordant_info)
            
            return discordant_analysis
            
        except Exception as e:
            self.logger.warning(f"Could not analyze discordant pairs: {e}")
            return []
    
    def _classify_discordant_type(self, read_pair: ReadPairInfo) -> str:
        """Classify the type of discordant read pair."""
        
        if read_pair.mate_unmapped:
            return "mate_unmapped"
        elif read_pair.mate_contig and read_pair.mate_contig != read_pair.contig:
            return "inter_contig"
        elif abs(read_pair.insert_size) > self.config.quality.max_insert_size:
            return "large_insert"
        else:
            return "orientation_anomaly"
    
    def _calculate_junction_confidence(self, evidence: BreakpointEvidence) -> float:
        """Calculate confidence in the junction based on multiple evidence types."""
        
        confidence_factors = []
        
        # Factor 1: Original breakpoint confidence
        confidence_factors.append(evidence.breakpoint.confidence)
        
        # Factor 2: Coverage profile (look for drops at breakpoint)
        if evidence.coverage_profile and evidence.coverage_profile.get("profile"):
            profile = evidence.coverage_profile["profile"]
            if len(profile) > 2:
                coverages = [cov for pos, cov in profile]
                bp_pos = evidence.breakpoint.position
                
                # Find coverage values closest to breakpoint
                bp_coverage = min(coverages, key=lambda x: abs(x - bp_pos), default=0)
                avg_coverage = np.mean(coverages) if coverages else 0
                
                if avg_coverage > 0:
                    coverage_factor = bp_coverage / avg_coverage
                    confidence_factors.append(1.0 - coverage_factor)  # Lower coverage = higher confidence
        
        # Factor 3: Discordant pair evidence
        if evidence.discordant_pairs:
            discordant_count = len(evidence.discordant_pairs)
            discordant_factor = min(1.0, discordant_count / 10.0)  # Normalize
            confidence_factors.append(discordant_factor)
        
        # Factor 4: Read alignment anomalies
        if evidence.read_alignments and evidence.read_alignments.get("anomalous_reads"):
            anomalous_count = len(evidence.read_alignments["anomalous_reads"])
            total_reads = evidence.read_alignments.get("total_reads", 1)
            anomaly_rate = anomalous_count / total_reads
            confidence_factors.append(anomaly_rate)
        
        # Calculate weighted average
        if confidence_factors:
            return np.mean(confidence_factors)
        else:
            return evidence.breakpoint.confidence
    
    def _attempt_local_assembly(self, breakpoint: Breakpoint, 
                              split_reads: List[Dict[str, Any]]) -> Optional[str]:
        """Attempt local assembly around the breakpoint (placeholder)."""
        
        # This would require implementing or calling external assembly tools
        # For now, return None as placeholder
        
        self.logger.debug(f"Local assembly attempted for {breakpoint.contig}:{breakpoint.position}")
        return None
    
    def generate_investigation_report(self, breakpoints: List[Breakpoint],
                                    output_dir: Path, 
                                    detailed: bool = True) -> Dict[str, Any]:
        """Generate comprehensive investigation report for multiple breakpoints."""
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger.info(f"Generating investigation report for {len(breakpoints)} breakpoints")
        
        investigations = []
        summary_stats = {
            "total_breakpoints": len(breakpoints),
            "high_confidence_breakpoints": 0,
            "contigs_affected": set(),
            "evidence_summary": defaultdict(int)
        }
        
        for i, bp in enumerate(breakpoints):
            self.logger.info(f"Processing breakpoint {i+1}/{len(breakpoints)}: {bp.contig}:{bp.position}")
            
            # Investigate breakpoint
            evidence = self.investigate_breakpoint(bp)
            investigations.append(evidence)
            
            # Update summary statistics
            if evidence.junction_confidence > 0.8:
                summary_stats["high_confidence_breakpoints"] += 1
            
            summary_stats["contigs_affected"].add(bp.contig)
            
            if evidence.discordant_pairs:
                summary_stats["evidence_summary"]["discordant_pairs"] += len(evidence.discordant_pairs)
            if evidence.split_reads:
                summary_stats["evidence_summary"]["split_reads"] += len(evidence.split_reads)
            
            # Write individual breakpoint report
            if detailed:
                bp_report_file = output_dir / f"breakpoint_{bp.contig}_{bp.position}.json"
                with open(bp_report_file, 'w') as f:
                    json.dump(self._serialize_evidence(evidence), f, indent=2)
        
        # Convert set to list for JSON serialization
        summary_stats["contigs_affected"] = list(summary_stats["contigs_affected"])
        
        # Write summary report
        summary_report = {
            "summary": summary_stats,
            "investigations": [self._serialize_evidence(inv) for inv in investigations]
        }
        
        summary_file = output_dir / "investigation_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary_report, f, indent=2)
        
        # Generate HTML report (optional)
        if detailed:
            self._generate_html_report(investigations, output_dir, summary_stats)
        
        self.logger.info(f"Investigation report written to {output_dir}")
        
        return summary_report
    
    def _serialize_evidence(self, evidence: BreakpointEvidence) -> Dict[str, Any]:
        """Convert BreakpointEvidence to JSON-serializable format."""
        
        return {
            "breakpoint": {
                "contig": evidence.breakpoint.contig,
                "position": evidence.breakpoint.position,
                "confidence": evidence.breakpoint.confidence,
                "evidence": evidence.breakpoint.evidence
            },
            "flanking_sequences": evidence.flanking_sequences,
            "coverage_profile": evidence.coverage_profile,
            "read_alignments": evidence.read_alignments,
            "split_reads": evidence.split_reads,
            "discordant_pairs": evidence.discordant_pairs,
            "local_assembly": evidence.local_assembly,
            "junction_confidence": evidence.junction_confidence
        }
    
    def _generate_html_report(self, investigations: List[BreakpointEvidence],
                            output_dir: Path, summary_stats: Dict[str, Any]) -> None:
        """Generate HTML report with visualizations."""
        
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Breakpoint Investigation Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                .summary {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .breakpoint {{ border: 1px solid #ccc; margin: 20px 0; padding: 15px; border-radius: 5px; }}
                .high-confidence {{ border-color: #ff6b6b; background-color: #ffe0e0; }}
                .medium-confidence {{ border-color: #ffa500; background-color: #fff5e0; }}
                .low-confidence {{ border-color: #cccccc; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
            </style>
        </head>
        <body>
            <h1>Breakpoint Investigation Report</h1>
            
            <div class="summary">
                <h2>Summary</h2>
                <p><strong>Total Breakpoints:</strong> {summary_stats['total_breakpoints']}</p>
                <p><strong>High Confidence:</strong> {summary_stats['high_confidence_breakpoints']}</p>
                <p><strong>Contigs Affected:</strong> {len(summary_stats['contigs_affected'])}</p>
            </div>
        """
        
        for inv in investigations:
            bp = inv.breakpoint
            confidence_class = (
                "high-confidence" if inv.junction_confidence > 0.8 
                else "medium-confidence" if inv.junction_confidence > 0.6 
                else "low-confidence"
            )
            
            html_content += f"""
            <div class="breakpoint {confidence_class}">
                <h3>{bp.contig}:{bp.position}</h3>
                <p><strong>Original Confidence:</strong> {bp.confidence:.3f}</p>
                <p><strong>Junction Confidence:</strong> {inv.junction_confidence:.3f}</p>
                <p><strong>Discordant Pairs:</strong> {len(inv.discordant_pairs)}</p>
                <p><strong>Split Reads:</strong> {len(inv.split_reads)}</p>
            </div>
            """
        
        html_content += """
        </body>
        </html>
        """
        
        html_file = output_dir / "investigation_report.html"
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info(f"HTML report generated: {html_file}")


# Command-line interface for standalone investigation
def run_investigation_cli():
    """Command-line interface for breakpoint investigation."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Investigate detected breakpoints in detail"
    )
    
    parser.add_argument("results_file", help="JSON results file from chimera detection")
    parser.add_argument("bam_file", help="BAM file used for detection")
    parser.add_argument("output_dir", help="Output directory for investigation results")
    parser.add_argument("--assembly", help="Assembly FASTA file (optional)")
    parser.add_argument("--min-confidence", type=float, default=0.7,
                       help="Minimum confidence threshold for investigation")
    parser.add_argument("--detailed", action="store_true",
                       help="Generate detailed individual reports")
    
    args = parser.parse_args()
    
    # Load results
    with open(args.results_file) as f:
        results = json.load(f)
    
    # Filter breakpoints by confidence
    breakpoints = []
    for bp_data in results.get("breakpoints", []):
        if bp_data["confidence"] >= args.min_confidence:
            bp = Breakpoint(
                contig=bp_data["contig"],
                position=bp_data["position"],
                confidence=bp_data["confidence"],
                evidence=bp_data.get("evidence", {}),
                left_metrics=None,
                right_metrics=None
            )
            breakpoints.append(bp)
    
    print(f"Investigating {len(breakpoints)} breakpoints with confidence >= {args.min_confidence}")
    
    # Create investigator
    from ..config.config import Config
    config = Config()
    investigator = BreakpointInvestigator(config, args.bam_file, args.assembly)
    
    # Run investigation
    report = investigator.generate_investigation_report(
        breakpoints, Path(args.output_dir), detailed=args.detailed
    )
    
    print(f"Investigation complete. Results in {args.output_dir}")
    print(f"High confidence breakpoints: {report['summary']['high_confidence_breakpoints']}")


if __name__ == "__main__":
    run_investigation_cli()