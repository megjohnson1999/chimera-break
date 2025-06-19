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
        """Generate interactive HTML breakpoint browser."""
        
        # Prepare data for JavaScript
        breakpoint_data = []
        contig_data = defaultdict(lambda: {"breakpoints": [], "total_length": 0, "is_circular": False})
        
        for inv in investigations:
            bp = inv.breakpoint
            
            # Extract topology information if available
            is_circular = False
            topology_confidence = 0.0
            if bp.evidence and "topology_analysis" in bp.evidence:
                is_circular = bp.evidence["topology_analysis"]
            if bp.evidence and "topology_adjustment" in bp.evidence:
                topology_confidence = abs(bp.evidence.get("topology_adjustment", 0))
            
            bp_data = {
                "contig": bp.contig,
                "position": bp.position,
                "confidence": bp.confidence,
                "junction_confidence": inv.junction_confidence,
                "discordant_pairs": len(inv.discordant_pairs),
                "split_reads": len(inv.split_reads),
                "is_circular": is_circular,
                "topology_confidence": topology_confidence,
                "evidence_summary": self._summarize_evidence(inv)
            }
            breakpoint_data.append(bp_data)
            
            # Update contig data
            contig_data[bp.contig]["breakpoints"].append({
                "position": bp.position,
                "confidence": bp.confidence,
                "junction_confidence": inv.junction_confidence
            })
            contig_data[bp.contig]["is_circular"] = is_circular
        
        # Estimate contig lengths from breakpoint positions
        for contig, data in contig_data.items():
            if data["breakpoints"]:
                positions = [bp["position"] for bp in data["breakpoints"]]
                data["total_length"] = max(positions) + 10000  # Add buffer
        
        # Count circular contigs
        circular_contigs = [c for c, data in contig_data.items() if data["is_circular"]]
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Chimera Breakpoint Browser</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .header h1 {{
            margin: 0 0 10px 0;
            font-size: 2.5em;
        }}
        .header .subtitle {{
            opacity: 0.9;
            font-size: 1.1em;
        }}
        .controls {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            display: flex;
            gap: 20px;
            align-items: center;
            flex-wrap: wrap;
        }}
        .control-group {{
            display: flex;
            flex-direction: column;
            gap: 5px;
        }}
        .control-group label {{
            font-weight: bold;
            color: #555;
            font-size: 0.9em;
        }}
        select, input[type="range"] {{
            padding: 8px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 0.9em;
        }}
        .summary-cards {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 20px;
        }}
        .summary-card {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .summary-card .number {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
        }}
        .summary-card .label {{
            color: #666;
            margin-top: 5px;
        }}
        .plot-container {{
            background: white;
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }}
        .details-panel {{
            background: white;
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            min-height: 200px;
        }}
        .details-panel h3 {{
            margin-top: 0;
            color: #333;
        }}
        .evidence-item {{
            background: #f8f9fa;
            padding: 10px;
            margin: 10px 0;
            border-left: 4px solid #667eea;
            border-radius: 4px;
        }}
        .circular-badge {{
            background: #ff6b6b;
            color: white;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 0.8em;
            margin-left: 10px;
        }}
        .confidence-high {{ color: #28a745; }}
        .confidence-medium {{ color: #ffc107; }}
        .confidence-low {{ color: #dc3545; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 Chimera Breakpoint Browser</h1>
        <div class="subtitle">Interactive analysis of chimeric contig breakpoints</div>
    </div>
    
    <div class="summary-cards">
        <div class="summary-card">
            <div class="number">{len(breakpoint_data)}</div>
            <div class="label">Total Breakpoints</div>
        </div>
        <div class="summary-card">
            <div class="number">{summary_stats['high_confidence_breakpoints']}</div>
            <div class="label">High Confidence</div>
        </div>
        <div class="summary-card">
            <div class="number">{len(summary_stats['contigs_affected'])}</div>
            <div class="label">Contigs Affected</div>
        </div>
        <div class="summary-card">
            <div class="number">{len(circular_contigs)}</div>
            <div class="label">Circular Contigs</div>
        </div>
    </div>
    
    <div class="controls">
        <div class="control-group">
            <label for="confidenceSlider">Minimum Confidence:</label>
            <input type="range" id="confidenceSlider" min="0" max="1" step="0.1" value="0.5">
            <span id="confidenceValue">0.5</span>
        </div>
        <div class="control-group">
            <label for="contigSelect">Filter by Contig:</label>
            <select id="contigSelect">
                <option value="">All Contigs</option>
                {chr(10).join([f'<option value="{contig}">{contig}{"" if not data["is_circular"] else " (Circular)"}</option>' for contig, data in contig_data.items()])}
            </select>
        </div>
        <div class="control-group">
            <label for="topologyFilter">Topology:</label>
            <select id="topologyFilter">
                <option value="">All</option>
                <option value="circular">Circular Only</option>
                <option value="linear">Linear Only</option>
            </select>
        </div>
    </div>
    
    <div class="plot-container">
        <div id="contigPlot" style="width:100%; height:600px;"></div>
    </div>
    
    <div class="details-panel">
        <h3>Breakpoint Details</h3>
        <div id="breakpointDetails">
            <p>Click on a breakpoint above to see detailed information.</p>
        </div>
    </div>

    <script>
        // Data from Python
        const breakpoints = {json.dumps(breakpoint_data, indent=8)};
        const contigs = {json.dumps(dict(contig_data), indent=8)};
        
        let filteredBreakpoints = [...breakpoints];
        
        // Create initial plot
        function createContigPlot() {{
            const traces = [];
            const contigNames = Object.keys(contigs);
            
            contigNames.forEach((contig, idx) => {{
                const contigBreakpoints = filteredBreakpoints.filter(bp => bp.contig === contig);
                const isCircular = contigs[contig].is_circular;
                
                if (contigBreakpoints.length === 0) return;
                
                const x = contigBreakpoints.map(bp => bp.position);
                const y = Array(contigBreakpoints.length).fill(idx);
                const colors = contigBreakpoints.map(bp => {{
                    if (bp.junction_confidence > 0.8) return '#28a745';
                    if (bp.junction_confidence > 0.6) return '#ffc107';
                    return '#dc3545';
                }});
                const sizes = contigBreakpoints.map(bp => Math.max(8, bp.junction_confidence * 15));
                
                const hovertext = contigBreakpoints.map(bp => 
                    `${{bp.contig}}:${{bp.position}}<br>` +
                    `Confidence: ${{bp.junction_confidence.toFixed(3)}}<br>` +
                    `Discordant pairs: ${{bp.discordant_pairs}}<br>` +
                    `Split reads: ${{bp.split_reads}}<br>` +
                    `Topology: ${{bp.is_circular ? 'Circular' : 'Linear'}}`
                );
                
                traces.push({{
                    x: x,
                    y: y,
                    mode: 'markers',
                    type: 'scatter',
                    name: contig + (isCircular ? ' (Circular)' : ''),
                    marker: {{
                        color: colors,
                        size: sizes,
                        symbol: isCircular ? 'circle' : 'diamond',
                        line: {{
                            width: isCircular ? 2 : 1,
                            color: isCircular ? '#ff6b6b' : '#666'
                        }}
                    }},
                    text: hovertext,
                    hovertemplate: '%{{text}}<extra></extra>',
                    customdata: contigBreakpoints
                }});
            }});
            
            const layout = {{
                title: 'Breakpoint Distribution Across Contigs',
                xaxis: {{
                    title: 'Position (bp)',
                    showgrid: true,
                    gridcolor: '#f0f0f0'
                }},
                yaxis: {{
                    title: 'Contigs',
                    tickvals: contigNames.map((_, idx) => idx),
                    ticktext: contigNames.map(name => contigs[name].is_circular ? name + ' ●' : name),
                    showgrid: true,
                    gridcolor: '#f0f0f0'
                }},
                plot_bgcolor: 'white',
                paper_bgcolor: 'white',
                hovermode: 'closest',
                showlegend: false,
                margin: {{ l: 100, r: 50, t: 50, b: 50 }}
            }};
            
            const config = {{
                responsive: true,
                displayModeBar: true,
                modeBarButtonsToRemove: ['select2d', 'lasso2d']
            }};
            
            Plotly.newPlot('contigPlot', traces, layout, config);
            
            // Add click handler
            document.getElementById('contigPlot').on('plotly_click', function(data) {{
                const point = data.points[0];
                if (point.customdata) {{
                    showBreakpointDetails(point.customdata[point.pointIndex]);
                }}
            }});
        }}
        
        function updateFilters() {{
            const minConfidence = parseFloat(document.getElementById('confidenceSlider').value);
            const selectedContig = document.getElementById('contigSelect').value;
            const topologyFilter = document.getElementById('topologyFilter').value;
            
            filteredBreakpoints = breakpoints.filter(bp => {{
                if (bp.junction_confidence < minConfidence) return false;
                if (selectedContig && bp.contig !== selectedContig) return false;
                if (topologyFilter === 'circular' && !bp.is_circular) return false;
                if (topologyFilter === 'linear' && bp.is_circular) return false;
                return true;
            }});
            
            createContigPlot();
        }}
        
        function showBreakpointDetails(breakpoint) {{
            const details = document.getElementById('breakpointDetails');
            const confidenceClass = breakpoint.junction_confidence > 0.8 ? 'confidence-high' : 
                                  breakpoint.junction_confidence > 0.6 ? 'confidence-medium' : 'confidence-low';
            
            details.innerHTML = `
                <h4>${{breakpoint.contig}}:${{breakpoint.position}} ${{breakpoint.is_circular ? '<span class="circular-badge">Circular</span>' : ''}}</h4>
                <div class="evidence-item">
                    <strong>Confidence:</strong> 
                    <span class="${{confidenceClass}}">${{breakpoint.junction_confidence.toFixed(3)}}</span>
                    ${{breakpoint.confidence !== breakpoint.junction_confidence ? 
                        ` (original: ${{breakpoint.confidence.toFixed(3)}})` : ''}}
                </div>
                <div class="evidence-item">
                    <strong>Evidence:</strong><br>
                    • Discordant read pairs: ${{breakpoint.discordant_pairs}}<br>
                    • Split reads: ${{breakpoint.split_reads}}<br>
                    • Topology analysis: ${{breakpoint.is_circular ? 'Likely circular genome' : 'Linear structure'}}
                    ${{breakpoint.topology_confidence > 0 ? 
                        `<br>• Confidence penalty: -${{breakpoint.topology_confidence.toFixed(2)}}` : ''}}
                </div>
                <div class="evidence-item">
                    <strong>Summary:</strong> ${{breakpoint.evidence_summary}}
                </div>
            `;
        }}
        
        // Event listeners
        document.getElementById('confidenceSlider').addEventListener('input', function() {{
            document.getElementById('confidenceValue').textContent = this.value;
            updateFilters();
        }});
        
        document.getElementById('contigSelect').addEventListener('change', updateFilters);
        document.getElementById('topologyFilter').addEventListener('change', updateFilters);
        
        // Initialize
        createContigPlot();
    </script>
</body>
</html>
        """
        
        html_file = output_dir / "breakpoint_browser.html"
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info(f"Interactive breakpoint browser generated: {html_file}")
    
    def _summarize_evidence(self, evidence: BreakpointEvidence) -> str:
        """Create a brief evidence summary for display."""
        summary_parts = []
        
        if evidence.discordant_pairs:
            summary_parts.append(f"{len(evidence.discordant_pairs)} discordant pairs")
        
        if evidence.split_reads:
            summary_parts.append(f"{len(evidence.split_reads)} split reads")
        
        if evidence.breakpoint.evidence.get("topology_analysis"):
            summary_parts.append("circular genome detected")
        
        if evidence.coverage_profile and evidence.coverage_profile.get("profile"):
            summary_parts.append("coverage anomaly")
        
        return "; ".join(summary_parts) if summary_parts else "limited evidence"


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