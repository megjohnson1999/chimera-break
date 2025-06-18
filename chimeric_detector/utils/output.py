"""Output formatting and reporting utilities."""

import json
import csv
from pathlib import Path
from typing import List, Dict, Any, Optional
import logging
from datetime import datetime

from ..core.window_analysis import Breakpoint
from ..config.config import OutputConfig


class OutputFormatter:
    """Handles output formatting for different formats."""
    
    def __init__(self, config: OutputConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
    def write_results(self, breakpoints: List[Breakpoint], output_path: Path, 
                     metadata: Optional[Dict[str, Any]] = None) -> None:
        """Write results in the configured format."""
        if self.config.format == "json":
            self._write_json(breakpoints, output_path, metadata)
        elif self.config.format == "tsv":
            self._write_tsv(breakpoints, output_path)
        elif self.config.format == "bed":
            self._write_bed(breakpoints, output_path)
        else:
            raise ValueError(f"Unknown output format: {self.config.format}")
            
    def _write_json(self, breakpoints: List[Breakpoint], output_path: Path,
                   metadata: Optional[Dict[str, Any]] = None) -> None:
        """Write results in JSON format."""
        results = {
            "metadata": metadata or {},
            "timestamp": datetime.now().isoformat(),
            "breakpoints": []
        }
        
        for bp in breakpoints:
            bp_data = {
                "contig": bp.contig,
                "position": bp.position,
                "confidence": round(bp.confidence, 4),
                "evidence": {k: round(v, 4) for k, v in bp.evidence.items()}
            }
            
            if self.config.include_read_evidence:
                bp_data["left_window"] = self._window_to_dict(bp.left_metrics)
                bp_data["right_window"] = self._window_to_dict(bp.right_metrics)
                
            results["breakpoints"].append(bp_data)
            
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
            
    def _write_tsv(self, breakpoints: List[Breakpoint], output_path: Path) -> None:
        """Write results in TSV format."""
        headers = [
            "contig", "position", "confidence",
            "proper_pair_drop", "insert_size_zscore", "discordant_rate"
        ]
        
        if self.config.include_read_evidence:
            headers.extend([
                "left_coverage", "left_proper_rate", "left_insert_median",
                "right_coverage", "right_proper_rate", "right_insert_median"
            ])
            
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(headers)
            
            for bp in breakpoints:
                row = [
                    bp.contig,
                    bp.position,
                    round(bp.confidence, 4),
                    round(bp.evidence.get("proper_pair_drop", 0), 4),
                    round(bp.evidence.get("insert_size_zscore", 0), 4),
                    round(bp.evidence.get("discordant_rate", 0), 4)
                ]
                
                if self.config.include_read_evidence:
                    row.extend([
                        round(bp.left_metrics.coverage, 2),
                        round(bp.left_metrics.proper_pair_rate, 4),
                        round(bp.left_metrics.insert_size_median, 0),
                        round(bp.right_metrics.coverage, 2),
                        round(bp.right_metrics.proper_pair_rate, 4),
                        round(bp.right_metrics.insert_size_median, 0)
                    ])
                    
                writer.writerow(row)
                
    def _write_bed(self, breakpoints: List[Breakpoint], output_path: Path) -> None:
        """Write results in BED format."""
        with open(output_path, 'w') as f:
            for bp in breakpoints:
                # BED format: chr start end name score
                start = max(0, bp.position - 50)
                end = bp.position + 50
                score = int(bp.confidence * 1000)
                
                f.write(f"{bp.contig}\t{start}\t{end}\tbreakpoint\t{score}\n")
                
    def _window_to_dict(self, window) -> Dict[str, Any]:
        """Convert window metrics to dictionary."""
        return {
            "start": window.start,
            "end": window.end,
            "coverage": round(window.coverage, 2),
            "read_count": window.read_count,
            "proper_pair_rate": round(window.proper_pair_rate, 4),
            "discordant_pair_rate": round(window.discordant_pair_rate, 4),
            "insert_size_median": round(window.insert_size_median, 0),
            "insert_size_mad": round(window.insert_size_mad, 0)
        }


class Logger:
    """Configurable logging utility."""
    
    @staticmethod
    def setup_logging(verbose: bool = False, debug: bool = False, 
                     log_file: Optional[Path] = None) -> None:
        """Set up logging configuration."""
        level = logging.WARNING
        if debug:
            level = logging.DEBUG
        elif verbose:
            level = logging.INFO
            
        handlers = [logging.StreamHandler()]
        if log_file:
            handlers.append(logging.FileHandler(log_file))
            
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=handlers
        )
        
        # Set specific loggers
        logging.getLogger("pysam").setLevel(logging.WARNING)