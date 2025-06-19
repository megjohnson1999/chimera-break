"""Validation modules for chimeric contig detection."""

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional
import logging
from pathlib import Path

from ..core.window_analysis import Breakpoint, WindowMetrics
from ..config.config import ValidationConfig
from .topology import TopologyValidator


class ValidationModule(ABC):
    """Base class for validation modules."""
    
    def __init__(self, config: ValidationConfig):
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
        
    @abstractmethod
    def validate(self, breakpoints: List[Breakpoint], **kwargs) -> List[Dict[str, Any]]:
        """Validate breakpoints and return validation results."""
        pass
        

class TaxonomyValidator(ValidationModule):
    """Validates breakpoints using taxonomic information."""
    
    def __init__(self, config: ValidationConfig):
        super().__init__(config)
        self.taxonomy_db = None
        self.assembly_sequences = None
        if config.taxonomy_db:
            self._load_taxonomy_db(config.taxonomy_db)
            
    def _load_taxonomy_db(self, db_path: str) -> None:
        """Load taxonomy database (implementation depends on format)."""
        self.logger.info(f"Loading taxonomy database from {db_path}")
        # For now, just mark as loaded - real implementation would load
        # a database like NCBI taxonomy, Kraken2 database, etc.
        self.taxonomy_db = {"loaded": True}
        
    def load_assembly(self, assembly_path: str) -> None:
        """Load assembly sequences for taxonomy analysis."""
        try:
            import pysam
            self.assembly_sequences = pysam.FastaFile(assembly_path)
            self.logger.info(f"Loaded assembly sequences from {assembly_path}")
        except ImportError:
            self.logger.error("pysam required for assembly loading")
            self.assembly_sequences = None
        except Exception as e:
            self.logger.error(f"Failed to load assembly: {e}")
            self.assembly_sequences = None
            
    def _classify_sequence(self, sequence: str) -> str:
        """Classify a sequence taxonomically (placeholder)."""
        # Real implementation would use:
        # - BLAST against NCBI nt database
        # - Kraken2/Kraken classification
        # - MMseqs2 taxonomy
        # - Custom HMM models
        
        # For demonstration, simple heuristic based on GC content
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
        
        if gc_content < 0.35:
            return "low_gc_virus"  # Many RNA viruses
        elif gc_content > 0.55:
            return "high_gc_bacteria"  # Some bacteria
        else:
            return "moderate_gc_organism"
            
    def validate(self, breakpoints: List[Breakpoint], **kwargs) -> List[Dict[str, Any]]:
        """Validate breakpoints based on taxonomic consistency."""
        if not self.assembly_sequences:
            self.logger.warning("No assembly sequences loaded - taxonomy validation disabled")
            return [{"error": "Assembly sequences required for taxonomy validation"}]
            
        results = []
        
        for bp in breakpoints:
            try:
                # Extract sequences flanking the breakpoint
                flank_size = 1000  # 1kb on each side
                left_start = max(0, bp.position - flank_size)
                left_end = bp.position
                right_start = bp.position
                right_end = bp.position + flank_size
                
                # Get sequences
                left_seq = self.assembly_sequences.fetch(bp.contig, left_start, left_end)
                right_seq = self.assembly_sequences.fetch(bp.contig, right_start, right_end)
                
                if len(left_seq) < 500 or len(right_seq) < 500:
                    # Not enough sequence for reliable classification
                    validation = {
                        "breakpoint": f"{bp.contig}:{bp.position}",
                        "taxonomy_consistent": None,
                        "confidence_adjustment": 0.0,
                        "notes": "Insufficient sequence for taxonomy analysis"
                    }
                else:
                    # Classify both flanking regions
                    left_taxon = self._classify_sequence(left_seq)
                    right_taxon = self._classify_sequence(right_seq)
                    
                    # Check consistency
                    consistent = left_taxon == right_taxon
                    confidence_penalty = 0.0 if consistent else 0.2
                    
                    validation = {
                        "breakpoint": f"{bp.contig}:{bp.position}",
                        "left_taxonomy": left_taxon,
                        "right_taxonomy": right_taxon,
                        "taxonomy_consistent": consistent,
                        "confidence_adjustment": -confidence_penalty,
                        "notes": f"Left: {left_taxon}, Right: {right_taxon}"
                    }
                    
            except Exception as e:
                validation = {
                    "breakpoint": f"{bp.contig}:{bp.position}",
                    "taxonomy_consistent": None,
                    "confidence_adjustment": 0.0,
                    "notes": f"Taxonomy analysis failed: {e}"
                }
                
            results.append(validation)
            
        return results
        

class SplitEvaluator(ValidationModule):
    """Evaluates the quality of potential contig splits."""
    
    def validate(self, breakpoints: List[Breakpoint], contig_lengths: Dict[str, int] = None) -> List[Dict[str, Any]]:
        """Evaluate if splitting at breakpoints produces viable contigs."""
        results = []
        
        # Group breakpoints by contig
        by_contig = {}
        for bp in breakpoints:
            if bp.contig not in by_contig:
                by_contig[bp.contig] = []
            by_contig[bp.contig].append(bp)
            
        for contig, contig_breakpoints in by_contig.items():
            # Sort by position
            contig_breakpoints.sort(key=lambda x: x.position)
            
            # Calculate split segments
            segments = []
            prev_pos = 0
            
            for bp in contig_breakpoints:
                segment_length = bp.position - prev_pos
                segments.append({
                    "start": prev_pos,
                    "end": bp.position,
                    "length": segment_length,
                    "viable": segment_length >= self.config.min_split_length
                })
                prev_pos = bp.position
                
            # Add final segment
            if contig_lengths and contig in contig_lengths:
                final_length = contig_lengths[contig] - prev_pos
                segments.append({
                    "start": prev_pos,
                    "end": contig_lengths[contig],
                    "length": final_length,
                    "viable": final_length >= self.config.min_split_length
                })
                
            # Evaluate overall split quality
            viable_segments = sum(1 for s in segments if s["viable"])
            total_segments = len(segments)
            
            evaluation = {
                "contig": contig,
                "breakpoint_count": len(contig_breakpoints),
                "segments": segments,
                "viable_segments": viable_segments,
                "total_segments": total_segments,
                "split_quality": viable_segments / total_segments if total_segments > 0 else 0
            }
            
            results.append(evaluation)
            
        return results
        

class ValidationPipeline:
    """Manages multiple validation modules including topology analysis."""
    
    def __init__(self, config: ValidationConfig, enable_topology: bool = True):
        self.config = config
        self.enable_topology = enable_topology
        self.logger = logging.getLogger(__name__)
        self.modules = []
        self.topology_validator = None
        
        if config.enabled and config.taxonomy_db:
            self.modules.append(TaxonomyValidator(config))
            
        if config.evaluate_splits:
            self.modules.append(SplitEvaluator(config))
        
        # Add topology validator if enabled
        if enable_topology:
            from ..config.config import Config
            # Create a full config object for topology analysis
            full_config = Config()
            full_config.validation = config
            self.topology_validator = TopologyValidator(full_config)
            
    def run_validation(self, breakpoints: List[Breakpoint], 
                      window_metrics_by_contig: Optional[Dict[str, List[WindowMetrics]]] = None,
                      bam_parser = None, **context) -> Dict[str, Any]:
        """Run all validation modules including topology analysis."""
        results = {}
        adjusted_breakpoints = breakpoints
        
        # Run topology validation first (it can adjust breakpoint confidences)
        if self.topology_validator and window_metrics_by_contig and bam_parser:
            self.logger.info("Running topology analysis")
            try:
                topology_results = self.topology_validator.validate_with_topology(
                    breakpoints, window_metrics_by_contig, bam_parser
                )
                results["TopologyValidator"] = topology_results
                
                # Use adjusted breakpoints for other validators
                adjusted_breakpoints = topology_results.get("adjusted_breakpoints", breakpoints)
                
            except Exception as e:
                self.logger.error(f"Error in topology validation: {e}")
                results["TopologyValidator"] = {"error": str(e)}
        
        # Run traditional validation modules
        for module in self.modules:
            module_name = module.__class__.__name__
            self.logger.info(f"Running {module_name}")
            
            try:
                module_results = module.validate(adjusted_breakpoints, **context)
                results[module_name] = module_results
            except Exception as e:
                self.logger.error(f"Error in {module_name}: {e}")
                results[module_name] = [{"error": str(e)}]
        
        # Add final breakpoint list to results
        results["final_breakpoints"] = [
            {
                "contig": bp.contig,
                "position": bp.position,
                "confidence": bp.confidence,
                "evidence": bp.evidence
            }
            for bp in adjusted_breakpoints
        ]
                
        return results