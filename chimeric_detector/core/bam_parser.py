"""BAM file parsing and read-pair extraction."""

import pysam
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Iterator
from collections import defaultdict
import logging

from ..config.config import QualityConfig


@dataclass
class ReadPairInfo:
    """Container for read pair information."""
    contig: str
    position: int
    insert_size: int
    is_proper_pair: bool
    is_discordant: bool
    mapping_quality: int
    mate_contig: Optional[str] = None
    mate_position: Optional[int] = None
    is_near_end: bool = False
    mate_unmapped: bool = False
    

class BAMParser:
    """Handles BAM file parsing with configurable quality filters."""
    
    def __init__(self, bam_path: str, config: QualityConfig):
        self.bam_path = bam_path
        self.config = config
        self.logger = logging.getLogger(__name__)
        self._bam = None
        self._insert_size_stats = {}
        
    def __enter__(self):
        self._bam = pysam.AlignmentFile(self.bam_path, "rb")
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._bam:
            self._bam.close()
            
    def get_contigs(self) -> List[str]:
        """Get list of reference contigs from BAM header."""
        return list(self._bam.references)
    
    def estimate_insert_size_distribution(self, contig: str, sample_size: int = 10000) -> Dict[str, float]:
        """Estimate insert size distribution from properly paired reads."""
        insert_sizes = []
        
        for read in self._bam.fetch(contig):
            if len(insert_sizes) >= sample_size:
                break
                
            if self._passes_quality_filters(read) and read.is_proper_pair:
                insert_sizes.append(abs(read.template_length))
        
        if len(insert_sizes) < 100:
            self.logger.warning(f"Only {len(insert_sizes)} proper pairs found for {contig}")
            return {"median": 500, "mad": 100, "mean": 500, "std": 100}
        
        insert_sizes = np.array(insert_sizes)
        
        # Use robust statistics
        median = np.median(insert_sizes)
        mad = np.median(np.abs(insert_sizes - median))
        
        # Also calculate mean/std for comparison
        mean = np.mean(insert_sizes)
        std = np.std(insert_sizes)
        
        stats = {
            "median": float(median),
            "mad": float(mad),
            "mean": float(mean),
            "std": float(std),
            "q1": float(np.percentile(insert_sizes, 25)),
            "q3": float(np.percentile(insert_sizes, 75))
        }
        
        self._insert_size_stats[contig] = stats
        return stats
    
    def get_read_pairs_in_window(self, contig: str, start: int, end: int) -> List[ReadPairInfo]:
        """Extract read pair information for a genomic window."""
        pairs = []
        
        for read in self._bam.fetch(contig, start, end):
            if not self._passes_quality_filters(read):
                continue
                
            if read.is_paired and not read.is_secondary and not read.is_supplementary:
                pair_info = self._extract_pair_info(read)
                if pair_info:
                    pairs.append(pair_info)
                    
        return pairs
    
    def iterate_windows(self, contig: str, window_size: int, step: int) -> Iterator[Tuple[int, int, List[ReadPairInfo]]]:
        """Iterate over windows with read pair data."""
        contig_length = self._bam.get_reference_length(contig)
        
        # Iterate through full windows
        for start in range(0, contig_length - window_size + 1, step):
            end = start + window_size
            pairs = self.get_read_pairs_in_window(contig, start, end)
            yield start, end, pairs
            
        # Handle final partial window if it exists and is significant
        final_start = ((contig_length - window_size) // step + 1) * step
        if final_start < contig_length and contig_length - final_start >= window_size // 2:
            pairs = self.get_read_pairs_in_window(contig, final_start, contig_length)
            yield final_start, contig_length, pairs
    
    def _passes_quality_filters(self, read) -> bool:
        """Check if read passes quality filters."""
        if read.mapping_quality < self.config.min_mapping_quality:
            return False
            
        if self.config.require_proper_pairs and not read.is_proper_pair:
            return False
            
        if read.is_paired and read.template_length != 0:
            insert_size = abs(read.template_length)
            if insert_size < self.config.min_insert_size or insert_size > self.config.max_insert_size:
                return False
                
        return True
    
    def _extract_pair_info(self, read) -> Optional[ReadPairInfo]:
        """Extract pair information from a read, handling contig end cases."""
        if not read.is_paired:
            return None
            
        contig_length = self._bam.get_reference_length(read.reference_name)
        read_end_pos = read.reference_end if read.reference_end else read.reference_start + len(read.query_sequence or "")
        
        # Determine if read is near contig ends (configurable buffer)
        end_buffer = self.config.end_read_buffer
        is_near_start = read.reference_start < end_buffer
        is_near_end = read_end_pos > contig_length - end_buffer
        is_near_contig_end = is_near_start or is_near_end
        
        # Handle insert size calculation
        insert_size = 0
        mate_unmapped = read.mate_is_unmapped
        
        if not mate_unmapped and read.template_length != 0:
            insert_size = abs(read.template_length)
        elif not mate_unmapped and read.reference_name == read.next_reference_name:
            # Calculate insert size manually for edge cases
            if read.next_reference_start is not None:
                if read.is_reverse:
                    insert_size = abs(read.reference_start - read.next_reference_start) + read_length
                else:
                    insert_size = abs(read.next_reference_start - read.reference_start) + read_length
        
        # Determine discordance - more nuanced for end reads
        is_discordant = False
        mate_contig = None
        mate_position = None
        
        if mate_unmapped:
            is_discordant = True
        elif read.reference_name != read.next_reference_name:
            # Inter-contig pairs
            is_discordant = True
            mate_contig = read.next_reference_name
            mate_position = read.next_reference_start
        else:
            # Same contig - check for unusual orientation or distance
            mate_position = read.next_reference_start
            
            # For reads near contig ends, be more lenient about proper pair classification
            if is_near_contig_end:
                # Check if mate maps outside contig bounds
                if mate_position is not None:
                    if mate_position < 0 or mate_position >= contig_length:
                        is_discordant = True
                    elif insert_size > self.config.max_insert_size * 1.5:  # More lenient for end reads
                        is_discordant = True
            else:
                # Standard discordance checks for internal reads
                if not read.is_proper_pair:
                    # Additional checks for why it's not proper
                    if insert_size > self.config.max_insert_size:
                        is_discordant = True
                    elif read.is_reverse == read.mate_is_reverse:  # Same orientation
                        is_discordant = True
                        
        return ReadPairInfo(
            contig=read.reference_name,
            position=read.reference_start,
            insert_size=insert_size,
            is_proper_pair=read.is_proper_pair,
            is_discordant=is_discordant,
            mapping_quality=read.mapping_quality,
            mate_contig=mate_contig,
            mate_position=mate_position,
            is_near_end=is_near_contig_end,
            mate_unmapped=mate_unmapped
        )
    
    def calculate_coverage_in_window(self, contig: str, start: int, end: int) -> float:
        """Calculate average coverage in a window."""
        coverage = np.zeros(end - start)
        
        for pileup_column in self._bam.pileup(contig, start, end, truncate=True):
            pos = pileup_column.reference_pos - start
            if 0 <= pos < len(coverage):
                coverage[pos] = pileup_column.nsegments
                
        return np.mean(coverage)