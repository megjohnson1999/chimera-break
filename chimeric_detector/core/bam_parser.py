"""BAM file parsing and read-pair extraction."""

import pysam
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Iterator, Generator
from collections import defaultdict
import logging
import mmap
import os
from pathlib import Path
from functools import lru_cache

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
    """Handles BAM file parsing with configurable quality filters and streaming support."""
    
    def __init__(self, bam_path: str, config: QualityConfig, use_streaming: bool = True, 
                 chunk_size: int = 10000, enable_mmap: bool = True):
        self.bam_path = bam_path
        self.config = config
        self.use_streaming = use_streaming
        self.chunk_size = chunk_size
        self.enable_mmap = enable_mmap
        self.logger = logging.getLogger(__name__)
        self._bam = None
        self._insert_size_stats = {}
        self._mmap_file = None
        self._file_size = 0
        self._coverage_cache = {}  # Simple instance-level cache
        
        # Check file size for memory management decisions
        if os.path.exists(bam_path):
            self._file_size = os.path.getsize(bam_path)
            if self._file_size > 1024 * 1024 * 1024:  # > 1GB
                self.logger.info(f"Large BAM file detected ({self._file_size / (1024**3):.1f}GB), enabling streaming mode")
                self.use_streaming = True
        
    def __enter__(self):
        """Open BAM file with comprehensive error handling."""
        try:
            self._bam = pysam.AlignmentFile(self.bam_path, "rb")
            
            # Validate BAM file has references
            if not self._bam.references:
                raise ValueError("BAM file contains no reference sequences")
                
            # Check if BAM file is indexed (required for random access)
            try:
                # Test if index exists by attempting to get stats
                self._bam.check_index()
            except (OSError, ValueError):
                logging.getLogger(__name__).warning(
                    f"BAM file {self.bam_path} may not be properly indexed. "
                    "Some operations may be slow or fail."
                )
                
        except FileNotFoundError:
            raise FileNotFoundError(f"BAM file not found: {self.bam_path}")
        except PermissionError:
            raise PermissionError(f"Permission denied reading BAM file: {self.bam_path}")
        except pysam.utils.SamtoolsError as e:
            raise ValueError(f"Invalid or corrupted BAM file {self.bam_path}: {e}")
        except Exception as e:
            raise RuntimeError(f"Unexpected error opening BAM file {self.bam_path}: {e}")
            
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close BAM file with error handling."""
        if self._bam:
            try:
                self._bam.close()
            except Exception as e:
                logging.getLogger(__name__).warning(f"Error closing BAM file: {e}")
        
        # Don't suppress exceptions from the with block
        return False
            
    def get_contigs(self) -> List[str]:
        """Get list of reference contigs from BAM header with error handling."""
        try:
            if not self._bam:
                raise RuntimeError("BAM file not opened")
            
            contigs = list(self._bam.references)
            if not contigs:
                raise ValueError("No contigs found in BAM file")
                
            return contigs
            
        except Exception as e:
            raise RuntimeError(f"Error reading contigs from BAM file: {e}")
    
    def estimate_insert_size_distribution(self, contig: str, sample_size: int = 10000) -> Dict[str, float]:
        """Estimate insert size distribution from properly paired reads with error handling."""
        try:
            if not self._bam:
                raise RuntimeError("BAM file not opened")
                
            # Validate contig exists
            if contig not in self._bam.references:
                raise ValueError(f"Contig '{contig}' not found in BAM file")
                
            insert_sizes = []
            
            try:
                for read in self._bam.fetch(contig):
                    if len(insert_sizes) >= sample_size:
                        break
                        
                    if self._passes_quality_filters(read) and read.is_proper_pair:
                        # Validate template length
                        if read.template_length is not None and read.template_length != 0:
                            insert_size = abs(read.template_length)
                            # Sanity check for realistic insert sizes
                            if 50 <= insert_size <= 50000:
                                insert_sizes.append(insert_size)
                                
            except (OSError, ValueError) as e:
                raise RuntimeError(f"Error reading from contig {contig}: {e}")
            
            if len(insert_sizes) < 100:
                self.logger.warning(f"Only {len(insert_sizes)} proper pairs found for {contig}, using defaults")
                return {"median": 500, "mad": 100, "mean": 500, "std": 100}
                
        except Exception as e:
            self.logger.error(f"Error estimating insert size for {contig}: {e}")
            # Return safe defaults to allow processing to continue
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
        """Extract read pair information for a genomic window with edge case handling."""
        if self.use_streaming:
            # Use streaming approach for large files
            return list(self.stream_read_pairs_in_window(contig, start, end))
        else:
            # Use traditional approach for smaller files
            return self._get_read_pairs_batch(contig, start, end)
    
    def stream_read_pairs_in_window(self, contig: str, start: int, end: int) -> Generator[ReadPairInfo, None, None]:
        """Stream read pair information for a genomic window to reduce memory usage."""
        try:
            if not self._bam:
                raise RuntimeError("BAM file not opened")
                
            # Validate parameters
            if start < 0:
                start = 0
            if end <= start:
                self.logger.warning(f"Invalid window coordinates: start={start}, end={end}")
                return
                
            # Validate contig exists
            if contig not in self._bam.references:
                raise ValueError(f"Contig '{contig}' not found in BAM file")
                
            contig_length = self._bam.get_reference_length(contig)
            if start >= contig_length:
                return  # Window is beyond contig end
                
            # Clamp end to contig length
            end = min(end, contig_length)
            
            try:
                # Process reads in chunks to control memory usage
                for read in self._bam.fetch(contig, start, end):
                    try:
                        if not self._passes_quality_filters(read):
                            continue
                            
                        if read.is_paired and not read.is_secondary and not read.is_supplementary:
                            pair_info = self._extract_pair_info(read)
                            if pair_info:
                                yield pair_info
                    except Exception as e:
                        self.logger.debug(f"Error processing read in {contig}:{start}-{end}: {e}")
                        continue
                        
            except (OSError, ValueError) as e:
                self.logger.error(f"Error fetching reads from {contig}:{start}-{end}: {e}")
                return
                
        except Exception as e:
            self.logger.error(f"Error streaming read pairs in window {contig}:{start}-{end}: {e}")
            return
    
    def _get_read_pairs_batch(self, contig: str, start: int, end: int) -> List[ReadPairInfo]:
        """Traditional batch processing for smaller windows."""
        try:
            if not self._bam:
                raise RuntimeError("BAM file not opened")
                
            # Validate parameters
            if start < 0:
                start = 0
            if end <= start:
                self.logger.warning(f"Invalid window coordinates: start={start}, end={end}")
                return []
                
            # Validate contig exists
            if contig not in self._bam.references:
                raise ValueError(f"Contig '{contig}' not found in BAM file")
                
            contig_length = self._bam.get_reference_length(contig)
            if start >= contig_length:
                return []  # Window is beyond contig end
                
            # Clamp end to contig length
            end = min(end, contig_length)
            
            pairs = []
            
            try:
                for read in self._bam.fetch(contig, start, end):
                    try:
                        if not self._passes_quality_filters(read):
                            continue
                            
                        if read.is_paired and not read.is_secondary and not read.is_supplementary:
                            pair_info = self._extract_pair_info(read)
                            if pair_info:
                                pairs.append(pair_info)
                    except Exception as e:
                        self.logger.debug(f"Error processing read in {contig}:{start}-{end}: {e}")
                        continue
                        
            except (OSError, ValueError) as e:
                self.logger.error(f"Error fetching reads from {contig}:{start}-{end}: {e}")
                return []
                
            return pairs
            
        except Exception as e:
            self.logger.error(f"Error getting read pairs in window {contig}:{start}-{end}: {e}")
            return []
    
    def iterate_windows(self, contig: str, window_size: int, step: int) -> Iterator[Tuple[int, int, List[ReadPairInfo]]]:
        """Iterate over windows with read pair data with comprehensive edge case handling."""
        try:
            if not self._bam:
                raise RuntimeError("BAM file not opened")
                
            # Validate parameters
            if window_size <= 0:
                raise ValueError(f"Window size must be positive, got {window_size}")
            if step <= 0:
                raise ValueError(f"Step size must be positive, got {step}")
                
            # Validate contig exists
            if contig not in self._bam.references:
                raise ValueError(f"Contig '{contig}' not found in BAM file")
                
            try:
                contig_length = self._bam.get_reference_length(contig)
            except Exception as e:
                raise RuntimeError(f"Error getting length for contig {contig}: {e}")
                
            if contig_length <= 0:
                self.logger.warning(f"Contig {contig} has zero or negative length: {contig_length}")
                return
                
            # Handle case where window is larger than contig
            if window_size > contig_length:
                self.logger.warning(f"Window size ({window_size}) larger than contig {contig} length ({contig_length})")
                # Still process the entire contig as one window
                pairs = self.get_read_pairs_in_window(contig, 0, contig_length)
                yield 0, contig_length, pairs
                return
                
            # Iterate through full windows
            for start in range(0, contig_length - window_size + 1, step):
                end = start + window_size
                try:
                    pairs = self.get_read_pairs_in_window(contig, start, end)
                    yield start, end, pairs
                except Exception as e:
                    self.logger.warning(f"Error processing window {contig}:{start}-{end}: {e}")
                    continue
                    
            # Handle final partial window if it exists and is significant
            if step < window_size:  # Only if windows overlap
                final_start = ((contig_length - window_size) // step + 1) * step
                if final_start < contig_length and contig_length - final_start >= window_size // 2:
                    try:
                        pairs = self.get_read_pairs_in_window(contig, final_start, contig_length)
                        yield final_start, contig_length, pairs
                    except Exception as e:
                        self.logger.warning(f"Error processing final window {contig}:{final_start}-{contig_length}: {e}")
                        
        except Exception as e:
            self.logger.error(f"Error iterating windows for {contig}: {e}")
            return
    
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
        
        # Calculate read length for manual insert size calculation
        read_length = len(read.query_sequence) if read.query_sequence else read.query_length
        
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
        """Calculate average coverage in a window with caching and memory optimization."""
        # Simple cache key
        cache_key = f"{contig}:{start}:{end}"
        if cache_key in self._coverage_cache:
            return self._coverage_cache[cache_key]
            
        try:
            if not self._bam:
                raise RuntimeError("BAM file not opened")
                
            window_size = end - start
            if window_size <= 0:
                return 0.0
                
            # For large windows, use sampling to reduce memory usage
            if window_size > 50000:  # 50kb threshold
                result = self._calculate_coverage_sampled(contig, start, end)
            else:
                result = self._calculate_coverage_full(contig, start, end)
            
            # Cache result (limit cache size)
            if len(self._coverage_cache) < 1000:
                self._coverage_cache[cache_key] = result
                
            return result
                
        except Exception as e:
            self.logger.error(f"Error calculating coverage in {contig}:{start}-{end}: {e}")
            return 0.0
    
    def _calculate_coverage_full(self, contig: str, start: int, end: int) -> float:
        """Calculate full coverage for smaller windows."""
        try:
            coverage = np.zeros(end - start, dtype=np.uint16)  # Use smaller data type
            
            for pileup_column in self._bam.pileup(contig, start, end, truncate=True):
                pos = pileup_column.reference_pos - start
                if 0 <= pos < len(coverage):
                    coverage[pos] = min(pileup_column.nsegments, 65535)  # Cap at uint16 max
                    
            return float(np.mean(coverage))
            
        except Exception as e:
            self.logger.warning(f"Error in full coverage calculation: {e}")
            return 0.0
    
    def _calculate_coverage_sampled(self, contig: str, start: int, end: int) -> float:
        """Calculate coverage using sampling for large windows."""
        try:
            # Sample every nth position to reduce memory usage
            window_size = end - start
            sample_rate = max(1, window_size // 10000)  # Sample ~10k positions max
            
            coverage_sum = 0
            sample_count = 0
            
            for pileup_column in self._bam.pileup(contig, start, end, truncate=True):
                pos = pileup_column.reference_pos - start
                if pos % sample_rate == 0:  # Sample every nth position
                    coverage_sum += pileup_column.nsegments
                    sample_count += 1
                    
            return float(coverage_sum / sample_count) if sample_count > 0 else 0.0
            
        except Exception as e:
            self.logger.warning(f"Error in sampled coverage calculation: {e}")
            return 0.0
    
    def get_memory_usage(self) -> Dict[str, float]:
        """Get current memory usage statistics."""
        import psutil
        process = psutil.Process()
        memory_info = process.memory_info()
        
        return {
            "rss_mb": memory_info.rss / (1024 * 1024),
            "vms_mb": memory_info.vms / (1024 * 1024),
            "file_size_mb": self._file_size / (1024 * 1024)
        }