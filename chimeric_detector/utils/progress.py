"""Progress tracking and checkpointing utilities."""

import json
import time
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List
from dataclasses import dataclass, asdict
from datetime import datetime

from .security import safe_file_open, create_secure_temp_file


@dataclass
class ProgressState:
    """Represents the current progress state."""
    total_contigs: int
    processed_contigs: int
    current_contig: Optional[str]
    start_time: float
    last_update: float
    completed_contigs: List[str]
    failed_contigs: List[str]
    partial_results: Dict[str, Any]
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ProgressState":
        """Create from dictionary."""
        return cls(**data)


class ProgressTracker:
    """Tracks processing progress and provides checkpointing."""
    
    def __init__(self, total_contigs: int, checkpoint_file: Optional[Path] = None):
        self.logger = logging.getLogger(__name__)
        self.checkpoint_file = checkpoint_file
        
        self.state = ProgressState(
            total_contigs=total_contigs,
            processed_contigs=0,
            current_contig=None,
            start_time=time.time(),
            last_update=time.time(),
            completed_contigs=[],
            failed_contigs=[],
            partial_results={}
        )
        
        # Load existing checkpoint if available
        if self.checkpoint_file and self.checkpoint_file.exists():
            self._load_checkpoint()
    
    def start_contig(self, contig: str) -> None:
        """Mark the start of processing a contig."""
        self.state.current_contig = contig
        self.state.last_update = time.time()
        
        self.logger.info(f"Starting contig {contig} ({self.state.processed_contigs + 1}/{self.state.total_contigs})")
        self._report_progress()
    
    def complete_contig(self, contig: str, results: Any = None) -> None:
        """Mark completion of a contig."""
        if contig != self.state.current_contig:
            self.logger.warning(f"Contig mismatch: expected {self.state.current_contig}, got {contig}")
        
        self.state.processed_contigs += 1
        self.state.completed_contigs.append(contig)
        self.state.current_contig = None
        self.state.last_update = time.time()
        
        if results is not None:
            self.state.partial_results[contig] = results
        
        self.logger.info(f"Completed contig {contig}")
        self._report_progress()
        
        # Save checkpoint periodically
        if self.checkpoint_file and self.state.processed_contigs % 10 == 0:
            self._save_checkpoint()
    
    def fail_contig(self, contig: str, error: str) -> None:
        """Mark failure of a contig."""
        if contig != self.state.current_contig:
            self.logger.warning(f"Contig mismatch: expected {self.state.current_contig}, got {contig}")
        
        self.state.processed_contigs += 1
        self.state.failed_contigs.append(contig)
        self.state.current_contig = None
        self.state.last_update = time.time()
        
        self.logger.error(f"Failed to process contig {contig}: {error}")
        self._report_progress()
        
        # Save checkpoint after failures
        if self.checkpoint_file:
            self._save_checkpoint()
    
    def is_complete(self) -> bool:
        """Check if all contigs have been processed."""
        return self.state.processed_contigs >= self.state.total_contigs
    
    def get_remaining_contigs(self, all_contigs: List[str]) -> List[str]:
        """Get list of contigs that still need processing."""
        processed = set(self.state.completed_contigs + self.state.failed_contigs)
        return [c for c in all_contigs if c not in processed]
    
    def get_results(self) -> Dict[str, Any]:
        """Get accumulated results."""
        return self.state.partial_results.copy()
    
    def get_summary(self) -> Dict[str, Any]:
        """Get progress summary."""
        elapsed = time.time() - self.state.start_time
        
        return {
            "total_contigs": self.state.total_contigs,
            "completed": len(self.state.completed_contigs),
            "failed": len(self.state.failed_contigs),
            "remaining": self.state.total_contigs - self.state.processed_contigs,
            "success_rate": len(self.state.completed_contigs) / max(self.state.processed_contigs, 1),
            "elapsed_time": elapsed,
            "estimated_remaining": self._estimate_remaining_time(),
            "current_contig": self.state.current_contig
        }
    
    def _report_progress(self) -> None:
        """Report current progress."""
        summary = self.get_summary()
        
        progress_pct = (self.state.processed_contigs / self.state.total_contigs) * 100
        
        self.logger.info(
            f"Progress: {self.state.processed_contigs}/{self.state.total_contigs} ({progress_pct:.1f}%) "
            f"- Success: {summary['completed']}, Failed: {summary['failed']}"
        )
        
        if summary['estimated_remaining'] > 0:
            remaining_str = self._format_time(summary['estimated_remaining'])
            self.logger.info(f"Estimated time remaining: {remaining_str}")
    
    def _estimate_remaining_time(self) -> float:
        """Estimate remaining processing time."""
        if self.state.processed_contigs == 0:
            return 0.0
        
        elapsed = time.time() - self.state.start_time
        rate = self.state.processed_contigs / elapsed
        remaining_contigs = self.state.total_contigs - self.state.processed_contigs
        
        return remaining_contigs / rate if rate > 0 else 0.0
    
    def _format_time(self, seconds: float) -> str:
        """Format time duration as human-readable string."""
        if seconds < 60:
            return f"{seconds:.0f}s"
        elif seconds < 3600:
            return f"{seconds/60:.1f}m"
        else:
            return f"{seconds/3600:.1f}h"
    
    def _save_checkpoint(self) -> None:
        """Save current state to checkpoint file."""
        if not self.checkpoint_file:
            return
            
        try:
            checkpoint_data = {
                "timestamp": datetime.now().isoformat(),
                "state": self.state.to_dict()
            }
            
            # Write atomically using temporary file
            temp_file = create_secure_temp_file(
                suffix=".json",
                prefix="checkpoint_",
                directory=self.checkpoint_file.parent
            )
            
            try:
                with safe_file_open(temp_file, 'w') as f:
                    json.dump(checkpoint_data, f, indent=2)
                
                # Move to final location
                temp_file.replace(self.checkpoint_file)
                self.logger.debug(f"Saved checkpoint to {self.checkpoint_file}")
                
            except Exception:
                # Clean up temp file on error
                if temp_file.exists():
                    temp_file.unlink()
                raise
                
        except Exception as e:
            self.logger.warning(f"Failed to save checkpoint: {e}")
    
    def _load_checkpoint(self) -> None:
        """Load state from checkpoint file."""
        try:
            with safe_file_open(self.checkpoint_file, 'r') as f:
                checkpoint_data = json.load(f)
            
            if "state" in checkpoint_data:
                self.state = ProgressState.from_dict(checkpoint_data["state"])
                self.logger.info(f"Loaded checkpoint from {self.checkpoint_file}")
                self.logger.info(f"Resuming from {self.state.processed_contigs}/{self.state.total_contigs} processed")
            else:
                self.logger.warning("Invalid checkpoint file format")
                
        except Exception as e:
            self.logger.warning(f"Failed to load checkpoint: {e}")
    
    def cleanup_checkpoint(self) -> None:
        """Remove checkpoint file after successful completion."""
        if self.checkpoint_file and self.checkpoint_file.exists():
            try:
                self.checkpoint_file.unlink()
                self.logger.info("Removed checkpoint file after successful completion")
            except Exception as e:
                self.logger.warning(f"Failed to remove checkpoint file: {e}")


class ProgressReporter:
    """Simple progress reporter for real-time updates."""
    
    def __init__(self, total: int, description: str = "Processing"):
        self.total = total
        self.current = 0
        self.description = description
        self.start_time = time.time()
        self.last_report = 0
        self.logger = logging.getLogger(__name__)
    
    def update(self, increment: int = 1) -> None:
        """Update progress."""
        self.current += increment
        
        # Report every 5% or every 10 items, whichever is less frequent
        report_interval = max(1, min(10, self.total // 20))
        
        if self.current - self.last_report >= report_interval or self.current >= self.total:
            self._report()
            self.last_report = self.current
    
    def _report(self) -> None:
        """Report current progress."""
        if self.total > 0:
            pct = (self.current / self.total) * 100
            self.logger.info(f"{self.description}: {self.current}/{self.total} ({pct:.1f}%)")