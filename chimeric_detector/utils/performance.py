"""Performance optimization utilities and monitoring."""

import os
import time
import psutil
import multiprocessing as mp
from typing import Dict, List, Optional, Callable, Any, Iterator, Tuple
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from contextlib import contextmanager
import logging
from pathlib import Path
import queue
import threading
from functools import wraps

from ..config.config import Config
from .progress import ProgressTracker, ProgressReporter


@dataclass
class PerformanceMetrics:
    """Container for performance metrics."""
    start_time: float = field(default_factory=time.time)
    end_time: Optional[float] = None
    peak_memory_mb: float = 0.0
    cpu_percent: float = 0.0
    items_processed: int = 0
    errors_encountered: int = 0
    
    @property
    def elapsed_time(self) -> float:
        """Get elapsed time in seconds."""
        end = self.end_time or time.time()
        return end - self.start_time
    
    @property
    def throughput(self) -> float:
        """Get items processed per second."""
        elapsed = self.elapsed_time
        return self.items_processed / elapsed if elapsed > 0 else 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "elapsed_time": self.elapsed_time,
            "peak_memory_mb": self.peak_memory_mb,
            "cpu_percent": self.cpu_percent,
            "items_processed": self.items_processed,
            "errors_encountered": self.errors_encountered,
            "throughput_per_sec": self.throughput
        }


class PerformanceMonitor:
    """Monitor system performance during processing."""
    
    def __init__(self, sample_interval: float = 1.0):
        self.sample_interval = sample_interval
        self.process = psutil.Process()
        self.metrics = PerformanceMetrics()
        self._monitoring = False
        self._monitor_thread = None
        self.logger = logging.getLogger(__name__)
    
    def start_monitoring(self):
        """Start performance monitoring in background thread."""
        if self._monitoring:
            return
            
        self._monitoring = True
        self._monitor_thread = threading.Thread(target=self._monitor_loop, daemon=True)
        self._monitor_thread.start()
        self.logger.debug("Performance monitoring started")
    
    def stop_monitoring(self):
        """Stop performance monitoring."""
        if not self._monitoring:
            return
            
        self._monitoring = False
        if self._monitor_thread:
            self._monitor_thread.join(timeout=2.0)
        self.metrics.end_time = time.time()
        self.logger.debug("Performance monitoring stopped")
    
    def _monitor_loop(self):
        """Background monitoring loop."""
        while self._monitoring:
            try:
                memory_info = self.process.memory_info()
                memory_mb = memory_info.rss / (1024 * 1024)
                self.metrics.peak_memory_mb = max(self.metrics.peak_memory_mb, memory_mb)
                
                # Sample CPU usage
                self.metrics.cpu_percent = self.process.cpu_percent()
                
                time.sleep(self.sample_interval)
                
            except Exception as e:
                self.logger.warning(f"Error in performance monitoring: {e}")
                time.sleep(self.sample_interval)
    
    def record_item_processed(self):
        """Record that an item was processed."""
        self.metrics.items_processed += 1
    
    def record_error(self):
        """Record that an error occurred."""
        self.metrics.errors_encountered += 1
    
    @contextmanager
    def monitoring_context(self):
        """Context manager for automatic monitoring."""
        self.start_monitoring()
        try:
            yield self
        finally:
            self.stop_monitoring()


@dataclass
class WorkItem:
    """Container for work items in parallel processing."""
    id: str
    contig: str
    args: Tuple[Any, ...]
    kwargs: Dict[str, Any] = field(default_factory=dict)


class ParallelProcessor:
    """Parallel processing manager for contig-level parallelism."""
    
    def __init__(self, config: Config, max_workers: Optional[int] = None):
        self.config = config
        self.max_workers = max_workers or min(mp.cpu_count(), 8)  # Cap at 8 for memory
        self.logger = logging.getLogger(__name__)
        
        # Adjust workers based on available memory
        available_memory_gb = psutil.virtual_memory().available / (1024**3)
        memory_per_worker = 2.0  # Estimate 2GB per worker
        max_memory_workers = int(available_memory_gb / memory_per_worker)
        
        if max_memory_workers < self.max_workers:
            self.logger.info(f"Reducing workers from {self.max_workers} to {max_memory_workers} due to memory constraints")
            self.max_workers = max(1, max_memory_workers)
    
    def process_contigs_parallel(self, work_items: List[WorkItem], 
                               process_func: Callable, 
                               progress_tracker: Optional[ProgressTracker] = None) -> Dict[str, Any]:
        """Process contigs in parallel using multiprocessing."""
        results = {}
        failed_items = []
        
        # Use process pool for CPU-bound work
        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all work items
            future_to_item = {}
            for item in work_items:
                future = executor.submit(process_func, *item.args, **item.kwargs)
                future_to_item[future] = item
            
            # Progress tracking
            if progress_tracker:
                reporter = ProgressReporter(
                    total=len(work_items),
                    description="Processing contigs in parallel"
                )
            
            # Collect results as they complete
            for future in as_completed(future_to_item):
                item = future_to_item[future]
                
                try:
                    result = future.result()
                    results[item.contig] = result
                    
                    if progress_tracker:
                        progress_tracker.complete_contig(item.contig, result)
                        reporter.update(1)
                        
                except Exception as e:
                    self.logger.error(f"Error processing {item.contig}: {e}")
                    failed_items.append(item)
                    
                    if progress_tracker:
                        progress_tracker.fail_contig(item.contig, str(e))
                        reporter.update(1)
        
        # Retry failed items sequentially
        if failed_items:
            self.logger.info(f"Retrying {len(failed_items)} failed items sequentially")
            for item in failed_items:
                try:
                    result = process_func(*item.args, **item.kwargs)
                    results[item.contig] = result
                    
                    if progress_tracker:
                        progress_tracker.complete_contig(item.contig, result)
                        
                except Exception as e:
                    self.logger.error(f"Sequential retry also failed for {item.contig}: {e}")
        
        return results
    
    def get_optimal_chunk_size(self, total_items: int) -> int:
        """Calculate optimal chunk size for parallel processing."""
        # Aim for 2-4 chunks per worker
        target_chunks = self.max_workers * 3
        chunk_size = max(1, total_items // target_chunks)
        
        # Cap chunk size to prevent memory issues
        max_chunk_size = 1000
        return min(chunk_size, max_chunk_size)


class WorkQueue:
    """Thread-safe work queue for parallel processing."""
    
    def __init__(self, maxsize: int = 0):
        self.queue = queue.Queue(maxsize=maxsize)
        self.results = queue.Queue()
        self.errors = queue.Queue()
        self._stop_event = threading.Event()
        self.logger = logging.getLogger(__name__)
    
    def add_work(self, item: WorkItem):
        """Add work item to queue."""
        self.queue.put(item)
    
    def add_sentinel(self):
        """Add sentinel to signal end of work."""
        self.queue.put(None)
    
    def get_work(self, timeout: Optional[float] = None) -> Optional[WorkItem]:
        """Get work item from queue."""
        try:
            return self.queue.get(timeout=timeout)
        except queue.Empty:
            return None
    
    def put_result(self, contig: str, result: Any):
        """Put result in results queue."""
        self.results.put((contig, result))
    
    def put_error(self, contig: str, error: Exception):
        """Put error in errors queue."""
        self.errors.put((contig, error))
    
    def get_results(self) -> List[Tuple[str, Any]]:
        """Get all results from queue."""
        results = []
        while not self.results.empty():
            try:
                results.append(self.results.get_nowait())
            except queue.Empty:
                break
        return results
    
    def get_errors(self) -> List[Tuple[str, Exception]]:
        """Get all errors from queue."""
        errors = []
        while not self.errors.empty():
            try:
                errors.append(self.errors.get_nowait())
            except queue.Empty:
                break
        return errors
    
    def stop(self):
        """Signal workers to stop."""
        self._stop_event.set()
    
    def is_stopped(self) -> bool:
        """Check if stop signal was sent."""
        return self._stop_event.is_set()


def profile_function(func: Callable) -> Callable:
    """Decorator to profile function execution."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        monitor = PerformanceMonitor()
        with monitor.monitoring_context():
            result = func(*args, **kwargs)
            monitor.record_item_processed()
        
        # Log performance metrics
        metrics = monitor.metrics
        logger = logging.getLogger(func.__module__)
        logger.info(f"{func.__name__} performance: "
                   f"{metrics.elapsed_time:.2f}s, "
                   f"{metrics.peak_memory_mb:.1f}MB peak, "
                   f"{metrics.throughput:.1f} items/sec")
        
        return result
    return wrapper


def batch_iterator(items: List[Any], batch_size: int) -> Iterator[List[Any]]:
    """Create batches from a list of items."""
    for i in range(0, len(items), batch_size):
        yield items[i:i + batch_size]


def memory_efficient_chunker(items: Iterator[Any], chunk_size: int) -> Iterator[List[Any]]:
    """Create memory-efficient chunks from an iterator."""
    chunk = []
    for item in items:
        chunk.append(item)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    
    if chunk:  # Yield remaining items
        yield chunk


class ResourceManager:
    """Manage system resources during processing."""
    
    def __init__(self, memory_limit_gb: Optional[float] = None):
        self.memory_limit_gb = memory_limit_gb
        self.logger = logging.getLogger(__name__)
        
        if not memory_limit_gb:
            # Set to 80% of available memory
            available_gb = psutil.virtual_memory().available / (1024**3)
            self.memory_limit_gb = available_gb * 0.8
    
    def check_memory_usage(self) -> bool:
        """Check if memory usage is within limits."""
        current_memory_gb = psutil.Process().memory_info().rss / (1024**3)
        
        if current_memory_gb > self.memory_limit_gb:
            self.logger.warning(f"Memory usage ({current_memory_gb:.1f}GB) exceeds limit ({self.memory_limit_gb:.1f}GB)")
            return False
        
        return True
    
    def suggest_gc(self) -> bool:
        """Suggest garbage collection if memory usage is high."""
        current_memory_gb = psutil.Process().memory_info().rss / (1024**3)
        threshold = self.memory_limit_gb * 0.7
        
        if current_memory_gb > threshold:
            self.logger.debug(f"Memory usage ({current_memory_gb:.1f}GB) above threshold, suggesting GC")
            return True
        
        return False
    
    def get_memory_stats(self) -> Dict[str, float]:
        """Get current memory statistics."""
        process = psutil.Process()
        memory_info = process.memory_info()
        system_memory = psutil.virtual_memory()
        
        return {
            "process_rss_gb": memory_info.rss / (1024**3),
            "process_vms_gb": memory_info.vms / (1024**3),
            "system_available_gb": system_memory.available / (1024**3),
            "system_used_percent": system_memory.percent,
            "memory_limit_gb": self.memory_limit_gb
        }


@contextmanager
def performance_context(description: str = "Operation"):
    """Context manager for performance monitoring."""
    monitor = PerformanceMonitor()
    logger = logging.getLogger(__name__)
    
    logger.info(f"Starting {description}")
    
    with monitor.monitoring_context():
        start_time = time.time()
        try:
            yield monitor
        finally:
            elapsed = time.time() - start_time
            metrics = monitor.metrics
            
            logger.info(f"Completed {description}: "
                       f"{elapsed:.2f}s elapsed, "
                       f"{metrics.peak_memory_mb:.1f}MB peak memory, "
                       f"{metrics.items_processed} items processed")
            
            if metrics.errors_encountered > 0:
                logger.warning(f"{metrics.errors_encountered} errors encountered during {description}")