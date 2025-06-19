"""Output formatting and reporting utilities."""

import json
import csv
import tempfile
import shutil
from pathlib import Path
from typing import List, Dict, Any, Optional
import logging
from datetime import datetime

from ..core.window_analysis import Breakpoint
from ..config.config import OutputConfig
from .security import safe_file_open, create_secure_temp_file, SecurityError


class OutputFormatter:
    """Handles output formatting for different formats."""
    
    def __init__(self, config: OutputConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def _atomic_write(self, content: str, output_path: Path) -> None:
        """Atomically write content to a file to prevent corruption on failure."""
        try:
            # Create temporary file in same directory as target
            temp_file = create_secure_temp_file(
                suffix=output_path.suffix,
                prefix=f"{output_path.stem}_",
                directory=output_path.parent
            )
            
            # Write to temporary file
            try:
                with safe_file_open(temp_file, 'w') as f:
                    f.write(content)
                    f.flush()  # Ensure data is written to disk
                
                # Atomically move temporary file to final location
                shutil.move(str(temp_file), str(output_path))
                self.logger.debug(f"Successfully wrote output to {output_path}")
                
            except Exception as e:
                # Clean up temporary file on error
                if temp_file.exists():
                    temp_file.unlink()
                raise
                
        except SecurityError as e:
            raise RuntimeError(f"Security error writing output: {e}")
        except PermissionError as e:
            raise RuntimeError(f"Permission denied writing to {output_path}: {e}")
        except OSError as e:
            raise RuntimeError(f"OS error writing to {output_path}: {e}")
        except Exception as e:
            raise RuntimeError(f"Unexpected error writing output: {e}")
        
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
            
        # Convert to JSON string and write atomically
        try:
            json_content = json.dumps(results, indent=2, ensure_ascii=False)
            self._atomic_write(json_content, output_path)
        except (TypeError, ValueError) as e:
            raise RuntimeError(f"Error serializing results to JSON: {e}")
            
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
            
        # Build TSV content in memory first
        try:
            from io import StringIO
            
            output = StringIO()
            writer = csv.writer(output, delimiter='\t')
            writer.writerow(headers)
            
            for bp in breakpoints:
                try:
                    row = [
                        bp.contig,
                        bp.position,
                        round(bp.confidence, 4),
                        round(bp.evidence.get("proper_pair_drop", 0), 4),
                        round(bp.evidence.get("insert_size_zscore", 0), 4),
                        round(bp.evidence.get("discordant_rate", 0), 4)
                    ]
                    
                    if self.config.include_read_evidence and bp.left_metrics and bp.right_metrics:
                        row.extend([
                            round(bp.left_metrics.coverage, 2),
                            round(bp.left_metrics.proper_pair_rate, 4),
                            round(bp.left_metrics.insert_size_median, 0),
                            round(bp.right_metrics.coverage, 2),
                            round(bp.right_metrics.proper_pair_rate, 4),
                            round(bp.right_metrics.insert_size_median, 0)
                        ])
                        
                    writer.writerow(row)
                except Exception as e:
                    self.logger.warning(f"Error formatting breakpoint {bp.contig}:{bp.position}: {e}")
                    continue
                    
            # Write atomically
            tsv_content = output.getvalue()
            self._atomic_write(tsv_content, output_path)
            
        except Exception as e:
            raise RuntimeError(f"Error generating TSV output: {e}")
                
    def _write_bed(self, breakpoints: List[Breakpoint], output_path: Path) -> None:
        """Write results in BED format."""
        try:
            bed_lines = []
            for bp in breakpoints:
                try:
                    # BED format: chr start end name score
                    start = max(0, bp.position - 50)
                    end = bp.position + 50
                    score = min(1000, max(0, int(bp.confidence * 1000)))  # Clamp score to valid range
                    
                    bed_lines.append(f"{bp.contig}\t{start}\t{end}\tbreakpoint\t{score}\n")
                except Exception as e:
                    self.logger.warning(f"Error formatting BED entry for {bp.contig}:{bp.position}: {e}")
                    continue
                    
            bed_content = "".join(bed_lines)
            self._atomic_write(bed_content, output_path)
            
        except Exception as e:
            raise RuntimeError(f"Error generating BED output: {e}")
                
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


class StructuredFormatter(logging.Formatter):
    """Custom formatter for structured logging with sanitization."""
    
    def __init__(self, sanitize_paths: bool = True):
        super().__init__()
        self.sanitize_paths = sanitize_paths
    
    def format(self, record):
        # Create structured log entry
        log_data = {
            "timestamp": datetime.fromtimestamp(record.created).isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
            "module": record.module,
            "function": record.funcName,
            "line": record.lineno
        }
        
        # Add thread/process info if available
        if hasattr(record, 'thread'):
            log_data["thread"] = record.thread
        if hasattr(record, 'process'):
            log_data["process"] = record.process
            
        # Add extra fields if present
        extra_fields = {}
        for key, value in record.__dict__.items():
            if key not in ['name', 'msg', 'args', 'levelname', 'levelno', 'pathname',
                          'filename', 'module', 'lineno', 'funcName', 'created', 
                          'msecs', 'relativeCreated', 'thread', 'threadName',
                          'processName', 'process', 'message']:
                extra_fields[key] = value
                
        if extra_fields:
            log_data["extra"] = extra_fields
            
        # Sanitize sensitive information
        if self.sanitize_paths:
            log_data["message"] = self._sanitize_message(log_data["message"])
            
        # Format as JSON for structured logging or simple format for console
        if hasattr(record, 'structured') and record.structured:
            return json.dumps(log_data, default=str)
        else:
            return f"{log_data['timestamp']} - {log_data['logger']} - {log_data['level']} - {log_data['message']}"
    
    def _sanitize_message(self, message: str) -> str:
        """Sanitize sensitive information from log messages."""
        import re
        
        # Remove or mask file paths that might contain sensitive info
        message = re.sub(r'/Users/[^/]+/', '/Users/***/', message)
        message = re.sub(r'/home/[^/]+/', '/home/***/', message)
        
        # Remove potential passwords or tokens (basic patterns)
        message = re.sub(r'password[=:]\s*\S+', 'password=***', message, flags=re.IGNORECASE)
        message = re.sub(r'token[=:]\s*\S+', 'token=***', message, flags=re.IGNORECASE)
        message = re.sub(r'key[=:]\s*\S+', 'key=***', message, flags=re.IGNORECASE)
        
        return message


class Logger:
    """Enhanced logging utility with structured logging and security features."""
    
    @staticmethod
    def setup_logging(verbose: bool = False, debug: bool = False, 
                     log_file: Optional[Path] = None, structured: bool = False,
                     sanitize_logs: bool = True) -> None:
        """Set up enhanced logging configuration."""
        # Determine log level
        level = logging.WARNING
        if debug:
            level = logging.DEBUG
        elif verbose:
            level = logging.INFO
        
        # Create formatters
        console_formatter = StructuredFormatter(sanitize_paths=sanitize_logs)
        file_formatter = StructuredFormatter(sanitize_paths=sanitize_logs)
        
        # Create handlers
        handlers = []
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        console_handler.setFormatter(console_formatter)
        handlers.append(console_handler)
        
        # File handler if specified
        if log_file:
            try:
                # Ensure log directory exists
                log_file.parent.mkdir(parents=True, exist_ok=True)
                
                file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
                file_handler.setLevel(logging.DEBUG)  # Always capture everything to file
                file_handler.setFormatter(file_formatter)
                handlers.append(file_handler)
                
                # Add rotation if file gets too large
                from logging.handlers import RotatingFileHandler
                rotating_handler = RotatingFileHandler(
                    log_file, maxBytes=10*1024*1024, backupCount=5
                )
                rotating_handler.setLevel(logging.DEBUG)
                rotating_handler.setFormatter(file_formatter)
                handlers.append(rotating_handler)
                
            except Exception as e:
                # Fall back to console only if file logging fails
                print(f"Warning: Could not set up file logging: {e}")
        
        # Configure root logger
        logging.basicConfig(
            level=level,
            handlers=handlers,
            force=True  # Override any existing configuration
        )
        
        # Configure specific loggers with appropriate levels
        logging.getLogger("pysam").setLevel(logging.WARNING)
        logging.getLogger("urllib3").setLevel(logging.WARNING)
        logging.getLogger("requests").setLevel(logging.WARNING)
        
        # Set up application loggers
        app_logger = logging.getLogger("chimeric_detector")
        app_logger.setLevel(level)
        
        # Log startup information
        startup_logger = logging.getLogger("chimeric_detector.startup")
        startup_logger.info(f"Logging initialized - Level: {logging.getLevelName(level)}")
        if log_file:
            startup_logger.info(f"Log file: {log_file}")
        if sanitize_logs:
            startup_logger.debug("Log sanitization enabled")
    
    @staticmethod
    def log_performance(logger: logging.Logger, operation: str, duration: float, 
                       count: Optional[int] = None, **kwargs) -> None:
        """Log performance metrics in a structured way."""
        perf_data = {
            "operation": operation,
            "duration_seconds": round(duration, 3),
            "performance": True
        }
        
        if count is not None:
            perf_data["count"] = count
            perf_data["rate_per_second"] = round(count / duration, 2) if duration > 0 else 0
        
        perf_data.update(kwargs)
        
        # Create log record with extra data
        logger.info(f"Performance: {operation} completed in {duration:.3f}s", extra=perf_data)
    
    @staticmethod
    def log_security_event(logger: logging.Logger, event_type: str, details: str, 
                          severity: str = "WARNING") -> None:
        """Log security-related events in a structured way."""
        security_data = {
            "security_event": True,
            "event_type": event_type,
            "severity": severity,
            "details": details
        }
        
        log_level = getattr(logging, severity.upper(), logging.WARNING)
        logger.log(log_level, f"Security event: {event_type} - {details}", extra=security_data)
    
    @staticmethod
    def log_error_with_context(logger: logging.Logger, error: Exception, context: Dict[str, Any]) -> None:
        """Log errors with additional context information."""
        error_data = {
            "error": True,
            "error_type": type(error).__name__,
            "error_message": str(error),
            "context": context
        }
        
        logger.error(f"Error in {context.get('operation', 'unknown')}: {error}", extra=error_data)