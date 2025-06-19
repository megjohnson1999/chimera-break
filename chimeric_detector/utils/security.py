"""Security utilities for input validation and safe file operations."""

import os
import re
from pathlib import Path
from typing import List, Optional, Union
import logging

logger = logging.getLogger(__name__)

# Maximum file sizes (in bytes)
MAX_BAM_FILE_SIZE = 50 * 1024 * 1024 * 1024  # 50GB
MAX_CONFIG_FILE_SIZE = 10 * 1024 * 1024  # 10MB
MAX_ASSEMBLY_FILE_SIZE = 10 * 1024 * 1024 * 1024  # 10GB

# Allowed file extensions
ALLOWED_BAM_EXTENSIONS = {'.bam', '.sam'}
ALLOWED_CONFIG_EXTENSIONS = {'.yaml', '.yml', '.json'}
ALLOWED_ASSEMBLY_EXTENSIONS = {'.fasta', '.fa', '.fna', '.fas'}
ALLOWED_OUTPUT_EXTENSIONS = {'.json', '.tsv', '.bed', '.txt'}

class SecurityError(Exception):
    """Raised when a security violation is detected."""
    pass

class PathTraversalError(SecurityError):
    """Raised when path traversal is detected."""
    pass

class FileSizeError(SecurityError):
    """Raised when file size exceeds limits."""
    pass

class FileTypeError(SecurityError):
    """Raised when file type is not allowed."""
    pass

def validate_path(path: Union[str, Path], 
                 allowed_dirs: Optional[List[Path]] = None,
                 must_exist: bool = False,
                 max_size: Optional[int] = None) -> Path:
    """
    Validate that a file path is safe to use.
    
    Args:
        path: The path to validate
        allowed_dirs: List of directories the path must be within
        must_exist: Whether the path must exist
        max_size: Maximum file size in bytes
        
    Returns:
        Resolved Path object
        
    Raises:
        PathTraversalError: If path traversal is detected
        FileNotFoundError: If must_exist=True and path doesn't exist
        FileSizeError: If file exceeds size limit
    """
    if isinstance(path, str):
        path = Path(path)
    
    # Check for path traversal attempts
    if ".." in path.parts:
        raise PathTraversalError(f"Path traversal detected in: {path}")
    
    # Resolve the path to handle symlinks and relative paths
    try:
        resolved_path = path.resolve()
    except (OSError, RuntimeError) as e:
        raise SecurityError(f"Failed to resolve path {path}: {e}")
    
    # Check if path is within allowed directories
    if allowed_dirs:
        allowed_dirs = [Path(d).resolve() for d in allowed_dirs]
        if not any(resolved_path.is_relative_to(allowed) for allowed in allowed_dirs):
            raise SecurityError(f"Path {resolved_path} is outside allowed directories")
    
    # Check if file exists when required
    if must_exist and not resolved_path.exists():
        raise FileNotFoundError(f"Required file does not exist: {resolved_path}")
    
    # Check file size if it exists and max_size is specified
    if max_size and resolved_path.exists() and resolved_path.is_file():
        file_size = resolved_path.stat().st_size
        if file_size > max_size:
            raise FileSizeError(f"File {resolved_path} ({file_size} bytes) exceeds maximum size ({max_size} bytes)")
    
    # Check for symbolic links (optional security measure)
    if path.is_symlink():
        logger.warning(f"Symbolic link detected: {path}")
    
    return resolved_path

def validate_bam_file(path: Union[str, Path]) -> Path:
    """Validate a BAM file path."""
    path = Path(path)
    
    # Check file extension
    if path.suffix.lower() not in ALLOWED_BAM_EXTENSIONS:
        raise FileTypeError(f"Invalid BAM file extension: {path.suffix}")
    
    # Validate path with BAM-specific constraints
    validated_path = validate_path(
        path,
        must_exist=True,
        max_size=MAX_BAM_FILE_SIZE
    )
    
    # Check for BAM index file
    index_paths = [
        validated_path.with_suffix('.bam.bai'),
        validated_path.with_suffix('.bai'),
        validated_path.parent / f"{validated_path.name}.bai"
    ]
    
    has_index = any(idx.exists() for idx in index_paths)
    if not has_index:
        logger.warning(f"No index file found for BAM file: {validated_path}")
    
    return validated_path

def validate_config_file(path: Union[str, Path]) -> Path:
    """Validate a configuration file path."""
    path = Path(path)
    
    # Check file extension
    if path.suffix.lower() not in ALLOWED_CONFIG_EXTENSIONS:
        raise FileTypeError(f"Invalid config file extension: {path.suffix}")
    
    return validate_path(
        path,
        must_exist=True,
        max_size=MAX_CONFIG_FILE_SIZE
    )

def validate_assembly_file(path: Union[str, Path]) -> Path:
    """Validate an assembly FASTA file path."""
    path = Path(path)
    
    # Check file extension
    if path.suffix.lower() not in ALLOWED_ASSEMBLY_EXTENSIONS:
        raise FileTypeError(f"Invalid assembly file extension: {path.suffix}")
    
    return validate_path(
        path,
        must_exist=True,
        max_size=MAX_ASSEMBLY_FILE_SIZE
    )

def validate_output_path(path: Union[str, Path], 
                        allowed_dirs: Optional[List[Path]] = None) -> Path:
    """Validate an output file path."""
    path = Path(path)
    
    # Check file extension
    if path.suffix.lower() not in ALLOWED_OUTPUT_EXTENSIONS:
        raise FileTypeError(f"Invalid output file extension: {path.suffix}")
    
    return validate_path(
        path,
        allowed_dirs=allowed_dirs,
        must_exist=False
    )

def sanitize_contig_name(name: str) -> str:
    """
    Sanitize a contig name to ensure it contains only safe characters.
    
    Args:
        name: The contig name to sanitize
        
    Returns:
        Sanitized contig name
        
    Raises:
        SecurityError: If name contains unsafe characters
    """
    if not name:
        raise SecurityError("Contig name cannot be empty")
    
    # Allow alphanumeric, underscore, dash, and dot
    if not re.match(r'^[a-zA-Z0-9_.-]+$', name):
        raise SecurityError(f"Contig name contains unsafe characters: {name}")
    
    # Prevent names that could be confused with file system elements
    if name in {'.', '..', 'CON', 'PRN', 'AUX', 'NUL'} or name.startswith('COM') or name.startswith('LPT'):
        raise SecurityError(f"Contig name is reserved: {name}")
    
    return name

def validate_numeric_parameter(value: Union[int, float], 
                             min_val: Union[int, float], 
                             max_val: Union[int, float], 
                             param_name: str) -> Union[int, float]:
    """
    Validate that a numeric parameter is within acceptable bounds.
    
    Args:
        value: The value to validate
        min_val: Minimum allowed value
        max_val: Maximum allowed value
        param_name: Name of the parameter for error messages
        
    Returns:
        The validated value
        
    Raises:
        SecurityError: If value is outside bounds
    """
    if not isinstance(value, (int, float)):
        raise SecurityError(f"{param_name} must be numeric, got {type(value)}")
    
    if not (min_val <= value <= max_val):
        raise SecurityError(f"{param_name} must be between {min_val} and {max_val}, got {value}")
    
    return value

def safe_file_open(path: Path, mode: str = 'r', **kwargs):
    """
    Safely open a file with additional security checks.
    
    Args:
        path: Path to the file
        mode: File open mode
        **kwargs: Additional arguments for open()
        
    Returns:
        File handle
        
    Raises:
        SecurityError: If security checks fail
    """
    # Additional check for write modes
    if any(m in mode for m in ['w', 'a', 'x']):
        # Ensure parent directory exists
        path.parent.mkdir(parents=True, exist_ok=True)
        
        # Check if we're about to overwrite an existing file
        if 'w' in mode and path.exists():
            logger.warning(f"Overwriting existing file: {path}")
    
    try:
        return open(path, mode, **kwargs)
    except (OSError, IOError) as e:
        raise SecurityError(f"Failed to open file {path}: {e}")

def create_secure_temp_file(suffix: str = '', prefix: str = 'chimera_', 
                           directory: Optional[Path] = None) -> Path:
    """
    Create a secure temporary file.
    
    Args:
        suffix: File suffix
        prefix: File prefix
        directory: Directory to create file in
        
    Returns:
        Path to the temporary file
    """
    import tempfile
    
    if directory:
        directory = validate_path(directory)
    
    fd, temp_path = tempfile.mkstemp(
        suffix=suffix,
        prefix=prefix,
        dir=directory
    )
    
    # Close the file descriptor and return the path
    os.close(fd)
    
    return Path(temp_path)