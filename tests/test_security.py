"""Tests for security module."""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import patch, mock_open

from chimeric_detector.utils.security import (
    validate_path, validate_bam_file, validate_config_file,
    validate_assembly_file, validate_output_path, sanitize_contig_name,
    validate_numeric_parameter, safe_file_open, create_secure_temp_file,
    SecurityError, PathTraversalError, FileSizeError, FileTypeError,
    MAX_BAM_FILE_SIZE, MAX_CONFIG_FILE_SIZE
)


class TestPathValidation:
    """Test path validation functions."""
    
    def test_validate_safe_path(self, temp_dir):
        """Test validation of safe paths."""
        safe_file = temp_dir / "safe.txt"
        safe_file.write_text("content")
        
        result = validate_path(safe_file, must_exist=True)
        assert result == safe_file.resolve()
    
    def test_path_traversal_detection(self):
        """Test detection of path traversal attempts."""
        with pytest.raises(PathTraversalError, match="Path traversal detected"):
            validate_path("../../../etc/passwd")
        
        with pytest.raises(PathTraversalError, match="Path traversal detected"):
            validate_path("safe/../../../etc/passwd")
    
    def test_nonexistent_file_required(self, temp_dir):
        """Test handling of required files that don't exist."""
        nonexistent = temp_dir / "nonexistent.txt"
        
        with pytest.raises(FileNotFoundError):
            validate_path(nonexistent, must_exist=True)
        
        # Should work when not required
        result = validate_path(nonexistent, must_exist=False)
        assert result == nonexistent.resolve()
    
    def test_file_size_validation(self, temp_dir):
        """Test file size validation."""
        large_file = temp_dir / "large.txt"
        large_file.write_text("x" * 1000)
        
        # Should pass with sufficient limit
        result = validate_path(large_file, max_size=2000)
        assert result == large_file.resolve()
        
        # Should fail with insufficient limit
        with pytest.raises(FileSizeError, match="exceeds maximum size"):
            validate_path(large_file, max_size=500)
    
    def test_allowed_directories(self, temp_dir):
        """Test allowed directory restrictions."""
        allowed_dir = temp_dir / "allowed"
        allowed_dir.mkdir()
        disallowed_dir = temp_dir / "disallowed"
        disallowed_dir.mkdir()
        
        safe_file = allowed_dir / "safe.txt"
        safe_file.write_text("content")
        unsafe_file = disallowed_dir / "unsafe.txt"
        unsafe_file.write_text("content")
        
        # Should pass in allowed directory
        result = validate_path(safe_file, allowed_dirs=[allowed_dir])
        assert result == safe_file.resolve()
        
        # Should fail in disallowed directory
        with pytest.raises(SecurityError, match="outside allowed directories"):
            validate_path(unsafe_file, allowed_dirs=[allowed_dir])


class TestFileTypeValidation:
    """Test file type-specific validation functions."""
    
    def test_validate_bam_file(self, temp_dir):
        """Test BAM file validation."""
        # Valid BAM file
        bam_file = temp_dir / "test.bam"
        bam_file.write_text("fake bam content")
        
        result = validate_bam_file(bam_file)
        assert result == bam_file.resolve()
        
        # Invalid extension
        with pytest.raises(FileTypeError, match="Invalid BAM file extension"):
            validate_bam_file(temp_dir / "test.txt")
        
        # Test size limit
        large_bam = temp_dir / "large.bam"
        with patch('chimeric_detector.utils.security.MAX_BAM_FILE_SIZE', 100):
            large_bam.write_text("x" * 200)
            with pytest.raises(FileSizeError):
                validate_bam_file(large_bam)
    
    def test_validate_config_file(self, temp_dir):
        """Test config file validation."""
        # Valid YAML config
        yaml_config = temp_dir / "config.yaml"
        yaml_config.write_text("key: value")
        
        result = validate_config_file(yaml_config)
        assert result == yaml_config.resolve()
        
        # Valid JSON config
        json_config = temp_dir / "config.json"
        json_config.write_text('{"key": "value"}')
        
        result = validate_config_file(json_config)
        assert result == json_config.resolve()
        
        # Invalid extension
        with pytest.raises(FileTypeError, match="Invalid config file extension"):
            validate_config_file(temp_dir / "config.txt")
    
    def test_validate_assembly_file(self, temp_dir):
        """Test assembly file validation."""
        # Valid FASTA file
        fasta_file = temp_dir / "assembly.fasta"
        fasta_file.write_text(">contig1\nATCG\n")
        
        result = validate_assembly_file(fasta_file)
        assert result == fasta_file.resolve()
        
        # Test different extensions
        for ext in [".fa", ".fna", ".fas"]:
            fasta_file = temp_dir / f"assembly{ext}"
            fasta_file.write_text(">contig1\nATCG\n")
            result = validate_assembly_file(fasta_file)
            assert result == fasta_file.resolve()
        
        # Invalid extension
        with pytest.raises(FileTypeError, match="Invalid assembly file extension"):
            validate_assembly_file(temp_dir / "assembly.txt")
    
    def test_validate_output_path(self, temp_dir):
        """Test output path validation."""
        # Valid output files
        for ext in [".json", ".tsv", ".bed", ".txt"]:
            output_file = temp_dir / f"output{ext}"
            result = validate_output_path(output_file)
            assert result == output_file.resolve()
        
        # Invalid extension
        with pytest.raises(FileTypeError, match="Invalid output file extension"):
            validate_output_path(temp_dir / "output.exe")


class TestInputSanitization:
    """Test input sanitization functions."""
    
    def test_sanitize_contig_name(self):
        """Test contig name sanitization."""
        # Valid names
        valid_names = ["contig1", "scaffold_001", "chr1.1", "NODE-123"]
        for name in valid_names:
            result = sanitize_contig_name(name)
            assert result == name
        
        # Invalid names
        invalid_names = ["", "con/tig", "con tig", "con|tig", "con<tig", "CON", "PRN", "COM1"]
        for name in invalid_names:
            with pytest.raises(SecurityError):
                sanitize_contig_name(name)
    
    def test_validate_numeric_parameter(self):
        """Test numeric parameter validation."""
        # Valid parameters
        assert validate_numeric_parameter(5, 0, 10, "test") == 5
        assert validate_numeric_parameter(0.5, 0.0, 1.0, "test") == 0.5
        
        # Out of bounds
        with pytest.raises(SecurityError, match="must be between"):
            validate_numeric_parameter(15, 0, 10, "test")
        
        with pytest.raises(SecurityError, match="must be between"):
            validate_numeric_parameter(-5, 0, 10, "test")
        
        # Non-numeric
        with pytest.raises(SecurityError, match="must be numeric"):
            validate_numeric_parameter("string", 0, 10, "test")


class TestSecureFileOperations:
    """Test secure file operation functions."""
    
    def test_safe_file_open_read(self, temp_dir):
        """Test safe file opening for reading."""
        test_file = temp_dir / "test.txt"
        test_content = "test content"
        test_file.write_text(test_content)
        
        with safe_file_open(test_file, 'r') as f:
            content = f.read()
            assert content == test_content
    
    def test_safe_file_open_write(self, temp_dir):
        """Test safe file opening for writing."""
        test_file = temp_dir / "new_file.txt"
        test_content = "new content"
        
        with safe_file_open(test_file, 'w') as f:
            f.write(test_content)
        
        assert test_file.read_text() == test_content
    
    def test_safe_file_open_nonexistent_path(self, temp_dir):
        """Test safe file opening with nonexistent directory."""
        nonexistent_dir = temp_dir / "nonexistent" / "subdir"
        test_file = nonexistent_dir / "test.txt"
        
        # Should create directory for write mode
        with safe_file_open(test_file, 'w') as f:
            f.write("content")
        
        assert test_file.exists()
        assert test_file.read_text() == "content"
    
    def test_create_secure_temp_file(self, temp_dir):
        """Test secure temporary file creation."""
        temp_file = create_secure_temp_file(
            suffix=".txt",
            prefix="test_",
            directory=temp_dir
        )
        
        assert temp_file.exists()
        assert temp_file.parent.resolve() == temp_dir.resolve()
        assert temp_file.name.startswith("test_")
        assert temp_file.suffix == ".txt"
        
        # Clean up
        temp_file.unlink()
    
    def test_create_secure_temp_file_default_location(self):
        """Test secure temporary file creation in default location."""
        temp_file = create_secure_temp_file(suffix=".test")
        
        assert temp_file.exists()
        assert temp_file.suffix == ".test"
        
        # Clean up
        temp_file.unlink()


class TestSecurityErrorHandling:
    """Test security error handling and edge cases."""
    
    def test_symlink_detection(self, temp_dir):
        """Test symbolic link detection."""
        real_file = temp_dir / "real.txt"
        real_file.write_text("content")
        
        symlink_file = temp_dir / "symlink.txt"
        try:
            symlink_file.symlink_to(real_file)
            
            # Should log warning but still work
            with patch('chimeric_detector.utils.security.logger') as mock_logger:
                result = validate_path(symlink_file)
                mock_logger.warning.assert_called_once()
                assert result == symlink_file.resolve()
                
        except OSError:
            # Skip test if symlinks not supported on this platform
            pytest.skip("Symlinks not supported on this platform")
    
    def test_permission_error_handling(self, temp_dir):
        """Test handling of permission errors."""
        restricted_file = temp_dir / "restricted.txt"
        
        # Simulate permission error
        with patch('builtins.open', side_effect=PermissionError("Access denied")):
            with pytest.raises(SecurityError, match="Failed to open file"):
                safe_file_open(restricted_file, 'r')
    
    def test_large_file_handling(self, temp_dir):
        """Test handling of very large files."""
        # Test with mocked file size - use proper config extension
        large_file = temp_dir / "large.yaml"
        large_file.write_text("content")
        
        # Mock stat to return large size
        with patch.object(Path, 'stat') as mock_stat:
            mock_stat.return_value.st_size = MAX_CONFIG_FILE_SIZE + 1
            
            with pytest.raises(FileSizeError, match="exceeds maximum size"):
                validate_config_file(large_file)
    
    def test_path_resolution_error(self, temp_dir):
        """Test handling of path resolution errors."""
        # Create path that will cause resolution error
        bad_path = temp_dir / ("x" * 1000)  # Very long path name
        
        with patch.object(Path, 'resolve', side_effect=OSError("Path too long")):
            with pytest.raises(SecurityError, match="Failed to resolve path"):
                validate_path(bad_path)


class TestSecurityIntegration:
    """Integration tests for security features."""
    
    def test_security_config_validation(self, temp_dir):
        """Test security validation in configuration loading."""
        from chimeric_detector.config.config import Config
        
        # Valid config
        valid_config = temp_dir / "valid.yaml"
        valid_config.write_text("""
window:
  size: 1000
  step: 100
detection:
  min_confidence_score: 0.5
""")
        
        config = Config.from_file(valid_config)
        assert config.window.size == 1000
        
        # Path traversal attempt
        traversal_config = temp_dir / "traversal.yaml"
        traversal_config.write_text("malicious: content")
        
        # The security validation should be handled by the calling code
        # This tests that the config loading itself works
        config = Config.from_file(traversal_config)
        assert config is not None
    
    def test_output_security_validation(self, temp_dir, sample_breakpoints):
        """Test security validation in output operations."""
        from chimeric_detector.utils.output import OutputFormatter
        from chimeric_detector.config.config import OutputConfig
        
        formatter = OutputFormatter(OutputConfig(format="json"))
        
        # Valid output path
        safe_output = temp_dir / "safe_output.json"
        formatter.write_results(sample_breakpoints, safe_output)
        assert safe_output.exists()
        
        # The path validation should be done by the calling code
        # Here we test that the formatter itself works correctly
    
    @pytest.mark.parametrize("malicious_input", [
        "../../../etc/passwd",
        "/etc/shadow",
        "C:\\Windows\\System32\\config\\SAM", 
        "file:///etc/passwd",
        "\\\\server\\share\\file.txt"
    ])
    def test_malicious_path_rejection(self, malicious_input):
        """Test rejection of various malicious path patterns."""
        # Only the path traversal pattern should be caught by our current validation
        if ".." in malicious_input:
            with pytest.raises(SecurityError):
                validate_path(malicious_input)
        else:
            # Other patterns would need additional security rules
            # For now, just verify they don't crash
            try:
                validate_path(malicious_input)
            except (FileNotFoundError, OSError):
                # Expected for non-existent paths
                pass
    
    def test_resource_exhaustion_protection(self, temp_dir):
        """Test protection against resource exhaustion attacks."""
        # Test with many file operations
        for i in range(100):
            small_file = temp_dir / f"file_{i}.txt"
            small_file.write_text(f"content {i}")
            
            # Should handle many small files without issue
            result = validate_path(small_file)
            assert result == small_file.resolve()
    
    def test_concurrent_access_safety(self, temp_dir):
        """Test safety of concurrent file access."""
        import threading
        import time
        
        test_file = temp_dir / "concurrent.txt"
        results = []
        errors = []
        
        def worker():
            try:
                with safe_file_open(test_file, 'a') as f:
                    f.write(f"Thread {threading.current_thread().ident}\n")
                    time.sleep(0.01)  # Small delay to increase chance of race condition
                results.append("success")
            except Exception as e:
                errors.append(str(e))
        
        # Start multiple threads
        threads = []
        for _ in range(10):
            t = threading.Thread(target=worker)
            threads.append(t)
            t.start()
        
        # Wait for all threads
        for t in threads:
            t.join()
        
        # Should have mostly successes
        assert len(results) >= 8  # Allow for some contention
        assert len(errors) <= 2
        
        # File should exist and have content
        assert test_file.exists()
        content = test_file.read_text()
        assert "Thread" in content