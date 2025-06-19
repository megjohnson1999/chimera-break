"""Tests for BAM parser with error handling."""

import pytest
import numpy as np
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path

from chimeric_detector.core.bam_parser import BAMParser, ReadPairInfo
from chimeric_detector.config.config import QualityConfig


class TestBAMParserInitialization:
    """Test BAM parser initialization and context management."""
    
    def test_bam_parser_creation(self, sample_config):
        """Test BAM parser creation."""
        parser = BAMParser("test.bam", sample_config.quality)
        assert parser.bam_path == "test.bam"
        assert parser.config == sample_config.quality
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_context_manager_success(self, mock_alignment_file, sample_config):
        """Test successful context manager operation."""
        mock_bam = MagicMock()
        mock_bam.references = ["contig1", "contig2"]
        mock_bam.check_index.return_value = True
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser as p:
            assert p._bam == mock_bam
            assert p == parser
        
        mock_bam.close.assert_called_once()
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_context_manager_file_not_found(self, mock_alignment_file, sample_config):
        """Test context manager with file not found error."""
        mock_alignment_file.side_effect = FileNotFoundError("File not found")
        
        parser = BAMParser("nonexistent.bam", sample_config.quality)
        
        with pytest.raises(FileNotFoundError, match="BAM file not found"):
            with parser:
                pass
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_context_manager_permission_error(self, mock_alignment_file, sample_config):
        """Test context manager with permission error."""
        mock_alignment_file.side_effect = PermissionError("Access denied")
        
        parser = BAMParser("restricted.bam", sample_config.quality)
        
        with pytest.raises(PermissionError, match="Permission denied reading BAM file"):
            with parser:
                pass
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_context_manager_corrupted_bam(self, mock_alignment_file, sample_config):
        """Test context manager with corrupted BAM file."""
        import pysam
        mock_alignment_file.side_effect = pysam.utils.SamtoolsError("Invalid BAM")
        
        parser = BAMParser("corrupted.bam", sample_config.quality)
        
        with pytest.raises(ValueError, match="Invalid or corrupted BAM file"):
            with parser:
                pass
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_context_manager_no_references(self, mock_alignment_file, sample_config):
        """Test context manager with BAM file containing no references."""
        mock_bam = MagicMock()
        mock_bam.references = []
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("empty.bam", sample_config.quality)
        
        with pytest.raises(RuntimeError, match="Unexpected error opening BAM file"):
            with parser:
                pass
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_context_manager_close_error(self, mock_alignment_file, sample_config):
        """Test context manager with error during close."""
        mock_bam = MagicMock()
        mock_bam.references = ["contig1"]
        mock_bam.check_index.return_value = True
        mock_bam.close.side_effect = Exception("Close error")
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        # Should not raise exception even if close fails
        with patch('chimeric_detector.core.bam_parser.logging.getLogger') as mock_logger:
            with parser:
                pass
            # Should log warning about close error
            mock_logger.return_value.warning.assert_called()


class TestBAMParserMethods:
    """Test BAM parser methods with error handling."""
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_get_contigs_success(self, mock_alignment_file, sample_config):
        """Test successful contig retrieval."""
        mock_bam = MagicMock()
        mock_bam.references = ["contig1", "contig2", "contig3"]
        mock_bam.check_index.return_value = True
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            contigs = parser.get_contigs()
            assert contigs == ["contig1", "contig2", "contig3"]
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_get_contigs_not_opened(self, mock_alignment_file, sample_config):
        """Test get_contigs when BAM file not opened."""
        parser = BAMParser("test.bam", sample_config.quality)
        
        with pytest.raises(RuntimeError, match="BAM file not opened"):
            parser.get_contigs()
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_get_contigs_no_contigs(self, mock_alignment_file, sample_config):
        """Test get_contigs with empty reference list."""
        mock_bam = MagicMock()
        mock_bam.references = []
        mock_bam.check_index.return_value = True
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with pytest.raises(RuntimeError, match="Unexpected error opening BAM file"):
            with parser:
                pass


class TestInsertSizeEstimation:
    """Test insert size estimation with error handling."""
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_estimate_insert_size_success(self, mock_alignment_file, sample_config, mock_bam_file):
        """Test successful insert size estimation."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        
        # Create mock reads with proper pairs
        mock_reads = []
        for i in range(200):
            read = MagicMock()
            read.is_paired = True
            read.is_proper_pair = True
            read.mapping_quality = 30
            read.template_length = 300 + i  # Varying insert sizes
            read.is_secondary = False
            read.is_supplementary = False
            mock_reads.append(read)
        
        mock_bam.fetch.return_value = mock_reads
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            stats = parser.estimate_insert_size_distribution("test_contig")
            
            assert "median" in stats
            assert "mad" in stats
            assert "mean" in stats
            assert "std" in stats
            assert stats["median"] > 0
            assert stats["mad"] >= 0
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_estimate_insert_size_contig_not_found(self, mock_alignment_file, sample_config):
        """Test insert size estimation with nonexistent contig."""
        mock_bam = MagicMock()
        mock_bam.references = ["contig1", "contig2"]
        mock_bam.check_index.return_value = True
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            stats = parser.estimate_insert_size_distribution("nonexistent_contig")
            
            # Should return default values
            assert stats["median"] == 500
            assert stats["mad"] == 100
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_estimate_insert_size_few_pairs(self, mock_alignment_file, sample_config):
        """Test insert size estimation with very few proper pairs."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        
        # Only a few reads
        mock_reads = []
        for i in range(10):
            read = MagicMock()
            read.is_paired = True
            read.is_proper_pair = True
            read.mapping_quality = 30
            read.template_length = 300
            read.is_secondary = False
            read.is_supplementary = False
            mock_reads.append(read)
        
        mock_bam.fetch.return_value = mock_reads
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            with patch.object(parser.logger, 'warning') as mock_warning:
                stats = parser.estimate_insert_size_distribution("test_contig")
                
                # Should warn about few pairs and return defaults
                mock_warning.assert_called()
                assert stats["median"] == 500
                assert stats["mad"] == 100
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_estimate_insert_size_fetch_error(self, mock_alignment_file, sample_config):
        """Test insert size estimation with fetch error."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        mock_bam.fetch.side_effect = OSError("Read error")
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            with patch.object(parser.logger, 'error') as mock_error:
                stats = parser.estimate_insert_size_distribution("test_contig")
                
                # Should log error and return defaults
                mock_error.assert_called()
                assert stats["median"] == 500
                assert stats["mad"] == 100
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_estimate_insert_size_with_invalid_template_lengths(self, mock_alignment_file, sample_config):
        """Test insert size estimation with invalid template lengths."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        
        # Mix of valid and invalid template lengths
        mock_reads = []
        valid_count = 0
        for i in range(200):
            read = MagicMock()
            read.is_paired = True
            read.is_proper_pair = True
            read.mapping_quality = 30
            read.is_secondary = False
            read.is_supplementary = False
            
            if i % 5 == 0:  # Some invalid values
                read.template_length = None
            elif i % 7 == 0:
                read.template_length = 0
            elif i % 11 == 0:
                read.template_length = 100000  # Too large
            else:
                read.template_length = 300 + i
                valid_count += 1
            
            mock_reads.append(read)
        
        mock_bam.fetch.return_value = mock_reads
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            stats = parser.estimate_insert_size_distribution("test_contig")
            
            # Should still work with valid reads
            if valid_count >= 100:
                assert stats["median"] > 250
                assert stats["median"] < 500
            else:
                # Fallback to defaults
                assert stats["median"] == 500


class TestWindowIteration:
    """Test window iteration with error handling."""
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_iterate_windows_success(self, mock_alignment_file, sample_config):
        """Test successful window iteration."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        mock_bam.get_reference_length.return_value = 5000
        mock_bam.fetch.return_value = []  # No reads for simplicity
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            windows = list(parser.iterate_windows("test_contig", 1000, 500))
            
            # Should have multiple windows
            assert len(windows) > 1
            
            # Check first window
            start, end, pairs = windows[0]
            assert start == 0
            assert end == 1000
            assert isinstance(pairs, list)
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_iterate_windows_invalid_parameters(self, mock_alignment_file, sample_config):
        """Test window iteration with invalid parameters."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            # Invalid window size
            windows = list(parser.iterate_windows("test_contig", 0, 100))
            assert len(windows) == 0
            
            # Invalid step size
            windows = list(parser.iterate_windows("test_contig", 1000, -1))
            assert len(windows) == 0
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_iterate_windows_contig_not_found(self, mock_alignment_file, sample_config):
        """Test window iteration with nonexistent contig."""
        mock_bam = MagicMock()
        mock_bam.references = ["contig1", "contig2"]
        mock_bam.check_index.return_value = True
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            windows = list(parser.iterate_windows("nonexistent", 1000, 500))
            assert len(windows) == 0
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_iterate_windows_window_larger_than_contig(self, mock_alignment_file, sample_config):
        """Test window iteration when window is larger than contig."""
        mock_bam = MagicMock()
        mock_bam.references = ["small_contig"]
        mock_bam.check_index.return_value = True
        mock_bam.get_reference_length.return_value = 500  # Small contig
        mock_bam.fetch.return_value = []
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            with patch('chimeric_detector.core.bam_parser.logging.getLogger') as mock_logger:
                windows = list(parser.iterate_windows("small_contig", 1000, 500))
                
                # Should still return one window covering entire contig
                assert len(windows) == 1
                start, end, pairs = windows[0]
                assert start == 0
                assert end == 500
                
                # Should log warning
                mock_logger.return_value.warning.assert_called()
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_iterate_windows_zero_length_contig(self, mock_alignment_file, sample_config):
        """Test window iteration with zero-length contig."""
        mock_bam = MagicMock()
        mock_bam.references = ["empty_contig"]
        mock_bam.check_index.return_value = True
        mock_bam.get_reference_length.return_value = 0
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            with patch('chimeric_detector.core.bam_parser.logging.getLogger') as mock_logger:
                windows = list(parser.iterate_windows("empty_contig", 1000, 500))
                
                # Should return no windows
                assert len(windows) == 0
                
                # Should log warning
                mock_logger.return_value.warning.assert_called()


class TestReadPairExtraction:
    """Test read pair extraction with error handling."""
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_get_read_pairs_success(self, mock_alignment_file, sample_config):
        """Test successful read pair extraction."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        mock_bam.get_reference_length.return_value = 10000
        
        # Create mock reads
        mock_reads = []
        for i in range(10):
            read = MagicMock()
            read.is_paired = True
            read.is_secondary = False
            read.is_supplementary = False
            read.mapping_quality = 30
            read.reference_name = "test_contig"
            read.reference_start = i * 100
            read.reference_end = read.reference_start + 100
            read.next_reference_name = "test_contig"
            read.next_reference_start = read.reference_start + 200
            read.template_length = 300
            read.is_proper_pair = True
            read.is_reverse = False
            read.mate_is_unmapped = False
            read.query_sequence = "A" * 100
            read.query_length = 100
            mock_reads.append(read)
        
        mock_bam.fetch.return_value = mock_reads
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            pairs = parser.get_read_pairs_in_window("test_contig", 0, 1000)
            
            assert len(pairs) > 0
            for pair in pairs:
                assert isinstance(pair, ReadPairInfo)
                assert pair.contig == "test_contig"
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_get_read_pairs_invalid_coordinates(self, mock_alignment_file, sample_config):
        """Test read pair extraction with invalid coordinates."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        mock_bam.get_reference_length.return_value = 10000
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            # Start >= end
            pairs = parser.get_read_pairs_in_window("test_contig", 1000, 500)
            assert len(pairs) == 0
            
            # Start beyond contig
            pairs = parser.get_read_pairs_in_window("test_contig", 15000, 16000)
            assert len(pairs) == 0
            
            # Negative start (should be corrected)
            pairs = parser.get_read_pairs_in_window("test_contig", -100, 1000)
            # Should not crash and should work with corrected coordinates
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_get_read_pairs_fetch_error(self, mock_alignment_file, sample_config):
        """Test read pair extraction with fetch error."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        mock_bam.get_reference_length.return_value = 10000
        mock_bam.fetch.side_effect = OSError("Fetch error")
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            with patch('chimeric_detector.core.bam_parser.logging.getLogger') as mock_logger:
                pairs = parser.get_read_pairs_in_window("test_contig", 0, 1000)
                
                # Should return empty list and log error
                assert len(pairs) == 0
                mock_logger.return_value.error.assert_called()


class TestReadPairInfoExtraction:
    """Test ReadPairInfo extraction with edge cases."""
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_extract_pair_info_success(self, mock_alignment_file, sample_config):
        """Test successful ReadPairInfo extraction."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        mock_bam.get_reference_length.return_value = 10000
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        # Create mock read
        read = MagicMock()
        read.is_paired = True
        read.reference_name = "test_contig"
        read.reference_start = 100
        read.reference_end = 200
        read.next_reference_name = "test_contig"
        read.next_reference_start = 300
        read.template_length = 300
        read.is_proper_pair = True
        read.is_reverse = False
        read.mate_is_unmapped = False
        read.query_sequence = "A" * 100
        read.query_length = 100
        
        with parser:
            pair_info = parser._extract_pair_info(read)
            
            assert pair_info is not None
            assert isinstance(pair_info, ReadPairInfo)
            assert pair_info.position == 100
            assert pair_info.mate_position == 300
            assert pair_info.is_proper_pair
            assert not pair_info.is_discordant
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_extract_pair_info_unpaired_read(self, mock_alignment_file, sample_config):
        """Test ReadPairInfo extraction with unpaired read."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        # Create unpaired read
        read = MagicMock()
        read.is_paired = False
        
        with parser:
            pair_info = parser._extract_pair_info(read)
            assert pair_info is None
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_extract_pair_info_missing_query_sequence(self, mock_alignment_file, sample_config):
        """Test ReadPairInfo extraction with missing query sequence."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        mock_bam.get_reference_length.return_value = 10000
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        # Create read with missing query sequence
        read = MagicMock()
        read.is_paired = True
        read.reference_name = "test_contig"
        read.reference_start = 100
        read.reference_end = None  # Missing
        read.next_reference_name = "test_contig"
        read.next_reference_start = 300
        read.template_length = 0  # Will trigger manual calculation
        read.is_proper_pair = True
        read.is_reverse = False
        read.mate_is_unmapped = False
        read.query_sequence = None  # Missing
        read.query_length = 100
        
        with parser:
            pair_info = parser._extract_pair_info(read)
            
            # Should still work with fallback logic
            assert pair_info is not None
            assert pair_info.position == 100


class TestErrorRecovery:
    """Test error recovery and graceful degradation."""
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_partial_failure_recovery(self, mock_alignment_file, sample_config):
        """Test recovery from partial failures during processing."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        mock_bam.get_reference_length.return_value = 10000
        
        # Mix of good and bad reads
        mock_reads = []
        for i in range(10):
            read = MagicMock()
            if i % 3 == 0:  # Every third read will cause an error
                read.reference_start = None  # This will cause an AttributeError
            else:
                read.is_paired = True
                read.is_secondary = False
                read.is_supplementary = False
                read.mapping_quality = 30
                read.reference_name = "test_contig"
                read.reference_start = i * 100
                read.reference_end = read.reference_start + 100
                read.next_reference_name = "test_contig"
                read.next_reference_start = read.reference_start + 200
                read.template_length = 300
                read.is_proper_pair = True
                read.is_reverse = False
                read.mate_is_unmapped = False
                read.query_sequence = "A" * 100
                read.query_length = 100
            mock_reads.append(read)
        
        mock_bam.fetch.return_value = mock_reads
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            with patch('chimeric_detector.core.bam_parser.logging.getLogger') as mock_logger:
                pairs = parser.get_read_pairs_in_window("test_contig", 0, 1000)
                
                # Should get pairs from successful reads
                assert len(pairs) > 0
                # Should log debug messages for failed reads
                mock_logger.return_value.debug.assert_called()
    
    @patch('chimeric_detector.core.bam_parser.pysam.AlignmentFile')
    def test_complete_failure_recovery(self, mock_alignment_file, sample_config):
        """Test recovery from complete method failure."""
        mock_bam = MagicMock()
        mock_bam.references = ["test_contig"]
        mock_bam.check_index.return_value = True
        mock_bam.get_reference_length.side_effect = Exception("Catastrophic failure")
        mock_alignment_file.return_value = mock_bam
        
        parser = BAMParser("test.bam", sample_config.quality)
        
        with parser:
            with patch('chimeric_detector.core.bam_parser.logging.getLogger') as mock_logger:
                windows = list(parser.iterate_windows("test_contig", 1000, 500))
                
                # Should return empty list and log error
                assert len(windows) == 0
                mock_logger.return_value.error.assert_called()