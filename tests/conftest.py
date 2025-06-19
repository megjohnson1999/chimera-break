"""Shared pytest fixtures and test utilities."""

import os
import pytest
import tempfile
import shutil
import json
import yaml
from pathlib import Path
from typing import Dict, List, Any, Generator
from unittest.mock import Mock, MagicMock
import pysam
import numpy as np

from chimeric_detector.config.config import Config, WindowConfig, DetectionConfig, QualityConfig, OutputConfig, ValidationConfig
from chimeric_detector.core.bam_parser import ReadPairInfo
from chimeric_detector.core.window_analysis import Breakpoint, WindowMetrics


@pytest.fixture(scope="session")
def test_data_dir():
    """Get test data directory."""
    test_dir = Path(__file__).parent / "data"
    test_dir.mkdir(exist_ok=True)
    return test_dir


@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


@pytest.fixture
def sample_config():
    """Create a sample configuration for testing."""
    return Config(
        window=WindowConfig(
            size=1000,
            step=100,
            min_coverage=5.0
        ),
        detection=DetectionConfig(
            min_proper_pair_drop=0.3,
            insert_size_zscore_threshold=2.0,
            min_confidence_score=0.6,
            statistical_method="robust"
        ),
        quality=QualityConfig(
            min_mapping_quality=20,
            min_insert_size=50,
            max_insert_size=2000,
            end_read_buffer=500
        ),
        output=OutputConfig(
            format="json",
            include_read_evidence=True
        ),
        validation=ValidationConfig(
            enabled=False
        )
    )


@pytest.fixture
def config_file(temp_dir, sample_config):
    """Create a temporary config file."""
    config_file = temp_dir / "test_config.yaml"
    config_data = sample_config.to_dict()
    
    with open(config_file, 'w') as f:
        yaml.dump(config_data, f)
    
    return config_file


@pytest.fixture
def invalid_config_file(temp_dir):
    """Create an invalid config file for testing error handling."""
    config_file = temp_dir / "invalid_config.yaml"
    
    with open(config_file, 'w') as f:
        f.write("invalid: yaml: content: [\n")  # Deliberately broken YAML
    
    return config_file


@pytest.fixture
def sample_read_pairs():
    """Create sample read pair data for testing."""
    return [
        ReadPairInfo(
            position=100,
            mate_position=200,
            insert_size=150,
            is_proper_pair=True,
            is_discordant=False,
            is_near_end=False,
            contig="test_contig",
            mate_contig="test_contig"
        ),
        ReadPairInfo(
            position=150,
            mate_position=250,
            insert_size=200,
            is_proper_pair=True,
            is_discordant=False,
            is_near_end=False,
            contig="test_contig",
            mate_contig="test_contig"
        ),
        ReadPairInfo(
            position=200,
            mate_position=300,
            insert_size=180,
            is_proper_pair=False,
            is_discordant=True,
            is_near_end=False,
            contig="test_contig",
            mate_contig="other_contig"
        )
    ]


@pytest.fixture
def sample_window_metrics():
    """Create sample window metrics for testing."""
    return [
        WindowMetrics(
            start=0,
            end=1000,
            coverage=10.5,
            read_count=25,
            proper_pair_rate=0.8,
            discordant_pair_rate=0.1,
            insert_size_median=180,
            insert_size_mad=15,
            insert_size_zscore=0.5
        ),
        WindowMetrics(
            start=100,
            end=1100,
            coverage=8.2,
            read_count=20,
            proper_pair_rate=0.4,  # Low proper pair rate (potential breakpoint)
            discordant_pair_rate=0.3,
            insert_size_median=250,
            insert_size_mad=25,
            insert_size_zscore=3.2  # High z-score (potential breakpoint)
        ),
        WindowMetrics(
            start=200,
            end=1200,
            coverage=12.1,
            read_count=30,
            proper_pair_rate=0.85,
            discordant_pair_rate=0.05,
            insert_size_median=175,
            insert_size_mad=12,
            insert_size_zscore=0.3
        )
    ]


@pytest.fixture
def sample_breakpoints():
    """Create sample breakpoints for testing."""
    return [
        Breakpoint(
            contig="test_contig",
            position=550,
            confidence=0.85,
            evidence={
                "proper_pair_drop": 0.4,
                "insert_size_zscore": 3.2,
                "discordant_rate": 0.25
            },
            left_metrics=WindowMetrics(
                start=0, end=1000, coverage=10.0, read_count=25,
                proper_pair_rate=0.8, discordant_pair_rate=0.1,
                insert_size_median=180, insert_size_mad=15, insert_size_zscore=0.5
            ),
            right_metrics=WindowMetrics(
                start=100, end=1100, coverage=8.0, read_count=20,
                proper_pair_rate=0.4, discordant_pair_rate=0.3,
                insert_size_median=250, insert_size_mad=25, insert_size_zscore=3.2
            )
        )
    ]


class MockAlignmentFile:
    """Mock pysam.AlignmentFile for testing."""
    
    def __init__(self, references=None, reference_lengths=None, reads=None):
        self.references = references or ["test_contig", "contig2"]
        self.reference_lengths = reference_lengths or [10000, 5000]
        self.reads = reads or []
        self.closed = False
        
        # Create reference name to length mapping
        self._ref_lengths = dict(zip(self.references, self.reference_lengths))
    
    def get_reference_length(self, reference):
        if reference not in self._ref_lengths:
            raise ValueError(f"Reference {reference} not found")
        return self._ref_lengths[reference]
    
    def fetch(self, reference=None, start=None, end=None):
        """Mock fetch method that returns filtered reads."""
        if reference and reference not in self.references:
            return []
        
        filtered_reads = []
        for read in self.reads:
            # Simple filtering logic
            if reference and hasattr(read, 'reference_name') and read.reference_name != reference:
                continue
            if start is not None and hasattr(read, 'reference_start') and read.reference_start < start:
                continue
            if end is not None and hasattr(read, 'reference_end') and read.reference_end > end:
                continue
            filtered_reads.append(read)
        
        return filtered_reads
    
    def check_index(self):
        """Mock index check - always passes for testing."""
        return True
    
    def close(self):
        self.closed = True
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class MockAlignedSegment:
    """Mock pysam.AlignedSegment for testing."""
    
    def __init__(self, **kwargs):
        # Set default values
        self.reference_name = kwargs.get('reference_name', 'test_contig')
        self.reference_start = kwargs.get('reference_start', 100)
        self.reference_end = kwargs.get('reference_end', self.reference_start + 100)
        self.next_reference_name = kwargs.get('next_reference_name', self.reference_name)
        self.next_reference_start = kwargs.get('next_reference_start', self.reference_start + 200)
        self.template_length = kwargs.get('template_length', 300)
        self.mapping_quality = kwargs.get('mapping_quality', 30)
        self.is_paired = kwargs.get('is_paired', True)
        self.is_proper_pair = kwargs.get('is_proper_pair', True)
        self.is_secondary = kwargs.get('is_secondary', False)
        self.is_supplementary = kwargs.get('is_supplementary', False)
        self.is_reverse = kwargs.get('is_reverse', False)
        self.mate_is_unmapped = kwargs.get('mate_is_unmapped', False)
        self.query_sequence = kwargs.get('query_sequence', 'A' * 100)
        self.query_length = kwargs.get('query_length', len(self.query_sequence) if self.query_sequence else 100)


@pytest.fixture
def mock_bam_file():
    """Create a mock BAM file for testing."""
    reads = [
        MockAlignedSegment(
            reference_start=100, reference_end=200, template_length=300,
            is_proper_pair=True, mapping_quality=30
        ),
        MockAlignedSegment(
            reference_start=150, reference_end=250, template_length=350,
            is_proper_pair=True, mapping_quality=25
        ),
        MockAlignedSegment(
            reference_start=200, reference_end=300, template_length=400,
            is_proper_pair=False, mapping_quality=15,
            next_reference_name='other_contig'  # Inter-contig pair
        )
    ]
    
    return MockAlignmentFile(
        references=["test_contig", "other_contig"],
        reference_lengths=[10000, 5000],
        reads=reads
    )


@pytest.fixture
def test_files_dir(test_data_dir):
    """Create directory for test files."""
    files_dir = test_data_dir / "files"
    files_dir.mkdir(exist_ok=True)
    return files_dir


@pytest.fixture
def sample_json_output(temp_dir, sample_breakpoints):
    """Create sample JSON output for testing."""
    output_file = temp_dir / "test_output.json"
    
    data = {
        "metadata": {
            "bam_file": "test.bam",
            "total_contigs": 2,
            "config": {}
        },
        "timestamp": "2024-01-01T00:00:00",
        "breakpoints": [
            {
                "contig": bp.contig,
                "position": bp.position,
                "confidence": bp.confidence,
                "evidence": bp.evidence
            } for bp in sample_breakpoints
        ]
    }
    
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    return output_file


@pytest.fixture
def security_test_paths(temp_dir):
    """Create various paths for security testing."""
    paths = {
        'safe_path': temp_dir / "safe_file.txt",
        'traversal_path': Path("../../../etc/passwd"),
        'nonexistent_path': temp_dir / "nonexistent.txt",
        'large_file': temp_dir / "large_file.txt"
    }
    
    # Create safe file
    with open(paths['safe_path'], 'w') as f:
        f.write("Safe content")
    
    # Create large file (for size testing)
    with open(paths['large_file'], 'w') as f:
        f.write("x" * 1000000)  # 1MB file
    
    return paths


def create_test_bam_file(output_path: Path, references: List[str] = None, 
                        reference_lengths: List[int] = None) -> Path:
    """Create a real BAM file for testing (requires pysam)."""
    if references is None:
        references = ["test_contig", "contig2"]
    if reference_lengths is None:
        reference_lengths = [10000, 5000]
    
    # Create header
    header = {
        'HD': {'VN': '1.6'},
        'SQ': [{'LN': length, 'SN': ref} for ref, length in zip(references, reference_lengths)]
    }
    
    # Create BAM file
    with pysam.AlignmentFile(str(output_path), "wb", header=header) as outf:
        # Add some sample reads
        for i in range(100):
            read = pysam.AlignedSegment()
            read.query_name = f"read_{i}"
            read.query_sequence = "A" * 100
            read.flag = 99 if i % 2 == 0 else 147  # Proper pair flags
            read.reference_id = 0
            read.reference_start = i * 50
            read.mapping_quality = 30
            read.cigar = [(0, 100)]  # 100M
            read.next_reference_id = 0
            read.next_reference_start = i * 50 + 200
            read.template_length = 300
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            
            outf.write(read)
    
    # Index the BAM file
    pysam.index(str(output_path))
    
    return output_path


@pytest.fixture
def real_test_bam(test_files_dir):
    """Create a real BAM file for integration testing."""
    bam_path = test_files_dir / "test.bam"
    if not bam_path.exists():
        try:
            return create_test_bam_file(bam_path)
        except Exception:
            pytest.skip("Cannot create test BAM file (pysam not available or permission issue)")
    return bam_path


# Performance testing utilities
class PerformanceTimer:
    """Simple performance timer for testing."""
    
    def __init__(self):
        import time
        self.start_time = None
        self.end_time = None
        self.time = time
    
    def __enter__(self):
        self.start_time = self.time.time()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = self.time.time()
    
    @property
    def elapsed(self):
        if self.start_time is None or self.end_time is None:
            return None
        return self.end_time - self.start_time


@pytest.fixture
def performance_timer():
    """Performance timer fixture."""
    return PerformanceTimer


# Test data generators
def generate_read_pairs(count: int, proper_pair_rate: float = 0.8) -> List[ReadPairInfo]:
    """Generate synthetic read pairs for testing."""
    pairs = []
    for i in range(count):
        is_proper = np.random.random() < proper_pair_rate
        pairs.append(ReadPairInfo(
            position=i * 10,
            mate_position=i * 10 + 200 + np.random.randint(-50, 50),
            insert_size=200 + np.random.randint(-100, 100),
            is_proper_pair=is_proper,
            is_discordant=not is_proper and np.random.random() < 0.5,
            mapping_quality=30 + np.random.randint(-10, 10),
            is_near_end=i < 10 or i > count - 10,
            contig="test_contig",
            mate_contig="test_contig" if is_proper else ("other_contig" if np.random.random() < 0.3 else "test_contig")
        ))
    return pairs


@pytest.fixture
def large_read_dataset():
    """Generate a large dataset for performance testing."""
    return generate_read_pairs(10000, proper_pair_rate=0.7)


# Parametrized test data
@pytest.fixture(params=[
    {"format": "json", "include_evidence": True},
    {"format": "tsv", "include_evidence": False},
    {"format": "bed", "include_evidence": False}
])
def output_format_config(request):
    """Parametrized output format configurations."""
    return OutputConfig(
        format=request.param["format"],
        include_read_evidence=request.param["include_evidence"]
    )


@pytest.fixture(params=["robust", "parametric", "nonparametric"])
def statistical_methods(request):
    """Parametrized statistical methods."""
    return request.param


# Error simulation fixtures
@pytest.fixture
def corrupted_files(temp_dir):
    """Create various corrupted files for error testing."""
    files = {}
    
    # Empty file
    empty_file = temp_dir / "empty.yaml"
    empty_file.touch()
    files['empty'] = empty_file
    
    # Binary file with YAML extension
    binary_file = temp_dir / "binary.yaml"
    with open(binary_file, 'wb') as f:
        f.write(b'\x00\x01\x02\x03\x04\x05')
    files['binary'] = binary_file
    
    # Very large file
    large_file = temp_dir / "large.yaml"
    with open(large_file, 'w') as f:
        f.write("key: " + "x" * 10000000)  # 10MB+ file
    files['large'] = large_file
    
    return files