"""Tests for configuration management."""

import pytest
import tempfile
import yaml
from pathlib import Path

from chimeric_detector.config.config import Config, WindowConfig, QualityConfig


def test_default_config():
    """Test default configuration creation."""
    config = Config()
    
    assert config.window.size == 1000
    assert config.quality.min_mapping_quality == 20
    assert config.detection.min_confidence_score == 0.7
    assert config.output.format == "json"


def test_config_from_dict():
    """Test configuration from dictionary."""
    data = {
        "window": {"size": 2000, "step": 200},
        "quality": {"min_mapping_quality": 30},
        "detection": {"min_confidence_score": 0.8}
    }
    
    config = Config.from_dict(data)
    
    assert config.window.size == 2000
    assert config.window.step == 200
    assert config.quality.min_mapping_quality == 30
    assert config.detection.min_confidence_score == 0.8


def test_config_from_yaml():
    """Test configuration from YAML file."""
    config_data = {
        "window": {"size": 1500},
        "output": {"format": "tsv"}
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(config_data, f)
        config_path = Path(f.name)
    
    try:
        config = Config.from_file(config_path)
        assert config.window.size == 1500
        assert config.output.format == "tsv"
    finally:
        config_path.unlink()


def test_config_to_dict():
    """Test configuration serialization."""
    config = Config()
    config.window.size = 2000
    
    data = config.to_dict()
    
    assert data["window"]["size"] == 2000
    assert "quality" in data
    assert "detection" in data