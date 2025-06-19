"""Configuration management for chimeric contig detector."""

import json
import yaml
import logging
from dataclasses import dataclass, field, asdict
from typing import Dict, Any, Optional, List
from pathlib import Path

from ..utils.security import safe_file_open, SecurityError


@dataclass
class WindowConfig:
    """Sliding window parameters."""
    size: int = 1000
    step: int = 100
    min_coverage: int = 10
    min_reads: int = 50
    
    
@dataclass
class QualityConfig:
    """Quality filtering parameters."""
    min_mapping_quality: int = 20
    min_base_quality: int = 20
    max_insert_size: int = 10000
    min_insert_size: int = 0
    require_proper_pairs: bool = False
    end_read_buffer: int = 300  # Distance from contig end to consider as "end read"
    lenient_end_reads: bool = True  # Be more lenient with end read classification
    
    
@dataclass
class DetectionConfig:
    """Breakpoint detection parameters."""
    min_proper_pair_drop: float = 0.3
    insert_size_zscore_threshold: float = 3.0
    discordant_pair_threshold: float = 0.2
    min_confidence_score: float = 0.7
    statistical_method: str = "robust"  # "robust", "parametric", "nonparametric"
    
    
@dataclass
class OutputConfig:
    """Output formatting options."""
    format: str = "json"  # "json", "tsv", "bed"
    include_read_evidence: bool = False
    verbose: bool = False
    debug: bool = False
    
    
@dataclass
class ValidationConfig:
    """Optional validation parameters."""
    enabled: bool = False
    taxonomy_db: Optional[str] = None
    evaluate_splits: bool = False
    min_split_length: int = 1000
    

@dataclass
class Config:
    """Main configuration container."""
    window: WindowConfig = field(default_factory=WindowConfig)
    quality: QualityConfig = field(default_factory=QualityConfig)
    detection: DetectionConfig = field(default_factory=DetectionConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    validation: ValidationConfig = field(default_factory=ValidationConfig)
    
    @classmethod
    def from_file(cls, config_path: Path) -> "Config":
        """Load configuration from YAML or JSON file with comprehensive error handling."""
        logger = logging.getLogger(__name__)
        
        try:
            # Use secure file opening
            with safe_file_open(config_path, 'r') as f:
                if config_path.suffix.lower() in ['.yaml', '.yml']:
                    try:
                        data = yaml.safe_load(f)
                    except yaml.YAMLError as e:
                        raise ValueError(f"Invalid YAML syntax in {config_path}: {e}")
                elif config_path.suffix.lower() == '.json':
                    try:
                        data = json.load(f)
                    except json.JSONDecodeError as e:
                        raise ValueError(f"Invalid JSON syntax in {config_path}: {e}")
                else:
                    raise ValueError(f"Unsupported config file format: {config_path.suffix}")
                
                # Validate that data was loaded
                if data is None:
                    logger.warning(f"Config file {config_path} is empty, using defaults")
                    data = {}
                elif not isinstance(data, dict):
                    raise ValueError(f"Config file must contain a dictionary, got {type(data)}")
                    
        except SecurityError as e:
            raise ValueError(f"Security error loading config: {e}")
        except FileNotFoundError:
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        except PermissionError:
            raise PermissionError(f"Permission denied reading config file: {config_path}")
        except Exception as e:
            raise RuntimeError(f"Unexpected error loading config from {config_path}: {e}")
        
        try:
            return cls.from_dict(data)
        except Exception as e:
            raise ValueError(f"Invalid configuration structure in {config_path}: {e}")
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Config":
        """Create config from dictionary."""
        config = cls()
        
        if 'window' in data:
            config.window = WindowConfig(**data['window'])
        if 'quality' in data:
            config.quality = QualityConfig(**data['quality'])
        if 'detection' in data:
            config.detection = DetectionConfig(**data['detection'])
        if 'output' in data:
            config.output = OutputConfig(**data['output'])
        if 'validation' in data:
            config.validation = ValidationConfig(**data['validation'])
            
        return config
    
    def update_from_args(self, args: Any) -> None:
        """Update config from command line arguments."""
        arg_dict = vars(args)
        
        # Map command line args to config sections
        mappings = {
            'window_size': ('window', 'size'),
            'window_step': ('window', 'step'),
            'min_coverage': ('window', 'min_coverage'),
            'min_mapping_quality': ('quality', 'min_mapping_quality'),
            'min_proper_pair_drop': ('detection', 'min_proper_pair_drop'),
            'output_format': ('output', 'format'),
            'verbose': ('output', 'verbose'),
            'debug': ('output', 'debug'),
            'validate_taxonomy': ('validation', 'enabled'),
            'taxonomy_db': ('validation', 'taxonomy_db'),
        }
        
        for arg_name, (section, field_name) in mappings.items():
            if arg_name in arg_dict and arg_dict[arg_name] is not None:
                setattr(getattr(self, section), field_name, arg_dict[arg_name])
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary."""
        return {
            'window': asdict(self.window),
            'quality': asdict(self.quality),
            'detection': asdict(self.detection),
            'output': asdict(self.output),
            'validation': asdict(self.validation),
        }
    
    def save(self, path: Path) -> None:
        """Save configuration to file."""
        data = self.to_dict()
        with open(path, 'w') as f:
            if path.suffix in ['.yaml', '.yml']:
                yaml.dump(data, f, default_flow_style=False)
            else:
                json.dump(data, f, indent=2)