"""Configuration management for chimeric contig detector."""

import json
import yaml
from dataclasses import dataclass, field, asdict
from typing import Dict, Any, Optional, List
from pathlib import Path


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
        """Load configuration from YAML or JSON file."""
        with open(config_path) as f:
            if config_path.suffix in ['.yaml', '.yml']:
                data = yaml.safe_load(f)
            else:
                data = json.load(f)
        
        return cls.from_dict(data)
    
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