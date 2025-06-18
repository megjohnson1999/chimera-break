# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0] - 2024-06-18

### Added
- Initial release of chimeric contig detector
- Read-pair based detection algorithm with sliding window analysis  
- Configurable parameters for all detection thresholds
- Robust statistical methods (median/MAD-based)
- Proper handling of reads near contig ends
- Multiple output formats (JSON, TSV, BED)
- Optional validation modules (taxonomy, split quality)
- Comprehensive logging with debug/verbose modes
- Full configuration system with YAML/JSON support
- Command-line interface with extensive options

### Features
- Insert size distribution estimation from data (no assumptions)
- Pluggable statistical methods (robust/parametric/nonparametric)
- Modular validation framework  
- Production-ready error handling and logging
- Designed for 340+ sample viral metagenomic co-assemblies