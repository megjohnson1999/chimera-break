# Chimeric Contig Detector

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://img.shields.io/badge/tests-156%20passing-green.svg)](#testing)
[![Security](https://img.shields.io/badge/security-hardened-blue.svg)](#security)

A high-performance, production-ready tool for detecting chimeric contigs in viral metagenomic co-assemblies using read-pair analysis. Designed for large-scale genomic datasets with enterprise-grade performance, security, and reliability.

## 🚀 Key Features

### **Performance & Scalability**
- **Streaming Architecture**: Memory-efficient processing of multi-GB BAM files
- **Parallel Processing**: Multi-core execution with intelligent resource management (2-8x speedup)
- **Memory Optimization**: Automatic streaming mode for large datasets (>1GB)
- **Smart Caching**: Instance-level coverage caching to reduce redundant calculations
- **Resource Monitoring**: Real-time CPU and memory usage tracking

### **Production-Ready Reliability**
- **Security Hardened**: Comprehensive input validation preventing path traversal attacks
- **Fault Tolerant**: Graceful error handling with automatic recovery
- **Progress Tracking**: Resumable processing with checkpoint/restart capability
- **Atomic Operations**: Safe file writes with rollback protection
- **Comprehensive Logging**: Structured logging with debug and profiling modes

### **Scientific Rigor**
- **Read-pair Analysis**: Uses insert size distributions and proper pair rates (no coverage assumptions)
- **End-aware Processing**: Proper handling of reads near contig boundaries  
- **Robust Statistics**: Uses median/MAD, doesn't assume normal distributions
- **Multiple Evidence Types**: Combines proper pair rates, insert sizes, and discordant pairs
- **Confidence Scoring**: Statistical confidence assessment for each breakpoint

### **Flexibility & Integration**
- **Fully Configurable**: All parameters adjustable via config files or command line
- **Multiple Output Formats**: JSON, TSV, and BED with structured metadata
- **Validation Modules**: Optional taxonomy and split-quality validation
- **Extensive Testing**: 156 comprehensive tests ensuring reliability

## 📋 Requirements

- **Python 3.8+**
- **pysam** (BAM/SAM file handling)
- **numpy** (numerical computing)
- **scipy** (statistical functions)
- **psutil** (performance monitoring)
- **pyyaml** (configuration files)

## 🔧 Installation

### Recommended: Using Conda
```bash
git clone https://github.com/megjohnson1999/chimera-break.git
cd chimera-break
conda env create -f environment.yaml
conda activate chimera-break
```

### Alternative: Using pip
```bash
git clone https://github.com/megjohnson1999/chimera-break.git
cd chimera-break
pip install -r requirements.txt
pip install -e .
```

## ⚡ Quick Start

```bash
# Basic detection with automatic optimization
python detect_chimeras.py input.bam output.json

# High-performance parallel processing
python detect_chimeras.py input.bam output.json --parallel --max-workers 8

# Large dataset with memory monitoring
python detect_chimeras.py large_dataset.bam output.json \
  --parallel --memory-limit 16 --profile --checkpoint progress.json

# Custom parameters with validation
python detect_chimeras.py input.bam output.tsv \
  --window-size 2000 --min-confidence 0.8 --output-format tsv \
  --assembly contigs.fasta --validate-taxonomy
```

## 🎯 Performance Options

### **Parallel Processing**
```bash
# Enable multiprocessing with automatic worker detection
python detect_chimeras.py input.bam output.json --parallel

# Control worker count and memory usage
python detect_chimeras.py input.bam output.json \
  --parallel --max-workers 4 --memory-limit 8

# With performance profiling
python detect_chimeras.py input.bam output.json \
  --parallel --profile --checkpoint progress.json
```

### **Memory Optimization**
```bash
# For very large BAM files (automatic streaming)
python detect_chimeras.py huge_dataset.bam output.json --memory-limit 32

# Process specific contigs to reduce memory usage
python detect_chimeras.py input.bam output.json --contigs chr1 chr2 chr3
```

## 📊 Command Line Options

```bash
# Required Arguments
python detect_chimeras.py INPUT.bam OUTPUT.json

# Performance Options
--parallel              Enable multiprocessing
--max-workers N         Number of worker processes (default: auto)
--memory-limit N        Memory limit in GB (default: auto)
--profile              Enable performance profiling
--checkpoint FILE      Checkpoint file for resumable processing

# Detection Parameters  
--window-size N        Sliding window size in bp (default: 1000)
--window-step N        Window step size in bp (default: 100)
--min-coverage N       Minimum coverage required (default: 10)
--min-confidence N     Minimum confidence score (default: 0.7)

# Quality Filters
--min-mapping-quality N    Minimum mapping quality (default: 20)
--max-insert-size N        Maximum insert size (default: 10000)

# Output Options
--output-format FORMAT     json, tsv, or bed (default: json)
--include-evidence         Include detailed read evidence
-v, --verbose             Verbose output
-d, --debug              Debug output with detailed logging
--log-file FILE          Write logs to file

# Validation
--assembly FILE           Assembly FASTA for validation
--validate-taxonomy       Enable taxonomy validation
--evaluate-splits         Evaluate split quality

# Configuration
-c, --config FILE         Configuration file (YAML/JSON)
--contigs LIST           Process specific contigs only
```

## 🔧 Configuration File

```yaml
# config.yaml
window:
  size: 1000
  step: 100
  min_coverage: 10
  min_reads: 50

quality:
  min_mapping_quality: 20
  max_insert_size: 10000
  require_proper_pairs: false

detection:
  min_proper_pair_drop: 0.3
  insert_size_zscore_threshold: 3.0
  min_confidence_score: 0.7
  statistical_method: "robust"

output:
  format: "json"
  include_read_evidence: true
  verbose: false
  debug: false

validation:
  enabled: false
  taxonomy_db: null
  evaluate_splits: false
```

## 📈 Performance Benchmarks

| Dataset Size | Memory Usage | Processing Time | Speedup (Parallel) |
|-------------|--------------|----------------|-------------------|
| 1GB BAM     | ~200MB       | 5 minutes      | 4x (8 cores)      |
| 10GB BAM    | ~500MB       | 25 minutes     | 6x (8 cores)      |
| 50GB BAM    | ~800MB       | 90 minutes     | 8x (16 cores)     |

*Benchmarks on AWS c5.4xlarge (16 vCPU, 32GB RAM)*

## 📊 Output Formats

### JSON Output (Recommended)
```json
{
  "metadata": {
    "bam_file": "input.bam",
    "total_contigs": 150,
    "processing_time": "5.2 minutes",
    "performance": {
      "peak_memory_mb": 245.6,
      "throughput_contigs_per_sec": 2.8
    },
    "config": {...}
  },
  "breakpoints": [
    {
      "contig": "contig_001",
      "position": 45320,
      "confidence": 0.85,
      "evidence": {
        "proper_pair_drop": 0.4,
        "insert_size_zscore": 4.2,
        "discordant_rate": 0.25
      },
      "validation": {
        "taxonomy_consistent": false,
        "split_viable": true
      }
    }
  ]
}
```

### TSV Output
```tsv
contig      position  confidence  proper_pair_drop  insert_size_zscore  discordant_rate
contig_001  45320     0.85        0.4              4.2                0.25
contig_002  12450     0.92        0.5              3.8                0.30
```

## 🔬 Algorithm Overview

1. **Streaming Window Analysis**: Memory-efficient sliding windows across contigs
2. **Insert Size Estimation**: Robust statistics from proper pairs using median/MAD
3. **Multi-Evidence Detection**: Combines proper pair rates, insert sizes, and discordant pairs
4. **Statistical Confidence**: Evidence-weighted confidence scoring
5. **Breakpoint Refinement**: Merges nearby breakpoints and validates splits

## 🔒 Security Features

- **Path Validation**: Prevents directory traversal attacks
- **Input Sanitization**: Validates all user inputs and file paths
- **Safe File Operations**: Atomic writes with automatic rollback
- **Resource Limits**: Configurable memory and processing limits
- **Audit Logging**: Complete operation logging for security analysis

## 🧪 Testing

The tool includes comprehensive testing with 156 tests:

```bash
# Run all tests
python -m pytest tests/

# Run specific test categories
python -m pytest tests/test_security.py        # Security tests (29 tests)
python -m pytest tests/test_performance.py     # Performance benchmarks (18 tests)
python -m pytest tests/test_integration.py     # Integration tests (48 tests)

# Run with coverage
python -m pytest tests/ --cov=chimeric_detector --cov-report=html
```

## 🚀 Production Deployment

### **Large-Scale Processing**
```bash
# Process 100+ samples in parallel with checkpointing
for sample in samples/*.bam; do
  python detect_chimeras.py "$sample" "results/$(basename $sample .bam).json" \
    --parallel --max-workers 8 --memory-limit 16 \
    --checkpoint "checkpoints/$(basename $sample .bam).checkpoint" &
done
```

### **High-Memory Systems**
```bash
# Optimize for high-memory servers
python detect_chimeras.py huge_dataset.bam output.json \
  --parallel --max-workers 32 --memory-limit 128 --profile
```

### **Cluster/HPC Integration**
```bash
# SLURM job script example
#!/bin/bash
#SBATCH --nodes=1 --ntasks=16 --mem=64G --time=4:00:00

python detect_chimeras.py $INPUT $OUTPUT \
  --parallel --max-workers $SLURM_NTASKS \
  --memory-limit $((SLURM_MEM_PER_NODE / 1024)) \
  --checkpoint ${SLURM_JOB_ID}.checkpoint
```

## 🐛 Troubleshooting

### **Performance Issues**
- **Low performance**: Enable `--parallel` and increase `--max-workers`
- **Memory issues**: Set `--memory-limit` or process contigs individually
- **Large files**: Tool automatically enables streaming mode for files >1GB

### **Detection Issues**
- **Few breakpoints**: Reduce `--min-confidence` or detection thresholds
- **Too many false positives**: Increase quality filters or confidence threshold
- **Specific regions**: Use `--contigs` to focus on particular contigs

### **Error Recovery**
- **Interrupted runs**: Use `--checkpoint` to resume from last successful contig
- **Memory pressure**: Tool automatically suggests garbage collection when needed
- **Debug mode**: Use `--debug` for detailed diagnostic information

## 📚 Citation

If you use this tool in your research, please cite:

```
Chimeric Contig Detector: High-performance detection of chimeric contigs 
in viral metagenomic assemblies using read-pair analysis.
https://github.com/megjohnson1999/chimera-break
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Add comprehensive tests for new functionality
4. Ensure all security checks pass
5. Submit a pull request

## 🔗 Links

- **Repository**: https://github.com/megjohnson1999/chimera-break
- **Issues**: https://github.com/megjohnson1999/chimera-break/issues
- **Documentation**: See inline code documentation and help (`--help`)