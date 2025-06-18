# Chimeric Contig Detector

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A robust tool for detecting chimeric contigs in viral metagenomic co-assemblies using read-pair analysis. Designed specifically for large-scale metagenomic datasets with 340+ samples.

## Features

- **Read-pair based detection**: Uses insert size distributions and proper pair rates (no coverage assumptions)
- **Fully configurable**: All parameters adjustable via config files or command line  
- **End-aware analysis**: Proper handling of reads near contig boundaries
- **Robust statistics**: Uses median/MAD, doesn't assume normal distributions
- **Modular validation**: Optional taxonomy and split-quality validation modules
- **Multiple output formats**: JSON, TSV, and BED formats with structured logging
- **Production ready**: Comprehensive error handling, logging, and debugging modes

## Installation

### Recommended: Using Conda (handles bioinformatics dependencies better)
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

### Dependencies
- **Python 3.8+**
- **pysam** (BAM/SAM file handling) - easier via conda
- **numpy** (numerical computing)
- **scipy** (statistical functions) 
- **pyyaml** (configuration files)

Note: `pysam` requires system libraries (htslib, zlib, etc.) that conda handles automatically.

## Quick Start

```bash
# Install with conda (recommended)
git clone https://github.com/megjohnson1999/chimera-break.git
cd chimera-break
conda env create -f environment.yaml
conda activate chimera-break

# Run detection (your BAM must be sorted and indexed)
python detect_chimeras.py input.bam output.json --debug

# With custom parameters
python detect_chimeras.py input.bam output.tsv \
  --window-size 2000 --min-confidence 0.8 --output-format tsv

# With taxonomy validation (requires assembly FASTA)
python detect_chimeras.py input.bam output.json \
  --assembly contigs.fasta --validate-taxonomy --taxonomy-db dummy
```

## Basic Usage

```bash
# Basic detection
python detect_chimeras.py input.bam output.json

# With custom configuration
python detect_chimeras.py input.bam output.json -c config.yaml

# Specific contigs only
python detect_chimeras.py input.bam output.json --contigs contig1 contig2

# With validation
python detect_chimeras.py input.bam output.json --validate-taxonomy --evaluate-splits
```

## Configuration

### What's Easily Configurable:

**Detection Parameters:**
- Window size and step size
- Quality filters (mapping quality, insert size limits)
- Statistical thresholds (proper pair drop, z-score cutoffs)
- Confidence scoring parameters

**Output Options:**
- Format (JSON/TSV/BED)
- Verbosity levels
- Evidence inclusion

**Validation Modules:**
- Taxonomy validation (with database)
- Split quality evaluation

### What Requires Code Changes:

**Statistical Methods:**
- New anomaly detection algorithms (implement `StatisticalMethod` class)
- Different confidence scoring approaches (modify `_calculate_confidence`)

**Validation Modules:**
- Custom validation logic (implement `ValidationModule` class)
- New evidence types (extend `Breakpoint` class)

**Output Formats:**
- Additional formats (extend `OutputFormatter`)

## Configuration File Example

```yaml
window:
  size: 1000
  step: 100
  min_coverage: 10

detection:
  min_proper_pair_drop: 0.3
  insert_size_zscore_threshold: 3.0
  min_confidence_score: 0.7
  statistical_method: "robust"

output:
  format: "json"
  include_read_evidence: true
```

## Output Format

### JSON Output
```json
{
  "metadata": {
    "bam_file": "input.bam",
    "total_contigs": 150,
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
      }
    }
  ]
}
```

### TSV Output
```
contig  position  confidence  proper_pair_drop  insert_size_zscore  discordant_rate
contig_001  45320  0.85  0.4  4.2  0.25
```

## Algorithm Overview

1. **Window-based Analysis**: Slides windows across contigs
2. **Insert Size Estimation**: Calculates robust statistics from proper pairs
3. **Anomaly Detection**: Identifies windows with:
   - Significant drops in proper pair rates
   - Insert size distribution anomalies
   - High discordant pair rates
4. **Confidence Scoring**: Combines multiple evidence types
5. **Breakpoint Merging**: Consolidates nearby breakpoints

## Validation

The tool includes optional validation modules:

### Taxonomy Validation
- **Purpose**: Checks if breakpoints separate regions with different taxonomic classifications
- **Requirements**: Assembly FASTA file (`--assembly contigs.fasta`)
- **Current Implementation**: Simple GC-content based classification (placeholder)
- **Future**: Can be extended with BLAST, Kraken2, or MMseqs2 classification
- **Usage**: 
  ```bash
  python detect_chimeras.py input.bam output.json \
    --assembly contigs.fasta --validate-taxonomy --taxonomy-db dummy
  ```

### Split Evaluation
- **Purpose**: Assesses viability of resulting contig segments after splitting
- **Requirements**: None (uses BAM header for contig lengths)
- **Function**: Checks if splits would create segments above minimum length threshold

## Command Line Options

```
Required:
  bam_file              Input BAM file (sorted and indexed)
  output                Output file path

Window parameters:
  --window-size         Sliding window size in bp
  --window-step         Window step size in bp
  --min-coverage        Minimum coverage required

Detection thresholds:
  --min-proper-pair-drop    Minimum drop in proper pair rate
  --insert-zscore-threshold Insert size z-score threshold
  --min-confidence          Minimum confidence score

Output options:
  --output-format       json, tsv, or bed
  --include-evidence    Include detailed evidence
  -v, --verbose         Verbose output
  -d, --debug          Debug output
```

## Testing on Real Data

The tool is designed for immediate testing:

1. Ensure BAM file is sorted and indexed
2. Start with default parameters
3. Use `--debug` mode to understand detection logic
4. Adjust thresholds based on your data characteristics
5. Enable validation modules as needed

## Troubleshooting

- **Low breakpoint counts**: Reduce `min_confidence_score` or detection thresholds
- **Too many false positives**: Increase quality filters or confidence threshold
- **Performance issues**: Increase window step size or limit to specific contigs
- **Memory issues**: Process contigs individually using `--contigs`