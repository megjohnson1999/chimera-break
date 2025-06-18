#!/usr/bin/env python3
"""
Chimeric contig detection tool for viral metagenomic co-assemblies.

This tool uses read-pair analysis to detect potential chimeric junctions
in assembled contigs without relying on coverage-based methods.
"""

import argparse
import sys
import logging
from pathlib import Path
from typing import List, Optional

from chimeric_detector.config.config import Config
from chimeric_detector.core.bam_parser import BAMParser
from chimeric_detector.core.window_analysis import WindowAnalyzer, WindowMetrics
from chimeric_detector.modules.validation import ValidationPipeline
from chimeric_detector.utils.output import OutputFormatter, Logger


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Detect chimeric contigs using read-pair analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("bam_file", type=str,
                       help="Input BAM file (must be sorted and indexed)")
    parser.add_argument("output", type=str,
                       help="Output file path")
    
    # Optional arguments
    parser.add_argument("-c", "--config", type=str,
                       help="Configuration file (YAML or JSON)")
    parser.add_argument("--contigs", type=str, nargs="+",
                       help="Specific contigs to analyze (default: all)")
    
    # Window parameters
    window_group = parser.add_argument_group("Window parameters")
    window_group.add_argument("--window-size", type=int,
                            help="Sliding window size in bp")
    window_group.add_argument("--window-step", type=int,
                            help="Window step size in bp")
    window_group.add_argument("--min-coverage", type=float,
                            help="Minimum coverage required in window")
    
    # Quality filters
    quality_group = parser.add_argument_group("Quality filters")
    quality_group.add_argument("--min-mapping-quality", type=int,
                             help="Minimum mapping quality")
    quality_group.add_argument("--max-insert-size", type=int,
                             help="Maximum insert size to consider")
    
    # Detection thresholds
    detection_group = parser.add_argument_group("Detection thresholds")
    detection_group.add_argument("--min-proper-pair-drop", type=float,
                               help="Minimum drop in proper pair rate")
    detection_group.add_argument("--insert-zscore-threshold", type=float,
                               help="Insert size z-score threshold")
    detection_group.add_argument("--min-confidence", type=float,
                               help="Minimum confidence score for reporting")
    
    # Validation options
    validation_group = parser.add_argument_group("Validation options")
    validation_group.add_argument("--validate-taxonomy", action="store_true",
                                help="Enable taxonomy validation")
    validation_group.add_argument("--taxonomy-db", type=str,
                                help="Path to taxonomy database")
    validation_group.add_argument("--evaluate-splits", action="store_true",
                                help="Evaluate quality of potential splits")
    
    # Output options
    output_group = parser.add_argument_group("Output options")
    output_group.add_argument("--output-format", choices=["json", "tsv", "bed"],
                            help="Output format")
    output_group.add_argument("--include-evidence", action="store_true",
                            help="Include detailed read evidence in output")
    output_group.add_argument("-v", "--verbose", action="store_true",
                            help="Enable verbose output")
    output_group.add_argument("-d", "--debug", action="store_true",
                            help="Enable debug output")
    output_group.add_argument("--log-file", type=str,
                            help="Write logs to file")
    
    return parser.parse_args()


def detect_chimeras_in_contig(bam_parser: BAMParser, analyzer: WindowAnalyzer,
                            contig: str, config: Config) -> List:
    """Detect chimeras in a single contig."""
    logger = logging.getLogger(__name__)
    logger.info(f"Processing contig: {contig}")
    
    # Estimate insert size distribution
    insert_stats = bam_parser.estimate_insert_size_distribution(contig)
    logger.debug(f"Insert size stats for {contig}: {insert_stats}")
    
    # Collect window metrics
    window_metrics = []
    
    for start, end, pairs in bam_parser.iterate_windows(
        contig, config.window.size, config.window.step
    ):
        coverage = bam_parser.calculate_coverage_in_window(contig, start, end)
        metrics = analyzer.calculate_window_metrics(pairs, coverage)
        
        if metrics:
            metrics.start = start
            metrics.end = end
            window_metrics.append(metrics)
    
    logger.info(f"Collected metrics for {len(window_metrics)} windows in {contig}")
    
    # Detect breakpoints
    breakpoints = analyzer.detect_breakpoints(window_metrics, contig, insert_stats)
    logger.info(f"Found {len(breakpoints)} potential breakpoints in {contig}")
    
    return breakpoints


def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Set up logging
    Logger.setup_logging(
        verbose=args.verbose,
        debug=args.debug,
        log_file=Path(args.log_file) if args.log_file else None
    )
    
    logger = logging.getLogger(__name__)
    logger.info("Starting chimeric contig detection")
    
    # Load configuration
    if args.config:
        config = Config.from_file(Path(args.config))
        logger.info(f"Loaded configuration from {args.config}")
    else:
        config = Config()
        
    # Update config with command line arguments
    config.update_from_args(args)
    
    # Initialize components
    analyzer = WindowAnalyzer(config.window, config.detection, config.quality)
    formatter = OutputFormatter(config.output)
    
    all_breakpoints = []
    
    # Process BAM file
    try:
        with BAMParser(args.bam_file, config.quality) as bam_parser:
            # Get contigs to process
            if args.contigs:
                contigs = args.contigs
            else:
                contigs = bam_parser.get_contigs()
                
            logger.info(f"Processing {len(contigs)} contigs")
            
            # Process each contig
            for contig in contigs:
                try:
                    breakpoints = detect_chimeras_in_contig(
                        bam_parser, analyzer, contig, config
                    )
                    all_breakpoints.extend(breakpoints)
                except Exception as e:
                    logger.error(f"Error processing contig {contig}: {e}")
                    if config.output.debug:
                        raise
                        
    except Exception as e:
        logger.error(f"Error processing BAM file: {e}")
        sys.exit(1)
        
    logger.info(f"Total breakpoints found: {len(all_breakpoints)}")
    
    # Run validation if enabled
    validation_results = None
    if config.validation.enabled:
        logger.info("Running validation modules")
        pipeline = ValidationPipeline(config.validation)
        
        # Get contig lengths for split evaluation
        contig_lengths = {}
        with BAMParser(args.bam_file, config.quality) as bam_parser:
            for contig in bam_parser.get_contigs():
                contig_lengths[contig] = bam_parser._bam.get_reference_length(contig)
                
        validation_results = pipeline.run_validation(
            all_breakpoints,
            contig_lengths=contig_lengths
        )
        
    # Write results
    metadata = {
        "bam_file": args.bam_file,
        "total_contigs": len(contigs),
        "config": config.to_dict()
    }
    
    if validation_results:
        metadata["validation"] = validation_results
        
    formatter.write_results(all_breakpoints, Path(args.output), metadata)
    logger.info(f"Results written to {args.output}")
    

if __name__ == "__main__":
    main()