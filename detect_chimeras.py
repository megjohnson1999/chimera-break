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
from typing import List, Optional, Tuple, Dict
from collections import defaultdict

from chimeric_detector.config.config import Config
from chimeric_detector.core.bam_parser import BAMParser
from chimeric_detector.core.window_analysis import WindowAnalyzer, WindowMetrics
from chimeric_detector.modules.validation import ValidationPipeline
from chimeric_detector.utils.output import OutputFormatter, Logger
from chimeric_detector.utils.progress import ProgressTracker, ProgressReporter
from chimeric_detector.utils.security import (
    validate_bam_file, validate_config_file, validate_assembly_file,
    validate_output_path, sanitize_contig_name, validate_numeric_parameter,
    SecurityError, PathTraversalError, FileSizeError, FileTypeError
)
from chimeric_detector.utils.performance import (
    ParallelProcessor, WorkItem, PerformanceMonitor, performance_context,
    ResourceManager
)


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
    
    # Optional assembly input for taxonomy validation
    parser.add_argument("--assembly", type=str,
                       help="Assembly FASTA file (required for taxonomy validation)")
    
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
    validation_group.add_argument("--disable-topology", action="store_true",
                                help="Disable topology-aware validation (DTR/circular genome detection)")
    validation_group.add_argument("--topology-threshold", type=float, default=0.6,
                                help="Confidence threshold for circular genome detection")
    
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
    output_group.add_argument("--checkpoint", type=str,
                            help="Checkpoint file for resumable processing")
    
    # Performance options
    perf_group = parser.add_argument_group("Performance options")
    perf_group.add_argument("--parallel", action="store_true",
                          help="Enable parallel processing")
    perf_group.add_argument("--max-workers", type=int,
                          help="Maximum number of worker processes")
    perf_group.add_argument("--memory-limit", type=float,
                          help="Memory limit in GB")
    perf_group.add_argument("--profile", action="store_true",
                          help="Enable performance profiling")
    
    return parser.parse_args()


def detect_chimeras_in_contig_parallel(bam_path: str, contig: str, config_dict: dict) -> Dict:
    """Detect chimeras in a single contig - parallel processing version.
    
    Returns:
        Dict with 'breakpoints' and 'window_metrics' keys
    """
    from chimeric_detector.config.config import Config
    from chimeric_detector.core.bam_parser import BAMParser
    from chimeric_detector.core.window_analysis import WindowAnalyzer
    
    # Reconstruct objects from serializable data
    config = Config.from_dict(config_dict)
    analyzer = WindowAnalyzer(config.window, config.detection, config.quality)
    
    logger = logging.getLogger(__name__)
    logger.info(f"Processing contig: {contig}")
    
    try:
        with BAMParser(bam_path, config.quality, use_streaming=True) as bam_parser:
            breakpoints, window_metrics = _detect_chimeras_core(bam_parser, analyzer, contig, config)
            return {
                'breakpoints': breakpoints,
                'window_metrics': window_metrics
            }
    except Exception as e:
        logger.error(f"Error in parallel processing of {contig}: {e}")
        return {'breakpoints': [], 'window_metrics': []}


def detect_chimeras_in_contig(bam_parser: BAMParser, analyzer: WindowAnalyzer,
                            contig: str, config: Config) -> Tuple[List, List]:
    """Detect chimeras in a single contig with comprehensive error handling.
    
    Returns:
        Tuple of (breakpoints, window_metrics)
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Processing contig: {contig}")
    
    return _detect_chimeras_core(bam_parser, analyzer, contig, config)


def _detect_chimeras_core(bam_parser: BAMParser, analyzer: WindowAnalyzer,
                         contig: str, config: Config) -> Tuple[List, List]:
    """Core chimera detection logic shared between serial and parallel versions.
    
    Returns:
        Tuple of (breakpoints, window_metrics)
    """
    logger = logging.getLogger(__name__)
    
    try:
        # Estimate insert size distribution
        try:
            insert_stats = bam_parser.estimate_insert_size_distribution(contig)
            logger.debug(f"Insert size stats for {contig}: {insert_stats}")
        except Exception as e:
            logger.error(f"Failed to estimate insert size for {contig}: {e}")
            # Use default stats to allow processing to continue
            insert_stats = {"median": 500, "mad": 100, "mean": 500, "std": 100}
        
        # Collect window metrics
        window_metrics = []
        windows_processed = 0
        windows_failed = 0
        
        try:
            for start, end, pairs in bam_parser.iterate_windows(
                contig, config.window.size, config.window.step
            ):
                windows_processed += 1
                try:
                    coverage = bam_parser.calculate_coverage_in_window(contig, start, end)
                    metrics = analyzer.calculate_window_metrics(pairs, coverage)
                    
                    if metrics:
                        metrics.start = start
                        metrics.end = end
                        window_metrics.append(metrics)
                        
                except Exception as e:
                    windows_failed += 1
                    logger.debug(f"Failed to process window {start}-{end} in {contig}: {e}")
                    continue
                    
        except Exception as e:
            logger.error(f"Fatal error during window iteration for {contig}: {e}")
            return [], []
        
        if windows_failed > 0:
            logger.warning(f"Failed to process {windows_failed}/{windows_processed} windows in {contig}")
            
        if not window_metrics:
            logger.warning(f"No valid window metrics collected for {contig}")
            return [], []
            
        logger.info(f"Collected metrics for {len(window_metrics)} windows in {contig}")
        
        # Detect breakpoints
        try:
            breakpoints = analyzer.detect_breakpoints(window_metrics, contig, insert_stats)
            logger.info(f"Found {len(breakpoints)} potential breakpoints in {contig}")
            return breakpoints, window_metrics
            
        except Exception as e:
            logger.error(f"Failed to detect breakpoints in {contig}: {e}")
            return [], window_metrics  # Still return window metrics even if breakpoint detection fails
            
    except Exception as e:
        logger.error(f"Unexpected error processing contig {contig}: {e}")
        return [], []


def main():
    """Main entry point."""
    args = parse_arguments()
    
    try:
        # Validate all input paths for security
        validated_bam = validate_bam_file(args.bam_file)
        validated_output = validate_output_path(args.output)
        
        # Validate optional inputs
        validated_config = None
        if args.config:
            validated_config = validate_config_file(args.config)
            
        validated_assembly = None
        if args.assembly:
            validated_assembly = validate_assembly_file(args.assembly)
            
        # Validate and sanitize contig names if provided
        validated_contigs = None
        if args.contigs:
            validated_contigs = [sanitize_contig_name(contig) for contig in args.contigs]
            
        # Validate numeric parameters if provided
        if args.window_size:
            validate_numeric_parameter(args.window_size, 100, 100000, "window_size")
        if args.window_step:
            validate_numeric_parameter(args.window_step, 10, 10000, "window_step")
        if args.min_coverage:
            validate_numeric_parameter(args.min_coverage, 0.1, 1000.0, "min_coverage")
        if args.min_confidence:
            validate_numeric_parameter(args.min_confidence, 0.0, 1.0, "min_confidence")
        if args.max_workers:
            validate_numeric_parameter(args.max_workers, 1, 32, "max_workers")
        if args.memory_limit:
            validate_numeric_parameter(args.memory_limit, 0.5, 1000.0, "memory_limit")
            
    except (SecurityError, PathTraversalError, FileSizeError, FileTypeError) as e:
        print(f"Security error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Input validation error: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Set up logging
    Logger.setup_logging(
        verbose=args.verbose,
        debug=args.debug,
        log_file=Path(args.log_file) if args.log_file else None
    )
    
    logger = logging.getLogger(__name__)
    logger.info("Starting chimeric contig detection")
    
    # Load configuration
    if validated_config:
        config = Config.from_file(validated_config)
        logger.info(f"Loaded configuration from {validated_config}")
    else:
        config = Config()
        
    # Update config with command line arguments
    config.update_from_args(args)
    
    # Initialize components
    analyzer = WindowAnalyzer(config.window, config.detection, config.quality)
    formatter = OutputFormatter(config.output)
    
    # Set up resource management
    resource_manager = ResourceManager(args.memory_limit)
    
    all_breakpoints = []
    all_window_metrics = defaultdict(list)  # Store window metrics by contig
    
    # Process BAM file
    try:
        with BAMParser(str(validated_bam), config.quality, use_streaming=True) as bam_parser:
            # Get contigs to process
            if validated_contigs:
                contigs = validated_contigs
            else:
                contigs = bam_parser.get_contigs()
                
            logger.info(f"Processing {len(contigs)} contigs")
            
            # Set up progress tracking
            checkpoint_file = Path(args.checkpoint) if args.checkpoint else None
            progress = ProgressTracker(len(contigs), checkpoint_file)
            
            # Resume from checkpoint if available
            remaining_contigs = progress.get_remaining_contigs(contigs)
            if len(remaining_contigs) < len(contigs):
                logger.info(f"Resuming from checkpoint: {len(remaining_contigs)} contigs remaining")
                # Load existing results
                existing_results = progress.get_results()
                for contig_results in existing_results.values():
                    if isinstance(contig_results, list):
                        all_breakpoints.extend(contig_results)
            
            # Choose processing method based on arguments and data size
            use_parallel = args.parallel and len(remaining_contigs) > 1
            
            if use_parallel:
                logger.info("Using parallel processing")
                
                # Set up parallel processor
                parallel_processor = ParallelProcessor(config, args.max_workers)
                
                # Create work items
                work_items = [
                    WorkItem(
                        id=f"contig_{contig}",
                        contig=contig,
                        args=(str(validated_bam), contig, config.to_dict())
                    )
                    for contig in remaining_contigs
                ]
                
                # Process in parallel with performance monitoring
                with performance_context(f"Parallel processing {len(work_items)} contigs") as monitor:
                    results = parallel_processor.process_contigs_parallel(
                        work_items, detect_chimeras_in_contig_parallel, progress
                    )
                    
                    # Collect results
                    for contig, result in results.items():
                        if isinstance(result, dict) and 'breakpoints' in result:
                            all_breakpoints.extend(result['breakpoints'])
                            all_window_metrics[contig] = result['window_metrics']
                            monitor.record_item_processed()
                        elif isinstance(result, list):  # Legacy format
                            all_breakpoints.extend(result)
                            monitor.record_item_processed()
                        else:
                            monitor.record_error()
                            
            else:
                logger.info("Using sequential processing")
                
                # Process each remaining contig sequentially
                with performance_context(f"Sequential processing {len(remaining_contigs)} contigs") as monitor:
                    for contig in remaining_contigs:
                        progress.start_contig(contig)
                        
                        try:
                            # Check memory usage
                            if resource_manager.suggest_gc():
                                import gc
                                gc.collect()
                                
                            breakpoints, window_metrics = detect_chimeras_in_contig(
                                bam_parser, analyzer, contig, config
                            )
                            all_breakpoints.extend(breakpoints)
                            all_window_metrics[contig] = window_metrics
                            progress.complete_contig(contig, breakpoints)
                            monitor.record_item_processed()
                            
                        except Exception as e:
                            error_msg = f"Error processing contig {contig}: {e}"
                            progress.fail_contig(contig, error_msg)
                            monitor.record_error()
                            
                            if config.output.debug:
                                raise
            
            # Clean up checkpoint on successful completion
            if progress.is_complete():
                progress.cleanup_checkpoint()
                
            # Log performance summary
            if args.profile:
                memory_stats = resource_manager.get_memory_stats()
                logger.info(f"Final memory usage: {memory_stats['process_rss_gb']:.1f}GB")
                        
    except Exception as e:
        logger.error(f"Error processing BAM file: {e}")
        sys.exit(1)
        
    logger.info(f"Total breakpoints found: {len(all_breakpoints)}")
    
    # Run validation if enabled
    validation_results = None
    if config.validation.enabled:
        logger.info("Running validation modules")
        
        # Check if topology validation is enabled
        enable_topology = not args.disable_topology
        pipeline = ValidationPipeline(config.validation, enable_topology=enable_topology)
        
        # Load assembly for taxonomy validation if provided
        if validated_assembly and config.validation.taxonomy_db:
            for module in pipeline.modules:
                if hasattr(module, 'load_assembly'):
                    module.load_assembly(str(validated_assembly))
        
        # Get contig lengths for split evaluation
        contig_lengths = {}
        with BAMParser(str(validated_bam), config.quality) as validation_bam_parser:
            for contig in validation_bam_parser.get_contigs():
                contig_lengths[contig] = validation_bam_parser._bam.get_reference_length(contig)
                
            validation_results = pipeline.run_validation(
                all_breakpoints,
                window_metrics_by_contig=dict(all_window_metrics) if enable_topology else None,
                bam_parser=validation_bam_parser if enable_topology else None,
                contig_lengths=contig_lengths
            )
        
    # Write results
    metadata = {
        "bam_file": str(validated_bam),
        "total_contigs": len(contigs),
        "config": config.to_dict()
    }
    
    if validation_results:
        metadata["validation"] = validation_results
        
    formatter.write_results(all_breakpoints, validated_output, metadata)
    logger.info(f"Results written to {validated_output}")
    

if __name__ == "__main__":
    main()