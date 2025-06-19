#!/usr/bin/env python3
"""
Standalone breakpoint investigation tool.

This tool provides detailed analysis of breakpoints detected by the chimera detector,
including sequence analysis, coverage profiles, read alignment patterns, and more.
"""

import argparse
import sys
import json
import logging
from pathlib import Path
from typing import List

from chimeric_detector.config.config import Config
from chimeric_detector.core.window_analysis import Breakpoint
from chimeric_detector.modules.investigation import BreakpointInvestigator
from chimeric_detector.utils.output import Logger
from chimeric_detector.utils.security import (
    validate_bam_file, validate_assembly_file, validate_output_path,
    SecurityError, PathTraversalError, FileSizeError, FileTypeError
)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Investigate detected breakpoints in detail",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("results_file", type=str,
                       help="JSON results file from chimera detection")
    parser.add_argument("bam_file", type=str,
                       help="BAM file used for detection (sorted and indexed)")
    parser.add_argument("output_dir", type=str,
                       help="Output directory for investigation results")
    
    # Optional arguments
    parser.add_argument("--assembly", type=str,
                       help="Assembly FASTA file for sequence analysis")
    parser.add_argument("--min-confidence", type=float, default=0.7,
                       help="Minimum confidence threshold for investigation")
    parser.add_argument("--flank-size", type=int, default=2000,
                       help="Size of flanking regions to analyze (bp)")
    parser.add_argument("--max-breakpoints", type=int, default=50,
                       help="Maximum number of breakpoints to investigate")
    
    # Analysis options
    parser.add_argument("--detailed", action="store_true",
                       help="Generate detailed individual breakpoint reports")
    parser.add_argument("--contigs", type=str, nargs="+",
                       help="Investigate only specific contigs")
    parser.add_argument("--positions", type=str, nargs="+",
                       help="Investigate specific positions (format: contig:position)")
    
    # Output options
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Enable verbose output")
    parser.add_argument("-d", "--debug", action="store_true",
                       help="Enable debug output")
    parser.add_argument("--log-file", type=str,
                       help="Write logs to file")
    
    return parser.parse_args()


def load_breakpoints_from_results(results_file: Path, min_confidence: float = 0.7,
                                 contigs_filter: List[str] = None,
                                 positions_filter: List[str] = None) -> List[Breakpoint]:
    """Load breakpoints from detection results file."""
    
    with open(results_file) as f:
        results = json.load(f)
    
    breakpoints = []
    
    for bp_data in results.get("breakpoints", []):
        # Apply confidence filter
        if bp_data["confidence"] < min_confidence:
            continue
        
        # Apply contig filter
        if contigs_filter and bp_data["contig"] not in contigs_filter:
            continue
        
        # Apply position filter
        if positions_filter:
            bp_key = f"{bp_data['contig']}:{bp_data['position']}"
            if bp_key not in positions_filter:
                continue
        
        # Create breakpoint object
        bp = Breakpoint(
            contig=bp_data["contig"],
            position=bp_data["position"],
            confidence=bp_data["confidence"],
            evidence=bp_data.get("evidence", {}),
            left_metrics=None,
            right_metrics=None
        )
        breakpoints.append(bp)
    
    return breakpoints


def main():
    """Main entry point."""
    args = parse_arguments()
    
    try:
        # Validate all input paths for security
        validated_results = Path(args.results_file)
        if not validated_results.exists():
            raise FileNotFoundError(f"Results file not found: {args.results_file}")
            
        validated_bam = validate_bam_file(args.bam_file)
        validated_output = validate_output_path(args.output_dir)
        
        # Validate optional assembly file
        validated_assembly = None
        if args.assembly:
            validated_assembly = validate_assembly_file(args.assembly)
            
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
    logger.info("Starting breakpoint investigation")
    
    try:
        # Load breakpoints from results
        breakpoints = load_breakpoints_from_results(
            validated_results,
            min_confidence=args.min_confidence,
            contigs_filter=args.contigs,
            positions_filter=args.positions
        )
        
        if not breakpoints:
            logger.warning("No breakpoints found matching criteria")
            print("No breakpoints found matching criteria")
            return
        
        # Limit number of breakpoints if requested
        if len(breakpoints) > args.max_breakpoints:
            logger.info(f"Limiting investigation to top {args.max_breakpoints} breakpoints by confidence")
            breakpoints = sorted(breakpoints, key=lambda x: x.confidence, reverse=True)[:args.max_breakpoints]
        
        logger.info(f"Investigating {len(breakpoints)} breakpoints")
        print(f"Investigating {len(breakpoints)} breakpoints...")
        
        # Create investigator
        config = Config()  # Use default config
        investigator = BreakpointInvestigator(
            config, 
            str(validated_bam), 
            str(validated_assembly) if validated_assembly else None
        )
        
        # Run investigation
        report = investigator.generate_investigation_report(
            breakpoints, 
            validated_output, 
            detailed=args.detailed
        )
        
        # Print summary
        summary = report["summary"]
        print(f"\nInvestigation Complete!")
        print(f"Output directory: {validated_output}")
        print(f"Total breakpoints analyzed: {summary['total_breakpoints']}")
        print(f"High confidence breakpoints: {summary['high_confidence_breakpoints']}")
        print(f"Contigs affected: {len(summary['contigs_affected'])}")
        
        if summary["evidence_summary"]:
            print(f"Evidence found:")
            for evidence_type, count in summary["evidence_summary"].items():
                print(f"  - {evidence_type}: {count}")
        
        # Recommend next steps
        high_conf_count = summary['high_confidence_breakpoints']
        if high_conf_count > 0:
            print(f"\n🎯 Recommended: Focus on {high_conf_count} high-confidence breakpoints")
            print(f"   Check: {validated_output}/investigation_report.html")
        else:
            print(f"\n💡 No high-confidence breakpoints found")
            print(f"   Consider lowering --min-confidence threshold")
        
        logger.info("Investigation completed successfully")
        
    except Exception as e:
        logger.error(f"Error during investigation: {e}")
        if args.debug:
            raise
        sys.exit(1)


if __name__ == "__main__":
    main()