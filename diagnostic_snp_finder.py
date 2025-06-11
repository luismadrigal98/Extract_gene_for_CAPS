#!/usr/bin/env python3
"""
Diagnostic SNP Finder for Parental Line Discrimination

This script identifies genic SNPs that can discriminate between parental lines,
specifically focusing on identifying variants diagnostic of 664c vs other parental lines.

Uses the original reference genome (where annotation is defined) for maximum robustness.

@author: Luis Javier Madrigal-Roca & John K. Kelly
@date: 2025/06/11
"""

import argparse
import sys
import os
import logging

# Add src to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from masking_vcf import mask_variants
from ancestry_inference import infer_ancestry_single
from screen_variants import screen_variants
from primer_design import design_primers

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("diagnostic_snp_finder.log"),
        logging.StreamHandler(sys.stdout)
    ]
)

def main():
    parser = argparse.ArgumentParser(
        description="Find diagnostic SNPs in genic regions for parental line discrimination"
    )
    
    # Input files
    parser.add_argument('--vcf', required=True, 
                       help='Input VCF file with F2 sequencing data')
    parser.add_argument('--reference', required=True,
                       help='Reference genome FASTA file (original, not remapped)')
    parser.add_argument('--gff3', required=True,
                       help='Gene annotation file (GFF3 format)')
    parser.add_argument('--roi_list', required=True,
                       help='ROI coordinates file (tab-separated with headers: ROI, Chr, Start, End)')
    parser.add_argument('--ancestry_map', required=True,
                       help='File mapping F2 samples to parental lines')
    
    # Output files
    parser.add_argument('--output_prefix', required=True,
                       help='Prefix for output files')
    
    # Analysis parameters
    parser.add_argument('--target_parent', default='664c',
                       help='Target parental line to discriminate (default: 664c)')
    parser.add_argument('--min_qual', type=float, default=20,
                       help='Minimum variant quality score')
    parser.add_argument('--max_variants', type=int, default=200,
                       help='Maximum number of variants for primer design')
    
    # Primer design parameters
    parser.add_argument('--primer3_settings', 
                       help='Primer3 settings file (optional)')
    parser.add_argument('--flanking_size', type=int, default=150,
                       help='Flanking sequence size for primer design')
    parser.add_argument('--parallel', action='store_true',
                       help='Use parallel processing for primer design')
    parser.add_argument('--num_workers', type=int,
                       help='Number of parallel workers (default: CPU count - 1)')
    
    # Optional features
    parser.add_argument('--keep_temp', action='store_true',
                       help='Keep temporary files for debugging')
    parser.add_argument('--temp_dir', 
                       help='Temporary directory (default: auto-generated)')
    
    args = parser.parse_args()
    
    # Define output files
    filtered_vcf = f"{args.output_prefix}_filtered_genic.vcf"
    ancestry_results = f"{args.output_prefix}_ancestry_inferred.tsv"
    screened_variants = f"{args.output_prefix}_screened_variants.tsv"
    primer_results = f"{args.output_prefix}_diagnostic_primers.txt"
    selected_primers = f"{args.output_prefix}_selected_primers.txt"
    
    logging.info("=== Diagnostic SNP Finder Pipeline ===")
    logging.info(f"Target parental line: {args.target_parent}")
    logging.info(f"Output prefix: {args.output_prefix}")
    
    try:
        # Step 1: Filter VCF for genic regions within ROI
        logging.info("Step 1: Filtering VCF for genic regions...")
        mask_variants(
            vcf=args.vcf,
            gff3=args.gff3,
            roi_list=args.roi_list,
            output=filtered_vcf,
            only_biallelic=True,
            min_quality=args.min_qual,
            filter_indels=True
        )
        
        # Step 2: Infer parental ancestry
        logging.info("Step 2: Inferring parental ancestry...")
        infer_ancestry_single(
            vcf=filtered_vcf,
            ROI_list=args.roi_list,
            ancestry_log=args.ancestry_map,
            output=ancestry_results,
            use_assembly_when_f2_missing=True,
            min_depth=3
        )
        
        # Step 3: Screen variants for primer compliance
        logging.info("Step 3: Screening variants for primer compliance...")
        screen_variants(
            tsv_list=[ancestry_results],
            output_dir=os.path.dirname(screened_variants) or ".",
            allele_col_pattern="allele",
            reliability_thr=0.8,
            diff_parent=args.target_parent,
            amplicon_size=args.flanking_size * 2,
            primer_size=20,
            displacement=True,
            displacement_steps=5
        )
        
        # Step 4: Design primers
        logging.info("Step 4: Designing diagnostic primers...")
        num_designed = design_primers(
            input_files=[screened_variants],
            reference_fasta=args.reference,
            output_file=primer_results,
            settings_file=args.primer3_settings,
            max_variants=args.max_variants,
            flanking_size=args.flanking_size,
            parallel=args.parallel,
            num_workers=args.num_workers,
            keep_temp=args.keep_temp,
            temp_dir=args.temp_dir,
            contrast=True,
            num_primers=50,
            selection_criteria="balanced",
            selected_output=selected_primers
        )
        
        logging.info("=== Pipeline Completed Successfully ===")
        logging.info(f"Designed primers for {num_designed} variants")
        logging.info(f"Results written to:")
        logging.info(f"  - All primers: {primer_results}")
        logging.info(f"  - Selected primers: {selected_primers}")
        
        # Optional: Add validation step
        validation_output = f"{args.output_prefix}_validation_report.md"
        logging.info(f"To validate primers, run:")
        logging.info(f"python -c \"from src.validation_utilities import validate_primers; "
                    f"validate_primers('{selected_primers}', ['genome1.fasta'], '{validation_output}')\"")
        
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
