#!/usr/bin/env python3
"""
Ultra-Fast Diagnostic SNP Finder - Highly Optimized Version

This version uses vectorized operations and parallel processing for maximum speed.
Expected to run in 15-30 minutes instead of hours.

@author: Luis Javier Madrigal-Roca & John K. Kelly  
@date: 2025/06/12
"""

import argparse
import sys
import os
import logging
import pandas as pd
import numpy as np
from tqdm import tqdm
import multiprocessing as mp

# Add src to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from masking_vcf import mask_variants
from ancestry_inference import infer_ancestry_single
from fast_screen_variants import (fast_screen_variants_parallel, 
                                 fast_filter_diagnostic_variants,
                                 fast_apply_spacing_filter)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("ultra_fast_snp_finder.log"),
        logging.StreamHandler(sys.stdout)
    ]
)

def ultra_fast_filter_snps(ancestry_file, output_file, target_parent='664c', 
                          min_reliability='medium', min_spacing=2000, max_snps=50,
                          n_workers=None):
    """
    Ultra-fast filtering using vectorized operations and parallel processing
    """
    logging.info(f"Ultra-fast filtering from {ancestry_file}")
    
    # Read data
    df = pd.read_csv(ancestry_file, sep='\t')
    logging.info(f"Loaded {len(df)} variants")
    
    # Quality filtering using vectorized operations
    reliability_values = {'low': 1, 'medium': 2, 'high': 3}
    min_rel_value = reliability_values[min_reliability]
    
    # Fast vectorized filtering
    quality_mask = (
        (df['complete_info'] == True) &
        (df['overall_reliability'].map(reliability_values) >= min_rel_value) &
        (df['has_f2_data'] == True)
    )
    
    filtered = df[quality_mask].copy()
    logging.info(f"After quality filtering: {len(filtered)} variants")
    
    if len(filtered) == 0:
        logging.warning("No variants passed quality filters!")
        return 0
    
    # Fast diagnostic filtering  
    diagnostic = fast_filter_diagnostic_variants(filtered, target_parent)
    logging.info(f"Found {len(diagnostic)} diagnostic variants")
    
    if len(diagnostic) == 0:
        logging.warning("No diagnostic variants found!")
        return 0
    
    # Fast parallel primer screening (optional - can be skipped for speed)
    logging.info("Performing fast primer compliance screening...")
    screened = fast_screen_variants_parallel(
        diagnostic, 
        primer_size=20, 
        amplicon_size=300, 
        displacement_steps=3,  # Reduced for speed
        n_workers=n_workers
    )
    
    # Keep only primer compliant variants
    primer_compliant = screened[screened['primer_compliant'] == True].copy()
    logging.info(f"Primer compliant variants: {len(primer_compliant)}")
    
    if len(primer_compliant) == 0:
        logging.warning("No primer compliant variants found!")
        # Fallback to diagnostic variants without primer checking
        primer_compliant = diagnostic.copy()
        logging.info(f"Using {len(primer_compliant)} diagnostic variants without primer screening")
    
    # Fast spacing filter
    spaced = fast_apply_spacing_filter(primer_compliant, min_spacing)
    logging.info(f"After spacing filter ({min_spacing}bp): {len(spaced)} variants")
    
    # Quality scoring and selection
    if len(spaced) > 0:
        # Simple quality score based on reliability and completeness
        def fast_quality_score(row):
            score = reliability_values[row['overall_reliability']]
            if row['complete_info']:
                score += 2
            if row['has_f2_data']:
                score += 1
            # Add QUAL score component
            if pd.notna(row.get('QUAL', 0)):
                score += min(2, row['QUAL'] / 100)  # Normalize QUAL
            return score
        
        spaced['quality_score'] = spaced.apply(fast_quality_score, axis=1)
        
        # Sort by quality and take top candidates
        final = spaced.sort_values(['quality_score', 'QUAL'], ascending=[False, False])
        
        if len(final) > max_snps:
            final = final.head(max_snps)
            logging.info(f"Selected top {max_snps} SNPs")
        
        # Write results
        final.to_csv(output_file, sep='\t', index=False)
        logging.info(f"Wrote {len(final)} ultra-fast diagnostic SNPs to {output_file}")
        
        # Summary statistics
        logging.info("=== ULTRA-FAST SNP SUMMARY ===")
        logging.info(f"Quality score range: {final['quality_score'].min():.1f} - {final['quality_score'].max():.1f}")
        logging.info(f"Reliability distribution:")
        for rel, count in final['overall_reliability'].value_counts().items():
            logging.info(f"  {rel}: {count}")
            
        # Show top candidates
        logging.info(f"Top 5 diagnostic SNPs:")
        for idx, row in final.head(5).iterrows():
            chrom, pos = row['CHROM'], row['POS']
            quality = row['quality_score']
            reliability = row['overall_reliability']
            target_allele = row.get(f'{target_parent}_allele', 'N/A')
            logging.info(f"  {chrom}:{pos} (Q:{quality:.1f}, R:{reliability}, {target_parent}:{target_allele})")
    
    return len(final) if len(spaced) > 0 else 0

def main():
    parser = argparse.ArgumentParser(
        description="Ultra-fast diagnostic SNP finder - optimized for maximum speed"
    )
    
    # Input files
    parser.add_argument('--vcf', required=True, 
                       help='Input VCF file with F2 sequencing data')
    parser.add_argument('--gff3', required=True,
                       help='Gene annotation file (GFF3 format)')
    parser.add_argument('--roi_list', required=True,
                       help='ROI coordinates file')
    parser.add_argument('--ancestry_map', required=True,
                       help='File mapping F2 samples to parental lines')
    
    # Output files
    parser.add_argument('--output_prefix', required=True,
                       help='Prefix for output files')
    
    # Analysis parameters
    parser.add_argument('--target_parent', default='664c',
                       help='Target parental line (default: 664c)')
    parser.add_argument('--min_qual', type=float, default=60,
                       help='Minimum variant quality score')
    parser.add_argument('--min_reliability', choices=['low', 'medium', 'high'], 
                       default='medium', help='Minimum reliability level')
    parser.add_argument('--min_spacing', type=int, default=2000,
                       help='Minimum distance between SNPs (bp)')
    parser.add_argument('--max_snps', type=int, default=50,
                       help='Maximum number of SNPs to output')
    
    # Performance options
    parser.add_argument('--n_workers', type=int, 
                       help='Number of parallel workers (default: auto)')
    parser.add_argument('--skip_primer_screen', action='store_true',
                       help='Skip primer screening for maximum speed')
    
    args = parser.parse_args()
    
    # Set up parallel workers
    if args.n_workers is None:
        args.n_workers = min(mp.cpu_count() - 1, 4)
    
    # Define output files
    filtered_vcf = f"{args.output_prefix}_filtered_genic.vcf"
    ancestry_results = f"{args.output_prefix}_ancestry_inferred.tsv"
    diagnostic_snps = f"{args.output_prefix}_ultra_fast_diagnostic_snps.tsv"
    
    logging.info("=== ULTRA-FAST Diagnostic SNP Finder ===")
    logging.info(f"Target parental line: {args.target_parent}")
    logging.info(f"Parallel workers: {args.n_workers}")
    logging.info(f"Min spacing: {args.min_spacing}bp")
    logging.info(f"Skip primer screening: {args.skip_primer_screen}")
    
    try:
        # Step 1: Filter VCF for genic regions
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
        
        # Step 3: Ultra-fast SNP filtering
        logging.info("Step 3: Ultra-fast diagnostic SNP filtering...")
        if args.skip_primer_screen:
            # Skip primer screening completely for maximum speed
            logging.info("Skipping primer screening for maximum speed...")
            # Use simplified version without primer compliance
            from fast_screen_variants import fast_filter_diagnostic_variants, fast_apply_spacing_filter
            
            df = pd.read_csv(ancestry_results, sep='\t')
            reliability_values = {'low': 1, 'medium': 2, 'high': 3}
            min_rel_value = reliability_values[args.min_reliability]
            
            # Fast quality filtering
            filtered = df[
                (df['complete_info'] == True) &
                (df['overall_reliability'].map(reliability_values) >= min_rel_value) &
                (df['has_f2_data'] == True)
            ].copy()
            
            # Fast diagnostic filtering
            diagnostic = fast_filter_diagnostic_variants(filtered, args.target_parent)
            
            # Fast spacing filter
            spaced = fast_apply_spacing_filter(diagnostic, args.min_spacing)
            
            if len(spaced) > args.max_snps:
                spaced = spaced.head(args.max_snps)
            
            spaced.to_csv(diagnostic_snps, sep='\t', index=False)
            num_snps = len(spaced)
        else:
            # Use full screening with primer compliance
            num_snps = ultra_fast_filter_snps(
                ancestry_file=ancestry_results,
                output_file=diagnostic_snps,
                target_parent=args.target_parent,
                min_reliability=args.min_reliability,
                min_spacing=args.min_spacing,
                max_snps=args.max_snps,
                n_workers=args.n_workers
            )
        
        logging.info("=== ULTRA-FAST SNP Discovery Completed ===")
        logging.info(f"Found {num_snps} diagnostic SNPs in record time!")
        logging.info(f"Results written to: {diagnostic_snps}")
        
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        import traceback
        logging.error(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()
