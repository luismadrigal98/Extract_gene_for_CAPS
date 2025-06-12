#!/usr/bin/env python3
"""
Fast Diagnostic SNP Finder - Streamlined for SNP Discovery Only

This script identifies high-quality genic SNPs that can discriminate between parental lines
without the expensive primer design step. Optimized for speed.

@author: Luis Javier Madrigal-Roca & John K. Kelly
@date: 2025/06/12
"""

import argparse
import sys
import os
import logging

# Add src to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from masking_vcf import mask_variants
from ancestry_inference import infer_ancestry_single

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("fast_snp_finder.log"),
        logging.StreamHandler(sys.stdout)
    ]
)

def filter_high_quality_snps(ancestry_file, output_file, target_parent='664c', 
                            min_reliability='medium', min_spacing=1000):
    """
    Filter ancestry results to find high-quality diagnostic SNPs with good spacing.
    
    Args:
        ancestry_file: Path to ancestry inference results
        output_file: Path to write filtered SNPs
        target_parent: Target parental line (e.g., '664c')
        min_reliability: Minimum reliability level ('low', 'medium', 'high')
        min_spacing: Minimum distance between selected SNPs (bp)
    """
    import pandas as pd
    
    logging.info(f"Filtering high-quality SNPs from {ancestry_file}")
    
    # Read ancestry results
    df = pd.read_csv(ancestry_file, sep='\t')
    logging.info(f"Loaded {len(df)} variants from ancestry inference")
    
    # Reliability hierarchy
    reliability_values = {'low': 1, 'medium': 2, 'high': 3}
    min_rel_value = reliability_values[min_reliability]
    
    # Filter for complete information and high reliability
    filtered = df[
        (df['complete_info'] == True) &
        (df['overall_reliability'].map(reliability_values) >= min_rel_value) &
        (df['has_f2_data'] == True)
    ].copy()
    
    logging.info(f"After quality filtering: {len(filtered)} variants")
    
    if len(filtered) == 0:
        logging.warning("No variants passed quality filters!")
        return
    
    # Find diagnostic SNPs (target parent differs from others)
    target_col = f"{target_parent}_allele"
    if target_col not in filtered.columns:
        logging.error(f"Target parent column '{target_col}' not found!")
        return
    
    # Get all parental allele columns
    allele_cols = [col for col in filtered.columns if col.endswith('_allele')]
    other_cols = [col for col in allele_cols if col != target_col]
    
    # Find variants where target parent differs from all others
    def is_diagnostic(row):
        target_allele = row[target_col]
        if target_allele == 'N':
            return False
        
        for other_col in other_cols:
            other_allele = row[other_col]
            if other_allele == 'N' or other_allele == target_allele:
                return False
        return True
    
    diagnostic = filtered[filtered.apply(is_diagnostic, axis=1)].copy()
    logging.info(f"Found {len(diagnostic)} diagnostic variants for {target_parent}")
    
    if len(diagnostic) == 0:
        logging.warning("No diagnostic variants found!")
        return
    
    # Sort by position for spacing calculation
    diagnostic = diagnostic.sort_values(['CHROM', 'POS'])
    
    # Apply spacing filter
    selected = []
    last_pos = -min_spacing
    last_chrom = None
    
    for _, variant in diagnostic.iterrows():
        chrom = variant['CHROM']
        pos = variant['POS']
        
        # If different chromosome or sufficient spacing
        if chrom != last_chrom or (pos - last_pos) >= min_spacing:
            selected.append(variant)
            last_pos = pos
            last_chrom = chrom
    
    selected_df = pd.DataFrame(selected)
    logging.info(f"After spacing filter ({min_spacing}bp): {len(selected_df)} variants")
    
    # Sort by quality metrics for final selection
    if len(selected_df) > 0:
        # Add quality score based on reliability and completeness
        selected_df['quality_score'] = selected_df.apply(
            lambda x: reliability_values[x['overall_reliability']] + 
                     (1 if x['complete_info'] else 0), axis=1
        )
        
        # Sort by quality score descending, then by position
        selected_df = selected_df.sort_values(['quality_score', 'CHROM', 'POS'], 
                                            ascending=[False, True, True])
        
        # Write results
        selected_df.to_csv(output_file, sep='\t', index=False)
        logging.info(f"Wrote {len(selected_df)} high-quality diagnostic SNPs to {output_file}")
        
        # Log summary
        logging.info("=== DIAGNOSTIC SNP SUMMARY ===")
        logging.info(f"Target parent: {target_parent}")
        logging.info(f"Total processed: {len(df)}")
        logging.info(f"Quality filtered: {len(filtered)}")
        logging.info(f"Diagnostic: {len(diagnostic)}")
        logging.info(f"Final selected: {len(selected_df)}")
        logging.info(f"Reliability distribution:")
        for rel in selected_df['overall_reliability'].value_counts().items():
            logging.info(f"  {rel[0]}: {rel[1]}")
    
    return len(selected_df) if len(selected_df) > 0 else 0

def main():
    parser = argparse.ArgumentParser(
        description="Fast diagnostic SNP finder - optimized for SNP discovery without primer design"
    )
    
    # Input files
    parser.add_argument('--vcf', required=True, 
                       help='Input VCF file with F2 sequencing data')
    parser.add_argument('--gff3', required=True,
                       help='Gene annotation file (GFF3 format)')
    parser.add_argument('--roi_list', required=True,
                       help='ROI coordinates file (tab-separated with headers: ROI/ROI_name, Chr/Chrom, Start, End)')
    parser.add_argument('--ancestry_map', required=True,
                       help='File mapping F2 samples to parental lines')
    
    # Output files
    parser.add_argument('--output_prefix', required=True,
                       help='Prefix for output files')
    
    # Analysis parameters
    parser.add_argument('--target_parent', default='664c',
                       help='Target parental line to discriminate (default: 664c)')
    parser.add_argument('--min_qual', type=float, default=60,
                       help='Minimum variant quality score')
    parser.add_argument('--min_reliability', choices=['low', 'medium', 'high'], 
                       default='medium', help='Minimum reliability level for SNPs')
    parser.add_argument('--min_spacing', type=int, default=1000,
                       help='Minimum distance between selected SNPs (bp)')
    parser.add_argument('--max_snps', type=int, default=50,
                       help='Maximum number of SNPs to output')
    
    # Speed optimizations
    parser.add_argument('--skip_primer_check', action='store_true', default=True,
                       help='Skip primer design validation (faster)')
    
    args = parser.parse_args()
    
    # Define output files
    filtered_vcf = f"{args.output_prefix}_filtered_genic.vcf"
    ancestry_results = f"{args.output_prefix}_ancestry_inferred.tsv"
    diagnostic_snps = f"{args.output_prefix}_diagnostic_snps.tsv"
    
    logging.info("=== Fast Diagnostic SNP Finder ===")
    logging.info(f"Target parental line: {args.target_parent}")
    logging.info(f"Output prefix: {args.output_prefix}")
    logging.info(f"Min reliability: {args.min_reliability}")
    logging.info(f"Min spacing: {args.min_spacing}bp")
    
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
        
        # Step 3: Filter for high-quality diagnostic SNPs
        logging.info("Step 3: Filtering for diagnostic SNPs...")
        num_snps = filter_high_quality_snps(
            ancestry_file=ancestry_results,
            output_file=diagnostic_snps,
            target_parent=args.target_parent,
            min_reliability=args.min_reliability,
            min_spacing=args.min_spacing
        )
        
        # Take top SNPs if we have too many
        if num_snps and num_snps > args.max_snps:
            logging.info(f"Selecting top {args.max_snps} SNPs from {num_snps} candidates...")
            import pandas as pd
            df = pd.read_csv(diagnostic_snps, sep='\t')
            top_snps = df.head(args.max_snps)
            final_output = f"{args.output_prefix}_top_{args.max_snps}_diagnostic_snps.tsv"
            top_snps.to_csv(final_output, sep='\t', index=False)
            logging.info(f"Final selection written to: {final_output}")
        
        logging.info("=== Fast SNP Discovery Completed Successfully ===")
        logging.info(f"Found {num_snps if num_snps else 0} diagnostic SNPs")
        logging.info(f"Results written to: {diagnostic_snps}")
        
        if not args.skip_primer_check and num_snps and num_snps > 0:
            logging.info("\nTo design primers for these SNPs, run:")
            logging.info(f"python diagnostic_snp_finder.py --vcf {args.vcf} --roi_list {args.roi_list} "
                        f"--ancestry_map {args.ancestry_map} --output_prefix primer_design_run "
                        f"--target_parent {args.target_parent} [other options]")
        
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        import traceback
        logging.error(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()
