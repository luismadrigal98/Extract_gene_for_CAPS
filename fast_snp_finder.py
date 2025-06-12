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
                            min_reliability='high', min_spacing=1000, require_complete_f2=True):
    """
    Filter ancestry results to find high-quality diagnostic SNPs with complete F2 information.
    
    Args:
        ancestry_file: Path to ancestry inference results
        output_file: Path to write filtered SNPs
        target_parent: Target parental line (e.g., '664c')
        min_reliability: Minimum reliability level ('low', 'medium', 'high')
        min_spacing: Minimum distance between selected SNPs (bp)
        require_complete_f2: If True, only keep variants with complete F2 data
    """
    import pandas as pd
    
    logging.info(f"Filtering high-quality SNPs from {ancestry_file}")
    logging.info(f"Requiring complete F2 information: {require_complete_f2}")
    
    # Read ancestry results
    df = pd.read_csv(ancestry_file, sep='\t')
    logging.info(f"Loaded {len(df)} variants from ancestry inference")
    
    # Reliability hierarchy
    reliability_values = {'low': 1, 'medium': 2, 'high': 3}
    min_rel_value = reliability_values[min_reliability]
    
    # Start with basic quality filters
    filtered = df[
        (df['complete_info'] == True) &
        (df['overall_reliability'].map(reliability_values) >= min_rel_value) &
        (df['has_f2_data'] == True)
    ].copy()
    
    logging.info(f"After basic quality filtering: {len(filtered)} variants")
    
    # Additional filter for complete F2 evidence if requested
    if require_complete_f2:
        # Look for variants with high-quality evidence from multiple sources
        # Priority: direct_parental > f2_inference > haplotype_inference > low_depth
        
        # Find parental source columns
        source_cols = [col for col in filtered.columns if col.endswith('_source')]
        data_status_cols = [col for col in filtered.columns if col.endswith('_data_status')]
        reliability_cols = [col for col in filtered.columns if col.endswith('_reliability')]
        context_cols = [col for col in filtered.columns if col.endswith('_context_agreement')]
        
        def has_complete_evidence(row):
            """Check if variant has complete F2 information from all evidence sources"""
            # Get all parent names
            parents = list(set([col.replace('_source', '').replace('_data_status', '').replace('_reliability', '').replace('_allele', '') 
                              for col in filtered.columns if any(suffix in col for suffix in ['_source', '_data_status', '_reliability', '_allele'])]))
            
            if not parents:
                return False
                
            # Count evidence quality for each parent
            high_quality_parents = 0
            total_parents = len(parents)
            
            for parent in parents:
                parent_source = row.get(f"{parent}_source", "")
                parent_status = row.get(f"{parent}_data_status", "")
                parent_reliability = row.get(f"{parent}_reliability", "")
                parent_context = row.get(f"{parent}_context_agreement", 1.0)  # Default to good agreement
                
                # Score this parent's evidence quality
                quality_score = 0
                
                # Source quality (highest priority)
                if parent_source == 'direct_parental':
                    quality_score += 3
                elif parent_source in ['f2_inference', 'f2_haplotype_inference']:
                    quality_score += 2
                elif parent_source == 'haplotype_block':
                    quality_score += 1
                
                # Data status bonus
                if parent_status == 'direct':
                    quality_score += 2
                elif parent_status == 'inferred':
                    quality_score += 1
                
                # Reliability bonus
                if parent_reliability == 'high':
                    quality_score += 2
                elif parent_reliability == 'medium':
                    quality_score += 1
                
                # Context agreement bonus (if available)
                if isinstance(parent_context, (int, float)) and parent_context >= 0.8:
                    quality_score += 1
                elif isinstance(parent_context, (int, float)) and parent_context >= 0.6:
                    quality_score += 0.5
                
                # Consider this parent as "high quality" if score >= 4 (combination of good source + status + reliability)
                if quality_score >= 4:
                    high_quality_parents += 1
            
            # Require that most parents (at least 75%) have high-quality evidence
            return high_quality_parents >= (total_parents * 0.75)
        
        if source_cols:  # Only apply if we have source information
            complete_f2_filtered = filtered[filtered.apply(has_complete_evidence, axis=1)].copy()
            logging.info(f"After complete F2 evidence filter: {len(complete_f2_filtered)} variants")
            
            if len(complete_f2_filtered) > 0:
                filtered = complete_f2_filtered
            else:
                logging.warning("Complete F2 filter removed all variants - using less stringent filter")
                # Fallback to simpler filter if the strict one removes everything
                simple_f2_filtered = filtered[
                    (filtered['complete_info'] == True) &
                    (filtered['has_f2_data'] == True)
                ].copy()
                logging.info(f"Using fallback filter: {len(simple_f2_filtered)} variants")
                filtered = simple_f2_filtered
        else:
            logging.warning("No source columns found - skipping F2 completeness filter")
    
    if len(filtered) == 0:
        logging.warning("No variants passed quality filters!")
        return 0
    
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
        # Add comprehensive quality score based on multiple factors
        def calculate_quality_score(row):
            """Calculate comprehensive quality score for variant selection"""
            score = 0
            
            # Base reliability score (3 points max)
            score += reliability_values[row['overall_reliability']]
            
            # Complete info bonus (2 points)
            if row['complete_info']:
                score += 2
            
            # F2 data availability bonus (1 point)  
            if row['has_f2_data']:
                score += 1
            
            # Context agreement bonuses (scan for context columns)
            context_bonus = 0
            context_count = 0
            for col in row.index:
                if col.endswith('_context_agreement') and pd.notna(row[col]):
                    try:
                        context_val = float(row[col])
                        if context_val >= 0.8:
                            context_bonus += 1
                        elif context_val >= 0.6:
                            context_bonus += 0.5
                        context_count += 1
                    except (ValueError, TypeError):
                        pass
            
            # Average context bonus (up to 2 points)
            if context_count > 0:
                score += min(2, context_bonus / context_count * 2)
            
            # Evidence source diversity bonus (up to 2 points)
            source_diversity = 0
            source_types = set()
            for col in row.index:
                if col.endswith('_source') and pd.notna(row[col]):
                    source_val = str(row[col])
                    if source_val in ['direct_parental', 'f2_inference', 'f2_haplotype_inference', 'haplotype_block']:
                        source_types.add(source_val)
            
            # Bonus for having multiple types of evidence
            if 'direct_parental' in source_types:
                source_diversity += 1
            if any(s in source_types for s in ['f2_inference', 'f2_haplotype_inference']):
                source_diversity += 1
            if len(source_types) >= 2:
                source_diversity += 0.5
            
            score += min(2, source_diversity)
            
            return round(score, 2)
        
        selected_df['quality_score'] = selected_df.apply(calculate_quality_score, axis=1)
        
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
        logging.info(f"Min spacing applied: {min_spacing}bp")
        
        # Log quality score distribution
        if 'quality_score' in selected_df.columns:
            logging.info(f"Quality score range: {selected_df['quality_score'].min():.1f} - {selected_df['quality_score'].max():.1f}")
            logging.info(f"Average quality score: {selected_df['quality_score'].mean():.1f}")
        
        # Log reliability distribution
        logging.info(f"Reliability distribution:")
        for rel in selected_df['overall_reliability'].value_counts().items():
            logging.info(f"  {rel[0]}: {rel[1]}")
        
        # Log evidence source summary
        logging.info(f"Evidence source summary:")
        source_summary = {}
        for col in selected_df.columns:
            if col.endswith('_source'):
                parent = col.replace('_source', '')
                sources = selected_df[col].value_counts()
                source_summary[parent] = sources.to_dict()
        
        for parent, sources in source_summary.items():
            logging.info(f"  {parent}: {sources}")
        
        # Show top candidates with their details
        logging.info(f"Top 5 diagnostic SNPs:")
        top_5 = selected_df.head(5)
        for idx, row in top_5.iterrows():
            chrom, pos = row['CHROM'], row['POS']
            quality = row.get('quality_score', 'N/A')
            reliability = row['overall_reliability']
            target_allele = row.get(f'{target_parent}_allele', 'N/A')
            logging.info(f"  {chrom}:{pos} (Q:{quality}, R:{reliability}, {target_parent}:{target_allele})")
    
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
    parser.add_argument('--min_depth', type=int, default=3,
                       help='Minimum read depth to consider a call reliable')
    parser.add_argument('--max_depth', type=int, default=200,
                       help='Maximum read depth to consider a call reliable (filters out high-coverage artifacts)')
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
            min_depth=args.min_depth,
            max_depth=args.max_depth
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
