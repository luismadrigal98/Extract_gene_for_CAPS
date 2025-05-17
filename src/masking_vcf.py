"""
This script masks a VCF file by preserving only the variants that are located in genic regions, as defined by a GFF3 file. It also filters the variants based 
on a list of regions of interest (ROI) provided by the user. The output is a VCF file containing only the relevant variants.

@author: Luis Javier Madrigal-Roca & John K. Kelly

@date 2025/05/16

"""

import pandas as pd
import logging
from tqdm import tqdm
from intervaltree import IntervalTree

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def read_vcf(path):
    """Read a VCF file into a pandas DataFrame."""
    logger.info(f"Reading VCF file: {path}")
    try:
        # First, peek at the file to identify header structure
        vcf_header_lines = []
        column_names = None
        
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('##'):
                    vcf_header_lines.append(line.strip())
                elif line.startswith('#CHROM'):
                    column_names = line.strip().split('\t')
                    # Replace #CHROM with CHROM for consistent access
                    column_names[0] = 'CHROM'  
                    break
        
        if not column_names:
            raise ValueError("Could not find column header line in VCF file")
        
        logger.info(f"Found {len(column_names)} columns in VCF header")
            
        # Now read the data portion
        vcf_df = pd.read_csv(
            path, 
            comment='#',
            sep='\t',
            names=column_names,
            dtype={'CHROM': str}
        )
        
        # Log column names for debugging
        logger.debug(f"DataFrame columns: {vcf_df.columns.tolist()}")
        
        # Store original headers for writing back
        vcf_df.attrs['headers'] = vcf_header_lines
        vcf_df.attrs['column_header'] = '\t'.join(column_names).replace('CHROM', '#CHROM') + '\n'
        
        logger.info(f"Successfully loaded {len(vcf_df)} variants from VCF file")
        return vcf_df
    
    except Exception as e:
        logger.error(f"Error reading VCF file: {e}", exc_info=True)
        raise

def mask_variants(vcf, gff3, roi_list, output):
    """
    Mask the variants in the VCF file based on the GFF3 file and ROI list.
    
    Parameters:
    vcf (str): Path to the input VCF file.
    gff3 (str): Path to the input GFF3 file.
    roi_list (str): Path to the ROI list file.
    output (str): Path to the output VCF file.
    
    Returns:
    None
    """
    logger.info("Starting variant masking process")
    
    try:
        # Read the VCF file
        vcf_df = read_vcf(vcf)
        initial_variant_count = len(vcf_df)
        
        # Read the GFF3 file
        logger.info(f"Reading GFF3 file: {gff3}")
        gff3_df = pd.read_csv(gff3, sep='\t', comment='#', header=None)
        logger.info(f"Loaded {len(gff3_df)} features from GFF3 file")

        # Filter the GFF3 file to keep only genic regions
        gff3_df = gff3_df[gff3_df[2].str.contains('gene')]
        logger.info(f"Filtered to {len(gff3_df)} gene features")
        
        # Read the ROI list
        logger.info(f"Reading ROI list: {roi_list}")
        try:
            roi_df = pd.read_csv(roi_list, sep='\t', header=0)
            logger.info(f"Loaded {len(roi_df)} regions of interest")
        except Exception as e:
            logger.warning(f"Error with tab delimiter, trying with flexible whitespace: {e}")
            roi_df = pd.read_csv(roi_list, delim_whitespace=True, header=0)
            logger.info(f"Successfully loaded {len(roi_df)} regions of interest with flexible whitespace")

        # Filter the VCF variants to keep only the chromosomes and positions in the ROI list
        logger.info("Filtering variants by regions of interest")
        roi_filtered_df = vcf_df[vcf_df.apply(
            lambda x: any((x['CHROM'] == r['Chrom']) and 
                            (r['Start'] <= x['POS'] <= r['End'])
                            for _, r in roi_df.iterrows()),
            axis=1)]
        
        roi_filtered_count = len(roi_filtered_df)
        logger.info(f"Filtered to {roi_filtered_count} variants within regions of interest " +
                    f"({roi_filtered_count/initial_variant_count:.2%} of original)")
        
        # Build interval trees per chromosome (much faster)
        logger.info("Building interval trees for gene regions")
        gene_trees = {}
        for _, row in tqdm(gff3_df.iterrows(), total=len(gff3_df), desc="Building gene interval trees"):
            chrom = row[0]
            if chrom not in gene_trees:
                gene_trees[chrom] = IntervalTree()
            # Create an Interval object or use the dictionary-style assignment
            gene_trees[chrom][row[3]:row[4]+1] = row[8]  # This is the correct way
        
        logger.info(f"Built interval trees for {len(gene_trees)} chromosomes")

        # Efficient lookup
        logger.info("Finding variants in genic regions")
        roi_filtered_df = roi_filtered_df.copy()  # Create an explicit copy
        roi_filtered_df['in_genic'] = roi_filtered_df.apply(
            lambda x: any(gene_trees.get(x['CHROM'], IntervalTree()).overlap(x['POS'], x['POS']+1)),
            axis=1
        )

        # Filter the VCF variants to keep only those in genic regions
        final_df = roi_filtered_df[roi_filtered_df['in_genic']]
        final_count = len(final_df)
        
        # Drop the 'in_genic' column
        final_df = final_df.drop(columns=['in_genic'])
        
        logger.info(f"Final result: {final_count} variants in genic regions " +
                    f"({final_count/roi_filtered_count:.2%} of ROI variants, " +
                    f"{final_count/initial_variant_count:.2%} of all variants)")

        # Add this to verify counts are correct:
        if 'headers' in vcf_df.attrs:
            header_count = len(vcf_df.attrs['headers']) + 1  # +1 for the column header
            logger.info(f"VCF contains {header_count} header lines and {final_count} variant lines")
        
        # Write the filtered variants to the output VCF file
        logger.info(f"Writing filtered variants to {output}")
        
        # Reconstruct proper VCF file with headers
        if 'headers' in vcf_df.attrs and 'column_header' in vcf_df.attrs:
            with open(output, 'w') as f:
                # Write metadata headers
                for header in vcf_df.attrs['headers']:
                    f.write(header + '\n')  # Add newline after each header
                # Write column header and data
                f.write(vcf_df.attrs['column_header'])
                final_df.to_csv(f, sep='\t', index=False, header=False)
        else:
            # Fallback if headers weren't captured
            final_df.to_csv(output, sep='\t', index=False)
        
        logger.info(f"Successfully wrote filtered VCF to {output}")
        
    except Exception as e:
        logger.error(f"Error in mask_variants: {e}", exc_info=True)
        raise