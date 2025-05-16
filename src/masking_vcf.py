"""
This script masks a VCF file by preserving only the variants that are located in genic regions, as defined by a GFF3 file. It also filters the variants based 
on a list of regions of interest (ROI) provided by the user. The output is a VCF file containing only the relevant variants.

@author: Luis Javier Madrigal-Roca & John K. Kelly

@date 2025/05/16

"""

import io
import pandas as pd
from tqdm import tqdm

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

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
    
    # Read the VCF file
    vcf_df = read_vcf(vcf)
    
    # Read the GFF3 file
    gff3_df = pd.read_csv(gff3, sep='\t', comment='#', header=None)

    # Filter the GFF3 file to keep only genic regions
    gff3_df = gff3_df[gff3_df[2].str.contains('gene')]
    
    # Read the ROI list
    roi_df = pd.read_csv(roi_list, sep='\t', header=0)

    # Filter the VCF variants to keep only the chromosomes and positions in the ROI list
    
    vcf_df = vcf_df[vcf_df['CHROM'].isin(roi_df[1]) & vcf_df['POS'] > roi_df[2].astype(int) & vcf_df['POS'] < roi_df[3].astype(int)]
    
    # Determine which variants are in genic regions

    vcf_df['in_genic'] = False

    for _, row in tqdm(gff3_df.iterrows(), total=gff3_df.shape[0], desc='Filtering variants'):
        chrom = row[0]
        start = row[3] + 1
        end = row[4] + 1

        # Check if the variant is in the genic region
        vcf_df.loc[(vcf_df['CHROM'] == chrom) & (vcf_df['POS'] >= start) & (vcf_df['POS'] <= end), 'in_genic'] = True

    # Filter the VCF variants to keep

    vcf_df = vcf_df[vcf_df['in_genic']]

    # Drop the 'in_genic' column
    vcf_df = vcf_df.drop(columns=['in_genic'])
    
    # Write the filtered variants to the output VCF file
    vcf_df.to_csv(output, sep='\t', index=False)