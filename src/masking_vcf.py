"""
This script masks a VCF file by preserving only the variants that are located in genic regions, as defined by a GFF3 file. It also filters the variants based 
on a list of regions of interest (ROI) provided by the user. The output is a VCF file containing only the relevant variants.

@author: Luis Javier Madrigal-Roca & John K. Kelly

@date 2025/05/16

"""

import io
import pandas as pd
from tqdm import tqdm
from intervaltree import IntervalTree

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
    vcf_df = vcf_df[vcf_df.apply(
        lambda x: any((x['CHROM'] == r['Chrom']) and 
                      (r['Start'] <= x['POS'] <= r['End'])
                      for _, r in roi_df.iterrows()),
        axis=1)]
    
    # Build interval trees per chromosome (much faster)
    gene_trees = {}
    for _, row in gff3_df.iterrows():
        chrom = row[0]
        if chrom not in gene_trees:
            gene_trees[chrom] = IntervalTree()
        gene_trees[chrom].add(row[3], row[4]+1, row[8])  # Start, end, gene info

    # Efficient lookup
    vcf_df['in_genic'] = vcf_df.apply(
        lambda x: any(gene_trees.get(x['CHROM'], IntervalTree()).overlap(x['POS'], x['POS']+1)),
        axis=1
    )

    # Filter the VCF variants to keep

    vcf_df = vcf_df[vcf_df['in_genic']]

    # Drop the 'in_genic' column
    vcf_df = vcf_df.drop(columns=['in_genic'])
    
    # Write the filtered variants to the output VCF file
    vcf_df.to_csv(output, sep='\t', index=False)