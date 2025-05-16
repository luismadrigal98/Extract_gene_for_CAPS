"""
This script is going to apply the set of rules we defined to detect the diagnostic markers.

They can be summarized as follows:

1. The variant must be in the ROI list. (previous scripts)
2. The variant must be in genic regions. (previous scripts)
3. 

"""

import os
import sys
from tqdm import tqdm
import pandas as pd

# Add the path to the src directory to the system path
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from src.masking_vcf import read_vcf

def screen_variants(vcf_file, ROI_list, output_dir, min_qual, min_dp,
                    distance_to_closest_marker, non_informative_thr_F2s,
                    heterozygous_thr_support_F2s):
    """
    Screen the variants for diagnostic markers.
    
    Parameters:
        vcf_file (str): Path to the input VCF file.
        ROI_list (str): Path to the list of regions of interest (ROI).
        output_dir (str): Directory where the output files will be saved.
        min_qual (float): Minimum quality score for the variants to be considered.
        min_dp (int): Minimum depth of coverage for the variants to be considered.
        distance_to_closest_marker (int): Distance to the closest marker.
        non_informative_thr_F2s (int): Non informative threshold for the F2s.
    """
    
    # Read the VCF file

    vcf_df = read_vcf(vcf_file)

    # Apply the global filtering first

    # Filter by quality score
    vcf_df = vcf_df[vcf_df['QUAL'] >= min_qual]

    # Define sample columns
    sample_columns = vcf_df.columns[9:]

    # Filter by depth (Allelic depth)
    for _, row in tqdm(vcf_df.iterrows(), total=len(vcf_df), desc="Filtering by depth"):
        # Check the AD field for each sample, and ensure that the total depth is above the threshold
        # GT:PL:AD        ./.:0,0,0:0,0   <<< This is the format of the field
        vcf_df['pass_depth'] = vcf_df[sample_columns].apply(
            
        )
        vcf_df = vcf_df[vcf_df['pass_depth']]
        vcf_df.drop(columns=['pass_depth'], inplace=True)