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
import logging

# Add the path to the src directory to the system path
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

setup_logging = lambda: logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("screen_variants.log"),
        logging.StreamHandler()
    ]
)

read_tsv = lambda file_path: pd.read_csv(file_path, sep='\t', header=0, low_memory=False)

def screen_variants(tsv_list, output_dir, allele_col_pattern, reliability_thr, diff_parent, amplicon_size, primer_size, displacement, displacement_steps):
    """
    Screen the variants for diagnostic markers.
    
    Parameters:
    tsv_list (list): List of TSV files to screen.
    output_dir (str): Directory to save the output files.
    allele_col_pattern (str): Pattern to match the allele columns.
    reliability_thr (float): Minimum reliability threshold.
    diff_parent (str): Identifier for the common parental allele. E. g. "664c"
    amplicon_size (int): Size of the amplicon.
    primer_size (int): Size of the primer.
    displacement (bool): Whether we can move the amplicon window to accommodate primers in a varian region.
    displacement_steps (int): Number of steps to move the amplicon window.
    """
    
    # Read the TSV files

    tsv_data = {}

    logging.info("Reading TSV files...")

    for tsv_file in tsv_list:
        
        name = os.path.basename(tsv_file).split('.')[0]
        logging.info(f"Processing file: {name}")

        # Read the TSV file
        tsv_data[tsv_file] = read_tsv(tsv_file)
        
        # Check if the file is empty
        if tsv_data[tsv_file].empty:
            logging.warning(f"File {tsv_file} is empty. Skipping.")
            continue
    
    # Process each TSV file

    for file in tqdm(tsv_data.keys(), desc="Processing TSV files", unit="file", total=len(tsv_data)):
        
        df = tsv_data[file]
        
        # Filter the columns based on the allele column pattern
        allele_cols = [col for col in df.columns if allele_col_pattern in col]
        
        if not allele_cols:
            logging.warning(f"No allele columns found in {file}. Skipping.")
            continue
        
        # Filter the rows based on the reliability threshold
        df_filtered = df[df['overall_reliability'] >= reliability_thr]

        logging.info(f"Filtered {len(df_filtered)} rows based on reliability threshold in {file}.")
        
        if df_filtered.empty:
            logging.warning(f"No rows found after filtering by reliability in {file}. Skipping.")
            continue
        
        # MAIN LOGIC HERE
        # First, in order to be useful the variant, the common parental allele must be different from the alternative alleles

        logging.info(f"Filtering variants based on common parental allele {diff_parent} in {file}.")

        # Filter out rows where the common parental allele is the same as the alternative alleles

        for col in allele_cols:
            diff_col = allele_cols[diff_parent in allele_cols]
            alt_cols = [col for col in allele_cols if col != diff_col]

            # Check if the common parental allele is different from the alternative alleles
            df_filtered = df_filtered[~df_filtered[diff_col].isin(df_filtered[alt_cols].values.flatten())]
            logging.info(f"Filtered {len(df_filtered)} rows based on common parental allele in {file}.")
        
        # Save the processed DataFrame to the output directory
        output_file = os.path.join(output_dir, os.path.basename(file))
        df_filtered.to_csv(output_file, sep='\t', index=False)
        
        logging.info(f"Processed file saved to {output_file}")