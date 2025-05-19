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
    
    # Check if the output directory exists, if not create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for tsv_file in tsv_list:
        name = os.path.basename(tsv_file).split('.')[0]
        logging.info(f"Processing file: {name}")
        
        # Read the TSV file
        tsv_data[tsv_file] = read_tsv(tsv_file)
        
        if tsv_data[tsv_file].empty:
            logging.warning(f"File {tsv_file} is empty. Skipping.")
            continue
    
    for file in tqdm(tsv_data.keys(), desc="Processing TSV files", unit="file", total=len(tsv_data)):
        df = tsv_data[file]
        
        # Get the min and max positions for this ROI
        roi_min_pos = df['POS'].min()
        roi_max_pos = df['POS'].max()
        
        # Filter the columns based on the allele column pattern
        allele_cols = [col for col in df.columns if allele_col_pattern in col]
        
        if not allele_cols:
            logging.warning(f"No allele columns found in {file}. Skipping.")
            continue
        
        # Filter the rows based on the reliability threshold
        if reliability_thr == 'high':
            threshold = 'high'
        elif reliability_thr == 'medium':
            threshold = ['high', 'medium']
        else:  # 'low' includes all levels
            threshold = ['high', 'medium', 'low']
            
        df_filtered = df[df['overall_reliability'].isin(threshold) if isinstance(threshold, list) else df['overall_reliability'] == threshold]
        logging.info(f"Filtered {len(df_filtered)} rows based on reliability threshold in {file}.")
        
        if df_filtered.empty:
            logging.warning(f"No rows found after filtering by reliability in {file}. Skipping.")
            continue
        
        # Find the common parental allele column
        diff_cols = [col for col in allele_cols if diff_parent in col]
        if not diff_cols:
            logging.warning(f"Common parental allele column not found for {diff_parent}. Skipping.")
            continue
            
        diff_col = diff_cols[0]  # Use the first matching column
        alt_cols = [col for col in allele_cols if col != diff_col]
        
        # Filter variants where ALL alt parents differ from common parent
        mask = pd.Series(True, index=df_filtered.index)  # Start with all True
        for alt_col in alt_cols:
            mask &= (df_filtered[diff_col] != df_filtered[alt_col])  # AND operation instead of OR
        
        df_filtered = df_filtered[mask]
        logging.info(f"Filtered to {len(df_filtered)} variants with different parental alleles.")
        
        if df_filtered.empty:
            logging.warning(f"No variants remaining after parental allele filtering. Skipping.")
            continue
        
        # Select only relevant columns (variant info, allele columns, and individual reliability, in this order)

        relevant_cols = ['POS', 'REF', 'ALT', 'overall_reliability'] + [diff_col] + alt_cols + [col for col in df_filtered.columns if col.endswith('_reliability')]

        df_filtered = df_filtered[relevant_cols]

        # Create new columns for primer compliance
        df_filtered['primer_compliant'] = False
        df_filtered['amplicon_start'] = None
        df_filtered['amplicon_end'] = None
        df_filtered['displacement'] = 0
        
        # THIS IS THE KEY FIX: Get all variant positions for primer region checking
        all_variant_positions = set(df_filtered['POS'].tolist())
        
        # Check primer placement for each variant
        for idx, row in tqdm(df_filtered.iterrows(), desc="Checking primers", unit="row", total=len(df_filtered)):
            variant_pos = row['POS']
            
            # Calculate amplicon boundaries
            amplicon_start = variant_pos - (amplicon_size // 2)
            amplicon_end = variant_pos + (amplicon_size // 2)
            
            # Skip if amplicon extends beyond ROI
            if amplicon_start < roi_min_pos or amplicon_end > roi_max_pos:
                continue
            
            # Define primer regions
            f_primer_start = amplicon_start
            f_primer_end = amplicon_start + primer_size
            r_primer_start = amplicon_end - primer_size
            r_primer_end = amplicon_end
            
            # Check if primer regions overlap with other variants
            f_primer_has_variant = any(pos != variant_pos and f_primer_start <= pos <= f_primer_end 
                                        for pos in all_variant_positions)
            r_primer_has_variant = any(pos != variant_pos and r_primer_start <= pos <= r_primer_end 
                                        for pos in all_variant_positions)
            
            # If no variants in primer regions, mark as compliant
            if not f_primer_has_variant and not r_primer_has_variant:
                df_filtered.at[idx, 'primer_compliant'] = True
                df_filtered.at[idx, 'amplicon_start'] = amplicon_start
                df_filtered.at[idx, 'amplicon_end'] = amplicon_end
                continue

            # If displacement is allowed, try displacing the amplicon
            is_compliant = False
            if displacement and not is_compliant:
                for step in range(1, displacement_steps + 1):
                    # Try shifting left
                    left_shift = -step
                    left_amplicon_start = amplicon_start + left_shift
                    left_amplicon_end = amplicon_end + left_shift
                    
                    left_f_primer_start = left_amplicon_start
                    left_f_primer_end = left_amplicon_start + primer_size
                    left_r_primer_start = left_amplicon_end - primer_size
                    left_r_primer_end = left_amplicon_end
                    
                    left_f_primer_has_variant = any(pos != variant_pos and left_f_primer_start <= pos <= left_f_primer_end 
                                                for pos in all_variant_positions)
                    left_r_primer_has_variant = any(pos != variant_pos and left_r_primer_start <= pos <= left_r_primer_end 
                                                for pos in all_variant_positions)
                    
                    # Check if variant is still in amplicon
                    variant_in_left_amplicon = left_amplicon_start <= variant_pos <= left_amplicon_end
                    
                    if not left_f_primer_has_variant and not left_r_primer_has_variant and variant_in_left_amplicon:
                        df_filtered.at[idx, 'primer_compliant'] = True
                        df_filtered.at[idx, 'amplicon_start'] = left_amplicon_start
                        df_filtered.at[idx, 'amplicon_end'] = left_amplicon_end
                        df_filtered.at[idx, 'displacement'] = left_shift
                        is_compliant = True
                        break
                    
                    # Try shifting right
                    right_shift = step
                    right_amplicon_start = amplicon_start + right_shift
                    right_amplicon_end = amplicon_end + right_shift
                    
                    right_f_primer_start = right_amplicon_start
                    right_f_primer_end = right_amplicon_start + primer_size
                    right_r_primer_start = right_amplicon_end - primer_size
                    right_r_primer_end = right_amplicon_end
                    
                    right_f_primer_has_variant = any(pos != variant_pos and right_f_primer_start <= pos <= right_f_primer_end 
                                                for pos in all_variant_positions)
                    right_r_primer_has_variant = any(pos != variant_pos and right_r_primer_start <= pos <= right_r_primer_end 
                                                for pos in all_variant_positions)
                    
                    # Check if variant is still in amplicon
                    variant_in_right_amplicon = right_amplicon_start <= variant_pos <= right_amplicon_end
                    
                    if not right_f_primer_has_variant and not right_r_primer_has_variant and variant_in_right_amplicon:
                        df_filtered.at[idx, 'primer_compliant'] = True
                        df_filtered.at[idx, 'amplicon_start'] = right_amplicon_start
                        df_filtered.at[idx, 'amplicon_end'] = right_amplicon_end
                        df_filtered.at[idx, 'displacement'] = right_shift
                        is_compliant = True
                        break

        # Save the processed DataFrame to the output directory
        output_file = os.path.join(output_dir, os.path.basename(file))
        df_filtered.to_csv(output_file, sep='\t', index=False)
        
        logging.info(f"Processed file saved to {output_file}")