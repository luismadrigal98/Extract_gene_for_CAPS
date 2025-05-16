"""
Utility functions for remapping ROI between two assemblies using minimap2.

@author: Luis J. Madrigal-Roca

@date: 2025/05/16

"""

import os
import subprocess
import pandas as pd
import logging
import shutil
from tempfile import mkdtemp
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("marker_wizard.log"),
        logging.StreamHandler()
    ]
)

def index_fasta(fasta_file, samtools_exe = "/home/l338m483/.conda/envs/PyR/bin/samtools"):
    """
    Index a fasta file using samtools faidx.
    
    :param fasta_file: Path to the fasta file to be indexed.
    :param samtools_exe: Path to the samtools executable. If not provided, it will use the default in the PATH.
    :return: None
    """
    
    # Check if the fasta file exists
    if not os.path.isfile(fasta_file):
        logging.error(f"Fasta file {fasta_file} does not exist.")
        raise FileNotFoundError(f"Fasta file {fasta_file} does not exist.")
    
    # Check if there is already an index file
    index_file = f"{fasta_file}.fai"
    if os.path.isfile(index_file):
        logging.info(f"Index file {index_file} already exists. Skipping indexing.")
        return

    # Build the command to index the fasta file
    command = [samtools_exe, "faidx", fasta_file]
    command = " ".join(command)
    logging.info(f"Indexing fasta file with command: {command}")
    subprocess.run(command, shell=True, check=True)

def build_fasta_for_ROI(ROI_list, current_ref, output_file):
    """
    Build a fasta file for the regions of interest (ROI) from the current reference genome.
    
    :param ROI_list: Path to the ROI list file. This is a tab-separated file with the format: <chromosome> <start> <end>.
    :param current_ref: Path to the current reference genome file in fasta format.
    :param output_file: Path to the output fasta file where the sequences will be saved.
    
    :return: Path to the output fasta file. 
    """
    
    # Read the ROI list - try with flexible whitespace parsing
    try:
        # First attempt with explicit tab
        roi_df = pd.read_csv(ROI_list, sep="\t", header=0, names=["ROI_name", "Chrom", "Start", "End"])
    except pandas.errors.ParserError:
        # If that fails, try with flexible whitespace
        logging.info("Tab delimiter failed, trying with flexible whitespace parsing")
        roi_df = pd.read_csv(ROI_list, delim_whitespace=True, header=0, names=["ROI_name", "Chrom", "Start", "End"])
        
    # Read the current reference genome
    ref_sequences = SeqIO.to_dict(SeqIO.parse(current_ref, "fasta"))
    
    # Extract sequences for each ROI and write to output fasta file
    with open(output_file, "w") as out_fasta:
        for _, row in tqdm(roi_df.iterrows(), total=roi_df.shape[0], desc="Extracting ROIs"):
            name = row["ROI_name"]
            chrom = row["Chrom"]
            start = row["Start"] + 1  # Convert to 1-based indexing
            end = row["End"] + 1
            
            if chrom in ref_sequences:
                seq = ref_sequences[chrom].seq[start:end]
                SeqIO.write(SeqIO.SeqRecord(seq, id=f"{name}_{chrom}:{start}-{end}", description=""), out_fasta, "fasta")
            else:
                logging.warning(f"Chromosome {chrom} not found in reference genome.")
    
    logging.info(f"Fasta file for ROI saved to {output_file}")

def run_minimap2(current_features, new_ref, output_file, minimap2, minimap2_opts = "-x asm5 -t 10 -p 0.9 -N 0"):
    """
    Run minimap2 to remap the current features to the new reference genome.
    
    :param current_features: Path to the current features file. This is a fasta file with the sequences to be remapped.
    :param new_ref: Path to the new reference genome file in fasta format.
    :param output_file: Path to the output file where the remapped features will be saved.
    :param minimap2: Path to the minimap2 executable.
    :param minimap2_opts: Options for minimap2.
    :return: Path to the output file with remapped features.
    """
        
    # Define the minimap2 command
    command = [minimap2, minimap2_opts, '-o', output_file, new_ref, current_features]
    command = " ".join(command)
    logging.info(f"Running minimap2 with command: {command}")
    subprocess.run(command, shell=True, check=True)
    
    return output_file

def process_minimap2_output(paf_file, output_file):
    """
    This script processes the output of minimap2 and generate a new coordinate file to work over the new reference genome.

    :param paf_file: Path to the minimap2 output file in PAF format.
    :param output_file: Path to the output file where the processed coordinates will be saved.

    :return: None
    """
    
    # Read the PAF file
    paf_df = pd.read_csv(paf_file, sep='\t', header=None)

    # Get the number of columns in the PAF file
    num_columns = paf_df.shape[1]

    # Define the mandatory column names (first 12 columns are standard in PAF)
    mandatory_columns = [
        'query_name', 'query_length', 'query_start', 'query_end', 'strand',
        'target_name', 'target_length', 'target_start', 'target_end',
        'num_matching_bases', 'alignment_block_length', 'mapping_quality'
        ]
    
    # Define the column names based on the number of columns
    # Add column names based on what's available
    if num_columns >= 12:
        # Add mandatory columns
        column_names = mandatory_columns.copy()
        
        # Add remaining columns with generic names
        for i in range(12, num_columns):
            column_names.append(f'col_{i+1}')
            
        paf_df.columns = column_names
    else:
        raise ValueError(f"PAF file has too few columns: {num_columns}, expected at least 12")
    
    # Filter duplicate entries (duplicated query names) and keep only the one with the highest mapping quality
    paf_df = paf_df.sort_values('mapping_quality', ascending=False).drop_duplicates('query_name')

    # Checking that all duplicates were removed
    assert paf_df['query_name'].duplicated().sum() == 0

    # Reordering by query name
    paf_df = paf_df.sort_values('query_name')

    # Export the processed PAF data into a new ROI file, tab delimited
    with open(output_file, 'w') as out_f:
        # Write the header
        out_f.write("ROI_name\tChrom\tStart\tEnd\n")
        
        # Iterate over the rows of the PAF dataframe
        for _, row in tqdm(paf_df.iterrows(), total=paf_df.shape[0], desc="Processing PAF output"):
            roi_name = row['query_name']
            chrom = row['target_name']
            start = int(row['target_start']) + 1
            end = int(row['target_end']) + 1
            out_f.write(f"{roi_name}\t{chrom}\t{start}\t{end}\n")
    logging.info(f"Processed PAF output saved to {output_file}")

def remap_variants(current_ref, new_ref, ROI_list, output_file, minimap2, minimap2_opts, samtools_exe, temp_dir, keep_temp = False):
    """
    Remap the variants from the current reference genome to the new reference genome.

    :param current_ref: Path to the current reference genome file in fasta format.
    :param new_ref: Path to the new reference genome file in fasta format.
    :param ROI_list: Path to the ROI list file. This is a tab-separated file with the format: <chromosome> <start> <end>.
    :param output_file: Path to the output file where the remapped variants will be saved.
    :param minimap2: Path to the minimap2 executable.
    :param minimap2_opts: Options for minimap2.
    :param temp_dir: Temporary directory for intermediate files.
    :param keep_temp: If True, keep the temporary files. If False, delete them after use.
    
    :return: None
    """
    
    # Index the new reference genome
    index_fasta(new_ref, samtools_exe = samtools_exe)

    # Create the temporary directory if it doesn't exist
    if temp_dir is None:
        temp_dir = mkdtemp(dir = os.getcwd())
        logging.info(f"Temporary directory created at {temp_dir}")
    elif not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    
    # Build the fasta file for the ROI
    roi_fasta = os.path.join(temp_dir, "ROI_sequences.fasta")
    build_fasta_for_ROI(ROI_list, current_ref, roi_fasta)
    
    # Create a temporary name for the paf file
    paf_file = os.path.join(temp_dir, "minimap2_output.paf")

    # Run minimap2 to remap the variants
    run_minimap2(roi_fasta, new_ref, paf_file, minimap2, minimap2_opts)

    # Process the minimap2 output and retrieve the new coordinates
    process_minimap2_output(paf_file, output_file)
    
    logging.info(f"Remapping completed. Output saved to {output_file}")

    # Clean up temporary files if specified
    if not keep_temp:
        shutil.rmtree(temp_dir)
        logging.info(f"Temporary files in {temp_dir} have been deleted.")
    else:
        logging.info(f"Temporary files in {temp_dir} have been kept.")