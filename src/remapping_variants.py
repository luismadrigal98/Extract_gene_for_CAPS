"""
Utility functions for remapping ROI between two assemblies using minimap2.

@author: Luis J. Madrigal-Roca

@date: 2025/05/16

"""

import os
import subprocess
import pandas as pd
import logging
from Bio import SeqIO
from Bio.Seq import Seq

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