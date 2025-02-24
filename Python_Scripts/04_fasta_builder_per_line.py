"""

This script is designed to construct FASTA files from given sequences, processing each line of input to generate a corresponding FASTA formatted 
output. It is capable of handling multiple sequences, ensuring each is properly formatted with a header and sequence data according to FASTA 
standards. The script is intended for use in bioinformatics applications where FASTA file creation is required from raw or pre-processed sequence 
data.

Usage:
    The script can be executed with command line arguments specifying input data sources and output file preferences. It requires input data 
    that includes sequence identifiers and nucleotide or amino acid sequences.

Features:
    - Reads sequence data from specified input sources.
    - Formats the sequence data into FASTA format, ensuring compliance with the standard.
    - Supports multiple sequences, processing each line of input as a separate sequence entry.
    - Allows customization of output file names and directories.

Requirements:
    This script requires Python 3.x and may depend on external libraries for handling specific file formats or data processing tasks.

Author: Luis Javier Madrigal Roca
Date: 2024-07-15
"""

import argparse
import subprocess
import os

def process_fasta(fasta_lines):
    """
    Process a list of FASTA file lines, merging sequence lines that belong to the same header.
    """
    processed_lines = []
    current_sequence = ""
    
    for line in fasta_lines:
        if line.startswith(">"):
            if current_sequence:  # If there's a sequence accumulated, append it before starting a new one
                processed_lines.append(current_sequence)
                current_sequence = ""
            processed_lines.append(line)  # Append the new header
        else:
            current_sequence += line.strip()  # Accumulate sequence lines
    if current_sequence:  # Don't forget to append the last sequence
        processed_lines.append(current_sequence)
    
    return processed_lines

def seq_extractor(gene_start, gene_end, chromosome, fasta):
    """
    Given a chromosome number and the position of the gene of interest, extract the gene sequence from the reference fasta file.
    """
    gene_sequence = ""
    found_chromosome = False  # Flag to indicate when the correct chromosome has been found

    for line in fasta:
        if found_chromosome:  # If the chromosome was found in the previous iteration, this line is the sequence
            gene_sequence = line[gene_start-1:gene_end]  # Extract the gene sequence using the provided start and end
            break  # Exit the loop after capturing the sequence
        if line.startswith(f">{chromosome}"):
            found_chromosome = True  # Set the flag to True when the correct chromosome header is found

    return gene_sequence

def main():
    parser = argparse.ArgumentParser(description = "Gene sequence extractor for different fasta files. It uses as input a multifasta file with all the gene sequences that are of interest")

    parser.add_argument("-i", "--input", help = "Gene fasta file(s)", required = True, nargs='+')
    parser.add_argument("-targ", "--target", help = "Target line reference directories from which you want to pull the corresponding gene sequences. The directories must correspond to the input files", required = True, nargs='+')
    parser.add_argument("-op", "--output_prefix", help = "Output fasta file", required = False, default = "output_seq")
    parser.add_argument("-od", "--output_directory", help = "Directory where all the outputs are going to be stored", required = False)

    args = parser.parse_args()

    args.input = [os.path.abspath(f) for f in args.input]
    args.target = [os.path.abspath(d) for d in args.target]
    if args.output_directory:
        args.output_directory = os.path.abspath(args.output_directory)
        os.makedirs(args.output_directory, exist_ok=True)

    fasta_files = args.input
    target_directories = args.target

    targ = [os.listdir(target) for target in target_directories]

    for subdir in range(len(target_directories)):
        for file in targ[subdir]:
            gene_sequence = {}

            temp_dir = os.path.abspath(f"Temporary_{file}")
            if not os.path.exists(temp_dir):
                try:
                    os.makedirs(temp_dir)
                except OSError as e:
                    print(f"Error creating temporary directory {temp_dir}: {e}")
                    continue

            output_file = f"out_{os.path.basename(fasta_files[subdir])}_vs_{file}"
            
            print(f"Temporary directory: {temp_dir}")
            print(f"Output file: {output_file}")
            print(f"Full output path: {os.path.join(temp_dir, output_file)}")

            subprocess.run(["minimap2", "-o", os.path.join(temp_dir, output_file), os.path.join(target_directories[subdir], file), fasta_files[subdir]], check=True)

            with open(os.path.join(temp_dir, output_file)) as f:
                lines = f.readlines()
                lines = [line.strip().split("\t") for line in lines]

            for line in lines:
                name = line[0]
                chromosome = line[5]
                st = line[7]
                end = line[8]

                with open(os.path.join(target_directories[subdir], file), "r") as f:
                    reference = f.read().splitlines()

                processed_reference = process_fasta(reference)

                gene_sequence[name] = seq_extractor(int(st), int(end), chromosome, processed_reference)

            output_filename = f"{args.output_prefix}_{os.path.basename(fasta_files[subdir])}_vs_{file}"
            output_path = os.path.join(args.output_directory, output_filename)

            print(f"Output filename: {output_filename}")
            print(f"Full output path: {output_path}")

            with open(output_path, "w") as f:
                for key, value in gene_sequence.items():
                    f.write(f">{key}\n{value}\n")

            subprocess.run(["rm", "-r", temp_dir], check=True)

if __name__ == "__main__":
    main()