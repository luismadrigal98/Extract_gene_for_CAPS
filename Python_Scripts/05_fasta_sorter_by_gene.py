"""
This script is used to remake the fasta files grouping in th same file the genes with the same name.

Usage:

@Author: Luis J. Madrigal Roca
@Date: 2024-07-15

"""

import os
import argparse

def fasta_sorter(gene_name, fasta):
    """
    Given a gene name and a fasta file, extract the gene sequences from the reference fasta files.
    """
    gene_sequence = ""
    found_gene = False  # Flag to indicate when the correct gene has been found

    with open(fasta, "r") as f:
        fasta_lines = f.readlines()
        for line in fasta_lines:
            if found_gene:  # If the gene was found in the previous iteration, this line is the sequence
                gene_sequence = line.strip()  # Use strip() to remove newline characters
                break  # Exit the loop after capturing the sequence
            if line.startswith(f">{gene_name}"):
                found_gene = True  # Set the flag to True when the correct gene header is found

    return gene_sequence

def main():
    parser = argparse.ArgumentParser(description = "Fasta sorter by gene name")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-id", "--input_directory", help="Directory with the gene fasta file(s)", nargs='+')
    group.add_argument("-i", "--input", help="Gene fasta file(s)", nargs='+')

    parser.add_argument("-g", "--gene_name", help="Gene name(s) to extract", required=True, nargs='+')

    args = parser.parse_args()

    if not (args.input_directory or args.input):
        parser.error("At least one of --input_directory or --input is required.")
    
    genes = args.gene_name

    if args.input_directory:
        fasta_files = []
        for directory in args.input_directory:
            for file in os.listdir(directory):
                if file.endswith(".fasta") or file.endswith(".fa"):
                    fasta_files.append(os.path.join(directory, file))

    else:
        fasta_files = args.input

    for gene in genes:
        
        record = {}
        
        for index, fasta in enumerate(fasta_files):
            gene_sequence = fasta_sorter(gene, fasta)
            record[fasta_files[index].strip().split('.')[2]] = gene_sequence
        
        output_file = f"{gene}_grouped.fasta"
        with open(output_file, "w") as f:
            for key in record:
                f.write(f">{key}\n")
                f.write(f"{record[key]}\n")
    
if __name__ == "__main__":
    main()