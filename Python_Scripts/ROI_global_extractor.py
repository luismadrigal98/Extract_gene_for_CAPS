"""
This script extracts an entire Region of Interest (ROI) from a reference genome
and maps it to corresponding regions in target sequences.

Usage:
    python roi_extractor_mapper.py -r "Chr_11:16838722-16843579" -f reference.fasta -t target1.fasta target2.fasta -o output_prefix

Author: Luis J. Madrigal Roca
"""

import argparse
import os
import subprocess
import tempfile
import shutil

def process_fasta(fasta_lines):
    """
    Process a list of FASTA file lines, merging sequence lines that belong to the same header.
    """
    processed_lines = []
    current_header = None
    current_sequence = ""
    
    for line in fasta_lines:
        line = line.strip()
        if line.startswith(">"):
            if current_header:
                processed_lines.append((current_header, current_sequence))
                current_sequence = ""
            current_header = line[1:]  # Remove '>' from header
        else:
            current_sequence += line
            
    if current_header and current_sequence:
        processed_lines.append((current_header, current_sequence))
    
    return processed_lines

def roi_extractor(chromosome, start, end, fasta_data):
    """
    Extract the sequence of a specific ROI from reference fasta data.
    
    Args:
        chromosome: Chromosome or contig name
        start: Start position (1-based)
        end: End position (1-based)
        fasta_data: Processed fasta data as list of (header, sequence) tuples
        
    Returns:
        The extracted ROI sequence or empty string if not found
    """
    for header, sequence in fasta_data:
        if header.split()[0] == chromosome:
            return sequence[start-1:end]
    return ""

def map_roi(roi_file, target_file, temp_dir):
    """
    Map ROI to target sequence using minimap2 and return mapped positions.
    """
    # First check if the target file exists
    if not os.path.exists(target_file):
        print(f"Warning: Target file not found: {target_file}")
        return []
    
    # Create a more descriptive output filename with basename only
    roi_basename = os.path.basename(roi_file)
    target_basename = os.path.basename(target_file)
    output_file = os.path.join(temp_dir, f"map_{roi_basename}_vs_{target_basename}.paf")
    
    try:
        # Add more verbose output for debugging
        print(f"Running minimap2 with target: {target_file}")
        print(f"Output will be written to: {output_file}")
        
        # Check if temp_dir exists and is writable
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir, exist_ok=True)
            print(f"Created temp directory: {temp_dir}")
        
        if not os.access(temp_dir, os.W_OK):
            print(f"Warning: Temp directory {temp_dir} is not writable")
            return []
        
        # Use the -a option for SAM output which may be more reliable
        cmd = ["minimap2", target_file, roi_file, "-o", output_file]
        print(f"Running command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            check=False,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True
        )
        
        # Handle error cases
        if result.returncode != 0:
            print(f"Error running minimap2 for {target_file}: {result.stderr}")
            print(f"Command output: {result.stdout}")
            return []
            
        # Check if output file was created
        if not os.path.exists(output_file):
            print(f"Warning: minimap2 output file not created: {output_file}")
            return []
        
        if os.path.getsize(output_file) == 0:
            print(f"Warning: minimap2 output file is empty: {output_file}")
            return []
            
        # Process mapping results
        mappings = []
        with open(output_file) as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) >= 10:  # PAF format has at least 12 fields
                    roi_name = fields[0]
                    target_name = fields[5]
                    target_start = int(fields[7])
                    target_end = int(fields[8])
                    mappings.append((roi_name, target_name, target_start, target_end))
        
        print(f"Found {len(mappings)} mappings in {output_file}")
        return mappings
        
    except Exception as e:
        print(f"Exception while processing {target_file}: {str(e)}")
        return []

def main():
    parser = argparse.ArgumentParser(description="ROI sequence extractor and mapper")
    
    parser.add_argument("-r", "--roi", help="ROI in format 'Chr:start-end' or a file with ROIs (one per line)", required=True)
    parser.add_argument("-f", "--reference", help="Reference fasta file", required=True)
    parser.add_argument("-t", "--targets", help="Target fasta file(s) to map the ROI against", required=True, nargs='+')
    parser.add_argument("-o", "--output_prefix", help="Prefix for output files", default="roi_sequences")
    parser.add_argument("--temp_dir", help="Directory to store temporary files", default="./temp_files")
    parser.add_argument("--remove_temp", help="Remove temporary directory after processing", action="store_true")

    args = parser.parse_args()
    
    # Create a custom temp directory in the current working directory
    temp_dir = args.temp_dir
    os.makedirs(temp_dir, exist_ok=True)
    print(f"Using temporary directory: {temp_dir}")
    
    # Parse ROI information
    roi_list = []
    if os.path.isfile(args.roi):
        with open(args.roi) as f:
            for line in f:
                if line.strip():
                    roi_parts = line.strip().split(":")
                    if len(roi_parts) == 2 and "-" in roi_parts[1]:
                        chrom = roi_parts[0]
                        start, end = map(int, roi_parts[1].split("-"))
                        roi_list.append((chrom, start, end))
    else:
        roi_parts = args.roi.split(":")
        if len(roi_parts) == 2 and "-" in roi_parts[1]:
            chrom = roi_parts[0]
            start, end = map(int, roi_parts[1].split("-"))
            roi_list.append((chrom, start, end))
    
    if not roi_list:
        print("Error: Invalid ROI format. Use 'Chr:start-end' or provide a file with ROIs.")
        return
    
    # Read reference fasta
    with open(args.reference, "r") as f:
        reference_lines = f.read().splitlines()
    
    processed_reference = process_fasta(reference_lines)
    
    try:
        for _, (chrom, start, end) in enumerate(roi_list):
            roi_name = f"{chrom}_{start}_{end}"
            roi_seq = roi_extractor(chrom, start, end, processed_reference)
            
            if not roi_seq:
                print(f"Warning: Could not extract ROI sequence for {roi_name}")
                continue
                
            # Create ROI fasta file in the user-specified temp directory
            roi_fasta = os.path.join(temp_dir, f"roi_{roi_name}.fasta")
            with open(roi_fasta, "w") as f:
                f.write(f">{roi_name}\n{roi_seq}\n")
            
            print(f"Created ROI file: {roi_fasta}")
            
            # Map ROI to each target
            mapped_sequences = {}
            
            # Filter out any invalid target files
            valid_targets = []
            for target_file in args.targets:
                if os.path.exists(target_file):
                    valid_targets.append(target_file)
                else:
                    print(f"Warning: Target file not found: {target_file}")
            
            if not valid_targets:
                print("Error: No valid target files provided.")
                return

            for target_file in valid_targets:
                # Read target fasta
                with open(target_file, "r") as f:
                    target_lines = f.read().splitlines()
                
                processed_target = process_fasta(target_lines)
                
                # Map ROI to target using minimap2
                mappings = map_roi(roi_fasta, target_file, temp_dir)
                
                # Extract mapped regions
                for roi_name, target_name, target_start, target_end in mappings:
                    target_seq = ""
                    for header, sequence in processed_target:
                        if header.split()[0] == target_name:
                            target_seq = sequence[target_start-1:target_end]
                            break
                    
                    if target_seq:
                        target_basename = os.path.basename(target_file).split('.')[0]
                        mapped_sequences[f"{target_basename}_{target_name}_{target_start}_{target_end}"] = target_seq
            
            # Write output file with all mapped sequences
            output_file = f"{args.output_prefix}_ROI_{roi_name}.fasta"
            with open(output_file, "w") as f:
                f.write(f">{chrom}_{start}_{end}_reference\n{roi_seq}\n")
                for header, sequence in mapped_sequences.items():
                    f.write(f">{header}\n{sequence}\n")
                
            print(f"ROI {roi_name} and {len(mapped_sequences)} mapped sequences saved to {output_file}")
    
    finally:
        if remove_temp:
            shutil.rmtree(temp_dir)
            print(f"Removed temporary directory: {temp_dir}")

if __name__ == "__main__":
    main()