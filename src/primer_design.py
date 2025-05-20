"""
This script is used to design primers for a given DNA sequence.

It leverages Primer3 to design optimal primers for variants identified in previous steps.
The script extracts sequences from a reference genome around target variants, sets up
appropriate Primer3 parameters, and calls Primer3 to generate primer candidates.

@author: Luis Javier Madrigal-Roca & John K. Kelly

@date 2025/05/19

"""

import os
import subprocess
import pandas as pd
import logging
from Bio import SeqIO
from Bio.Seq import Seq
import tempfile
from tqdm import tqdm

def extract_sequence(fasta_file, chrom, start, end):
    """
    Extract a sequence from a FASTA file based on chromosome and position.
    
    Parameters:
    fasta_file (str): Path to the FASTA file.
    chrom (str): Chromosome name.
    start (int): Start position (1-based).
    end (int): End position (1-based).
    
    Returns:
    str: Extracted sequence.
    """
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id == chrom or record.id.split()[0] == chrom:
                # Convert to 0-based indexing for Python slice
                return str(record.seq[start-1:end])
        
        logging.error(f"Chromosome {chrom} not found in reference file.")
        return None
    except Exception as e:
        logging.error(f"Error extracting sequence: {e}")
        return None

def create_primer3_input(name, sequence, target_pos, target_length, settings, output_file):
    """
    Create a Primer3 input file for a given sequence and target position.
    
    Parameters:
    name(str): Name of the sequence.
    sequence (str): The DNA sequence.
    target_pos (int): The position of the variant within the sequence (0-based). If variant is greater than 1 bp (non a SNP), it is the start of the target.
    settings (dict): Dictionary of Primer3 settings.
    output_file (str): Output file path.
    
    Returns:
    str: Path to the created input file.
    """
    # Calculate target region (1 bp at the variant position)
    target = f"{target_pos},{target_length}"
    
    with open(output_file, 'w') as f:
        f.write(f"SEQUENCE_ID={name}\n")
        f.write(f"SEQUENCE_TEMPLATE={sequence}\n")
        f.write(f"SEQUENCE_TARGET={target}\n")
        
        # Add all settings
        for key, value in settings.items():
            f.write(f"{key}={value}\n")
        
        f.write("=\n")
    
    return output_file

def run_primer3(input_file, primer3_exe='~/.conda/envs/salmon/bin/primer3_core', settings_file=None, primer3_args=""):
    """
    Run Primer3 on the input file.
    
    Parameters:
    input_file (str): Path to the Primer3 input file.
    primer3_exe (str): Path to the Primer3 executable.
    settings_file (str, optional): Path to the Primer3 settings file.
    primer3_args (str, optional): Additional Primer3 command line arguments.
    
    Returns:
    str: Primer3 output as string.
    """
    cmd = [primer3_exe]
    
    # Add command line arguments
    if primer3_args:
        cmd.extend(primer3_args.split())
    
    # Add settings file if provided
    if settings_file:
        cmd.extend(["--p3_settings_file", settings_file])
    
    try:
        # Run primer3 with input file redirected to stdin
        with open(input_file, 'r') as f:
            result = subprocess.run(
                cmd, 
                stdin=f,
                text=True,
                capture_output=True,
                check=True
            )
        return result.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running Primer3: {e}")
        logging.error(f"Stderr: {e.stderr}")
        return None
    except Exception as e:
        logging.error(f"Error: {e}")
        return None

def parse_primer3_output(output):
    """
    Parse Primer3 output into a structured format.
    
    Parameters:
    output (str): Primer3 output string.
    
    Returns:
    dict: Dictionary with parsed primer data.
    """
    result = {}
    
    # Parse the Boulder-IO format
    for line in output.split("\n"):
        if line.strip() == "=":
            break
            
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
    
    # Extract primer pairs
    primers = []
    num_returned = int(result.get('PRIMER_PAIR_NUM_RETURNED', '0'))
    
    for i in range(num_returned):
        primer_pair = {
            'pair_penalty': float(result.get(f'PRIMER_PAIR_{i}_PENALTY', '0')),
            'product_size': int(result.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', '0')),
            'left': {
                'sequence': result.get(f'PRIMER_LEFT_{i}_SEQUENCE', ''),
                'start': int(result.get(f'PRIMER_LEFT_{i}', '0,0').split(',')[0]),
                'length': int(result.get(f'PRIMER_LEFT_{i}', '0,0').split(',')[1]),
                'tm': float(result.get(f'PRIMER_LEFT_{i}_TM', '0')),
                'gc_percent': float(result.get(f'PRIMER_LEFT_{i}_GC_PERCENT', '0')),
                'self_any': result.get(f'PRIMER_LEFT_{i}_SELF_ANY', '0'),
                'self_end': result.get(f'PRIMER_LEFT_{i}_SELF_END', '0'),
                'penalty': float(result.get(f'PRIMER_LEFT_{i}_PENALTY', '0'))
            },
            'right': {
                'sequence': result.get(f'PRIMER_RIGHT_{i}_SEQUENCE', ''),
                'start': int(result.get(f'PRIMER_RIGHT_{i}', '0,0').split(',')[0]),
                'length': int(result.get(f'PRIMER_RIGHT_{i}', '0,0').split(',')[1]),
                'tm': float(result.get(f'PRIMER_RIGHT_{i}_TM', '0')),
                'gc_percent': float(result.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', '0')),
                'self_any': result.get(f'PRIMER_RIGHT_{i}_SELF_ANY', '0'),
                'self_end': result.get(f'PRIMER_RIGHT_{i}_SELF_END', '0'),
                'penalty': float(result.get(f'PRIMER_RIGHT_{i}_PENALTY', '0'))
            },
        }
        primers.append(primer_pair)
    
    return {
        'primers': primers,
        'num_returned': num_returned,
        'explain': {
            'left': result.get('PRIMER_LEFT_EXPLAIN', ''),
            'right': result.get('PRIMER_RIGHT_EXPLAIN', ''),
            'pair': result.get('PRIMER_PAIR_EXPLAIN', '')
        },
        'warnings': result.get('PRIMER_WARNING', ''),
        'errors': result.get('PRIMER_ERROR', '')
    }

def design_primers(input_files, reference_fasta, output_file, settings_file=None, 
                    primer3_exe="~/.conda/envs/salmon/bin/primer3_core", primer3_args="", quality_threshold="high", 
                    min_high=3, min_medium=2, max_low=0, flanking_size=150, target_length = 1,
                    max_variants=50,
                    error_log=None):
    """
    Design primers for variants using Primer3.
    
    Parameters:
    input_files (list): List of TSV files from screen_variants output.
    reference_fasta (str): Path to reference FASTA file.
    output_file (str): Path to output file for primer results.
    settings_file (str, optional): Path to Primer3 settings file.
    primer3_exe (str): Path to primer3_core executable.
    primer3_args (str): Additional command line arguments for Primer3.
    quality_threshold (str): Minimum quality threshold (high, medium, low).
    min_high (int): Minimum number of high quality variants.
    min_medium (int): Minimum number of medium quality variants.
    max_low (int): Maximum number of low quality variants.
    flanking_size (int): Size of flanking region on each side of variant.
    target_length (int): Length of the target region for primer design (default to 1 when SNPs).
    max_variants (int): Maximum number of variants to design primers for.
    error_log (str, optional): Path to error log file.
    
    Returns:
    int: Number of variants for which primers were successfully designed.
    """
    # Setup logging
    if error_log:
        error_handler = logging.FileHandler(error_log)
        error_handler.setLevel(logging.ERROR)
        logging.getLogger().addHandler(error_handler)
    
    logging.info(f"Starting primer design for variants from {len(input_files)} input files")
    
    # Default primer3 settings
    default_settings = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": "1",
        "PRIMER_PICK_RIGHT_PRIMER": "1",
        "PRIMER_PICK_INTERNAL_OLIGO": "0",
        "PRIMER_OPT_SIZE": "20",
        "PRIMER_MIN_SIZE": "18",
        "PRIMER_MAX_SIZE": "25",
        "PRIMER_OPT_TM": "60.0",
        "PRIMER_MIN_TM": "57.0",
        "PRIMER_MAX_TM": "63.0",
        "PRIMER_MIN_GC": "20.0",
        "PRIMER_MAX_GC": "80.0",
        "PRIMER_MAX_POLY_X": "5",
        "PRIMER_NUM_RETURN": "5",
        "PRIMER_PRODUCT_SIZE_RANGE": "100-300",
        "PRIMER_EXPLAIN_FLAG": "1"
    }
    
    # Load settings file if provided
    if settings_file:
        with open(settings_file, 'r') as f:
            for line in f:
                if '=' in line and line.startswith('PRIMER_'):
                    key, value = line.strip().split('=', 1)
                    default_settings[key] = value
    
    # List to store all primer results
    all_results = []
    
    # Count of successfully designed primers
    successful_designs = 0
    
    # Process each input file
    for input_file in input_files:
        logging.info(f"Processing input file: {input_file}")
        
        try:
            # Read input file
            df = pd.read_csv(input_file, sep='\t', low_memory=False)

            # Ensure POS is numeric for flanking calculations
            df['POS'] = pd.to_numeric(df['POS'], errors='coerce')

            if df.empty:
                logging.warning(f"Empty file: {input_file}")
                continue
            
            # Filter variants based on quality thresholds
            if quality_threshold == "high":
                df_filtered = df[df['overall_reliability'] == 'high']
            elif quality_threshold == "medium":
                df_filtered = df[df['overall_reliability'].isin(['high', 'medium'])]
            else:  # low or any other value
                df_filtered = df
            
            # Further filter based on counts
            high_count = (df_filtered['overall_reliability'] == 'high').sum()
            medium_count = (df_filtered['overall_reliability'] == 'medium').sum()
            low_count = (df_filtered['overall_reliability'] == 'low').sum()
            
            if high_count < min_high:
                logging.warning(f"Not enough high quality variants in {input_file}: {high_count} < {min_high}")
                continue
                
            if medium_count < min_medium:
                logging.warning(f"Not enough medium quality variants in {input_file}: {medium_count} < {min_medium}")
                continue
                
            if low_count > max_low:
                logging.warning(f"Too many low quality variants in {input_file}: {low_count} > {max_low}")
                continue
            
            # Filter only primer compliant variants
            df_primer_compliant = df_filtered[df_filtered['primer_compliant'] == True]
            if df_primer_compliant.empty:
                logging.warning(f"No primer compliant variants in {input_file}")
                continue
            
            # Order the variants by QUAL

            df_primer_compliant = df_primer_compliant.sort_values(by='QUAL', ascending=False)

            # Limit number of variants
            print("DEBUGGING: TILL HERE WITHOUT A PROBLEM")
            # if max_variants != -9:
            #     if len(df_primer_compliant) > max_variants and max_variants:
            #         logging.warning(f"Limiting from {len(df_primer_compliant)} to {max_variants} variants in {input_file}")
            #         df_primer_compliant = df_primer_compliant.head(max_variants)
            
            # Process each variant
            for _, variant in tqdm(df_primer_compliant.iterrows(), 
                                desc=f"Designing primers for variants in {os.path.basename(input_file)}",
                                total=len(df_primer_compliant)):
                
                chrom = variant['CHROM']
                pos = target_pos = variant['POS']
                ref = variant['REF']
                alt = variant['ALT']
                
                # Calculate region to extract
                start = max(1, pos - flanking_size)
                end = pos + flanking_size
                
                # Extract sequence
                sequence = extract_sequence(reference_fasta, chrom, start, end)
                if not sequence:
                    logging.error(f"Failed to extract sequence for {chrom}:{pos}")
                    continue
                
                # Create temporary input file
                with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as temp_input:
                    temp_input_path = temp_input.name
                
                # Create primer3 input
                try:
                    create_primer3_input(f"M_{chrom}_{pos}", sequence, target_pos, target_length, default_settings, temp_input_path)
                except Exception as e:
                    logging.error(f"Failed to create Primer3 input for {chrom}:{pos}: {e}")
                    os.unlink(temp_input_path)  # Clean up
                    continue
                
                # Run primer3
                primer3_output = run_primer3(temp_input_path, settings_file, primer3_args)
                os.unlink(temp_input_path)  # Clean up
                
                if not primer3_output:
                    logging.error(f"Failed to run Primer3 for {chrom}:{pos}")
                    continue
                
                # Parse primer3 output
                parsed_output = parse_primer3_output(primer3_output)
                if parsed_output['num_returned'] == 0:
                    logging.warning(f"No primers found for {chrom}:{pos}")
                    continue
                
                # Add variant info to results
                result = {
                    'chrom': chrom,
                    'position': pos,
                    'ref': ref,
                    'alt': alt,
                    'region_start': start,
                    'region_end': end,
                    'sequence': sequence,
                    'target_position': target_pos,
                    'reliability': variant['overall_reliability'],
                    'primer_results': parsed_output
                }
                
                all_results.append(result)
                successful_designs += 1
                
        except Exception as e:
            logging.error(f"Error processing file {input_file}: {e}")
            continue
    
    # Write results to output file
    try:
        with open(output_file, 'w') as f:
            # Write header
            f.write("CHROM\tPOS\tREF\tALT\tReliability\tPrimer_Rank\t")
            f.write("Left_Primer\tLeft_Start\tLeft_TM\tLeft_GC\t")
            f.write("Right_Primer\tRight_Start\tRight_TM\tRight_GC\t")
            f.write("Product_Size\tPenalty\n")
            
            # Write primer data
            for result in all_results:
                chrom = result['chrom']
                pos = result['position']
                ref = result['ref']
                alt = result['alt']
                reliability = result['reliability']
                
                for i, primer in enumerate(result['primer_results']['primers']):
                    f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{reliability}\t{i+1}\t")
                    f.write(f"{primer['left']['sequence']}\t{result['region_start'] + primer['left']['start']}\t{primer['left']['tm']}\t{primer['left']['gc_percent']}\t")
                    f.write(f"{primer['right']['sequence']}\t{result['region_start'] + primer['right']['start'] - primer['right']['length'] + 1}\t{primer['right']['tm']}\t{primer['right']['gc_percent']}\t")
                    f.write(f"{primer['product_size']}\t{primer['pair_penalty']}\n")
                    
        logging.info(f"Successfully designed primers for {successful_designs} variants out of {sum(len(pd.read_csv(file, sep='\t')) for file in input_files)} total variants")
        logging.info(f"Results written to {output_file}")
    except Exception as e:
        logging.error(f"Error writing output file: {e}")
    
    return successful_designs