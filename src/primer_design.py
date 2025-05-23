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
import sys
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import random
import string
import concurrent.futures
# Local imports
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from primer_contrast import select_best_primers
from primer_design_parallel import process_variant_for_parallel

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

def run_primer3(input_file, primer3_exe='~/.conda/envs/salmon/bin/primer3_core', settings_file=None, primer3_args="", also_get_formatted=False, timeout=300):
    """
    Run Primer3 on the input file.
    
    Parameters:
    input_file (str): Path to Primer3 input file
    primer3_exe (str): Path to primer3_core executable
    settings_file (str, optional): Path to Primer3 settings file
    primer3_args (str): Additional arguments for Primer3
    also_get_formatted (bool): If True, also run with --format_output and return both outputs
    timeout (int): Maximum time in seconds to wait for primer3 execution
    
    Returns:
    str or tuple: If also_get_formatted is False, returns Boulder IO output.
                    If True, returns (boulder_output, formatted_output) tuple.
    """
    # First run: Get Boulder IO output (for parsing)
    cmd = [primer3_exe]
    
    # Add command line arguments
    if primer3_args:
        # Make sure --format_output is not in the arguments
        args_list = primer3_args.split()
        args_list = [arg for arg in args_list if arg != '--format_output' and arg != '-format_output']
        cmd.extend(args_list)
    
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
                check=True,
                timeout=timeout  # Add timeout parameter
            )
        boulder_output = result.stdout
        
        # If formatted output is not requested, return just the Boulder IO output
        if not also_get_formatted:
            return boulder_output
            
        # Second run: Get formatted output (for human reading)
        format_cmd = cmd.copy()
        format_cmd.append("--format_output")
        
        with open(input_file, 'r') as f:
            format_result = subprocess.run(
                format_cmd,
                stdin=f,
                text=True,
                capture_output=True,
                check=True,
                timeout=timeout  # Add timeout to formatted run too
            )
        formatted_output = format_result.stdout
        
        # Return both outputs
        return (boulder_output, formatted_output)
        
    except subprocess.TimeoutExpired:
        variant_id = os.path.basename(input_file).replace("primer3_input.txt", "")
        logging.error(f"Primer3 process timed out after {timeout} seconds for variant {variant_id}")
        return None
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
                    min_high=3, min_medium=2, max_low=0, flanking_size=150, target_length=1,
                    max_variants=50, keep_temp=False, temp_dir=None, error_log=None,
                    contrast=False, num_primers=50, selection_criteria="balanced", selected_output=None,
                    parallel=False, num_workers=None, timeout=300):
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
    target_length (int): Length of the target region for primer design.
    max_variants (int): Maximum number of variants to design primers for.
    keep_temp (bool): Whether to keep temporary files for inspection.
    temp_dir (str, optional): Directory to store temporary files (created if not exists).
    error_log (str, optional): Path to error log file.
    contrast (bool): Whether to perform contrast selection of primers.
    num_primers (int): Number of primers to select.
    selection_criteria (str): Criteria for selecting primers (e.g., "balanced", "specificity").
    selected_output (str, optional): Path to output file for selected primers.
    parallel (bool): Whether to use parallel processing for primer design.
    num_workers (int, optional): Number of worker processes to use. Default: CPU count-1
    timeout (int): Maximum time in seconds to wait for primer3 execution. Default: 300
    
    Returns:
    int: Number of variants for which primers were successfully designed.
    """
    # Setup logging
    if error_log:
        error_handler = logging.FileHandler(error_log)
        error_handler.setLevel(logging.ERROR)
        logging.getLogger().addHandler(error_handler)
    
    logging.info(f"Starting primer design for variants from {len(input_files)} input files")
    
    # Create temp directory if keeping files
    if keep_temp:
        if temp_dir is None:
            temp_dir = os.path.join(os.getcwd(), "primer3_temp_files")
        
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
            logging.info(f"Created temporary directory: {temp_dir}")
        else:
            logging.info(f"Using existing temporary directory: {temp_dir}")
    
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
    
    # Create random ID generator for temp files in parallel mode
    def random_id(length=6):
        return ''.join(random.choices(string.ascii_lowercase + string.digits, k=length))
    
    # Configure parallel processing
    if parallel:
        if num_workers is None:
            num_workers = max(1, multiprocessing.cpu_count() - 1)
        logging.info(f"Parallel processing enabled with {num_workers} workers")
    
    # Define function for processing a single variant
    def process_variant(variant_data, process_id=None, primer3_exe=primer3_exe):
        # Unpack variant data
        idx, variant = variant_data
        
        # Use unique ID for temp files in parallel mode
        variant_id = f"{process_id}_{idx}" if process_id else str(idx)
        
        # Extract the variant info
        chrom = variant['CHROM']
        pos = variant['POS']
        ref = variant['REF']
        alt = variant['ALT']
        variant_id_str = f"{chrom}_{pos}_{ref}_{alt}"
        
        # Calculate region to extract
        start = max(1, pos - flanking_size)
        end = pos + flanking_size
        
        # Calculate the target position
        target_pos = pos - start

        # Extract sequence
        sequence = extract_sequence(reference_fasta, chrom, start, end)
        if not sequence:
            logging.error(f"Failed to extract sequence for {chrom}:{pos}")
            return None
        
        # Create a variant-specific directory for temp files if keeping
        if keep_temp:
            variant_dir = os.path.join(temp_dir, f"{chrom}_{pos}")
            if not os.path.exists(variant_dir):
                os.makedirs(variant_dir)
            input_file_path = os.path.join(variant_dir, "primer3_input.txt")
            output_file_path = os.path.join(variant_dir, "primer3_output.txt")
        else:
            # Use temporary files that will be deleted
            temp_input = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt')
            input_file_path = temp_input.name
            temp_input.close()
        
        # Create primer3 input
        try:
            create_primer3_input(f"M_{chrom}_{pos}", sequence, target_pos, target_length, default_settings, input_file_path)
        except Exception as e:
            logging.error(f"Failed to create Primer3 input for {chrom}:{pos}: {e}")
            if not keep_temp:
                os.unlink(input_file_path)
            return None
        
        # Run primer3
        # Expand the directory of the executable
        primer3_exe = os.path.expanduser(primer3_exe)
        
        # Run primer3 with both output formats
        primer3_result = run_primer3(
            input_file=input_file_path, 
            primer3_exe=primer3_exe, 
            settings_file=settings_file, 
            primer3_args=primer3_args,
            also_get_formatted=keep_temp,  # Use formatted output when keeping temp files
            timeout=timeout  # Pass the timeout parameter
        )

        # Handle the result based on output format
        if keep_temp and primer3_result and isinstance(primer3_result, tuple):
            boulder_output, formatted_output = primer3_result
            
            # Save both outputs
            with open(output_file_path, 'w') as f:
                f.write(formatted_output)  # Save human-readable output to file
                
            # Parse the boulder output for further processing
            parsed_output = parse_primer3_output(boulder_output)
            
            # Check if primer3 was successful - use boulder_output in this branch
            if not boulder_output:
                logging.error(f"Failed to run Primer3 for {chrom}:{pos}")
                return None
        else:
            # Just boulder output
            primer3_output = primer3_result
            
            # Save output if keeping temp files
            if keep_temp and primer3_output:
                with open(output_file_path, 'w') as f:
                    f.write(primer3_output)
            
            parsed_output = parse_primer3_output(primer3_output)
            
            # Check if primer3 was successful - use primer3_output in this branch
            if not primer3_output:
                logging.error(f"Failed to run Primer3 for {chrom}:{pos}")
                return None

        if parsed_output['num_returned'] == 0:
            logging.warning(f"No primers found for {chrom}:{pos}")
            return None
        
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
        
        return result
    
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
            if max_variants != -9:
                if len(df_primer_compliant) > max_variants and max_variants:
                    logging.warning(f"Limiting from {len(df_primer_compliant)} to {max_variants} variants in {input_file}")
                    df_primer_compliant = df_primer_compliant.head(max_variants)
            
            if parallel:
                # Process in parallel using batches with as_completed
                from concurrent.futures import as_completed
                
                results = []
                variant_items = [(i, row) for i, row in df_primer_compliant.iterrows()]
                process_ids = [random_id() for _ in range(len(variant_items))]
                
                # Use batch size equal to number of workers for optimal resource use
                batch_size = num_workers * 2  # Use 2x workers for better throughput
                
                with tqdm(total=len(variant_items), desc=f"Designing primers for {os.path.basename(input_file)} (parallel)") as pbar:
                    for i in range(0, len(variant_items), batch_size):
                        batch_variants = variant_items[i:i+batch_size]
                        batch_process_ids = process_ids[i:i+batch_size]
                        
                        # Create args for this batch
                        batch_args = []
                        for j, variant_item in enumerate(batch_variants):
                            args = (variant_item, reference_fasta, flanking_size, target_length, 
                                    default_settings, keep_temp, temp_dir, primer3_exe, settings_file, 
                                    primer3_args, batch_process_ids[j], timeout)  # Add timeout here
                            batch_args.append(args)
                        
                        # Process this batch
                        with ProcessPoolExecutor(max_workers=num_workers) as executor:
                            # Submit jobs for this batch
                            future_to_args = {executor.submit(process_variant_for_parallel, args): args 
                                                for args in batch_args}
                            
                            # Process results as they complete
                            for future in as_completed(future_to_args):
                                try:
                                    result = future.result(timeout=600)  # Add overall future timeout too
                                    if result:
                                        results.append(result)
                                        successful_designs += 1
                                except concurrent.futures.TimeoutError:
                                    # Handle case where entire future times out
                                    logging.error(f"Future execution timed out after 600 seconds")
                                except Exception as e:
                                    # Handle other exceptions
                                    logging.error(f"Error in parallel processing: {str(e)}")
                                finally:
                                    pbar.update(1)
        
                # After parallel processing, add results to the main all_results list
                all_results.extend(results)
                logging.info(f"Added {len(results)} results from parallel processing to main results list")
            else:
                # Original sequential processing with optimization
                if not parallel:
                    batch_size = 10  # Process in small batches to maintain progress visibility
                    variant_items = [(i, row) for i, row in df_primer_compliant.iterrows()]
                    
                    with tqdm(total=len(variant_items), desc=f"Designing primers for {os.path.basename(input_file)}") as pbar:
                        for i in range(0, len(variant_items), batch_size):
                            batch_variants = variant_items[i:i+batch_size]
                            
                            # Process this small batch sequentially
                            for variant_item in batch_variants:
                                result = process_variant(variant_item)
                                if result:
                                    all_results.append(result)
                                    successful_designs += 1
                                pbar.update(1)
    
        except Exception as e:
            logging.error(f"Error processing file {input_file}: {e}")
            continue
    
    # Write all results to output file
    try:
        # Add debug information
        logging.info(f"Writing {len(all_results)} primer designs to {output_file}")
        
        with open(output_file, 'w') as f:
            # Write header
            f.write("CHROM\tPOS\tREF\tALT\tReliability\tPrimer_Rank\t")
            f.write("Left_Primer\tLeft_Start\tLeft_TM\tLeft_GC\t")
            f.write("Right_Primer\tRight_Start\tRight_TM\tRight_GC\t")
            f.write("Product_Size\tPenalty\n")
            
            # Check if all_results contains any entries
            if not all_results:
                logging.error("all_results list is empty - no primer data to write")
                return 0
                
            # Write primer data for all designs
            for result_idx, result in enumerate(all_results):
                if result_idx % 50 == 0:  # Log progress for large files
                    logging.info(f"Writing primer {result_idx+1}/{len(all_results)}")
                
                chrom = result['chrom']
                pos = result['position']
                ref = result['ref']
                alt = result['alt']
                reliability = result['reliability']
                
                for primer_idx, primer in enumerate(result['primer_results']['primers']):
                    f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{reliability}\t{primer_idx+1}\t")
                    f.write(f"{primer['left']['sequence']}\t{result['region_start'] + primer['left']['start']}\t{primer['left']['tm']}\t{primer['left']['gc_percent']}\t")
                    f.write(f"{primer['right']['sequence']}\t{result['region_start'] + primer['right']['start'] - primer['right']['length'] + 1}\t{primer['right']['tm']}\t{primer['right']['gc_percent']}\t")
                    f.write(f"{primer['product_size']}\t{primer['pair_penalty']}\n")
                    
                    # Flush the file buffer periodically to ensure data is written to disk
                    if (result_idx * len(result['primer_results']['primers']) + primer_idx) % 100 == 0:
                        f.flush()
        
        logging.info(f"Successfully wrote all {len(all_results)} primer designs to file")
        
    except Exception as e:
        logging.error(f"Error writing to output file: {str(e)}", exc_info=True)
        # Try to write to an alternative file to recover data
        try:
            alt_output = output_file + ".recovery"
            logging.info(f"Attempting to write data to recovery file: {alt_output}")
            with open(alt_output, 'w') as f:
                import json
                json.dump([{k: str(v) if not isinstance(v, (str, int, float, bool, list, dict)) else v 
                            for k, v in r.items()} for r in all_results], f)
            logging.info(f"Recovery data written to {alt_output}")
        except Exception as e2:
            logging.error(f"Failed to write recovery data: {str(e2)}")
    
    # Post-process and select primers if contrast is enabled
    if contrast and all_results:
        logging.info(f"Selecting top {num_primers} primers using '{selection_criteria}' criteria")
        selected_primers = select_best_primers(all_results, num_primers, selection_criteria)
        
        # Write selected primers to specified output file if provided
        if selected_output and selected_primers:
            with open(selected_output, 'w') as f:
                # Write header with extra column for score
                f.write("Rank\tCHROM\tPOS\tREF\tALT\tReliability\tComposite_Score\t")
                f.write("Left_Primer\tLeft_Start\tLeft_TM\tLeft_GC\t")
                f.write("Right_Primer\tRight_Start\tRight_TM\tRight_GC\t")
                f.write("Product_Size\tPenalty\n")
                
                # Write data for selected primers
                for rank, result in enumerate(selected_primers, 1):
                    primer = result['selected_primer']
                    f.write(f"{rank}\t{result['chrom']}\t{result['position']}\t{result['ref']}\t{result['alt']}\t")
                    f.write(f"{result['reliability']}\t{result['composite_score']:.4f}\t")
                    f.write(f"{primer['left']['sequence']}\t{result['region_start'] + primer['left']['start']}\t{primer['left']['tm']}\t{primer['left']['gc_percent']}\t")
                    f.write(f"{primer['right']['sequence']}\t{result['region_start'] + primer['right']['start'] - primer['right']['length'] + 1}\t{primer['right']['tm']}\t{primer['right']['gc_percent']}\t")
                    f.write(f"{primer['product_size']}\t{primer['pair_penalty']}\n")
            
            logging.info(f"Selected {len(selected_primers)} primers written to {selected_output}")
    
    if successful_designs > 0:
        logging.info(f"Successfully designed primers for {successful_designs} variants out of {sum(len(pd.read_csv(file, sep='\t')) for file in input_files)} total variants")
        logging.info(f"Results written to {output_file}")
    else:
        logging.warning("No primers were successfully designed. Try again with different parameters.")
    
    return successful_designs