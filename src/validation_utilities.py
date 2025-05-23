"""
This script will contain utility functions for validating the action of primers in silico.

@author: Luis Javier Madrigal-Roca & John K. Kelly

@date 2025/05/16

"""

from Bio.Blast import NCBIXML
import pandas as pd
import logging
import os
import subprocess
import xml.parsers.expat
import tempfile
import shutil
from Bio import SeqIO

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def create_blast_db(fasta_file, db_name):
    """Create a BLAST database from a FASTA file"""
    logger.info(f"Creating BLAST database for {fasta_file}")
    cmd = ["makeblastdb", "-in", fasta_file, "-dbtype", "nucl", "-out", db_name]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Error creating BLAST database: {result.stderr}")
        raise RuntimeError(f"makeblastdb failed with exit code {result.returncode}")
        
    logger.info(f"Created BLAST database at {db_name}")

def write_primers_as_fasta(primers, output_file):
    """
    Write a list of primers to a FASTA file.
    
    Parameters:
    primers (list): List of primer sequences.
    output_file (str): Path to the output FASTA file.
    
    Returns:
    None
    """
    
    logger.info(f"Writing primers to {output_file}")
    
    with open(output_file, 'w') as f:
        for i, primer in enumerate(primers):
            f.write(f">primer_{i}\n{primer}\n")
    
    logger.info(f"Primers written to {output_file} successfully.")

def blast_primers(primer_fasta, db_path, output_file, evalue=0.1, task="blastn-short", word_size=7):
    """Run BLAST to find primer binding sites with more stringent parameters"""
    logger.info(f"BLASTing primers against {db_path}")
    
    cmd = [
        "blastn", 
        "-query", primer_fasta,
        "-db", db_path,
        "-out", output_file,
        "-outfmt", "5",  # XML output
        "-evalue", str(evalue),  # Ensure this is used correctly
        "-task", task,
        "-word_size", str(word_size),  # Word size appropriate for primers
        "-dust", "no",   # Don't filter low complexity regions
        "-max_target_seqs", "20"  # Limit maximum targets to reduce processing
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Error running BLAST: {result.stderr}")
        raise RuntimeError(f"blastn failed with exit code {result.returncode}")
        
    logger.info(f"BLAST results written to {output_file}")

def read_blast_results(blast_output, db_name=None, min_identity_pct=90, min_coverage=80, check_3prime=True):
    """
    Read and parse BLAST results, filtering for high-quality primer binding sites only.
    """
    logger.info(f"Reading BLAST results from {blast_output}")
    
    # Check if file exists and is not empty
    if not os.path.exists(blast_output) or os.path.getsize(blast_output) == 0:
        logger.warning(f"BLAST output file {blast_output} is empty or does not exist")
        return []
    
    hits = []
    try:
        with open(blast_output) as result_handle:
            try:
                blast_records = list(NCBIXML.parse(result_handle))
                
                if not blast_records:
                    logger.warning(f"No BLAST records found in {blast_output}")
                    return []
                
                total_hits = 0
                quality_hits = []
                
                for blast_record in blast_records:
                    query_length = blast_record.query_length
                    query_name = blast_record.query
                    
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            total_hits += 1
                            
                            # Calculate identity percentage and coverage
                            identity_pct = (hsp.identities / hsp.align_length) * 100
                            coverage = (hsp.align_length / query_length) * 100
                            
                            # Check if 3' end of primer matches (last 5bp)
                            three_prime_ok = True
                            if check_3prime:
                                if "left" in query_name.lower():
                                    three_prime_ok = (hsp.query[-5:] == hsp.sbjct[-5:])
                                else:
                                    three_prime_ok = (hsp.query[-5:] == hsp.sbjct[-5:])
                            
                            # Filter for quality hits only
                            if (identity_pct >= min_identity_pct and 
                                coverage >= min_coverage and
                                hsp.expect <= 0.1 and
                                three_prime_ok):
                                
                                is_forward = hsp.sbjct_start < hsp.sbjct_end
                                strand = 'plus' if is_forward else 'minus'
                                
                                quality_hits.append({
                                    'query': query_name,
                                    'hit_id': alignment.hit_id,
                                    'e_value': hsp.expect,
                                    'identity_pct': identity_pct,
                                    'coverage': coverage, 
                                    'strand': strand,
                                    'sbjct_start': hsp.sbjct_start,
                                    'sbjct_end': hsp.sbjct_end,
                                    'genome': os.path.basename(db_name) if db_name else "unknown"
                                })
                
                logger.info(f"Found {len(quality_hits)} high-quality binding sites out of {total_hits} total hits")
                return quality_hits
                
            except (xml.parsers.expat.ExpatError, ValueError) as e:
                logger.error(f"XML parsing error in {blast_output}: {e}")
                return []
    except Exception as e:
        logger.error(f"Error processing BLAST output file {blast_output}: {e}")
        return []

def extract_amplicon_sequences(amplicon_coordinates, fasta_file):
    """
    Extract amplicon sequences from the original FASTA file based on BLAST results.
    
    Parameters:
    amplicon_coordinates (list): List of amplicon coordinates.
    fasta_file (str): Path to the original FASTA file.
    
    Returns:
    list: Extracted amplicon sequences.
    """
    
    logger.info(f"Extracting amplicon sequences from {fasta_file}")
    
    # Read the original FASTA file
    from Bio import SeqIO
    
    amplicon_sequences = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        for coord in amplicon_coordinates:
            start = coord['amplicon_start']
            end = coord['amplicon_end']
            amplicon_seq = record.seq[start:end]
            amplicon_sequences.append(str(amplicon_seq))
    
    logger.info(f"Extracted {len(amplicon_sequences)} amplicon sequences")
    
    # Write the amplicon sequences to a new FASTA file
    output_file = os.path.splitext(fasta_file)[0] + "_amplicons.fasta"

    with open(output_file, 'w') as f:
        for i, seq in enumerate(amplicon_sequences):
            f.write(f">amplicon_{i}\n{seq}\n")
    
    return amplicon_sequences

def validate_primers(primers_file, genomes, output_file, temp_dir=None, keep_temp=False,
                        evalue=0.1, task="blastn-short", word_size=7,
                        min_identity_pct=90.0, min_coverage=80.0, check_3prime=True, 
                        specific_output=None):  # Added parameter for specific primers output
    """
    Validate primers by BLASTing against target genomes and analyzing potential amplicons.
    
    Parameters:
    primers_file (str): Path to the input file with primers.
    genomes (list): List of genome FASTA files to BLAST against.
    output_file (str): Path to the output file for validation results.
    temp_dir (str, optional): Temporary directory for intermediate files.
    keep_temp (bool, optional): Whether to keep the temporary files.
    evalue (float, optional): E-value threshold for BLAST.
    task (str, optional): BLAST task to use.
    word_size (int, optional): Word size for BLAST.
    min_identity_pct (float, optional): Minimum identity percentage for BLAST hits.
    min_coverage (float, optional): Minimum coverage percentage for BLAST hits.
    check_3prime (bool, optional): Whether to check 3' end matching.
    specific_output (str, optional): Path to output file for primers specific to all genomes.
    """
    
    # Create temp directory
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp(prefix="primer_validation_")
        created_temp = True
    else:
        os.makedirs(temp_dir, exist_ok=True)
        created_temp = False
        
    logger.info(f"Using temporary directory: {temp_dir}")
    
    # Load primers from file
    primers_df = pd.read_csv(primers_file, sep='\t')
    
    # Process each primer pair individually
    validation_results = []
    
    # Create BLAST databases for each genome
    genome_dbs = {}
    for genome_file in genomes:
        db_name = os.path.join(temp_dir, os.path.basename(genome_file))
        create_blast_db(genome_file, db_name)
        genome_dbs[genome_file] = db_name
    
    # Process each primer individually - no need to preload genomes since we're not extracting sequences
    for idx, row in primers_df.iterrows():
        # Extract primer information
        chrom = row['CHROM']
        pos = row['POS']
        ref = row['REF']
        alt = row['ALT']
        
        # Create a unique identifier with index
        primer_id = f"{chrom}_{pos}_{ref}_{alt}_id{idx}"
        
        logger.info(f"Processing primer pair: {primer_id}")
        
        # Use the first primer pair
        left_primer = row['Left_Primer']
        right_primer = row['Right_Primer']
        
        # Create primer-specific directory
        primer_dir = os.path.join(temp_dir, primer_id)
        os.makedirs(primer_dir, exist_ok=True)
        
        # Create FASTA files for primers
        left_fasta = os.path.join(primer_dir, "left.fasta")
        right_fasta = os.path.join(primer_dir, "right.fasta")
        
        with open(left_fasta, 'w') as f:
            f.write(f">left_{primer_id}\n{left_primer}\n")
        
        with open(right_fasta, 'w') as f:
            f.write(f">right_{primer_id}\n{right_primer}\n")
        
        primer_result = {
            'primer_id': primer_id,
            'chrom': chrom,
            'pos': pos,
            'ref': ref,
            'alt': alt,
            'left_primer': left_primer,
            'right_primer': right_primer,
            'genomes': {}
        }
        
        # Initialize results for all genomes first
        for genome_file in genomes:
            genome_name = os.path.basename(genome_file)
            primer_result['genomes'][genome_name] = {
                'specific': False,
                'reason': 'not_processed',
                'amplicon_length': 'N/A'
            }
        
        # Check each genome
        for genome_file in genomes:
            genome_name = os.path.basename(genome_file)
            db_path = genome_dbs[genome_file]
            
            # Define the output file paths for BLAST results
            left_out = os.path.join(primer_dir, f"{genome_name}_left_blast.xml")
            right_out = os.path.join(primer_dir, f"{genome_name}_right_blast.xml")
            
            # BLAST both primers
            blast_primers(left_fasta, db_path, left_out, evalue=evalue, task=task, word_size=word_size)
            blast_primers(right_fasta, db_path, right_out, evalue=evalue, task=task, word_size=word_size)
            
            # Parse results
            left_hits = read_blast_results(left_out, db_path, 
                                           min_identity_pct=min_identity_pct, 
                                           min_coverage=min_coverage, 
                                           check_3prime=check_3prime)
            right_hits = read_blast_results(right_out, db_path,
                                            min_identity_pct=min_identity_pct, 
                                            min_coverage=min_coverage, 
                                            check_3prime=check_3prime)
            
            # Detailed result reporting with clear diagnostics
            if len(left_hits) == 0 and len(right_hits) == 0:
                primer_result['genomes'][genome_name] = {
                    'specific': False,
                    'reason': 'no_hits_for_both_primers',
                    'amplicon_length': 'N/A'
                }
            elif len(left_hits) == 0:
                primer_result['genomes'][genome_name] = {
                    'specific': False,
                    'reason': 'no_hit_for_left_primer',
                    'amplicon_length': 'N/A'
                }
            elif len(right_hits) == 0:
                primer_result['genomes'][genome_name] = {
                    'specific': False,
                    'reason': 'no_hit_for_right_primer',
                    'amplicon_length': 'N/A'
                }
            elif len(left_hits) > 1:
                primer_result['genomes'][genome_name] = {
                    'specific': False,
                    'reason': 'multiple_hits_for_left_primer',
                    'amplicon_length': 'N/A'
                }
            elif len(right_hits) > 1:
                primer_result['genomes'][genome_name] = {
                    'specific': False,
                    'reason': 'multiple_hits_for_right_primer',
                    'amplicon_length': 'N/A'
                }
            else:
                # Process single hits case
                left_hit = left_hits[0]
                right_hit = right_hits[0]
                
                # Check if they hit same chromosome
                if left_hit['hit_id'] != right_hit['hit_id']:
                    primer_result['genomes'][genome_name] = {
                        'specific': False,
                        'reason': 'hits_on_different_chromosomes',
                        'amplicon_length': 'N/A'
                    }
                    continue
                    
                # Check orientation for proper amplification
                is_valid_orientation = False
                amplicon_size = 'N/A'
                
                if left_hit['strand'] == 'plus' and right_hit['strand'] == 'minus':
                    # Normal orientation: left → right
                    is_valid_orientation = True
                    amplicon_start = min(left_hit['sbjct_start'], right_hit['sbjct_start'])
                    amplicon_end = max(left_hit['sbjct_end'], right_hit['sbjct_end'])
                    amplicon_size = amplicon_end - amplicon_start + 1
                elif left_hit['strand'] == 'minus' and right_hit['strand'] == 'plus':
                    # Reverse orientation
                    is_valid_orientation = True
                    amplicon_start = min(left_hit['sbjct_start'], right_hit['sbjct_end']) 
                    amplicon_end = max(left_hit['sbjct_start'], right_hit['sbjct_end'])
                    amplicon_size = amplicon_end - amplicon_start + 1
                
                if is_valid_orientation:
                    # Just record the results without writing amplicon files
                    primer_result['genomes'][genome_name] = {
                        'specific': True,
                        'amplicon_length': amplicon_size,
                        'left_pos': f"{left_hit['sbjct_start']}-{left_hit['sbjct_end']} ({left_hit['strand']})",
                        'right_pos': f"{right_hit['sbjct_start']}-{right_hit['sbjct_end']} ({right_hit['strand']})"
                    }
                else:
                    primer_result['genomes'][genome_name] = {
                        'specific': False,
                        'reason': 'invalid_primer_orientation',
                        'amplicon_length': 'N/A'
                    }
        
        validation_results.append(primer_result)
    
    # Write validation summary
    with open(output_file, 'w') as f:
        f.write("# Primer Validation Results\n\n")
        
        for result in validation_results:
            f.write(f"## Primer {result['primer_id']}\n")
            f.write(f"- Left: {result['left_primer']}\n")
            f.write(f"- Right: {result['right_primer']}\n\n")
            
            f.write("| Genome | Specific | Amplicon Length | Notes |\n")
            f.write("|--------|----------|----------------|-------|\n")
            
            for genome_name, details in result['genomes'].items():
                specific = "Yes" if details.get('specific', False) else "No"
                amplicon_length = details.get('amplicon_length', "N/A")
                
                if details.get('specific', False):
                    notes = f"Valid amplicon {amplicon_length}bp"
                else:
                    notes = details.get('reason', "Unknown issue")
                
                f.write(f"| {genome_name} | {specific} | {amplicon_length} | {notes} |\n")
            
            f.write("\n")
    
    # NEW CODE: Output primers that are specific to all genomes
    if specific_output:
        # Find primers that are specific to all genomes
        specific_primers = []
        for result in validation_results:
            all_specific = True
            amplicon_lengths = {}
            
            # Check if primer is specific for all genomes
            for genome_name, details in result['genomes'].items():
                if not details.get('specific', False):
                    all_specific = False
                    break
                amplicon_lengths[genome_name] = details.get('amplicon_length', 'N/A')
            
            # If specific for all genomes, add to our list
            if all_specific:
                primer_data = {
                    'CHROM': result['chrom'],
                    'POS': result['pos'],
                    'REF': result['ref'],
                    'ALT': result['alt'],
                    'Left_Primer': result['left_primer'],
                    'Right_Primer': result['right_primer'],
                    'Primer_ID': result['primer_id']
                }
                
                # Add amplicon lengths for each genome
                for genome_name, length in amplicon_lengths.items():
                    primer_data[f'Amplicon_Length_{genome_name}'] = length
                
                specific_primers.append(primer_data)
        
        # Write specific primers to file if any found
        if specific_primers:
            logger.info(f"Writing {len(specific_primers)} specific primers to {specific_output}")
            
            # Convert to DataFrame and write to file
            specific_df = pd.DataFrame(specific_primers)
            
            # Order columns logically
            columns = ['Primer_ID', 'CHROM', 'POS', 'REF', 'ALT', 'Left_Primer', 'Right_Primer']
            # Add amplicon length columns
            for genome_name in genomes:
                genome_base = os.path.basename(genome_name)
                if f'Amplicon_Length_{genome_base}' in specific_df.columns:
                    columns.append(f'Amplicon_Length_{genome_base}')
            
            # Write to TSV file
            specific_df.to_csv(specific_output, sep='\t', index=False, columns=columns)
            logger.info(f"Specific primers ready for ordering saved to {specific_output}")
        else:
            logger.warning("No primers were specific to all genomes. No specific_output file created.")
    
    # Clean up temp directory if needed
    if not keep_temp and created_temp:
        shutil.rmtree(temp_dir)
        logger.info(f"Temporary directory {temp_dir} removed")
    else:
        logger.info(f"Temporary files kept in {temp_dir}")
    
    return validation_results