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
                     min_identity_pct=90.0, min_coverage=80.0, check_3prime=True):
    """
    Validate primers by BLASTing against target genomes and analyzing potential amplicons.
    
    Parameters:
    primers_file (str): File with designed primers
    genomes (list): List of target genome FASTA files
    output_file (str): Output file for validation results
    temp_dir (str, optional): Directory to store temporary files
    keep_temp (bool): Whether to keep temporary files
    evalue (float): E-value threshold for BLAST
    task (str): BLAST task
    word_size (int): Word size for BLAST
    min_identity_pct (float): Minimum percent identity for valid binding
    min_coverage (float): Minimum percent of primer covered by alignment
    check_3prime (bool): Whether to check 3' end matches specifically
    """
    import tempfile
    import shutil
    import pandas as pd
    from Bio import SeqIO
    
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
    
    # Group by variant to get primer pairs
    validation_results = []
    all_amplicons_fasta = os.path.join(os.path.dirname(output_file), "all_amplicons.fasta")
    all_amplicons_file = open(all_amplicons_fasta, 'w')
    
    # Create BLAST databases for each genome
    genome_dbs = {}
    for genome_file in genomes:
        db_name = os.path.join(temp_dir, os.path.basename(genome_file))
        create_blast_db(genome_file, db_name)
        genome_dbs[genome_file] = db_name
    
    # Process each primer pair
    for (chrom, pos, ref, alt), group in primers_df.groupby(['CHROM', 'POS', 'REF', 'ALT']):
        if len(group) == 0:
            continue
            
        primer_id = f"{chrom}_{pos}_{ref}_{alt}"
        logger.info(f"Processing primer pair for {primer_id}")
        
        # Use the first primer pair
        first_row = group.iloc[0]
        left_primer = first_row['Left_Primer']
        right_primer = first_row['Right_Primer']
        
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
        
        amplicon_fasta = os.path.join(primer_dir, "amplicons.fasta")
        with open(amplicon_fasta, 'w') as primer_amplicons:
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
                
                # Check if primers are specific (one hit per primer)
                if len(left_hits) == 1 and len(right_hits) == 1:
                    left_hit = left_hits[0]
                    right_hit = right_hits[0]
                    
                    # Check if they're on the same chromosome
                    if left_hit['hit_id'] == right_hit['hit_id']:
                        # Check orientation - for proper PCR, we need:
                        # left primer on + strand AND right primer on - strand
                        # OR left primer on - strand AND right primer on + strand
                        if (left_hit['strand'] == 'plus' and right_hit['strand'] == 'minus') or \
                           (left_hit['strand'] == 'minus' and right_hit['strand'] == 'plus'):
                            
                            # Determine amplicon coordinates based on strand
                            if left_hit['strand'] == 'plus':
                                amplicon_start = left_hit['sbjct_end']  # End of left primer binding
                                amplicon_end = right_hit['sbjct_start'] # Start of right primer binding
                            else:
                                # Reversed orientation
                                amplicon_start = right_hit['sbjct_end']  # End of right primer binding
                                amplicon_end = left_hit['sbjct_start']   # Start of left primer binding
                            
                            # Make sure amplicon_start < amplicon_end
                            if amplicon_start > amplicon_end:
                                amplicon_start, amplicon_end = amplicon_end, amplicon_start
                            
                            # Extract the amplicon sequence
                            with open(genome_file, 'r') as genome_handle:
                                for record in SeqIO.parse(genome_handle, "fasta"):
                                    if record.id == left_hit['hit_id'].split(' ')[0]:
                                        amplicon_seq = record.seq[amplicon_start-1:amplicon_end-1]
                                        amplicon_len = len(amplicon_seq)
                                        
                                        # Write to primer-specific amplicon file
                                        primer_amplicons.write(f">{primer_id}_{genome_name}\n{amplicon_seq}\n")
                                        
                                        # Also write to master amplicon file
                                        all_amplicons_file.write(f">{primer_id}_{genome_name}\n{amplicon_seq}\n")
                                        
                                        primer_result['genomes'][genome_name] = {
                                            'specific': True,
                                            'amplicon_length': amplicon_len,
                                            'amplicon_seq': str(amplicon_seq)
                                        }
                                        continue
                        else:
                            logger.warning(f"Unexpected primer orientations for {primer_id} in {genome_name}")
                            primer_result['genomes'][genome_name] = {
                                'specific': False,
                                'reason': 'unexpected_orientation'
                            }
                    else:
                        logger.warning(f"Primers hit different chromosomes for {primer_id} in {genome_name}")
                        primer_result['genomes'][genome_name] = {
                            'specific': False,
                            'reason': 'different_chromosomes'
                        }
                else:
                    primer_result['genomes'][genome_name] = {
                        'specific': False,
                        'reason': 'multiple_hits',
                        'left_hits': len(left_hits),
                        'right_hits': len(right_hits)
                    }
        
        validation_results.append(primer_result)
    
    # Close master amplicon file
    all_amplicons_file.close()
    
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
                notes = details.get('reason', "") if not details.get('specific', False) else ""
                
                f.write(f"| {genome_name} | {specific} | {amplicon_length} | {notes} |\n")
            
            f.write("\n")
    
    logger.info(f"Validation results written to {output_file}")
    logger.info(f"All amplicon sequences written to {all_amplicons_fasta}")
    
    # Clean up temp directory if needed
    if not keep_temp and created_temp:
        shutil.rmtree(temp_dir)
        logger.info(f"Temporary directory {temp_dir} removed")
    else:
        logger.info(f"Temporary files kept in {temp_dir}")
    
    return validation_results