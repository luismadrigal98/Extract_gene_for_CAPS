"""
This script will contain utility functions for validating the action of primers in silico.

@author: Luis Javier Madrigal-Roca & John K. Kelly

@date 2025/05/16

"""

from Bio.Blast import NCBICommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import pandas as pd
import logging
import os
import subprocess

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

def blast_primers(primer_fasta, db_path, output_file, evalue=10, task="blastn-short"):
    """Run BLAST to find primer binding sites"""
    logger.info(f"BLASTing primers against {db_path}")
    
    cmd = [
        "blastn", 
        "-query", primer_fasta,
        "-db", db_path,
        "-out", output_file,
        "-outfmt", "5",  # XML output
        "-evalue", str(evalue),
        "-task", task
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Error running BLAST: {result.stderr}")
        raise RuntimeError(f"blastn failed with exit code {result.returncode}")
        
    logger.info(f"BLAST results written to {output_file}")

def read_blast_results(blast_output, db_name):
    """
    Read and parse BLAST results from an XML file.
    
    Parameters:
    blast_output (str): Path to the BLAST output XML file.
    db_name (str): Name of the BLAST database.
    
    Returns:
    list: Parsed BLAST results.
    """
    
    from Bio.Blast import NCBIXML
    
    logger.info(f"Reading BLAST results from {blast_output}")
    
    with open(blast_output) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        results = {}
        results['raw_out'] = []
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    results['raw_out'].append({
                        'query': blast_record.query,
                        'hit_id': alignment.hit_id,
                        'hit_def': alignment.hit_def,
                        'e_value': hsp.expect,
                        'identity': hsp.identities,
                        'length': hsp.align_length
                    })
    
    # Extract the start and end positions of the primers for future extraction of in silico amplicon
    for i, result in enumerate(results['raw_out']):
        if 'primer' in result['query']:
            start = int(result['query'].split('_')[1])
            end = start + result['length']
            results['raw_out'][i]['start'] = start
            results['raw_out'][i]['end'] = end
    
    # Infer the amplicon coordinates
    # The amplicon is located between the end of the left primer and the start of the right primer

    for i, result in enumerate(results['raw_out']):
        if 'left' in result['query']:
            results['raw_out'][i]['amplicon_start'] = result['end']
        elif 'right' in result['query']:
            results['raw_out'][i]['amplicon_end'] = result['start']

    amplicon_coordinates = []
    for i, result in enumerate(results['raw_out']):
        if 'amplicon_start' in result and 'amplicon_end' in result:
            amplicon_coordinates.append({
                'amplicon_start': result['amplicon_start'],
                'amplicon_end': result['amplicon_end']
            })

    logger.info(f"Parsed {len(results)} BLAST results")
    return amplicon_coordinates

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

def validate_primers(primers_file, genomes, output_file, temp_dir=None, keep_temp=False):
    """
    Validate primers by BLASTing against target genomes and analyzing potential amplicons.
    
    Parameters:
    primers_file (str): File with designed primers
    genomes (list): List of target genome FASTA files
    output_file (str): Output file for validation results
    temp_dir (str, optional): Directory to store temporary files
    keep_temp (bool): Whether to keep temporary files
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
                
                # BLAST both primers
                left_out = os.path.join(primer_dir, f"{genome_name}_left_blast.xml")
                right_out = os.path.join(primer_dir, f"{genome_name}_right_blast.xml")
                
                blast_primers(left_fasta, db_path, left_out)
                blast_primers(right_fasta, db_path, right_out)
                
                # Parse results
                left_hits = read_blast_results(left_out)
                right_hits = read_blast_results(right_out)
                
                # Check if primers are specific (one hit per primer)
                if len(left_hits) == 1 and len(right_hits) == 1:
                    left_hit = left_hits[0]
                    right_hit = right_hits[0]
                    
                    # Check if they're on the same chromosome
                    if left_hit['hit_id'] == right_hit['hit_id']:
                        # Determine amplicon coordinates based on primer orientations
                        # For simplicity, let's assume left primer binds to forward strand
                        # and right primer binds to reverse strand
                        if left_hit['strand'] == 'plus' and right_hit['strand'] == 'minus':
                            amplicon_start = left_hit['sbjct_end'] + 1  # After left primer
                            amplicon_end = right_hit['sbjct_end'] + 1   # Before right primer
                            
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