"""
This script is going to determine, based on the genotyping data of the F2s, the allele identity of the parental lines.

@author: Luis Javier Madrigal-Roca & John K. Kelly

@date 2025/05/17

"""

import pandas as pd
import sys

# Add local directory to system path

sys.path.append('.')
from src.masking_vcf import read_vcf

def calculate_likelihood(genotype_counts, p1_allele, p2_allele, error_rate=0.05):
    """
    Calculate likelihood of parental genotype combination using F2 segregation patterns
    
    Parameters:
    -----------
    genotype_counts: dict
        Counts of genotypes observed in F2s {'0/0': n1, '0/1': n2, '1/1': n3}
    p1_allele: str
        Inferred allele for parent 1 ('0' or '1')
    p2_allele: str
        Inferred allele for parent 2 ('0' or '1')
    error_rate: float
        Estimated genotyping error rate (default: 0.05)
        
    Returns:
    --------
    float: Log-likelihood score
    """
    import math
    
    # Expected F2 proportions based on parental genotypes
    expected_proportions = {
        # Parent1/Parent2 genotype: {expected F2 genotype proportions}
        ('0', '0'): {'0/0': 1.0, '0/1': 0.0, '1/1': 0.0},  # Both parents homozygous ref
        ('0', '1'): {'0/0': 0.25, '0/1': 0.5, '1/1': 0.25},  # One parent homozygous ref, one homozygous alt
        ('1', '0'): {'0/0': 0.25, '0/1': 0.5, '1/1': 0.25},  # One parent homozygous alt, one homozygous ref
        ('1', '1'): {'0/0': 0.0, '0/1': 0.0, '1/1': 1.0}   # Both parents homozygous alt
    }
    
    # Get expected proportions for this parental combination
    key = (p1_allele, p2_allele)
    if key not in expected_proportions:
        return float('-inf')  # Invalid combination
        
    expected = expected_proportions[key]
    
    # Calculate log-likelihood with error model
    total_count = sum(count for gt, count in genotype_counts.items() if gt != './.')
    if total_count == 0:
        return float('-inf')
    
    log_likelihood = 0
    for genotype, count in genotype_counts.items():
        if genotype == './.':
            continue
        
        # Apply error model
        prob = 0
        for true_gt, true_prob in expected.items():
            if true_gt == genotype:
                # Correct genotype called
                prob += true_prob * (1 - error_rate)
            else:
                # Error in genotype calling
                prob += true_prob * (error_rate / 2)  # Divide by 2 as there are 2 possible error states
                
        if prob > 0:
            log_likelihood += count * math.log(prob)
        else:
            log_likelihood += count * -100  # Very low probability for zero cases
    
    return log_likelihood

def infer_parental_genotypes(genotype_counts, error_rate=0.05):
    """
    Infer most likely parental genotype combination using maximum likelihood
    """
    # All possible parental allele combinations for biallelic site
    parent_combinations = [
        ('0', '0'),  # Both parents homozygous reference
        ('0', '1'),  # Parent 1 reference, Parent 2 alternate
        ('1', '0'),  # Parent 1 alternate, Parent 2 reference
        ('1', '1')   # Both parents homozygous alternate
    ]
    
    # Calculate likelihood for each combination
    likelihoods = {}
    for p1, p2 in parent_combinations:
        likelihood = calculate_likelihood(genotype_counts, p1, p2, error_rate)
        likelihoods[(p1, p2)] = likelihood
    
    # Find combination with maximum likelihood
    max_key = max(likelihoods, key=likelihoods.get)
    max_likelihood = likelihoods[max_key]
    
    # Calculate relative likelihood (how much better is the best model vs others)
    total_likelihood = sum(math.exp(l - max_likelihood) for l in likelihoods.values())
    confidence = 1.0 / total_likelihood
    
    return {
        'p1_allele': max_key[0],
        'p2_allele': max_key[1],
        'log_likelihood': max_likelihood,
        'confidence': confidence
    }

def infer_block_ancestry(vcf_df, samples, chrom, start, end, error_rate=0.05):
    """
    Infer parental ancestry for a block of variants in a region like an inversion
    """
    # Extract variants in the region
    region_variants = vcf_df[(vcf_df['CHROM'] == chrom) & 
                            (vcf_df['POS'] >= start) & 
                            (vcf_df['POS'] <= end)]
    
    # Collect genotype counts for all variants in region
    all_counts = []
    for _, variant in region_variants.iterrows():
        genotype_counts = {'0/0': 0, '0/1': 0, '1/1': 0, './.': 0}
        
        for sample in samples:
            if sample in vcf_df.columns:
                gt = extract_genotype(variant[sample])
                genotype_counts[gt] = genotype_counts.get(gt, 0) + 1
        
        all_counts.append(genotype_counts)
    
    # Calculate joint likelihood across all variants
    parent_combinations = [('0', '0'), ('0', '1'), ('1', '0'), ('1', '1')]
    joint_likelihoods = {combo: 0 for combo in parent_combinations}
    
    for counts in all_counts:
        for combo in parent_combinations:
            p1, p2 = combo
            likelihood = calculate_likelihood(counts, p1, p2, error_rate)
            joint_likelihoods[combo] += likelihood
    
    # Find maximum likelihood combination
    max_combo = max(joint_likelihoods, key=joint_likelihoods.get)
    max_likelihood = joint_likelihoods[max_combo]
    
    return {
        'p1_allele': max_combo[0],
        'p2_allele': max_combo[1],
        'log_likelihood': max_likelihood,
        'variant_count': len(all_counts)
    }

def infer_ancestry(vcf, ROI_list, ancestry_log, output, estimate_likelihood=False):
    """
    Infer ancestry and parental alleles for genetic variants based on F2 segregation patterns.
    This function analyzes genotype patterns in F2 populations to infer the most likely
    alleles present in their parental lines. It processes a VCF file, groups samples by
    their parental crosses, and estimates which allele each parent contributed based on
    segregation patterns in the F2 generation.
    Parameters
    ----------
    vcf : str
        Path to the VCF file containing genotype data for all samples
    ROI_list : list or str
        List of regions of interest or path to file containing regions of interest
    ancestry_log : str
        Path to tab-separated file containing sample relationship data (must include
        columns: 'ID', 'FC', 'Common', 'Alt')
    output : str
        Path where the output file will be saved (tab-separated format)
    estimate_likelihood : bool, optional
        Whether to calculate likelihood of the ancestry inference (default: False)
    Returns
    -------
    None
        Results are written to the specified output file containing:
        - Variant information (chromosome, position, reference, alternative alleles)
        - Inferred alleles for each parent
        - Likelihood of inference (if requested)
    Notes
    -----
    The ancestry log file must contain at least four columns:
    - 'ID': Sample identifier
    - 'FC': Generation code (F1, F2, etc.)
    - 'Common': Common parent identifier
    - 'Alt': Alternative parent identifier
    The inference works best with adequate sample sizes for each F2 group.

    """
    
    # 1. Load relationship data
    ancestry_df = pd.read_csv(ancestry_log, sep='\t')
    
    # 2. Group F2 samples by their parental lines
    f2_groups = {}
    for _, row in ancestry_df.iterrows():
        if row['FC'] == 'F2':
            key = f"{row['Common']}_{row['Alt']}"
            if key not in f2_groups:
                f2_groups[key] = []
            f2_groups[key].append(row['ID'])
    
    # 3. Load VCF and extract variants in ROI
    vcf_df = read_vcf(vcf)
    
    # 4. For each variant, calculate allele frequencies per F2 group
    results = []
    for _, variant in vcf_df.iterrows():
        variant_results = {'CHROM': variant['CHROM'], 'POS': variant['POS'], 
                          'REF': variant['REF'], 'ALT': variant['ALT']}
        
        for cross, samples in f2_groups.items():
            common_parent, alt_parent = cross.split('_')
            
            # Count genotypes (0/0, 0/1, 1/1) across F2s in this group
            genotype_counts = {'0/0': 0, '0/1': 0, '1/1': 0, './.': 0}
            
            for sample in samples:
                if sample in vcf_df.columns:
                    gt = extract_genotype(variant[sample])  # Helper function needed
                    genotype_counts[gt] = genotype_counts.get(gt, 0) + 1
            
            # Simple inference: 
            # - If most F2s are homozygous ref: common parent is likely homozygous ref
            # - If most F2s are homozygous alt: common parent is likely homozygous alt
            # - If many F2s are het: parents likely have different genotypes
            
            ratio_het = genotype_counts['0/1'] / sum(v for k,v in genotype_counts.items() if k != './.')
            ratio_hom_alt = genotype_counts['1/1'] / sum(v for k,v in genotype_counts.items() if k != './.')
            
            # Infer parental genotypes and calculate likelihood if requested
            # (Simplified inference shown here)
            if ratio_het > 0.4:  # Expected ~0.5 for het
                common_allele = '0'
                alt_allele = '1'
                likelihood = calculate_likelihood(genotype_counts, '0', '1') if estimate_likelihood else None
            else:
                # More complex inference needed
                pass
                
            variant_results[f"{common_parent}_allele"] = common_allele
            variant_results[f"{alt_parent}_allele"] = alt_allele
            if estimate_likelihood:
                variant_results[f"likelihood"] = likelihood
                
        results.append(variant_results)
    
    # 5. Output results
    pd.DataFrame(results).to_csv(output, sep='\t', index=False)