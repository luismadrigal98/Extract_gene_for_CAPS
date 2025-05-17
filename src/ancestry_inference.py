"""
This script is going to determine, based on the genotyping data of the F2s, the allele identity of the parental lines.

@author: Luis Javier Madrigal-Roca & John K. Kelly

@date 2025/05/17

"""

import pandas as pd
import sys
import math

# Add local directory to system path

sys.path.append('.')
from src.masking_vcf import read_vcf

def calculate_likelihood(genotype_counts, p1_allele, p2_allele, quality_data=None, error_rate=0.05):
    """
    Calculate likelihood with quality-weighted error model
    
    Args:
        genotype_counts: Dict of counts for each genotype
        p1_allele, p2_allele: Parental alleles to test
        quality_data: Optional dict with quality info {genotype: (count, depth, qual)}
        error_rate: Base error rate to adjust with quality metrics
    """
    import math
    
    # Expected proportions (same as before)
    expected_proportions = {
        ('0', '0'): {'0/0': 1.0, '0/1': 0.0, '1/1': 0.0},
        ('0', '1'): {'0/0': 0.25, '0/1': 0.5, '1/1': 0.25},
        ('1', '0'): {'0/0': 0.25, '0/1': 0.5, '1/1': 0.25},
        ('1', '1'): {'0/0': 0.0, '0/1': 0.0, '1/1': 1.0}
    }
    
    key = (p1_allele, p2_allele)
    if key not in expected_proportions:
        return float('-inf')
        
    expected = expected_proportions[key]
    
    # Calculate log-likelihood with quality-adjusted error model
    log_likelihood = 0
    
    for genotype, count in genotype_counts.items():
        if genotype == './.':
            continue
            
        # Adjust error rate based on quality if available
        adj_error_rate = error_rate
        if quality_data and genotype in quality_data:
            # Higher depth and qual score = lower error rate
            _, depth, qual = quality_data[genotype]
            depth_factor = min(depth / 20, 1.0)  # Cap at depth of 20
            adj_error_rate = error_rate * (1 - depth_factor) * (1 - qual)
        
        # Calculate probability with adjusted error rate
        prob = 0
        for true_gt, true_prob in expected.items():
            if true_gt == genotype:
                prob += true_prob * (1 - adj_error_rate)
            else:
                prob += true_prob * (adj_error_rate / 2)
                
        if prob > 0:
            log_likelihood += count * math.log(prob)
        else:
            log_likelihood += count * -100
    
    return log_likelihood

def infer_parental_genotypes(genotype_counts, error_rate=0.05, quality_data=None):
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
        likelihood = calculate_likelihood(genotype_counts, p1, p2, quality_data, error_rate)
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

def extract_genotype(genotype_field, return_quality=False):
    """
    Extract genotype and quality metrics from VCF format field.
    
    Input formats can be:
    - Simple genotype: "0/0"
    - Complex field: "1/1:90,9,0:0,3" (GT:PL:AD format)
    - Missing data: "./."
    
    Args:
        genotype_field: The genotype field from VCF
        return_quality: If True, return quality metrics as well
    
    Returns:
        If return_quality=False: standardized genotype only
        If return_quality=True: tuple of (genotype, depth, qual_score)
    """
    if pd.isna(genotype_field):
        return "./." if not return_quality else ("./.", 0, 0)
        
    if isinstance(genotype_field, str):
        # Default values
        gt = genotype_field
        depth = 0
        qual = 0
        
        if ':' in genotype_field:
            parts = genotype_field.split(':')
            gt = parts[0]
            
            # Extract depth from AD field (typically 3rd field)
            if len(parts) >= 3 and ',' in parts[2]:
                try:
                    depths = list(map(int, parts[2].split(',')))
                    depth = sum(depths)  # Total read depth
                except ValueError:
                    pass
                    
            # Extract quality from PL field (typically 2nd field)
            if len(parts) >= 2 and ',' in parts[1]:
                try:
                    phred_scores = list(map(int, parts[1].split(',')))
                    # Lower PL value = higher confidence
                    if gt == '0/0' and len(phred_scores) > 0:
                        qual = 1 / (1 + phred_scores[0])
                    elif gt in ('0/1', '1/0') and len(phred_scores) > 1:
                        qual = 1 / (1 + phred_scores[1])
                    elif gt == '1/1' and len(phred_scores) > 2:
                        qual = 1 / (1 + phred_scores[2])
                except ValueError:
                    pass
            
        # Standardize genotype format
        std_gt = "./."
        if gt in ['0/0', '0|0']:
            std_gt = '0/0'
        elif gt in ['0/1', '1/0', '0|1', '1|0']:
            std_gt = '0/1'
        elif gt in ['1/1', '1|1']:
            std_gt = '1/1'
        
        return std_gt if not return_quality else (std_gt, depth, qual)
    
    return "./." if not return_quality else ("./.", 0, 0)

def infer_ancestry(vcf, ROI_list, ancestry_log, output, context_window=20):
    """
    Infer ancestry and parental alleles for genetic variants based on F2 segregation patterns.
    Process each ROI separately and calculate both individual variant and contextual group likelihoods.
    """
    # 1. Load relationship data
    ancestry_df = pd.read_csv(ancestry_log, header=0, sep='\t')
    
    # 2. Group F2 samples by their parental lines
    f2_groups = {}
    for _, row in ancestry_df.iterrows():
        if row['FC'] == 'F2':
            key = f"{row['Common']}_{row['Alt']}"
            if key not in f2_groups:
                f2_groups[key] = []
            f2_groups[key].append(row['ID'])
    
    # 3. Load VCF data
    vcf_df = read_vcf(vcf)
    
    # 4. Load ROI list
    try:
        roi_df = pd.read_csv(ROI_list, sep='\t')
    except:
        roi_df = pd.read_csv(ROI_list, delim_whitespace=True)
    
    # 5. Process each ROI separately
    for _, roi in roi_df.iterrows():
        roi_name = roi['ROI_name']
        chrom = roi['Chrom']
        start = roi['Start']
        end = roi['End']
        
        print(f"Processing ROI: {roi_name} ({chrom}:{start}-{end})")
        
        # Filter variants for this ROI
        roi_variants = vcf_df[(vcf_df['CHROM'] == chrom) & 
                                (vcf_df['POS'] >= start) & 
                                (vcf_df['POS'] <= end)].copy()
        
        if len(roi_variants) == 0:
            print(f"No variants found in ROI {roi_name}")
            continue
        
        # Process each variant individually
        results = []
        for idx, variant in roi_variants.iterrows():
            variant_record = {'CHROM': variant['CHROM'], 'POS': variant['POS'],
                            'REF': variant['REF'], 'ALT': variant['ALT']}
            
            # For each cross (e.g., 664c_767c)
            for cross, samples in f2_groups.items():
                common_parent, alt_parent = cross.split('_')
                
                # Extract genotypes for this variant
                genotype_counts = {'0/0': 0, '0/1': 0, '1/1': 0, './.': 0}
                quality_data = {'0/0': [], '0/1': [], '1/1': []}

                # Add to the genotype extraction section
                for sample in samples:
                    if sample in roi_variants.columns:
                        # Get genotype and quality metrics
                        gt, depth, qual = extract_genotype(variant[sample], return_quality=True)
                        genotype_counts[gt] += 1
                        if gt != './.':
                            quality_data[gt].append((1, depth, qual))  # Count, depth, qual

                # After collecting all genotype counts
                total_samples = sum(genotype_counts.values())
                missing_samples = genotype_counts['./.']
                missing_ratio = missing_samples / total_samples if total_samples > 0 else 1.0

                # Skip likelihood calculation if too much data is missing
                if missing_ratio > 0.8:  # You can adjust this threshold
                    print(f"Skipping variant at {variant['CHROM']}:{variant['POS']} for cross {cross}: {missing_ratio:.2f} missing data")
                    variant_record[f"{common_parent}_allele"] = "N"  # Indicate insufficient data
                    variant_record[f"{alt_parent}_allele"] = "N"
                    variant_record['likelihood'] = float('nan')
                    variant_record['confidence'] = 0
                else:
                    # Calculate average quality metrics for each genotype
                    avg_quality = {}
                    for gt, data in quality_data.items():
                        if data:  # Only if we have data for this genotype
                            # Calculate averages from collected (count, depth, qual) tuples
                            total_depth = sum(d for _, d, _ in data)
                            total_qual = sum(q for _, _, q in data)
                            count = len(data)
                            avg_quality[gt] = (count, total_depth/count if count > 0 else 0, total_qual/count if count > 0 else 0)
                    
                    # Only calculate inference if we have enough data
                    inference = infer_parental_genotypes(genotype_counts, quality_data=avg_quality)
                    
                    # Store results
                    variant_record[f"{common_parent}_allele"] = inference['p1_allele']
                    variant_record[f"{alt_parent}_allele"] = inference['p2_allele']
                    variant_record['likelihood'] = inference['log_likelihood']
                    variant_record['confidence'] = inference['confidence']

                # Always store genotype counts for reference
                variant_record[f'hom_ref_count'] = genotype_counts['0/0']
                variant_record[f'het_count'] = genotype_counts['0/1']
                variant_record[f'hom_alt_count'] = genotype_counts['1/1']
                variant_record[f'missing_count'] = genotype_counts['./.']
                variant_record[f'missing_ratio'] = missing_ratio
                
            results.append(variant_record)
        
        # Now add contextual analysis
        results_with_context = []
        for i, result in enumerate(results):
            # Define the context window (looking at neighboring variants)
            start_idx = max(0, i - context_window//2)
            end_idx = min(len(results) - 1, i + context_window//2)
            context_variants = results[start_idx:end_idx+1]
            
            # Skip the current variant from context
            context_variants = [v for v in context_variants if v['POS'] != result['POS']]
            
            # Calculate group likelihood
            for parent in set([k.split('_')[0] for k in result.keys() if k.endswith('_allele')]):
                # Count how many neighboring variants agree with this one
                parent_allele = result.get(f"{parent}_allele")
                if parent_allele:
                    matching = sum(1 for v in context_variants if v.get(f"{parent}_allele") == parent_allele)
                    total = sum(1 for v in context_variants if f"{parent}_allele" in v)
                    result[f"{parent}_context_agreement"] = matching/total if total > 0 else None
                    
                    # Flag potential errors
                    if matching/total < 0.25 and total >= 5:
                        result[f"{parent}_potential_error"] = True
            
            results_with_context.append(result)
        
        # Write output for this ROI with organized columns
        roi_output = f"{output}_{roi_name}.tsv"

        # Define column ordering for better readability
        variant_cols = ['CHROM', 'POS', 'REF', 'ALT']
        parental_cols = [col for col in results_with_context[0].keys() if col.endswith('_allele') and not col.endswith('context_agreement')]
        likelihood_cols = ['likelihood', 'confidence']
        count_cols = [col for col in results_with_context[0].keys() if any(x in col for x in ['count', 'context_agreement', 'potential_error'])]

        # Order columns in a logical way
        ordered_cols = variant_cols + parental_cols + likelihood_cols + count_cols

        # Write the ordered DataFrame
        pd.DataFrame(results_with_context)[ordered_cols].to_csv(roi_output, sep='\t', index=False)
        print(f"Results for ROI {roi_name} saved to {roi_output}")

        # Create simplified results by excluding certain columns or selecting only essential ones
        simplified_results = []
        for result in results_with_context:
            # Create a simplified version of each result with only the essential information
            simplified_result = {
                'CHROM': result['CHROM'],
                'POS': result['POS'],
                'REF': result['REF'],
                'ALT': result['ALT'],
                'confidence': result['confidence']
            }
            
            # Add all parental allele columns
            for col in result.keys():
                if col.endswith('_allele') and not col.endswith('context_agreement'):
                    simplified_result[col] = result[col]
            
            simplified_results.append(simplified_result)

        # For the simplified output, maintain consistent ordering
        simplified_output = f"{output}_simplified_{roi_name}.tsv"
        simplified_cols = ['CHROM', 'POS', 'REF', 'ALT'] + [col for col in simplified_results[0].keys() 
                                                          if col not in ['CHROM', 'POS', 'REF', 'ALT', 'confidence']] + ['confidence']
        pd.DataFrame(simplified_results)[simplified_cols].to_csv(simplified_output, sep='\t', index=False)
        print(f"Simplified results for ROI {roi_name} saved to {simplified_output}")