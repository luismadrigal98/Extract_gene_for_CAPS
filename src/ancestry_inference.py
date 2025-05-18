"""
This script is going to determine, based on the genotyping data of the F2s, the allele identity of the parental lines.

@author: Luis Javier Madrigal-Roca & John K. Kelly

@date 2025/05/17

"""

import pandas as pd
import sys
import math
from tqdm import tqdm
import logging

# Set up logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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
    - Missing data: "./." or ./.:0,0,0:0,0
    
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

def infer_ancestry_multiple(vcf, ROI_list, ancestry_log, output, context_window=20, approach='single', check_parentals = False):
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
        
        # Before processing variants, identify all common parents across crosses
        all_parents = set()
        common_parents = set()
        for cross in f2_groups.keys():
            parents = cross.split('_')
            for p in parents:
                if p in all_parents:
                    common_parents.add(p)
                all_parents.add(p)

        # For each variant
        for idx, variant in roi_variants.iterrows():
            variant_record = {'CHROM': variant['CHROM'], 'POS': variant['POS'],
                            'REF': variant['REF'], 'ALT': variant['ALT']}
            
            # Create dictionary to store cross-specific data
            cross_data = {}
            
            # Store cross-specific allele predictions temporarily
            cross_predictions = {}
            
            # Store confidence values for later weighted averaging
            parent_confidences = {p: [] for p in all_parents}
            parent_alleles = {p: [] for p in all_parents}
            
            # Process each cross separately first
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
                if missing_ratio > 0.95:  # You can adjust this threshold
                    print(f"Skipping variant at {variant['CHROM']}:{variant['POS']} for cross {cross}: {missing_ratio:.2f} missing data")
                    cross_predictions[cross] = {common_parent: "N", alt_parent: "N"}
                    cross_data[cross] = {"confidence": 0, "missing_ratio": missing_ratio}
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
                    cross_predictions[cross] = {
                        common_parent: inference['p1_allele'], 
                        alt_parent: inference['p2_allele']
                    }
                    cross_data[cross] = {
                        "confidence": inference['confidence'],
                        "likelihood": inference['log_likelihood']
                    }
                    
                    # Record predictions with their confidence for each parent
                    parent_confidences[common_parent].append((inference['p1_allele'], inference['confidence']))
                    parent_confidences[alt_parent].append((inference['p2_allele'], inference['confidence']))
                    parent_alleles[common_parent].append(inference['p1_allele'])
                    parent_alleles[alt_parent].append(inference['p2_allele'])
            
            # Reconcile predictions for common parents (like 664c)
            final_alleles = {}
            for parent in all_parents:
                if parent in common_parents:
                    # For common parents, use confidence-weighted consensus
                    if not parent_confidences[parent]:
                        final_alleles[parent] = "N"
                    else:
                        # Use most confident prediction that isn't "N"
                        valid_predictions = [(allele, conf) for allele, conf in parent_confidences[parent] if allele != "N"]
                        if valid_predictions:
                            final_alleles[parent] = max(valid_predictions, key=lambda x: x[1])[0]
                        else:
                            final_alleles[parent] = "N"
                else:
                    # For unique parents, use their single prediction
                    if not parent_confidences[parent]:
                        final_alleles[parent] = "N" 
                    else:
                        # Use most confident prediction
                        final_alleles[parent] = max(parent_confidences[parent], key=lambda x: x[1])[0]
            
            # Store final reconciled alleles
            for parent, allele in final_alleles.items():
                variant_record[f"{parent}_allele"] = allele
                
            # Store cross-specific data
            for cross, data in cross_data.items():
                for key, value in data.items():
                    variant_record[f"{cross}_{key}"] = value

            # Before the cross loop, add cumulative counters
            total_genotype_counts = {'0/0': 0, '0/1': 0, '1/1': 0, './.': 0}
            total_samples_checked = 0

            # Before processing crosses, create a dict to track alleles per parent
            parental_alleles = {}

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
                if missing_ratio > 0.95:  # You can adjust this threshold
                    print(f"Skipping variant at {variant['CHROM']}:{variant['POS']} for cross {cross}: {missing_ratio:.2f} missing data")
                    # Store in the parent-specific dictionary instead of overwriting
                    parental_alleles[common_parent] = "N"
                    parental_alleles[alt_parent] = "N"
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
                    parental_alleles[common_parent] = inference['p1_allele']
                    parental_alleles[alt_parent] = inference['p2_allele']
                    variant_record['likelihood'] = inference['log_likelihood']
                    variant_record['confidence'] = inference['confidence']

                # Always store genotype counts for reference
                variant_record[f'hom_ref_count'] = genotype_counts['0/0']
                variant_record[f'het_count'] = genotype_counts['0/1']
                variant_record[f'hom_alt_count'] = genotype_counts['1/1']
                variant_record[f'missing_count'] = genotype_counts['./.']
                variant_record[f'missing_ratio'] = missing_ratio

                # Update total counts
                total_genotype_counts['0/0'] += genotype_counts['0/0']
                total_genotype_counts['0/1'] += genotype_counts['0/1']
                total_genotype_counts['1/1'] += genotype_counts['1/1'] 
                total_genotype_counts['./.'] += genotype_counts['./.']
                total_samples_checked += sum(genotype_counts.values())
                
                # Store per-cross counts if needed
                variant_record[f'{cross}_hom_ref'] = genotype_counts['0/0']
                # ... other cross-specific fields ...
                
            # After the cross loop, transfer all parental alleles to variant_record
            for parent, allele in parental_alleles.items():
                variant_record[f"{parent}_allele"] = allele

            # After the cross loop, store total counts
            variant_record[f'hom_ref_count'] = total_genotype_counts['0/0']
            variant_record[f'het_count'] = total_genotype_counts['0/1']
            variant_record[f'hom_alt_count'] = total_genotype_counts['1/1']
            variant_record[f'missing_count'] = total_genotype_counts['./.']
            variant_record[f'missing_ratio'] = total_genotype_counts['./.'] / total_samples_checked if total_samples_checked > 0 else 1.0
                
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

def infer_ancestry_single(vcf, ROI_list, ancestry_log, output):
    """
    Infer parental alleles for single F2 individuals.
    
    Args:
        vcf: Path to VCF file
        ROI_list: Path to regions of interest file
        ancestry_log: Path to relationship map file
        output: Base name for output files
    """
    # Read input files
    vcf_df = read_vcf(vcf)
    ancestry_df = pd.read_csv(ancestry_log, header=0, sep='\t')
    
    try:
        roi_df = pd.read_csv(ROI_list, sep='\t')
    except:
        roi_df = pd.read_csv(ROI_list, delim_whitespace=True)
    
    # Get F2 samples and their parents
    f2_samples = {}
    for _, row in ancestry_df.iterrows():
        if row['FC'] == 'F2':
            f2_samples[str(row['ID'])] = (row['Common'], row['Alt'])
    
    # Process each ROI
    for _, roi in roi_df.iterrows():
        roi_name = roi['ROI_name'].split('_')[0].strip()
        chrom = roi['Chrom']
        start = roi['Start']
        end = roi['End']
        
        logging.info(f"Processing ROI: {roi_name} ({chrom}:{start}-{end})")
        
        # Filter variants in this ROI
        roi_variants = vcf_df[(vcf_df['CHROM'] == chrom) & 
                             (vcf_df['POS'] >= start) & 
                             (vcf_df['POS'] <= end)].copy()
        
        if len(roi_variants) == 0:
            logging.info(f"No variants found in ROI {roi_name}")
            continue
        
        # Create DataFrame to store results
        results = []
        
        # Process each variant position
        for _, variant in tqdm(roi_variants.iterrows(), desc=f"Processing variants in {roi_name}"):
            variant_record = {
                'CHROM': variant['CHROM'],
                'POS': variant['POS'],
                'REF': variant['REF'],
                'ALT': variant['ALT']
            }
            
            # Store evidence from each F2 sample
            sample_evidence = {}
            
            # Process each F2 sample
            for sample_id, (common_parent, alt_parent) in f2_samples.items():
                if sample_id not in roi_variants.columns:
                    continue
                    
                # Get F2 and parental genotypes
                f2_gt = extract_genotype(variant[sample_id])
                common_gt = extract_genotype(variant[common_parent]) if common_parent in roi_variants.columns else './.'
                alt_gt = extract_genotype(variant[alt_parent]) if alt_parent in roi_variants.columns else './.'
                
                # Skip if F2 genotype is missing
                if f2_gt == './.':
                    continue
                
                # Initialize evidence record
                evidence = {
                    'common_allele': None,
                    'alt_allele': None,
                    'confidence': 0.0,
                    'support_type': 'None'
                }
                
                # Infer parental alleles based on F2 genotype
                if f2_gt == '0/0':  # F2 is homozygous reference
                    if common_gt == '0/0' and alt_gt == '0/0':
                        # Both parents are reference - strong evidence
                        evidence['common_allele'] = '0'
                        evidence['alt_allele'] = '0'
                        evidence['confidence'] = 0.9
                        evidence['support_type'] = 'Assembly/Complete'
                    elif common_gt == '0/0' and alt_gt == './.':
                        # Only common parent genotyped - moderate evidence
                        evidence['common_allele'] = '0'
                        evidence['alt_allele'] = '0'  # Inferred
                        evidence['confidence'] = 0.7
                        evidence['support_type'] = 'Assembly/Partial'
                    elif common_gt == './.' and alt_gt == '0/0':
                        # Only alt parent genotyped - moderate evidence
                        evidence['common_allele'] = '0'  # Inferred
                        evidence['alt_allele'] = '0'
                        evidence['confidence'] = 0.7
                        evidence['support_type'] = 'Assembly/Partial'
                    else:
                        # F2 doesn't match parental genotypes - weak inference
                        evidence['common_allele'] = '0'
                        evidence['alt_allele'] = '0'
                        evidence['confidence'] = 0.5
                        evidence['support_type'] = 'F2_evidence'
                
                elif f2_gt == '0/1':  # F2 is heterozygous
                    # In an inversion, heterozygous F2s typically have one allele from each parent
                    evidence['common_allele'] = '0'
                    evidence['alt_allele'] = '1'
                    evidence['confidence'] = 0.8
                    evidence['support_type'] = 'F2_segregation'
                    
                    # Check if assemblies support this
                    if common_gt == '0/0' and alt_gt == '1/1':
                        evidence['confidence'] = 0.95
                        evidence['support_type'] = 'Assembly/Complete'
                
                elif f2_gt == '1/1':  # F2 is homozygous alternate
                    if common_gt == '1/1' and alt_gt == '1/1':
                        # Both parents are alternate - strong evidence
                        evidence['common_allele'] = '1'
                        evidence['alt_allele'] = '1'
                        evidence['confidence'] = 0.9
                        evidence['support_type'] = 'Assembly/Complete'
                    elif common_gt == '1/1' and alt_gt == './.':
                        # Only common parent genotyped - moderate evidence
                        evidence['common_allele'] = '1'
                        evidence['alt_allele'] = '1'  # Inferred
                        evidence['confidence'] = 0.7
                        evidence['support_type'] = 'Assembly/Partial'
                    elif common_gt == './.' and alt_gt == '1/1':
                        # Only alt parent genotyped - moderate evidence
                        evidence['common_allele'] = '1'  # Inferred
                        evidence['alt_allele'] = '1'
                        evidence['confidence'] = 0.7
                        evidence['support_type'] = 'Assembly/Partial'
                    else:
                        # F2 doesn't match parental genotypes - weak inference
                        evidence['common_allele'] = '1'
                        evidence['alt_allele'] = '1'
                        evidence['confidence'] = 0.5
                        evidence['support_type'] = 'F2_evidence'
                
                # Store evidence from this sample
                if evidence['common_allele'] is not None:
                    key = f"{sample_id}_{common_parent}_{alt_parent}"
                    sample_evidence[key] = evidence
            
            # Combine evidence across samples for this variant
            if sample_evidence:
                # Get weighted consensus for common parent
                common_votes = {'0': 0.0, '1': 0.0}
                alt_votes = {'0': 0.0, '1': 0.0}
                
                for evidence in sample_evidence.values():
                    if evidence['common_allele'] in common_votes:
                        common_votes[evidence['common_allele']] += evidence['confidence']
                    if evidence['alt_allele'] in alt_votes:
                        alt_votes[evidence['alt_allele']] += evidence['confidence']
                
                # Determine final alleles
                common_allele = max(common_votes, key=common_votes.get) if sum(common_votes.values()) > 0 else 'N'
                alt_allele = max(alt_votes, key=alt_votes.get) if sum(alt_votes.values()) > 0 else 'N'
                confidence = max(sum(common_votes.values()), sum(alt_votes.values())) / len(sample_evidence)
                
                # Store in variant record
                for parent_id in set([p for s, (p, _) in f2_samples.items()] + [p for s, (_, p) in f2_samples.items()]):
                    if parent_id in [p for s, (p, _) in f2_samples.items()]:
                        # This is a common parent
                        variant_record[f"{parent_id}_allele"] = common_allele
                    else:
                        # This is an alt parent
                        variant_record[f"{parent_id}_allele"] = alt_allele
                
                variant_record['confidence'] = confidence
                variant_record['sample_count'] = len(sample_evidence)
            else:
                # No evidence for this variant
                for parent_id in set([p for s, (p, _) in f2_samples.items()] + [p for s, (_, p) in f2_samples.items()]):
                    variant_record[f"{parent_id}_allele"] = 'N'
                variant_record['confidence'] = 0.0
                variant_record['sample_count'] = 0
            
            results.append(variant_record)
        
        # Save results for this ROI
        if results:
            output_file = f"{output}_{roi_name}.tsv"
            pd.DataFrame(results).to_csv(output_file, sep='\t', index=False)
            logging.info(f"Saved results for ROI {roi_name} to {output_file}")
        else:
            logging.info(f"No results for ROI {roi_name}")
