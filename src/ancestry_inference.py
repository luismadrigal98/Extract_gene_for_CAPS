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
    Infer parental alleles for single F2 individuals with reliability indicators.
    
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
            
        # STEP 1: Calculate region-wide patterns for each F2 (for context)
        f2_regional_patterns = {}
        
        for sample_id in f2_samples:
            if sample_id not in roi_variants.columns:
                continue
                
            # Count genotypes across the region
            genotype_counts = {'0/0': 0, '0/1': 0, '1/1': 0, './.': 0}
            depth_data = {'0/0': [], '0/1': [], '1/1': []}
            
            for _, variant in roi_variants.iterrows():
                if sample_id in variant:
                    gt, depth, _ = extract_genotype(variant[sample_id], return_quality=True)
                    genotype_counts[gt] += 1
                    if gt != './.':
                        depth_data[gt].append(depth)
            
            # Calculate zygosity ratios and average depths
            total = sum([genotype_counts['0/0'], genotype_counts['0/1'], genotype_counts['1/1']])
            if total > 0:
                # Get ratios
                hom_ref_ratio = genotype_counts['0/0'] / total
                het_ratio = genotype_counts['0/1'] / total
                hom_alt_ratio = genotype_counts['1/1'] / total
                
                # Get average depths
                avg_depths = {}
                for gt in depth_data:
                    if depth_data[gt]:
                        avg_depths[gt] = sum(depth_data[gt]) / len(depth_data[gt])
                    else:
                        avg_depths[gt] = 0
                
                # Determine predominant pattern
                if max(hom_ref_ratio, het_ratio, hom_alt_ratio) < 0.6:
                    pattern = "mixed"
                elif hom_ref_ratio > het_ratio and hom_ref_ratio > hom_alt_ratio:
                    pattern = "homozygous_ref"
                elif het_ratio > hom_ref_ratio and het_ratio > hom_alt_ratio:
                    pattern = "heterozygous"
                else:
                    pattern = "homozygous_alt"
                
                # Store pattern info
                f2_regional_patterns[sample_id] = {
                    'pattern': pattern,
                    'hom_ref_ratio': hom_ref_ratio,
                    'het_ratio': het_ratio,
                    'hom_alt_ratio': hom_alt_ratio,
                    'avg_depths': avg_depths
                }
                
                logging.info(f"Sample {sample_id} shows {pattern} pattern across region "
                             f"(ref:{hom_ref_ratio:.2f} het:{het_ratio:.2f} alt:{hom_alt_ratio:.2f})")
            else:
                f2_regional_patterns[sample_id] = {'pattern': 'unknown'}
        
        # Track which parents are common across families
        common_parents = {}
        alt_parents = {}
        
        # Identify common parents and alternative parents
        for sample_id, (common, alt) in f2_samples.items():
            if common in common_parents:
                common_parents[common].append(sample_id)
            else:
                common_parents[common] = [sample_id]
                
            if alt in alt_parents:
                alt_parents[alt].append(sample_id)
            else:
                alt_parents[alt] = [sample_id]
        
        # Create list to store results
        results = []
        
        # STEP 2: Process each variant position
        for _, variant in tqdm(roi_variants.iterrows(), desc=f"Processing variants in {roi_name}"):
            variant_record = {
                'CHROM': variant['CHROM'],
                'POS': variant['POS'],
                'REF': variant['REF'],
                'ALT': variant['ALT']
            }
            
            # Store raw F2 genotypes for reference
            raw_f2_data = {}
            for sample_id in f2_samples:
                if sample_id in variant:
                    gt, depth, _ = extract_genotype(variant[sample_id], return_quality=True)
                    # Store original genotype and depth
                    raw_f2_data[f"{sample_id}_genotype"] = gt
                    raw_f2_data[f"{sample_id}_depth"] = depth
                    
                    # Store allelic depths if available
                    if ':' in variant[sample_id]:
                        fields = variant[sample_id].split(':')
                        if len(fields) >= 3 and ',' in fields[2]:
                            try:
                                allelic_depths = [int(d) for d in fields[2].split(',')]
                                if len(allelic_depths) >= 2:
                                    raw_f2_data[f"{sample_id}_ref_depth"] = allelic_depths[0]
                                    raw_f2_data[f"{sample_id}_alt_depth"] = allelic_depths[1]
                                    
                                    # Flag potential issues with heterozygous calls
                                    if gt == '0/1':
                                        if depth <= 2:
                                            raw_f2_data[f"{sample_id}_quality"] = "very_low_depth"
                                        elif allelic_depths[0] > 3*allelic_depths[1]:
                                            raw_f2_data[f"{sample_id}_quality"] = "skewed_ref"
                                        elif allelic_depths[1] > 3*allelic_depths[0]:
                                            raw_f2_data[f"{sample_id}_quality"] = "skewed_alt"
                                        else:
                                            raw_f2_data[f"{sample_id}_quality"] = "balanced"
                                    else:
                                        if depth <= 2:
                                            raw_f2_data[f"{sample_id}_quality"] = "low_depth"
                                        else:
                                            raw_f2_data[f"{sample_id}_quality"] = "good"
                            except ValueError:
                                raw_f2_data[f"{sample_id}_quality"] = "unknown"
            
            # Process common parents - these need consensus across families
            for common_parent, samples in common_parents.items():
                evidence_for_parent = []
                raw_evidence = []
                
                for sample_id in samples:
                    if sample_id not in roi_variants.columns:
                        continue
                    
                    # Get sample family information
                    common_p, alt_p = f2_samples[sample_id]
                    
                    # Skip if this sample doesn't involve this common parent
                    if common_p != common_parent:
                        continue
                    
                    # Extract genotype and quality metrics
                    f2_gt, f2_depth, f2_qual = extract_genotype(variant[sample_id], return_quality=True)
                    
                    # Skip if F2 genotype is missing
                    if f2_gt == './.':
                        continue
                    
                    # Get regional pattern
                    region_pattern = f2_regional_patterns.get(sample_id, {}).get('pattern', 'unknown')
                    
                    # Calculate confidence level based on depth and consistency with regional pattern
                    confidence_level = "high"
                    if f2_depth == 1:
                        confidence_level = "very_low"
                    elif f2_depth == 2:
                        confidence_level = "low"
                    elif f2_depth <= 4:
                        confidence_level = "medium"
                        
                    # Check allelic depths for heterozygous calls
                    allelic_imbalance = False
                    if f2_gt == '0/1' and ':' in variant[sample_id]:
                        fields = variant[sample_id].split(':')
                        if len(fields) >= 3 and ',' in fields[2]:
                            try:
                                allelic_depths = [int(d) for d in fields[2].split(',')]
                                if len(allelic_depths) >= 2:
                                    if allelic_depths[0] > 3*allelic_depths[1] or allelic_depths[1] > 3*allelic_depths[0]:
                                        allelic_imbalance = True
                                        if confidence_level != "very_low":
                                            confidence_level = "low"  # Downgrade confidence for imbalanced het calls
                            except ValueError:
                                pass
                    
                    # Check regional consistency
                    if region_pattern != 'unknown' and region_pattern != 'mixed':
                        if (region_pattern == 'homozygous_ref' and f2_gt != '0/0') or \
                           (region_pattern == 'homozygous_alt' and f2_gt != '1/1') or \
                           (region_pattern == 'heterozygous' and f2_gt != '0/1'):
                            # Call conflicts with regional pattern
                            if confidence_level not in ["very_low", "low"]:
                                confidence_level = "low"  # Downgrade confidence for pattern conflicts
                    
                    # Store raw evidence
                    raw_evidence.append({
                        'sample': sample_id,
                        'genotype': f2_gt, 
                        'depth': f2_depth,
                        'allelic_imbalance': allelic_imbalance,
                        'confidence': confidence_level,
                        'region_pattern': region_pattern
                    })
                    
                    # Calculate confidence weight (0.1-1.0)
                    confidence_weight = {
                        "very_low": 0.1,
                        "low": 0.3,
                        "medium": 0.7,
                        "high": 1.0
                    }[confidence_level]
                    
                    # Infer parental allele from F2 genotype with weighted confidence
                    if f2_gt == '0/0':
                        evidence_for_parent.append(('0', 0.8 * confidence_weight))
                    
                    elif f2_gt == '0/1':
                        evidence_for_parent.append(('0', 0.7 * confidence_weight))
                    
                    elif f2_gt == '1/1':
                        evidence_for_parent.append(('1', 0.8 * confidence_weight))
                
                # Determine consensus for this common parent
                variant_record[f"{common_parent}_evidence"] = str(raw_evidence)
                
                if evidence_for_parent:
                    # Count weighted votes for each allele
                    allele_votes = {'0': 0.0, '1': 0.0}
                    for allele, confidence in evidence_for_parent:
                        allele_votes[allele] += confidence
                    
                    # Choose the allele with the most votes
                    common_parent_allele = max(allele_votes, key=allele_votes.get)
                    common_parent_confidence = max(allele_votes.values()) / sum(allele_votes.values())
                    variant_record[f"{common_parent}_allele"] = common_parent_allele
                    variant_record[f"{common_parent}_confidence"] = round(common_parent_confidence, 2)
                    
                    # Add reliability flag
                    if common_parent_confidence < 0.6:
                        variant_record[f"{common_parent}_reliability"] = "low"
                    elif common_parent_confidence < 0.8:
                        variant_record[f"{common_parent}_reliability"] = "medium"
                    else:
                        variant_record[f"{common_parent}_reliability"] = "high"
                else:
                    # No evidence for this parent at this position
                    variant_record[f"{common_parent}_allele"] = 'N'
                    variant_record[f"{common_parent}_confidence"] = 0.0
                    variant_record[f"{common_parent}_reliability"] = "none"
            
            # Process alternative parents - these are handled independently
            for alt_parent, samples in alt_parents.items():
                evidence_for_parent = []
                raw_evidence = []
                
                for sample_id in samples:
                    if sample_id not in roi_variants.columns:
                        continue
                    
                    # Get sample family information
                    common_p, alt_p = f2_samples[sample_id]
                    
                    # Skip if this sample doesn't involve this alternative parent
                    if alt_p != alt_parent:
                        continue
                    
                    # Extract genotype and quality metrics
                    f2_gt, f2_depth, f2_qual = extract_genotype(variant[sample_id], return_quality=True)
                    
                    # Skip if F2 genotype is missing
                    if f2_gt == './.':
                        continue
                    
                    # Get regional pattern
                    region_pattern = f2_regional_patterns.get(sample_id, {}).get('pattern', 'unknown')
                    
                    # Calculate confidence level based on depth and consistency
                    confidence_level = "high"
                    if f2_depth == 1:
                        confidence_level = "very_low"
                    elif f2_depth == 2:
                        confidence_level = "low"
                    elif f2_depth <= 4:
                        confidence_level = "medium"
                        
                    # Check allelic depths for heterozygous calls
                    allelic_imbalance = False
                    if f2_gt == '0/1' and ':' in variant[sample_id]:
                        fields = variant[sample_id].split(':')
                        if len(fields) >= 3 and ',' in fields[2]:
                            try:
                                allelic_depths = [int(d) for d in fields[2].split(',')]
                                if len(allelic_depths) >= 2:
                                    if allelic_depths[0] > 3*allelic_depths[1] or allelic_depths[1] > 3*allelic_depths[0]:
                                        allelic_imbalance = True
                                        if confidence_level != "very_low":
                                            confidence_level = "low"
                            except ValueError:
                                pass
                    
                    # Store raw evidence
                    raw_evidence.append({
                        'sample': sample_id,
                        'genotype': f2_gt, 
                        'depth': f2_depth,
                        'allelic_imbalance': allelic_imbalance,
                        'confidence': confidence_level,
                        'region_pattern': region_pattern
                    })
                    
                    # Get common parent's inferred allele for this position
                    common_p_allele = variant_record.get(f"{common_p}_allele", 'N')
                    
                    # Calculate confidence weight
                    confidence_weight = {
                        "very_low": 0.1,
                        "low": 0.3,
                        "medium": 0.7,
                        "high": 1.0
                    }[confidence_level]
                    
                    # Infer alternative parent's allele based on F2 genotype and common parent's allele
                    if f2_gt == '0/0':
                        if common_p_allele == '0':
                            # Both parents have reference allele
                            evidence_for_parent.append(('0', 0.8 * confidence_weight))
                        elif common_p_allele == '1':
                            # Inconsistent pattern - don't add evidence
                            pass
                            
                    elif f2_gt == '0/1':
                        if common_p_allele == '0':
                            # Common has reference, alt has alternate
                            evidence_for_parent.append(('1', 0.7 * confidence_weight))
                        elif common_p_allele == '1':
                            # Common has alternate, alt has reference
                            evidence_for_parent.append(('0', 0.7 * confidence_weight))
                        else:
                            # Common parent allele unknown, don't make strong inference
                            pass
                            
                    elif f2_gt == '1/1':
                        if common_p_allele == '1':
                            # Both parents have alternate allele
                            evidence_for_parent.append(('1', 0.8 * confidence_weight))
                        elif common_p_allele == '0':
                            # Inconsistent pattern - don't add evidence
                            pass
                
                # Determine allele for this alternative parent
                variant_record[f"{alt_parent}_evidence"] = str(raw_evidence)
                
                if evidence_for_parent:
                    # Count weighted votes for each allele
                    allele_votes = {'0': 0.0, '1': 0.0}
                    for allele, confidence in evidence_for_parent:
                        allele_votes[allele] += confidence
                    
                    # Choose the allele with the most votes
                    alt_parent_allele = max(allele_votes, key=allele_votes.get)
                    alt_parent_confidence = max(allele_votes.values()) / sum(allele_votes.values()) if sum(allele_votes.values()) > 0 else 0
                    variant_record[f"{alt_parent}_allele"] = alt_parent_allele
                    variant_record[f"{alt_parent}_confidence"] = round(alt_parent_confidence, 2)
                    
                    # Add reliability flag
                    if alt_parent_confidence < 0.6:
                        variant_record[f"{alt_parent}_reliability"] = "low"
                    elif alt_parent_confidence < 0.8:
                        variant_record[f"{alt_parent}_reliability"] = "medium"
                    else:
                        variant_record[f"{alt_parent}_reliability"] = "high"
                else:
                    # No evidence for this parent at this position
                    variant_record[f"{alt_parent}_allele"] = 'N'
                    variant_record[f"{alt_parent}_confidence"] = 0.0
                    variant_record[f"{alt_parent}_reliability"] = "none"
            
            # Store raw F2 data for reference
            variant_record.update(raw_f2_data)
            
            # Calculate overall metrics for this variant
            allele_counts = sum(1 for k, v in variant_record.items() 
                                if k.endswith('_allele') and v != 'N')
            variant_record['allele_count'] = allele_counts
            
            results.append(variant_record)
        
        # STEP 3: After processing all variants, add contextual validation
        results_with_context = []
        
        for i, result in enumerate(results):
            # Define context window (surrounding variants)
            window_size = 20
            start_idx = max(0, i - window_size//2)
            end_idx = min(len(results) - 1, i + window_size//2)
            context_variants = results[start_idx:end_idx+1]
            
            # Skip current variant
            context_variants = [v for v in context_variants if v['POS'] != result['POS']]
            
            # Check for agreement with surrounding context for each parent
            for parent in list(common_parents.keys()) + list(alt_parents.keys()):
                if f"{parent}_allele" in result and result[f"{parent}_allele"] != 'N':
                    parent_allele = result[f"{parent}_allele"]
                    
                    # Count matching context variants
                    matches = 0
                    total = 0
                    
                    for variant in context_variants:
                        if f"{parent}_allele" in variant and variant[f"{parent}_allele"] != 'N':
                            total += 1
                            if variant[f"{parent}_allele"] == parent_allele:
                                matches += 1
                    
                    if total >= 3:  # Only consider if we have enough context
                        result[f"{parent}_context_agreement"] = round(matches/total, 2)
                        
                        # Flag discordant variants
                        if matches/total < 0.3:
                            result[f"{parent}_context_conflict"] = True
            
            results_with_context.append(result)
        
        # Save complete results with all information
        if results_with_context:
            output_file = f"{output}_{roi_name}.tsv"
            pd.DataFrame(results_with_context).to_csv(output_file, sep='\t', index=False)
            logging.info(f"Saved complete results for ROI {roi_name} to {output_file}")
            
            # Create simplified results file (just key information)
            simplified_results = []
            for result in results_with_context:
                simplified_result = {
                    'CHROM': result['CHROM'],
                    'POS': result['POS'],
                    'REF': result['REF'],
                    'ALT': result['ALT'],
                }
                
                # Add all parental allele information
                for col in result.keys():
                    if col.endswith('_allele') or col.endswith('_reliability') or \
                       col.endswith('_context_agreement') or col.endswith('_confidence'):
                        simplified_result[col] = result[col]
                
                simplified_results.append(simplified_result)
            
            simplified_output = f"{output}_simplified_{roi_name}.tsv"
            pd.DataFrame(simplified_results).to_csv(simplified_output, sep='\t', index=False)
            logging.info(f"Saved simplified results for ROI {roi_name} to {simplified_output}")
        else:
            logging.info(f"No results for ROI {roi_name}")
