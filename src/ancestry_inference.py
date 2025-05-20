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

def infer_ancestry_multiple(vcf, ROI_list, ancestry_log, output, context_window=20, approach='single'):
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

        # Initialize results list
        results = []  # Fix: Initialize the results list before using it
        
        # For each variant
        for idx, variant in roi_variants.iterrows():
            variant_record = {'CHROM': variant['CHROM'], 'POS': variant['POS'],
                            'REF': variant['REF'], 'ALT': variant['ALT']}
            
            # Store cross-specific data
            cross_data = {}
            cross_predictions = {}
            
            # Store confidence values for later weighted averaging
            parent_confidences = {p: [] for p in all_parents}
            parent_alleles = {p: [] for p in all_parents}
            total_genotype_counts = {'0/0': 0, '0/1': 0, '1/1': 0, './.': 0}
            total_samples_checked = 0
            
            # Process each cross - consolidated genotype processing loop
            for cross, samples in f2_groups.items():
                common_parent, alt_parent = cross.split('_')
                
                # Extract genotypes for this variant
                genotype_counts = {'0/0': 0, '0/1': 0, '1/1': 0, './.': 0}
                quality_data = {'0/0': [], '0/1': [], '1/1': []}

                # Extract genotype data
                for sample in samples:
                    if sample in roi_variants.columns:
                        # Get genotype and quality metrics
                        gt, depth, qual = extract_genotype(variant[sample], return_quality=True)
                        genotype_counts[gt] += 1
                        if gt != './.':
                            quality_data[gt].append((1, depth, qual))
                
                # Update total counts
                for gt, count in genotype_counts.items():
                    total_genotype_counts[gt] += count
                total_samples_checked += sum(genotype_counts.values())
                
                # Calculate missing ratio
                total_samples = sum(genotype_counts.values())
                missing_samples = genotype_counts['./.']
                missing_ratio = missing_samples / total_samples if total_samples > 0 else 1.0

                # Skip likelihood calculation if too much data is missing
                if missing_ratio > 0.95:
                    print(f"Skipping variant at {variant['CHROM']}:{variant['POS']} for cross {cross}: {missing_ratio:.2f} missing data")
                    cross_predictions[cross] = {common_parent: "N", alt_parent: "N"}
                    cross_data[cross] = {"confidence": 0, "missing_ratio": missing_ratio}
                else:
                    # Calculate average quality metrics for each genotype
                    avg_quality = {}
                    for gt, data in quality_data.items():
                        if data:
                            total_depth = sum(d for _, d, _ in data)
                            total_qual = sum(q for _, _, q in data)
                            count = len(data)
                            avg_quality[gt] = (count, total_depth/count if count > 0 else 0, total_qual/count if count > 0 else 0)
                    
                    # Calculate inference
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
                
                # Store per-cross counts
                variant_record[f'{cross}_hom_ref'] = genotype_counts['0/0']
                variant_record[f'{cross}_het'] = genotype_counts['0/1']
                variant_record[f'{cross}_hom_alt'] = genotype_counts['1/1']
                variant_record[f'{cross}_missing'] = genotype_counts['./.']
            
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
            
            # Store total counts
            variant_record[f'hom_ref_count'] = total_genotype_counts['0/0']
            variant_record[f'het_count'] = total_genotype_counts['0/1']
            variant_record[f'hom_alt_count'] = total_genotype_counts['1/1']
            variant_record[f'missing_count'] = total_genotype_counts['./.']
            variant_record[f'missing_ratio'] = total_genotype_counts['./.'] / total_samples_checked if total_samples_checked > 0 else 1.0
            
            # Append to results
            results.append(variant_record)
        
        # STEP 4: Add context validation 
        results_with_context = []
        
        for i, result in enumerate(results):
            # Make a copy of the result
            result_with_context = result.copy()
            
            # Define context window - more efficient approach
            window_size = context_window  # Use the parameter
            start_idx = max(0, i - window_size//2)
            end_idx = min(len(results) - 1, i + window_size//2)
            
            # Check for agreement with surrounding context for each parent
            for parent in all_parents:
                parent_allele = result.get(f"{parent}_allele", 'N')
                
                if parent_allele != 'N':
                    # Count matching context variants more efficiently
                    matches = 0
                    total = 0
                    
                    for j in range(start_idx, end_idx+1):
                        if j == i:  # Skip current variant
                            continue
                            
                        v = results[j]
                        if v.get(f"{parent}_allele", 'N') != 'N':
                            total += 1
                            if v.get(f"{parent}_allele") == parent_allele:
                                matches += 1
                    
                    if total >= 3:  # Only calculate if enough context
                        result_with_context[f"{parent}_context_agreement"] = round(matches/total, 2)
                        
                        # Flag potential errors
                        if matches/total < 0.25 and total >= 5:
                            result_with_context[f"{parent}_potential_error"] = True
            
            results_with_context.append(result_with_context)
        
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
                'confidence': result.get('confidence', None)
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

def infer_ancestry_single(vcf, ROI_list, ancestry_log, output, use_assembly_when_f2_missing=False, min_depth=3):
    """
    Infer parental alleles for single F2 individuals with integrated assembly data handling.
    
    Args:
        vcf: Path to VCF file containing both F2 and parental/assembly data
        ROI_list: Path to regions of interest file
        ancestry_log: Path to relationship map file
        output: Base name for output files
        use_assembly_when_f2_missing: If True, use assembly data for positions without F2 data
        min_depth: Minimum read depth to consider a call reliable
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
    parental_samples = set()
    
    for _, row in ancestry_df.iterrows():
        if row['FC'] == 'F2':
            f2_samples[str(row['ID'])] = (row['Common'], row['Alt'])
            parental_samples.add(row['Common'])
            parental_samples.add(row['Alt'])
        elif row['FC'] in ['P', 'Assembly']:  # Add any other identifiers for parental samples
            parental_samples.add(row['ID'])
    
    logging.info(f"Found {len(f2_samples)} F2 samples and {len(parental_samples)} parental samples")
    
    # Get all parents that need allele assignments
    all_parents = set()
    for _, (common, alt) in f2_samples.items():
        all_parents.add(common)
        all_parents.add(alt)
    
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
                                (vcf_df['POS'] <= end)]
        
        if len(roi_variants) == 0:
            logging.info(f"No variants found in ROI {roi_name}")
            continue
        
        # STEP 1: Extract parental haplotypes directly from VCF if available
        parental_haplotypes = {}
        
        for parent in all_parents:
            if parent in roi_variants.columns:
                # Extract haplotype from parental data
                haplotype = []
                
                for _, variant in roi_variants.iterrows():
                    if parent in variant:
                        gt, depth, _ = extract_genotype(variant[parent], return_quality=True)
                        pos = variant['POS']
                        
                        # Only use reliable homozygous calls
                        if (gt == '0/0' or gt == '1/1') and depth >= min_depth:
                            allele = '0' if gt == '0/0' else '1'
                            haplotype.append((pos, allele))
                
                if haplotype:
                    parental_haplotypes[parent] = haplotype
                    logging.info(f"Extracted {len(haplotype)} reliable positions for {parent} directly from VCF")
        
        # STEP 2: Extract regional patterns from F2s
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
        
        # Identify common and alternative parents
        common_parents = {}
        alt_parents = {}
        
        # Group F2s by their parents
        for sample_id, (common, alt) in f2_samples.items():
            if common in common_parents:
                common_parents[common].append(sample_id)
            else:
                common_parents[common] = [sample_id]
                
            if alt in alt_parents:
                alt_parents[alt].append(sample_id)
            else:
                alt_parents[alt] = [sample_id]
        
        results = []
        
        # STEP 3: Process each variant position
        for _, variant in tqdm(roi_variants.iterrows(), desc=f"Processing variants in {roi_name}"):
            # Check if this position has F2 data
            has_f2_data = any(sample in variant for sample in f2_samples)
            
            # Create the variant record for ALL positions
            variant_record = {
                'CHROM': variant['CHROM'],
                'POS': variant['POS'],
                'REF': variant['REF'],
                'ALT': variant['ALT'],
                'QUAL': variant['QUAL'],
                'has_f2_data': has_f2_data
            }
            
            # Process each parent
            for parent in all_parents:
                parent_allele = 'N'  # Changed from 'NA' to 'N' for consistency
                parent_confidence = 0.0
                parent_reliability = 'none'
                parent_source = 'no_data_available'  # Better description than just 'none'
                
                # TIER 1: Direct parental genotype (if available)
                if parent in variant:
                    gt, depth, _ = extract_genotype(variant[parent], return_quality=True)
                    
                    if gt == '0/0' and depth >= min_depth:
                        parent_allele = '0'
                        parent_confidence = min(0.5 + (depth / 20), 0.95)
                        parent_reliability = 'high' if depth >= 10 else 'medium'
                        parent_source = 'direct_parental'
                    
                    elif gt == '1/1' and depth >= min_depth:
                        parent_allele = '1'
                        parent_confidence = min(0.5 + (depth / 20), 0.95)
                        parent_reliability = 'high' if depth >= 10 else 'medium'
                        parent_source = 'direct_parental'
                    elif gt != './.':  # Low depth but still has data
                        parent_allele = '0' if gt == '0/0' else '1'
                        parent_confidence = 0.3  # Low confidence due to insufficient depth
                        parent_reliability = 'low'
                        parent_source = 'low_depth_parental'
                
                # TIER 2: F2-contingent haplotype evidence
                if parent_allele == 'N' and has_f2_data:
                    # First gather F2 samples relevant to this parent
                    relevant_f2s = []
                    if parent in common_parents:
                        relevant_f2s.extend(common_parents[parent])
                    if parent in alt_parents:
                        relevant_f2s.extend(alt_parents[parent])
                    
                    if relevant_f2s:
                        # Find F2s with data at this position
                        f2s_with_data = [f2 for f2 in relevant_f2s if f2 in variant and extract_genotype(variant[f2]) != './.']
                        
                        if f2s_with_data:
                            # Look for nearby variants where F2s show consistent segregation pattern
                            current_pos = variant['POS']
                            window_variants = roi_variants[(roi_variants['POS'] > current_pos - 50000) & 
                                                            (roi_variants['POS'] < current_pos + 50000)]
                            
                            # For each F2 with data at current position, check consistency across nearby positions
                            f2_haplotype_evidence = []
                            
                            for f2 in f2s_with_data:
                                # Get this F2's genotype at current position
                                current_gt = extract_genotype(variant[f2])
                                
                                # Check consistency in nearby variants
                                consistent_positions = 0
                                total_positions = 0
                                
                                for _, nearby_var in window_variants.iterrows():
                                    if nearby_var['POS'] == current_pos:
                                        continue  # Skip current position
                                    
                                    if f2 in nearby_var and extract_genotype(nearby_var[f2]) == current_gt:
                                        consistent_positions += 1
                                    total_positions += 1
                                
                                if total_positions > 0 and consistent_positions / total_positions > 0.7:
                                    # This F2 shows a consistent haplotype pattern
                                    # Use it to infer parent's allele
                                    
                                    # For common parent:
                                    if parent in common_parents and f2 in common_parents[parent]:
                                        if current_gt == '0/0':
                                            f2_haplotype_evidence.append(('0', 0.6))
                                        elif current_gt == '1/1':
                                            f2_haplotype_evidence.append(('1', 0.6))
                                        # Heterozygous F2s are less useful for haplotype inference
                                    
                                    # For alternative parent:
                                    elif parent in alt_parents and f2 in alt_parents[parent]:
                                        # Get common parent for this F2
                                        common_p = next(cp for cp, ap in f2_samples.items() if ap == parent)
                                        common_p_allele = variant_record.get(f"{common_p}_allele")
                                        
                                        if common_p_allele != 'N':
                                            if current_gt == '0/0' and common_p_allele == '0':
                                                f2_haplotype_evidence.append(('0', 0.6))
                                            elif current_gt == '1/1' and common_p_allele == '1':
                                                f2_haplotype_evidence.append(('1', 0.6))
                                            elif current_gt == '0/0' and common_p_allele == '1':
                                                f2_haplotype_evidence.append(('0', 0.7))
                                            elif current_gt == '1/1' and common_p_allele == '0':
                                                f2_haplotype_evidence.append(('1', 0.7))
                            
                            # Combine all haplotype evidence
                            if f2_haplotype_evidence:
                                allele_votes = {'0': 0.0, '1': 0.0}
                                for allele, confidence in f2_haplotype_evidence:
                                    allele_votes[allele] += confidence
                                
                                if sum(allele_votes.values()) > 0:
                                    haplotype_allele = max(allele_votes, key=allele_votes.get)
                                    haplotype_confidence = max(allele_votes.values()) / sum(allele_votes.values())
                                    
                                    parent_allele = haplotype_allele
                                    parent_confidence = 0.6 * haplotype_confidence  # Still moderate confidence
                                    parent_reliability = 'medium'
                                    parent_source = 'f2_haplotype_inference'

                # TIER 3: Standard F2 evidence (keep as is but as Tier 3)
                if parent_allele == 'N' and has_f2_data:
                    # Process parent as common parent
                    if parent in common_parents:
                        evidence_for_parent = []
                        
                        for sample_id in common_parents[parent]:
                            if sample_id not in variant:
                                continue
                            
                            # Get F2 genotype
                            f2_gt, f2_depth, _ = extract_genotype(variant[sample_id], return_quality=True)
                            if f2_gt == './.':
                                continue
                            
                            # Get region pattern
                            region_pattern = f2_regional_patterns.get(sample_id, {}).get('pattern', 'unknown')
                            
                            # Calculate confidence level based on depth
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
                            
                            # Calculate confidence weight
                            confidence_weight = {
                                "very_low": 0.1,
                                "low": 0.3,
                                "medium": 0.7,
                                "high": 1.0
                            }[confidence_level]
                            
                            # Infer common parent's allele
                            if f2_gt == '0/0':
                                evidence_for_parent.append(('0', 0.8 * confidence_weight))
                            elif f2_gt == '0/1':
                                evidence_for_parent.append(('0', 0.7 * confidence_weight))
                            elif f2_gt == '1/1':
                                evidence_for_parent.append(('1', 0.8 * confidence_weight))
                        
                        # Determine consensus for common parent
                        if evidence_for_parent:
                            # Count weighted votes
                            allele_votes = {'0': 0.0, '1': 0.0}
                            for allele, confidence in evidence_for_parent:
                                allele_votes[allele] += confidence
                            
                            # Choose allele with most votes
                            f2_allele = max(allele_votes, key=allele_votes.get)
                            f2_confidence = max(allele_votes.values()) / sum(allele_votes.values())
                            
                            # Only update if we don't have a better source or if F2 confidence is high
                            if parent_allele == 'N' or (parent_source != 'direct_parental' and f2_confidence > parent_confidence):
                                parent_allele = f2_allele
                                parent_confidence = f2_confidence
                                parent_source = 'f2_inference'
                                
                                # Set reliability based on confidence
                                if parent_confidence < 0.6:
                                    parent_reliability = 'low'
                                elif parent_confidence < 0.8:
                                    parent_reliability = 'medium'
                                else:
                                    parent_reliability = 'high'
                    
                    # Process parent as alternative parent
                    if parent in alt_parents:
                        evidence_for_parent = []
                        
                        for sample_id in alt_parents[parent]:
                            if sample_id not in variant:
                                continue
                            
                            # Get sample family info
                            common_p, alt_p = f2_samples[sample_id]
                            
                            # Skip if this sample doesn't involve this alt parent
                            if alt_p != parent:
                                continue
                            
                            # Get F2 genotype
                            f2_gt, f2_depth, _ = extract_genotype(variant[sample_id], return_quality=True)
                            if f2_gt == './.':
                                continue
                            
                            # Calculate confidence level based on depth
                            confidence_level = "high"
                            if f2_depth == 1:
                                confidence_level = "very_low"
                            elif f2_depth == 2:
                                confidence_level = "low"
                            elif f2_depth <= 4:
                                confidence_level = "medium"
                            
                            # Calculate confidence weight
                            confidence_weight = {
                                "very_low": 0.1,
                                "low": 0.3,
                                "medium": 0.7,
                                "high": 1.0
                            }[confidence_level]
                            
                            # Get common parent's inferred allele
                            common_p_allele = variant_record.get(f"{common_p}_allele", None)
                            
                            # Infer alternative parent's allele based on F2 genotype and common parent's allele
                            if f2_gt == '0/0':
                                if common_p_allele == '0' or common_p_allele is None:
                                    evidence_for_parent.append(('0', 0.8 * confidence_weight))
                            elif f2_gt == '0/1':
                                if common_p_allele == '0':
                                    evidence_for_parent.append(('1', 0.7 * confidence_weight))
                                elif common_p_allele == '1':
                                    evidence_for_parent.append(('0', 0.7 * confidence_weight))
                            elif f2_gt == '1/1':
                                if common_p_allele == '1' or common_p_allele is None:
                                    evidence_for_parent.append(('1', 0.8 * confidence_weight))
                        
                        # Determine consensus for alternative parent
                        if evidence_for_parent:
                            # Count weighted votes
                            allele_votes = {'0': 0.0, '1': 0.0}
                            for allele, confidence in evidence_for_parent:
                                allele_votes[allele] += confidence
                            
                            # Choose allele with most votes
                            f2_allele = max(allele_votes, key=allele_votes.get)
                            f2_confidence = max(allele_votes.values()) / sum(allele_votes.values())
                            
                            # Only update if we don't have a better source or if F2 confidence is high
                            if parent_allele == 'N' or (parent_source != 'direct_parental' and f2_confidence > parent_confidence):
                                parent_allele = f2_allele
                                parent_confidence = f2_confidence
                                parent_source = 'f2_inference'
                                
                                # Set reliability based on confidence
                                if parent_confidence < 0.6:
                                    parent_reliability = 'low'
                                elif parent_confidence < 0.8:
                                    parent_reliability = 'medium'
                                else:
                                    parent_reliability = 'high'
                
                # Store results for this parent
                variant_record[f"{parent}_allele"] = parent_allele
                variant_record[f"{parent}_confidence"] = round(parent_confidence, 2)
                variant_record[f"{parent}_reliability"] = parent_reliability
                variant_record[f"{parent}_source"] = parent_source

                # After all tiers are processed, categorize the data availability status
                if parent_source == 'no_data_available':
                    variant_record[f"{parent}_data_status"] = "missing"
                elif parent_source in ['direct_parental', 'haplotype_block']:
                    variant_record[f"{parent}_data_status"] = "direct"
                elif parent_source == 'f2_inference':
                    variant_record[f"{parent}_data_status"] = "inferred"
                elif parent_source == 'low_depth_parental':
                    variant_record[f"{parent}_data_status"] = "low_quality"
            
            # Add summary columns:
            
            # 1. Check if we have complete information for all parents
            all_parents_have_alleles = all(variant_record.get(f"{parent}_allele", 'N') != 'N' 
                                        for parent in all_parents)
            variant_record['complete_info'] = all_parents_have_alleles
            
            # 2. Calculate overall reliability score
            reliability_values = {
                "none": 0,
                "low": 1,
                "medium": 2,
                "high": 3
            }
            
            # Get reliability scores for each parent
            reliability_scores = [reliability_values[variant_record.get(f"{parent}_reliability", "none")] 
                                for parent in all_parents]
            
            # Calculate average reliability and map back to categories
            if avg_reliability >= 2.5:
                overall_rel = "high"
            elif avg_reliability >= 1.5:
                overall_rel = "medium"
            elif avg_reliability > 0:
                overall_rel = "low"
            else:
                overall_rel = "none"
                
            variant_record['overall_reliability'] = overall_rel
            
            # Record the data source distribution
            source_counts = {}
            for parent in all_parents:
                source = variant_record.get(f"{parent}_source", "unknown")
                if source not in source_counts:
                    source_counts[source] = 1
                else:
                    source_counts[source] += 1
            
            sources_str = "; ".join(f"{source}: {count}" for source, count in source_counts.items() if source != 'none')
            variant_record['evidence_summary'] = sources_str if sources_str else "none"
            
            results.append(variant_record)
        
        # STEP 4: Add context validation 
        results_with_context = []
        
        for i, result in enumerate(results):
            # Make a copy of the result
            result_with_context = result.copy()
            
            # Define context window
            window_size = 20
            start_idx = max(0, i - window_size//2)
            end_idx = min(len(results) - 1, i + window_size//2)
            context_variants = results[start_idx:end_idx+1]
            
            # Skip current variant
            context_variants = [v for v in context_variants if v['POS'] != result['POS']]
            
            # Check for agreement with surrounding context for each parent
            for parent in all_parents:
                parent_allele = result.get(f"{parent}_allele", 'N')
                
                if parent_allele != 'N':
                    # Count matching context variants
                    matches = 0
                    total = 0
                    
                    for v in context_variants:
                        if v.get(f"{parent}_allele", 'N') != 'N':
                            total += 1
                            if v.get(f"{parent}_allele") == parent_allele:
                                matches += 1
                    
                    if total >= 3:  # Only calculate if enough context
                        result_with_context[f"{parent}_context_agreement"] = round(matches/total, 2)
                        
                        # Flag potential errors
                        if matches/total < 0.3:
                            result_with_context[f"{parent}_context_conflict"] = True
                            
                            # Reduce confidence for variants that conflict with context
                            if f"{parent}_confidence" in result:
                                result_with_context[f"{parent}_confidence"] *= 0.5
                                
                                # Update reliability if confidence is reduced
                                new_conf = result_with_context[f"{parent}_confidence"]
                                if new_conf < 0.6:
                                    result_with_context[f"{parent}_reliability"] = "low"
                                elif new_conf < 0.8:
                                    result_with_context[f"{parent}_reliability"] = "medium"
            
            results_with_context.append(result_with_context)
        
        # Save results
        if results_with_context:
            output_file = f"{output}_{roi_name}.tsv"
            pd.DataFrame(results_with_context).to_csv(output_file, sep='\t', index=False)
            logging.info(f"Saved complete results for ROI {roi_name} to {output_file}")
            
            # Create simplified results
            simplified_results = []
            for result in results_with_context:
                simplified_result = {
                    'CHROM': result['CHROM'],
                    'POS': result['POS'],
                    'REF': result['REF'],
                    'ALT': result['ALT'],
                    'QUAL': result['QUAL'],
                    'has_f2_data': result['has_f2_data'],
                    'complete_info': result['complete_info'],
                    'overall_reliability': result['overall_reliability']
                }
                
                # Add all parental allele info
                for parent in all_parents:
                    if f"{parent}_allele" in result:
                        simplified_result[f"{parent}_allele"] = result[f"{parent}_allele"]
                    if f"{parent}_reliability" in result:
                        simplified_result[f"{parent}_reliability"] = result[f"{parent}_reliability"]
                    if f"{parent}_context_agreement" in result:
                        simplified_result[f"{parent}_context_agreement"] = result[f"{parent}_context_agreement"]
                
                simplified_results.append(simplified_result)
            
            simplified_output = f"{output}_simplified_{roi_name}.tsv"
            pd.DataFrame(simplified_results).to_csv(simplified_output, sep='\t', index=False)
            logging.info(f"Saved simplified results to {simplified_output}")
        else:
            logging.info(f"No results for ROI {roi_name}")