"""Primer Contrast Module"""

import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def select_best_primers(all_results, num_to_select=50, selection_criteria="balanced"):
    """
    Select the best primers from all designed primers.
    
    Parameters:
    all_results (list): List of all primer design results.
    num_to_select (int): Number of top primers to select.
    selection_criteria (str): Strategy for selection: 'balanced', 'tm_stability', 'size', 'specificity'
    
    Returns:
    list: Selected top results.
    """
    # Flatten the list of primers from all variants
    all_primers = []
    
    for result in all_results:
        for i, primer in enumerate(result['primer_results']['primers']):
            # Create a copy of the result for each primer
            primer_result = result.copy()
            primer_result['selected_primer'] = primer
            primer_result['primer_index'] = i
            
            # Calculate composite scores based on selection criteria
            if selection_criteria == "balanced":
                # Lower is better for this score
                score = (
                    primer['pair_penalty'] * 5.0 +  # Primer3's own penalty
                    abs(primer['left']['tm'] - primer['right']['tm']) * 3.0 +  # Tm difference
                    abs(primer['left']['gc_percent'] - 50) * 0.5 +  # GC optimal deviation
                    abs(primer['right']['gc_percent'] - 50) * 0.5 +  # GC optimal deviation
                    float(primer['left']['self_any']) * 1.0 +  # Self complementarity
                    float(primer['right']['self_any']) * 1.0 +  # Self complementarity
                    float(primer['left']['self_end']) * 2.0 +  # 3' self complementarity
                    float(primer['right']['self_end']) * 2.0    # 3' self complementarity
                )
                
                # Apply reliability factor (high=1.0, medium=1.2, low=1.5)
                reliability_factor = 1.0
                if result['reliability'] == 'medium':
                    reliability_factor = 1.2
                elif result['reliability'] == 'low':
                    reliability_factor = 1.5
                
                primer_result['composite_score'] = score * reliability_factor
            
            elif selection_criteria == "tm_stability":
                # Focus on Tm and stability
                score = (
                    abs(primer['left']['tm'] - primer['right']['tm']) * 5.0 +  # Tm difference
                    abs(primer['left']['tm'] - 60.0) * 2.0 +  # Deviation from optimal Tm
                    abs(primer['right']['tm'] - 60.0) * 2.0 +  # Deviation from optimal Tm
                    float(primer['left']['self_any']) * 1.0 +  # Self complementarity
                    float(primer['right']['self_any']) * 1.0    # Self complementarity
                )
                
                reliability_factor = 1.0
                if result['reliability'] == 'medium':
                    reliability_factor = 1.2
                elif result['reliability'] == 'low':
                    reliability_factor = 1.5
                
                primer_result['composite_score'] = score * reliability_factor
                
            elif selection_criteria == "size":
                # Focus on product size - prefer certain size ranges
                size = primer['product_size']
                
                # Calculate how close we are to ideal ranges
                if 80 <= size <= 120:
                    size_penalty = 0  # Ideal small range
                elif 150 <= size <= 200:
                    size_penalty = 0  # Ideal medium range
                elif 220 <= size <= 280:
                    size_penalty = 0  # Ideal large range
                else:
                    # Calculate distance to nearest ideal range
                    dist_to_small = min(abs(size - 80), abs(size - 120)) if size < 80 or size > 120 else 0
                    dist_to_medium = min(abs(size - 150), abs(size - 200)) if size < 150 or size > 200 else 0
                    dist_to_large = min(abs(size - 220), abs(size - 280)) if size < 220 or size > 280 else 0
                    size_penalty = min(dist_to_small, dist_to_medium, dist_to_large) * 2.0
                
                score = (
                    size_penalty +
                    primer['pair_penalty'] * 3.0 +  # Still consider Primer3's penalty
                    abs(primer['left']['tm'] - primer['right']['tm']) * 1.0  # Reduced weight for Tm difference
                )
                
                reliability_factor = 1.0
                if result['reliability'] == 'medium':
                    reliability_factor = 1.2
                elif result['reliability'] == 'low':
                    reliability_factor = 1.5
                
                primer_result['composite_score'] = score * reliability_factor
                
            elif selection_criteria == "specificity":
                # Focus on minimizing primer dimers and self-complementarity
                score = (
                    float(primer['left']['self_any']) * 3.0 +  # Self complementarity
                    float(primer['right']['self_any']) * 3.0 +  # Self complementarity
                    float(primer['left']['self_end']) * 4.0 +  # 3' self complementarity
                    float(primer['right']['self_end']) * 4.0 +  # 3' self complementarity
                    primer['pair_penalty'] * 2.0 +  # Still consider Primer3's penalty
                    abs(primer['left']['tm'] - primer['right']['tm'])  # Reduced weight for Tm difference
                )
                
                reliability_factor = 1.0
                if result['reliability'] == 'medium':
                    reliability_factor = 1.2
                elif result['reliability'] == 'low':
                    reliability_factor = 1.5
                
                primer_result['composite_score'] = score * reliability_factor
            
            all_primers.append(primer_result)
    
    # Sort by composite score (lower is better)
    all_primers.sort(key=lambda x: x['composite_score'])
    
    # Check that you have enough primers to select
    if len(all_primers) < num_to_select:
        logging.warning(f"Not enough primers to select. Found {len(all_primers)}, but need {num_to_select}.")
    num_to_select = min(num_to_select, len(all_primers))
    logging.info(f"Selecting top {num_to_select} primers based on '{selection_criteria}' criteria")

    # Take the top primers
    return all_primers[:num_to_select]