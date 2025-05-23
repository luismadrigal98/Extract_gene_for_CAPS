"""
Parallelized primer design using Primer3 for a list of variants.

"""

def process_variant_for_parallel(args):
    """Process a single variant for primer design (parallelizable version)"""
    # Import everything needed in this process
    import os
    import tempfile
    import subprocess
    import logging
    from Bio import SeqIO
    
    # Define extract_sequence locally since it's not available in worker process
    def extract_sequence_local(fasta_file, chrom, start, end):
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id == chrom or record.id.split()[0] == chrom:
                    return str(record.seq[start-1:end])
            return None
        except Exception:
            return None
    
    # Define create_primer3_input locally
    def create_primer3_input_local(name, sequence, target_pos, target_length, settings, output_file):
        target = f"{target_pos},{target_length}"
        with open(output_file, 'w') as f:
            f.write(f"SEQUENCE_ID={name}\n")
            f.write(f"SEQUENCE_TEMPLATE={sequence}\n")
            f.write(f"SEQUENCE_TARGET={target}\n")
            for key, value in settings.items():
                f.write(f"{key}={value}\n")
            f.write("=\n")
        return output_file
    
    # Define run_primer3 locally
    def run_primer3_local(input_file, primer3_exe, settings_file=None, primer3_args="", also_get_formatted=False, timeout=300):
        cmd = [primer3_exe]
        if primer3_args:
            args_list = primer3_args.split()
            args_list = [arg for arg in args_list if arg != '--format_output' and arg != '-format_output']
            cmd.extend(args_list)
        if settings_file:
            cmd.extend(["--p3_settings_file", settings_file])
        
        try:
            with open(input_file, 'r') as f:
                result = subprocess.run(cmd, stdin=f, text=True, capture_output=True, check=True, timeout=timeout)
            boulder_output = result.stdout
            
            if not also_get_formatted:
                return boulder_output
                
            format_cmd = cmd.copy()
            format_cmd.append("--format_output")
            with open(input_file, 'r') as f:
                format_result = subprocess.run(format_cmd, stdin=f, text=True, capture_output=True, check=True, timeout=timeout)
            return (boulder_output, format_result.stdout)
        except Exception:
            return None
    
    # Define parse_primer3_output locally
    def parse_primer3_output_local(output):
        result = {}
        for line in output.split("\n"):
            if line.strip() == "=":
                break
            if "=" in line:
                key, value = line.split("=", 1)
                result[key] = value
        
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
    
    # Now process the variant
    try:
        # Unpack arguments
        variant_data, reference_fasta, flanking_size, target_length, default_settings, \
        keep_temp, temp_dir, primer3_exe, settings_file, primer3_args, process_id, timeout = args
        
        # Unpack variant data
        idx, variant = variant_data
        
        # Extract the variant info
        chrom = variant['CHROM']
        pos = variant['POS']
        ref = variant['REF']
        alt = variant['ALT']
        
        # Calculate region to extract
        start = max(1, pos - flanking_size)
        end = pos + flanking_size
        target_pos = pos - start

        # Extract sequence
        sequence = extract_sequence_local(reference_fasta, chrom, start, end)
        if not sequence:
            return None
        
        # Create temp files
        if keep_temp:
            variant_dir = os.path.join(temp_dir, f"{chrom}_{pos}")
            if not os.path.exists(variant_dir):
                os.makedirs(variant_dir)
            input_file_path = os.path.join(variant_dir, "primer3_input.txt")
            output_file_path = os.path.join(variant_dir, "primer3_output.txt")
        else:
            temp_input = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt')
            input_file_path = temp_input.name
            temp_input.close()
        
        # Create primer3 input
        try:
            create_primer3_input_local(f"M_{chrom}_{pos}", sequence, target_pos, target_length, 
                                     default_settings, input_file_path)
        except Exception:
            if not keep_temp:
                try:
                    os.unlink(input_file_path)
                except:
                    pass
            return None
        
        # Run primer3
        primer3_exe_path = os.path.expanduser(primer3_exe)
        primer3_result = run_primer3_local(
            input_file=input_file_path, 
            primer3_exe=primer3_exe_path, 
            settings_file=settings_file, 
            primer3_args=primer3_args,
            also_get_formatted=keep_temp,
            timeout=timeout
        )

        # Handle the result
        if keep_temp and primer3_result and isinstance(primer3_result, tuple):
            boulder_output, formatted_output = primer3_result
            with open(output_file_path, 'w') as f:
                f.write(formatted_output)
            parsed_output = parse_primer3_output_local(boulder_output)
            if not boulder_output:
                return None
        else:
            primer3_output = primer3_result
            if keep_temp and primer3_output:
                with open(output_file_path, 'w') as f:
                    f.write(primer3_output)
            parsed_output = parse_primer3_output_local(primer3_output)
            if not primer3_output:
                return None

        if parsed_output['num_returned'] == 0:
            return None
        
        # Clean up temp files if not keeping them
        if not keep_temp:
            try:
                os.unlink(input_file_path)
            except:
                pass
        
        # Return the result
        return {
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
        
    except Exception as e:
        # Print to stderr since logging isn't available in worker process
        import sys
        print(f"Error processing variant: {str(e)}", file=sys.stderr)
        return None