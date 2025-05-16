"""
This is the main program of the repository. It will borrow ideas from individual script files but it is thought as a way to diminish the errors obtained when
based solely in the assembly of the parent lines for the design and identification of the markers.

It will be a pipeline that will take as input the assembly of the parents, the annotation of the genome, and the VCF file with the F2 data. It will output a list of markers
that can be used for the genotyping of the F2 population. The markers will be designed to be diagnostic of the variant of interest (originally thought to diagnose inversions).
The markers will be restricted to genic regions.

The pipeline will be divided into several steps:
1. **Annotation Processing**: This step will take the GFF3/GTF annotation file and the reference genome as input and will output a BED file with the genic regions. This will be used to create a mask for restricting the analysis to genes.
2. **Variant Processing**: This step will take the VCF file from the F2 sequencing and the BED file with the genes as input and will output a filtered VCF with candidate heterozygous sites in genes. This will be used to identify segregating variants with appropriate frequencies.
                        Also, at this level, we are going to filter or ROI, that is, the inversion coordinates.
3. **Assembly Validation**: This step will take the filtered VCF and the focal assembly as input and will output a validated variant list with confirmed differences. This will be used to check if the variants from the F2 data exist between the parental assemblies.
4. **Primer Design**: This step will take the validated variant list and the surrounding sequence context as input and will output primer pairs with quality metrics and predicted products. This will be used to design optimal primers for validated polymorphic sites.
The pipeline will be modular, with each step implemented in a separate script. The scripts will be designed to be run independently, but they will also be able to be run as part of the pipeline. The pipeline will be designed to be flexible and extensible, so that new steps can be added easily in the future.

The pipeline will be designed to be user-friendly, with clear documentation and examples. The scripts will be designed to be run from the command line, with options for input and output files. The pipeline will be designed to be run on a local machine or on a cluster, with options for parallel processing.

@author: Luis Javier Madrigal-Roca & John K. Kelly

@date 2025/05/16

"""

import os
import argparse
import sys
import logging

# Local dependencies

from src.remapping_variants import *

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("marker_wizard.log"),
        logging.StreamHandler(sys.stdout)
    ]
)

def main():

    # Define the parser

    parser = argparse.ArgumentParser(description="Marker Wizard: A pipeline for designing markers for genotyping F2 populations.")
    subparsers = parser.add_subparsers(dest='command', help='Subcommands')

    # Remapping the variant to new reference genome

    remap_parser = subparsers.add_parser('Remap', help='Remap the variant to the new reference genome')

    remap_parser.add_argument('--current_ref', type=str, required=True, help='Current reference genome file')
    remap_parser.add_argument('--new_ref', type=str, required=True, help='New reference genome file')
    remap_parser.add_argument('--ROI_list', type=str, required=True, help='List of regions of interest (ROI) to be remapped')
    remap_parser.add_argument('--output', type=str, required=True, help='Output file for the remapped variants')
    remap_parser.add_argument('--minimap2', type=str, required=False, help='Path to minimap2 executable',
                            default='/home/l338m483/.conda/envs/PyR/bin/minimap2')
    remap_parser.add_argument('--minimap2_opts', type=str, required=False, help='Options for minimap2',
                            default='-x asm10 -t 10 -p 0.9 -N 50 --secondary=no')
    remap_parser.add_argument('--samtools_exe', type=str, required=False, help='Path to samtools executable',
                            default='/home/l338m483/.conda/envs/PyR/bin/samtools')
    remap_parser.add_argument('--temp_dir', type=str, required=False, help='Temporary directory for intermediate files',
                            default=None)
    remap_parser.add_argument('--keep', action='store_true', help='Keep intermediate files')
    
    # Execute the right command
    args = parser.parse_args()
    if args.command == 'Remap':
        remap_variants(args.current_ref, args.new_ref, args.ROI_list, args.output, args.minimap2, args.minimap2_opts, args.samtools_exe,
                        args.temp_dir, args.keep)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

"""

## 4. Implementation Design for Multi-Evidence Marker Pipeline

```markdown
## Implementation Design: Multi-Evidence Marker Pipeline

The following modular design implements the multi-evidence approach:

### Components

1. **Annotation Processor** (`annotation_processor.py`)
   - Input: GFF3/GTF annotation file, reference genome
   - Output: BED file with genic regions
   - Function: Creates a mask for restricting analysis to genes

2. **Variant Processor** (`variant_processor.py`)
   - Input: VCF file from F2 sequencing, BED file with genes
   - Output: Filtered VCF with candidate heterozygous sites in genes
   - Function: Identifies segregating variants with appropriate frequencies

3. **Assembly Validator** (`assembly_validator.py`)
   - Input: Filtered VCF, focal assembly, alternative assembly
   - Output: Validated variant list with confirmed differences
   - Function: Checks if variants from F2 data exist between parental assemblies

4. **Primer Designer** (`primer_designer.py`)
   - Input: Validated variant list, surrounding sequence context
   - Output: Primer pairs with quality metrics and predicted products
   - Function: Designs optimal primers for validated polymorphic sites

### Execution Flow

```bash
# 1. Process annotation to create genic mask
python annotation_processor.py -a annotation.gff -r reference.fa -o genes.bed

# 2. Filter variants using F2 data and genic mask
python variant_processor.py -v f2_data.vcf -b genes.bed -o candidate_markers.vcf

# 3. Validate variants across parental assemblies
python assembly_validator.py -v candidate_markers.vcf -f focal.fa -a alternative.fa -o validated_markers.tsv

# 4. Design primers for validated markers
python primer_designer.py -m validated_markers.tsv -r reference.fa -o primer_pairs.tsv
"""