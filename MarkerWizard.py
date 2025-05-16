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
from src.masking_vcf import *

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

    # Diagnostic markers searching

    masking_parser = subparsers.add_parser('Mask', help='Using the vcf file, process genic regions (inferred from the gff3 file) and filter the variants to be used for the design of the markers based on the ROI list (which must be in concordance with the reference used for variant calling).')
    
    masking_parser.add_argument('--vcf', type=str, required=True, help='Input VCF file')
    masking_parser.add_argument('--gff3', type=str, required=True, help='Input GFF3 file')
    masking_parser.add_argument('--ROI_list', type=str, required=True, help='List of regions of interest (ROI) to be screened')
    masking_parser.add_argument('--output', type=str, required=True, help='Output file for the screened variants. This is a vcf file with the variants that are in the ROI list and are in genic regions.')

    # Execute the right command
    args = parser.parse_args()
    if args.command == 'Remap':
        remap_variants(args.current_ref, args.new_ref, args.ROI_list, args.output, args.minimap2, args.minimap2_opts, args.samtools_exe,
                        args.temp_dir, args.keep)
    elif args.command == 'Mask':
        mask_variants(args.vcf, args.gff3, args.ROI_list, args.output)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()