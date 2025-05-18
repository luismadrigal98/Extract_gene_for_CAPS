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
from src.screen_variants import *
from src.ancestry_inference import *

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

    # VCF initial filtering

    masking_parser = subparsers.add_parser('Mask', help='Using the vcf file, process genic regions (inferred from the gff3 file) and filter the variants to be used for the design of the markers based on the ROI list (which must be in concordance with the reference used for variant calling).')
    
    masking_parser.add_argument('--vcf', type=str, required=True, help='Input VCF file')
    masking_parser.add_argument('--gff3', type=str, required=True, help='Input GFF3 file')
    masking_parser.add_argument('--ROI_list', type=str, required=True, help='List of regions of interest (ROI) to be screened')
    masking_parser.add_argument('--output', type=str, required=True, help='Output file for the screened variants. This is a vcf file with the variants that are in the ROI list and are in genic regions.')
    masking_parser.add_argument('--only_biallelic', action='store_true', help='Only keep biallelic variants. This is going to be the default option, but you can change it if you want to keep all the variants.')
    masking_parser.add_argument('--min_qual', type=float, required=False, help='Minimum quality score for the variants to be considered. Default is 20.', default=20)

    # Determine ancestry from the F2 data
    # Given that those F2 are the closest ones to the actual plants employed in the experiment, we are going to try to infer the alleles from each of the parental lines.
    # The ROI are the ones that expand 664 inversions, so, there is not recombination in the F1, and in consequence, F2 can be either homozygous or heterozygous,
    # but no recombination should be evinced (signal of each allele origin should be clear then).

    inference_parser = subparsers.add_parser('Infer', help='Infer the ancestry of the F2 data')
    inference_parser.add_argument('--vcf', type=str, required=True, help='Input VCF file')
    inference_parser.add_argument('--ROI_list', type=str, required=True, help='List of regions of interest (ROI) to be screened')
    inference_parser.add_argument('--ancestry_log', type=str, required=True, help='This should be a tab delimited file will the sample ID for the F2s and the parental lines involved. This will allow to retrieve the identity of the alleles present in the original parental lines.') 
    inference_parser.add_argument('--output', type=str, required = True, help="Ouptput name. This is a file with the variant information and the inferred allele from each parental line. So, instead of having the actual samples, we will see the alleles for each parental. Also, the likelihood of the allele will be reported if specified.")
    inference_parser.add_argument('--context', type=int, default=20, help='How many variants to consider in the contextual analysis.')

    # Searching for diagnostic markers

    screen_parser = subparsers.add_parser('Screen', help='Screen the variants for diagnostic markers')

    screen_parser_global = screen_parser.add_argument_group("Global arguments")

    screen_parser_global.add_argument('--vcf', type=str, required=True, help='Input VCF file. This is the filtered vcf file with the variants that are in the ROI list and are in genic regions.')
    screen_parser_global.add_argument('--ROI_list', type=str, required=True, help='List of regions of interest (ROI) to be screened')
    screen_parser_global.add_argument('--output_dir', type=str, required=True, help='Output file for the screened variants. This is a directory where individual tab-delimited files per region of interest are going to be created. Each file contains the ID of the marker, and the coordinates for further extraction when primers are going to be design.')

    screen_parser_filtering = screen_parser.add_argument_group("Filtering by conditions to meet. High level filtering before starting the application of the rules to detect diagnostic markers.")

    screen_parser_filtering.add_argument('--min_qual', type=float, required=False, help='Minimum quality score for the variants to be considered. Default is 20.', default=20)
    screen_parser_filtering.add_argument('--min_dp', type=int, required=False, help='Minimum depth of coverage for the variants to be considered. Default is 5.', default=10) # Remember, you are working with shallow sequencing data, so the depth of coverage is going to be low. The default is 5, but you can change it to 10 if you want.

    screen_parser_rules = screen_parser.add_argument_group("Rules to detect diagnostic markers. These rules are going to be applied to the filtered variants.")

    screen_parser_rules.add_argument('--distance_to_closest_marker', type=int, required=False, help='Distance to the closest marker. Default is 1000 bp.', default=1000) # This is the distance to the closest marker. If the distance is too small, the markers are going to be too close to each other and they are going to be difficult to amplify.
    screen_parser_rules.add_argument('--non_informative_thr_F2s', type=int, required=False, help='Non informative threshold for the F2s. In other words, how many missing genotypes are we willing to accept. Default is 2.', default=2)
    screen_parser_rules.add_argument('--heterozygous_thr_support_F2s', type=int, required=False, help='Heterozygous threshold for the F2s. In other words, how many heterozygous genotypes are we willing to accept. Default is 2.', default=2)

    # Execute the right command
    args = parser.parse_args()
    if args.command == 'Remap':
        remap_variants(args.current_ref, args.new_ref, args.ROI_list, args.output, args.minimap2, args.minimap2_opts, args.samtools_exe,
                        args.temp_dir, args.keep)
    elif args.command == 'Mask':
        mask_variants(args.vcf, args.gff3, args.ROI_list, args.output, args.only_biallelic, args.min_qual)
    elif args.command == 'Infer':
        infer_ancestry(args.vcf, args.ROI_list, args.ancestry_log, args.output, args.context)
    elif args.command == 'Screen':
        screen_variants(args.vcf, args.ROI_list, args.output_dir, args.min_qual, args.min_dp,
                        args.distance_to_closest_marker, args.non_informative_thr_F2s, args.heterozygous_thr_support_F2s)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()