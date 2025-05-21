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
from src.primer_design import *

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
    inference_parser.add_argument('--ancestry_log', type=str, required=True, help='This should be a tab delimited file will the sample ID for the F2s and the parental lines involved.') 
    inference_parser.add_argument('--output', type=str, required=True, help="Output name for parental allele inference results")
    inference_parser.add_argument('--context', type=int, default=20, help='How many variants to consider in the contextual analysis.')
    inference_parser.add_argument('--approach', type=str, default='single', choices=['multiple', 'single'], help='Approach to use for the ancestry inference.')
    inference_parser.add_argument('--use_assembly_when_f2_missing', action='store_true', 
                            help='Use assembly data for positions without F2 data (if False, positions without F2 data will be skipped)')
    inference_parser.add_argument('--min_depth', type=int, default=3, help='Minimum read depth to consider a call reliable')

    # Searching for diagnostic markers  <<< Pending to be implemented

    screen_parser = subparsers.add_parser('Screen', help='Screen the variants for diagnostic markers')

    screen_parser_global = screen_parser.add_argument_group("Global arguments")

    screen_parser_global.add_argument('--inferred_alleles_tsv', nargs='+', type=str, required=True, help='This should be a tab delimited file with the inferred alleles for parental lines involved in the mapping population.')
    screen_parser_global.add_argument('--output_dir', type=str, required=True, help='Output file for the screened variants. This is a directory where individual tab-delimited files per region of interest are going to be created. Each file contains the ID of the marker, and the coordinates for further extraction when primers are going to be design.')

    screen_parser_filters = screen_parser.add_argument_group("Filtering arguments")

    screen_parser_filters.add_argument('--allele_col_pattern', type=str, required=False, default='_allele', help='Pattern to identify the allele columns in the inferred alleles file. Default is "_allele".')
    screen_parser_filters.add_argument('--overall_reliability_to_retain', type=str, required=False, default='high', choices=['high', 'medium', 'low'], help='Overall reliability to retain the markers. Default is "high".')
    screen_parser_filters.add_argument('--diff_parental', type=str, required=True, help='Common parental or parental for which we want to find a unique variant relative to alternative parentals.')
    screen_parser_filters.add_argument('--potential_size_of_amplicon', type=int, required=False, default=300, help="This is the desired size of the amplicon. This value is going to be used, assuming the variant will be in the middle of the amplicon, to search half of this size to the right and to the left, and see if a typical primer could be set up to amplify the region without being put in another variant region (which could disrupt amplification)")
    screen_parser_filters.add_argument('--potential_size_of_primers', type=int, required=False, default=20, help="This is the desired size of the primers. This value is going to be used to search for primers in the region of interest. Default is 20.")
    screen_parser_filters.add_argument('--displace_amplicon_window', action='store_true', help='Displace the amplicon window to the right and left of the variant. This is going to be used to search to try to accommodate the primers if they originally fall in a variant region. Default is False.')
    screen_parser_filters.add_argument('--displacement_tol', type=int, required=False, default=10, help='Displacement tolerance for the amplicon window. This is going to be used to search to try to accommodate the primers if they originally fall in a variant region. Default is 10. In other words, how many steps can we move the window to the right or to the left?')

    # Primer design step based on primer3
    design_parser = subparsers.add_parser('Design', help='Design primers for the variants')

    design_parser_inout = design_parser.add_argument_group("Input and output files")
    design_parser_inout.add_argument('--input_files', '-i', nargs='+', type=str, required=True, 
                                help='Input TSV files from screen_variants step')
    design_parser_inout.add_argument('--reference_fasta', '-r', type=str, required=True,
                                help='Reference FASTA file containing the sequences')
    design_parser_inout.add_argument('--output', '-o',type=str, required=True, 
                                help='Output file for the designed primers')
    design_parser_inout.add_argument('--error_log', '-e', type=str, required=False, 
                                help='Error log file for the designed primers')

    design_parser_filtering = design_parser.add_argument_group("Filtering options")
    design_parser_filtering.add_argument('--quality_threshold', type=str, default='high', choices=['high', 'medium', 'low'],
                                    help='Minimum quality threshold for variants (high, medium, low)')
    design_parser_filtering.add_argument('--min_high', type=int, default=3, 
                                    help='Minimum number of high reliability variants required')
    design_parser_filtering.add_argument('--min_medium', type=int, default=2, 
                                    help='Minimum number of medium reliability variants required')
    design_parser_filtering.add_argument('--max_low', type=int, default=0, 
                                    help='Maximum number of low reliability variants allowed')
    design_parser_filtering.add_argument('--max_variants', type=int, default=-9, 
                                    help='Maximum number of variants to process per file. Default -9 is meant to indicate no limit')

    design_parser_sequence = design_parser.add_argument_group("Sequence options")
    design_parser_sequence.add_argument('--flanking_size', type=int, default=150, 
                                    help='Size of flanking region on each side of variant')

    design_parser_primer3 = design_parser.add_argument_group("Primer3 arguments")
    design_parser_primer3.add_argument('--primer3_exe', type=str, required=False, 
                                help='Path to primer3 executable', 
                                default='~/.conda/envs/salmon/bin/primer3_core')
    design_parser_primer3.add_argument('--primer3_args', type=str, required=False, 
                                help='Command Line Options for primer3', 
                                default='--default_version=2 --format_output --strict_tags')
    design_parser_primer3.add_argument('--settings_file', type=str, required=False, 
                                help='Settings file for primer3')

    design_parser_misc = design_parser.add_argument_group("Miscellaneous options")
    design_parser_misc.add_argument('--keep_temp', action='store_true', 
                            help='Keep temporary Primer3 input/output files for inspection')
    design_parser_misc.add_argument('--temp_dir', type=str, required=False, 
                            help='Directory to store temporary files (created if not exists)')

    # >>>> COMMANDS ARE MANAGED HERE <<<<< #
    # Execute the right command
    args = parser.parse_args()
    if args.command == 'Remap':
        remap_variants(args.current_ref, args.new_ref, args.ROI_list, args.output, args.minimap2, args.minimap2_opts, args.samtools_exe,
                        args.temp_dir, args.keep)
    elif args.command == 'Mask':
        mask_variants(args.vcf, args.gff3, args.ROI_list, args.output, args.only_biallelic, args.min_qual)
    elif args.command == 'Infer':
        if args.approach == 'multiple':
            infer_ancestry_multiple(args.vcf, args.ROI_list, args.ancestry_log, args.output, args.context)
        elif args.approach == 'single':
            infer_ancestry_single(args.vcf, args.ROI_list, args.ancestry_log, args.output, 
                                    use_assembly_when_f2_missing=args.use_assembly_when_f2_missing,
                                    min_depth=args.min_depth)
    elif args.command == 'Screen':
        screen_variants(args.inferred_alleles_tsv, args.output_dir, args.allele_col_pattern, args.overall_reliability_to_retain,
                        args.diff_parental, args.potential_size_of_amplicon, args.potential_size_of_primers, args.displace_amplicon_window,
                        args.displacement_tol)
    elif args.command == 'Design':
        design_primers(args.input_files, args.reference_fasta, args.output, 
                        args.settings_file, args.primer3_exe, args.primer3_args,
                        args.quality_threshold, args.min_high, args.min_medium, 
                        args.max_low, args.flanking_size, 1,  # target_length parameter
                        args.max_variants, args.keep_temp, args.temp_dir, args.error_log)

    else:
        parser.print_help()

if __name__ == "__main__":
    main()