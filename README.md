# MarkerWizard: Multi-Evidence Diagnostic Marker Development Tool

MarkerWizard is a comprehensive toolkit for developing robust diagnostic markers for genotyping, particularly focused on identifying reliable polymorphic sites using multiple lines of evidence.

## Overview

This pipeline integrates various data sources to identify high-confidence markers:
1. **Remaps variants** between different reference assemblies
2. **Filters variants** to focus on genic regions
3. **Infers ancestry** using F2 segregation data
4. **Screens for diagnostic markers** based on defined criteria
5. **Designs primers** for the identified markers

## Requirements

- Python 3.x
- Dependencies:
  - pandas
  - tqdm
  - logging
  - intervaltree
  - Bio (BioPython)
- External tools:
  - minimap2 (for sequence alignment)
  - samtools
  - primer3 (for primer design)

## Installation

```bash
# Clone the repository
git clone https://github.com/luismadrigal98/Extract_gene_for_CAPS.git

# Install required Python packages
pip install pandas tqdm intervaltree biopython

# Ensure external tools are installed on your system
# For Ubuntu/Debian:
sudo apt-get install minimap2 samtools
```

## Usage

MarkerWizard is organized as a command-line tool with several subcommands:

### 1. Remap Variants

Remaps regions of interest (ROI) between two different genome assemblies.

```bash
python MarkerWizard.py Remap --current_ref <CURRENT_REF_FASTA> --new_ref <NEW_REF_FASTA> \
  --ROI_list <ROI_LIST> --output <OUTPUT_FILE> [options]
```

Options:
- `--minimap2`: Path to minimap2 executable
- `--minimap2_opts`: Options for minimap2
- `--samtools_exe`: Path to samtools executable
- `--temp_dir`: Temporary directory
- `--keep`: Keep intermediate files

### 2. Mask Variants

Filters VCF files to focus on variants in genic regions within specified ROIs.

```bash
python MarkerWizard.py Mask --vcf <VCF_FILE> --gff3 <GFF3_FILE> \
  --ROI_list <ROI_LIST> --output <OUTPUT_FILE> [options]
```

Options:
- `--only_biallelic`: Only keep biallelic variants
- `--min_qual`: Minimum quality score (default: 20)

### 3. Infer Ancestry

Infers parental alleles from F2 segregation data.

```bash
python MarkerWizard.py Infer --vcf <VCF_FILE> --ROI_list <ROI_LIST> \
  --ancestry_log <ANCESTRY_LOG> --output <OUTPUT_PREFIX> [options]
```

Options:
- `--context`: Number of variants to consider in contextual analysis
- `--approach`: Approach to use ('multiple' or 'single')
- `--use_assembly_when_f2_missing`: Use assembly data for positions without F2 data
- `--min_depth`: Minimum read depth to consider a call reliable

### 4. Screen Variants

Screens variants for diagnostic markers based on defined criteria.

```bash
python MarkerWizard.py Screen --inferred_alleles_tsv <ALLELES_TSV> [<ALLELES_TSV>...] \
  --output_dir <OUTPUT_DIRECTORY> --diff_parental <DIFF_PARENT> [options]
```

Options:
- `--allele_col_pattern`: Pattern to identify allele columns
- `--overall_reliability_to_retain`: Reliability level ('high', 'medium', 'low')
- `--potential_size_of_amplicon`: Desired amplicon size
- `--potential_size_of_primers`: Desired primer size
- `--displace_amplicon_window`: Allow displacement of amplicon window
- `--displacement_tol`: Displacement tolerance steps

### 5. Design Primers

Designs primers for the identified markers.

```bash
python MarkerWizard.py Design --output <OUTPUT_FILE> --settings_file <SETTINGS_FILE> [options]
```

Options:
- `--primer3_exe`: Path to primer3 executable
- `--primer3_clo`: Command line options for primer3

## Pipeline Workflow

MarkerWizard employs a multi-evidence approach that integrates:

1. **Remapping**: Translating regions of interest between genome assemblies
2. **Genic Focus**: Restricting analysis to gene regions using GFF3 annotations
3. **Segregation Evidence**: Using F2 data to identify sites with expected Mendelian patterns
4. **Parental Validation**: Confirming polymorphisms between parental assemblies
5. **Primer Design**: Designing optimal primers for validated polymorphic sites

This integrated approach significantly improves the success rate of diagnostic marker development by filtering out artifacts and focusing on biologically validated polymorphisms.

## Output Files

- **Remapping**: Coordinate translation file between assemblies
- **Masking**: Filtered VCF with variants in genic regions
- **Inference**: TSV files with inferred parental alleles
- **Screening**: TSV files with candidate diagnostic markers
- **Design**: Primer design output files

## Author

Luis J. Madrigal Roca & John K. Kelly

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.