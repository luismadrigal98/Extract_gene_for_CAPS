# MarkerWizard: Multi-Evidence Diagnostic Marker Development Tool

MarkerWizard is a comprehensive toolkit for developing robust diagnostic markers for genotyping, particularly focused on identifying reliable polymorphic sites using multiple lines of evidence. Now includes **ultra-fast optimized tools** for high-throughput SNP discovery.

## Overview

This pipeline integrates various data sources to identify high-confidence markers:
1. **Remaps variants** between different reference assemblies
2. **Filters variants** to focus on genic regions
3. **Infers ancestry** using F2 segregation data
4. **Screens for diagnostic markers** based on defined criteria
5. **Designs primers** for the identified markers

### 🚀 NEW: High-Performance Tools

- **ultra_fast_snp_finder.py**: Ultra-optimized SNP discovery (15-30 minutes vs hours)
- **fast_snp_finder.py**: Fast SNP discovery without primer design (1-2 hours)
- **fast_screen_variants.py**: Vectorized variant screening with parallel processing
- **diagnostic_snp_finder.py**: Complete pipeline with primer design

## Quick Start - Ultra-Fast SNP Discovery

For rapid diagnostic SNP identification, use the optimized tools:

### Option 1: Maximum Speed (15-30 minutes)
```bash
# Ultra-fast SNP discovery with vectorized operations and parallel processing
python ultra_fast_snp_finder.py \
  --vcf your_variants.vcf \
  --gff3 your_genes.gff3 \
  --roi_list your_regions.txt \
  --ancestry_map your_relationship_map.txt \
  --target_parent 664c \
  --min_qual 60 \
  --min_reliability medium \
  --min_spacing 2000 \
  --max_snps 20 \
  --n_workers 6 \
  --output_prefix UltraFast_Results

# For absolute maximum speed, add:
# --skip_primer_screen
```

### Option 2: Fast with Primer Validation (1-2 hours)
```bash
# Fast SNP discovery with primer compliance checking
python fast_snp_finder.py \
  --vcf your_variants.vcf \
  --gff3 your_genes.gff3 \
  --roi_list your_regions.txt \
  --ancestry_map your_relationship_map.txt \
  --target_parent 664c \
  --min_qual 60 \
  --min_reliability medium \
  --min_spacing 2000 \
  --max_snps 20 \
  --output_prefix Fast_Results
```

### Option 3: Complete Pipeline with Primer Design (2-4 hours)
```bash
# Full pipeline including primer design and validation
python diagnostic_snp_finder.py \
  --vcf your_variants.vcf \
  --reference your_reference.fasta \
  --gff3 your_genes.gff3 \
  --roi_list your_regions.txt \
  --ancestry_map your_relationship_map.txt \
  --target_parent 664c \
  --min_qual 60 \
  --max_variants 200 \
  --parallel \
  --num_workers 4 \
  --output_prefix Complete_Results
```

## Performance Comparison

| Tool | Speed | Features | Best For |
|------|--------|----------|----------|
| `ultra_fast_snp_finder.py` | ⚡ 15-30 min | Vectorized ops, parallel processing | Initial discovery |
| `fast_snp_finder.py` | 🚀 1-2 hours | Primer compliance, spacing filters | Quality candidates |
| `diagnostic_snp_finder.py` | 🔬 2-4 hours | Full primer design, validation | Production markers |

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
- `--context`: Number of neighboring variants to consider for F2 genotype consistency analysis (detects recombination events and sequencing errors)
- `--approach`: Approach to use ('multiple' or 'single')
- `--use_assembly_when_f2_missing`: Use assembly data for positions without F2 data
- `--min_depth`: Minimum read depth to consider a call reliable (default: 3)
- `--max_depth`: Maximum read depth to consider a call reliable - filters out high-coverage artifacts (default: 200)

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

### Fast SNP Discovery Output
- **`*_ultra_fast_diagnostic_snps.tsv`**: Ultra-fast SNP candidates with quality scores
- **`*_diagnostic_snps.tsv`**: Fast SNP discovery results with primer compliance
- **`*_filtered_genic.vcf`**: Quality-filtered VCF with genic variants only
- **`*_ancestry_inferred.tsv`**: Parental ancestry inference results

### Complete Pipeline Output
- **Remapping**: Coordinate translation file between assemblies
- **Masking**: Filtered VCF with variants in genic regions
- **Inference**: TSV files with inferred parental alleles
- **Screening**: TSV files with candidate diagnostic markers
- **Design**: Primer design output files

### Key Output Columns
- `CHROM`, `POS`: Genomic coordinates
- `REF`, `ALT`: Reference and alternative alleles
- `{parent}_allele`: Inferred alleles for each parent
- `overall_reliability`: Quality assessment (high/medium/low)
- `complete_info`: Boolean indicating complete evidence
- `has_f2_data`: Boolean indicating F2 data availability
- `quality_score`: Comprehensive quality metric
- `primer_compliant`: Primer design feasibility (if checked)

## Advanced Features

### Optimization Options
- **Parallel Processing**: Multi-core support for CPU-intensive steps
- **Vectorized Operations**: NumPy/pandas optimizations for large datasets
- **Memory Efficiency**: Optimized data structures and chunked processing
- **Configurable Filters**: Adjustable quality thresholds and spacing requirements

### Primer Design Integration
- **Primer3 Integration**: Automated primer design with customizable parameters
- **Conflict Detection**: Checks for adjacent variants that interfere with primers
- **Displacement Strategy**: Automatic amplicon repositioning to avoid conflicts
- **Quality Scoring**: Multi-factor primer quality assessment

### Input File Compatibility
The tools support flexible input formats:
- **ROI Files**: Both `ROI,Chr,Start,End` and `ROI_name,Chrom,Start,End` formats
- **VCF Files**: Standard VCF format with quality scores
- **GFF3 Files**: Gene annotation files for genic region filtering
- **Ancestry Maps**: Tab-separated files mapping samples to parental lines

## Troubleshooting & Tips

### Performance Optimization
1. **Choose the right tool for your needs**:
   - Use `ultra_fast_snp_finder.py` for initial exploration
   - Use `fast_snp_finder.py` for quality candidates
   - Use `diagnostic_snp_finder.py` for production markers

2. **Optimize parameters for speed**:
   - Increase `--min_spacing` to reduce candidates
   - Lower `--max_snps` for faster processing
   - Use `--skip_primer_screen` for maximum speed
   - Adjust `--n_workers` based on available CPU cores

3. **Memory management**:
   - For large VCF files, consider splitting by chromosome
   - Reduce `--max_variants` if running out of memory
   - Use SSD storage for temporary files when possible

### Common Issues
- **Column format errors**: Ensure ROI files have consistent headers
- **Missing dependencies**: Check that all required Python packages are installed
- **Primer3 not found**: Verify primer3 is installed and in PATH
- **Memory errors**: Reduce batch sizes or use more restrictive quality filters

### Quality Control
- **Check reliability distribution**: Aim for majority high/medium reliability variants
- **Validate spacing**: Ensure selected SNPs meet minimum distance requirements
- **Review diagnostic criteria**: Confirm target parent differs from all others
- **Primer compliance**: Verify no adjacent variants interfere with primer binding

### Example Cluster Job Scripts
See included job scripts:
- `ultra_fast_job.sh`: Ultra-fast SNP discovery (30 minutes, 8 cores)
- `fast_snp_job.sh`: Fast SNP discovery (2 hours, 4 cores)
- `remote_design.sh`: Complete pipeline (4+ hours, multiple cores)

## Author

Luis J. Madrigal Roca & John K. Kelly

## Version History

- **v2.0** (2025-06-12): Added ultra-fast optimized tools with vectorized operations and parallel processing
- **v1.5** (2025-06-11): Enhanced primer design with conflict detection and displacement strategies
- **v1.0** (2025): Initial release with basic pipeline functionality

## Performance Benchmarks

Tested on chromosome 6 inversion region with ~50,000 variants:

| Tool | Runtime | Memory | Output Quality |
|------|---------|--------|----------------|
| Ultra-fast (skip primer) | 15 min | 2GB | High-confidence candidates |
| Ultra-fast (with primer) | 30 min | 4GB | Primer-validated candidates |
| Fast SNP finder | 90 min | 6GB | Complete quality assessment |
| Full pipeline | 3+ hours | 8GB | Production-ready markers |

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

For bug reports or feature requests, please open an issue on the project repository.