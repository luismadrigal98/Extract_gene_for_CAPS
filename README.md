# Extract_gene_for_CAPS

A collection of Python scripts for processing and analyzing genomic data, particularly focused on gene extraction and sequence analysis from GFF3 and FASTA files.

## Overview

This pipeline consists of five scripts that work together to:
1. Extract genes from regions of interest in GFF3 files
2. Screen and sort genes by length
3. Extract gene sequences from reference genomes
4. Build FASTA files from sequence data
5. Sort and group genes across multiple FASTA files
6. Reverse complement sequences

## Requirements

- Python 3.x
- pandas
- gffpandas
- minimap2 (for sequence alignment)

## Installation

```bash
# Clone the repository
git clone https://github.com/luismadrigal98/Extract_gene_for_CAPS.git

# Install required Python packages
pip install pandas gffpandas

# Ensure minimap2 is installed on your system
# For Ubuntu/Debian:
sudo apt-get install minimap2
```

## Scripts Description (For step-by-step execution)

### 1. Genes by ROI (`01_genes_by_ROI.py`)
Extracts genes from specified regions of interest in a GFF3 file.

```bash
python 01_genes_by_ROI.py <GFF3_file> [-roi_file ROI_FILE] [-roi ROI...] [-of {csv,tsv}]
```

- `GFF3_file`: Input GFF3 file with gene annotations
- `-roi_file`: File containing regions of interest (format: ROI_name\tChromosome\tstart\tend)
- `-roi`: Direct command line input of regions (format: ROI_name_Chromosome_start_end)
- `-of`: Output format (csv/tsv)

### 2. Gene Screener (`02_gene_screener_by_length.py`)
Screens and sorts genes by length from the output of script 1.

```bash
python 02_gene_screener_by_length.py -i INPUT_FILE [-s SEPARATOR] [-of {csv,tsv}]
```

- `-i`: Input file (CSV/TSV from script 1)
- `-s`: Separator (',' or '\t')
- `-of`: Output format

### 3. Gene Sequence Extractor (`03_gene_sequence_extractor.py`)
Extracts gene sequences from a reference FASTA file.

```bash
python 03_gene_sequence_extractor.py -i INPUT -gd GENE_DICT -ch CHROMOSOME -f FASTA [-s SEP] [-on OUTPUT_NAME]
```

- `-i`: Gene of interest or input file
- `-gd`: Gene dictionary file (CSV/TSV)
- `-ch`: Chromosome location
- `-f`: Reference FASTA file
- `-s`: Separator
- `-on`: Output filename

### 4. FASTA Builder (`04_fasta_builder_per_line.py`)
Constructs FASTA files from sequence data, processing multiple sequences.

```bash
python 04_fasta_builder_per_line.py -i INPUT [INPUT ...] -targ TARGET [TARGET ...] [-op OUTPUT_PREFIX] [-od OUTPUT_DIRECTORY]
```

- `-i`: Input FASTA file(s)
- `-targ`: Target reference directories
- `-op`: Output prefix
- `-od`: Output directory

### 5. FASTA Sorter (`05_fasta_sorter_by_gene.py`)
Groups and sorts FASTA files by gene name.

```bash
python 05_fasta_sorter_by_gene.py (-id INPUT_DIRECTORY [INPUT_DIRECTORY ...] | -i INPUT [INPUT ...]) -g GENE_NAME [GENE_NAME ...]
```

- `-id`: Input directory containing FASTA files
- `-i`: Input FASTA files
- `-g`: Gene names to extract

## ROI Global Extractor

The Global Extractor (`ROI_global_extractor.py`) automates the extraction and mapping of entire regions of interest (ROIs) across different genome assemblies.

```bash
python ROI_global_extractor.py -r ROI [-r ROI...] -f REFERENCE_FASTA -t TARGET [TARGET...] -o OUTPUT_PREFIX [--temp_dir TEMP_DIR] [--remove_temp]
```

## Pipeline Workflow

1. Start with a GFF3 file and identify regions of interest
2. Extract and sort genes from these regions
3. Extract gene sequences from reference genome
4. Build and organize FASTA files
5. Sort and group genes across files

## Output Files

- Script 1: `{ROI_name}_{chromosome}_{start}_{end}_genes.{csv/tsv}`
- Script 2: `{input_base_name}_sorted.{csv/tsv}`
- Script 3: `{output_name}.fasta`
- Script 4: `{output_prefix}_{input_file}_vs_{target_file}.fasta`
- Script 5: `{gene_name}_grouped.fasta`

## Steps for the global extractor (automated extraction of the entire ROI of interest)

## Steps for the search of potential diagnostic markers using transferred annotation file for the focal line (with the variant or variants of interest), vcf file with F2-shallow sequenced plants and parental assemblies. This last step will leverage different levels of evidence, given that the previous approaches are sensitive to errors. This approach is also prone to error, but might be more robust in the long run. The idea is to use the annotation file to mask whatever that is not a gene, search for sites where we see heterozygous individuals (that could be due to differences between the common parental and each alternative line inside each family), see if in truth that difference also exist between the focal assembly and the alternative parental, and report those positions and variants for guiding the primer design.

## 2. Enhanced Multi-Evidence Approach Description

```markdown
## Multi-Evidence Marker Development Pipeline

This advanced approach addresses the limitations of single-evidence marker design by integrating multiple layers of genomic evidence to identify high-confidence polymorphic sites suitable for diagnostic marker development.

### Rationale
Traditional approaches often target assembly differences that may reflect technical artifacts rather than biological variation. This pipeline combines four independent lines of evidence to identify reliable polymorphic sites:

1. **Genic Focus**: Utilizes gene annotation to restrict analysis to coding regions, which generally have higher assembly quality and more reliable mapping
2. **Segregation Evidence**: Leverages F2 sequencing data to identify sites showing expected Mendelian segregation patterns
3. **Parental Confirmation**: Validates polymorphisms by confirming their presence between parental assemblies
4. **Functional Relevance**: Prioritizes variants with potential functional impact (missense, splice site) when possible

### Workflow
1. Use annotation files to mask non-genic regions in the focal genome
2. Process VCF files from F2 shallow-sequenced plants to identify heterozygous sites
3. For each heterozygous site, confirm the polymorphism exists between the focal and alternative parental assemblies
4. Score and rank sites based on multiple quality metrics
5. Generate primer design recommendations for the highest-ranking sites

This integrated approach significantly improves the success rate of diagnostic marker development by filtering out artifacts and focusing on biologically validated polymorphisms.
```

## Author

Luis J. Madrigal Roca

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.