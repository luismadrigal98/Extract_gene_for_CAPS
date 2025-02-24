"""
This script will pull the genes contained in a given region of interest (ROI) from a GFF3 file.

Input:
    - GFF3 file: The file containing the gene annotations.
    - Region of interest (string or list): The region of the genome for which the genes will be extracted.

Output:
    - A list of genes contained in the ROI in csv format or tsv.    

@Author: Luis J. Madrigal Roca
@Date: 2024-07-11

"""

import sys
import gffpandas.gffpandas as gffpd  # type: ignore
import argparse as ap

# Create the parser
parser = ap.ArgumentParser(description="Extract genes from a region of interest in a GFF3 file.")

# Add the arguments
parser.add_argument("GFF3_file", help="The GFF3 file containing the gene annotations.")
parser.add_argument("-roi_file", "--Region_of_interest_file", help="A file containing regions of interest. The format should be: '<ROI_name>\\t<Chromosome>\\t<start>\\t<end>'.", default=None)
parser.add_argument("-roi", "--Region_of_interest", nargs="*", help="The region of the genome for which the genes will be extracted directly from command line. The format should be: '<ROI_name>_<Chromosome>_<start>_<end>'.", default=None)
parser.add_argument("-of", "--output_format", help="Format of the output.", default='csv', choices=['csv', 'tsv'])

args = parser.parse_args()

# Read the GFF3 file
gff = gffpd.read_gff3(args.GFF3_file)

ROI_name = []
Chr = []
p0 = []
pf = []
ROI = []

if args.Region_of_interest_file:
    with open(args.Region_of_interest_file, "r") as file:
        for line in file.readlines()[1:]:
            chain = line.strip().split("\t")
            ROI_name.append(chain[0])
            Chr.append(chain[1])
            p0.append(chain[2])
            pf.append(chain[3])
            ROI.append(f"{chain[1]}:{chain[2]}-{chain[3]}")
elif args.Region_of_interest:
    if len(args.Region_of_interest) == 1:
        chain = args.Region_of_interest[0].split("_")
        if len(chain) == 4:
            ROI_name, Chr, p0, pf = [chain[0]], [chain[1]], [chain[2]], [chain[3]]
            ROI = [f"{Chr[0]}:{p0[0]}-{pf[0]}"]
        else:
            print("Invalid ROI format. Please use '<ROI_name>_<Chromosome>_<start>_<end>'.")
            sys.exit(1)
    else:
        for roi_str in args.Region_of_interest:
            chain = roi_str.split("_")
            if len(chain) == 4:
                ROI_name.append(chain[0])
                Chr.append(chain[1])
                p0.append(chain[2])
                pf.append(chain[3])
                ROI.append(f"{chain[1]}:{chain[2]}-{chain[3]}")
            else:
                print("Invalid ROI format in list. Please use '<ROI_name>_<Chromosome>_<start>_<end>' for each ROI.")
                sys.exit(1)
else:
    print("Please provide a region of interest either via file or directly as an argument.")
    sys.exit(1)

for region in range(len(ROI)):
    genes = gff.filter_feature_of_type(['gene'])
    mask = genes.df.apply(lambda row: row['seq_id'] == Chr[region] and int(row['start']) >= int(p0[region]) and int(row['end']) <= int(pf[region]), axis=1)
    df = genes.df[mask]

    if df.empty:
        print(f"No genes found in {ROI[region]}.")
    else:
        df.to_csv(f"{ROI_name[region]}_{ROI[region].replace(':', '_')}_genes.{args.output_format}", sep='\t' if args.output_format == 'tsv' else ',')
        print(f"Genes extracted for {ROI[region]} and saved to file.")