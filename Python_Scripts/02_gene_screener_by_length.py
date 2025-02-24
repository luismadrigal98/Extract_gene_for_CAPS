"""
This script will screen a csv or tsv file containing genes annotations of a ROI, derived from the script 01_genes_by_ROI.py.

Input:
    - csv or tsv files containing the genes of interest.

Output:
    - A clean csv or tsv file, with the genes sorted by length.   

@Author: Luis J. Madrigal Roca
@Date: 2024-07-11

"""

import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description = "Candidate gene screener by length")

parser.add_argument("-i", "--input", help = "Input file", required = True)
parser.add_argument("-s", "--sep", help = "Separator. Only ',' or '\t' currently admitted.", default = ",")
parser.add_argument("-of", "--output_format", help = "Format of the output file.", default = 'tsv', choices = ['csv', 'tsv'])

args = parser.parse_args()

if not ("csv" in args.input or "tsv" in args.input):
    print("Invalid input file format. Please use a csv or tsv file.")
    raise SystemExit

df = pd.read_csv(args.input, sep = args.sep)

## Example of the structure of the input file:
#,seq_id,source,type,start,end,score,strand,phase,attributes
#284058,Chr_11,phytozomev13,gene,16838722,16843579,.,-,.,ID=MgIM767.11G103200.v2.1;Name=MgIM767.11G103200;ancestorIdentifier=MgIM767.11G103200.v1.1

# Ensure the DataFrame has the required columns
if 'start' not in df.columns or 'end' not in df.columns:
    print("Error: Input file is missing 'start' or 'end' columns.")
    raise SystemExit

df["length"] = df["end"] - df["start"]
df["Gene_name"] = df["attributes"].str.split(";", expand = True)[1].str.split("=", expand = True)[1]

# Select the columns that are going to be in the output file
df = df[["seq_id", "type", "Gene_name", "start", "end", "length"]]
df = df.sort_values(by = "length", ascending = False)

# Extract the base name and extension
base_name, _ = os.path.splitext(args.input)

# Append "_sorted" to the base name and add the appropriate extension
if args.output_format == "csv":
    output_file = f"{base_name}_sorted.csv"
else:  # Assuming the only other option is "tsv"
    output_file = f"{base_name}_sorted.tsv"

# Save the DataFrame to the output file
df.to_csv(output_file, sep=',' if args.output_format == "csv" else '\t', index=False)