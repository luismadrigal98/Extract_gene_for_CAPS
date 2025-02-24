""""
This script will extract the inversions, the lines where they are present, and the coordinates of the inversions

"""

import pandas as pd
import os

GENO_DIR = "/home/l338m483/scratch/Corrected_inv_genotype_files"

files = []
results = {}

for filename in os.listdir(GENO_DIR):
    if filename.endswith(".genotype.txt"):
        files.append(filename)

for file in files:
    with open(os.path.join(GENO_DIR, file), 'r') as f:
        
        key = file.split(".")[0]
        
        # Read only the first line of the file (the info we need is only there)
        line = f.readline().strip()

        # Split the line by the tab character
        line = line.split("\t")

        # Extract the inversion information
        INV_ID = line[0]
        CHROM = line[1].split(":")[0]
        INV_START = line[1].split(":")[1].split(",")[0]
        INV_END = line[1].split(":")[1].split(",")[1]
        GENES_IN_INV = line[2]

        # Store the information in a dictionary
        results[key] = {
            "INV_ID": int(INV_ID),
            "CHROM": CHROM,
            "INV_START": int(INV_START),
            "INV_END": int(INV_END),
            "GENES_IN_INV": int(GENES_IN_INV)
        }

# Create a DataFrame from the dictionary
df = pd.DataFrame.from_dict(results, orient='index')
df.to_csv("inversions_info_map.csv")