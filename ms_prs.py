import sys

# PRS calculation for multiple sclerosis.

# Parsing the base data (GWAS summary statistics) vcf file.

def parse_base_data():
    base = open(sys.argv[1],"r")
    lines = base.readlines()
    print(lines[0])
    base.close()

parse_base_data()
