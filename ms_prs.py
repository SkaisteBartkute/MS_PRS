import sys

# PRS calculation for multiple sclerosis.

# Parsing the base data (GWAS summary statistics) vcf file.

def parse_base_data():
    # First argument of the file is a base data vcf file.
    base = open(sys.argv[1],"r")
    lines = base.readlines()
    tmp_list = []
    for line in lines:
        # Ignoring metadata.
        if not line.startswith("#"):
            tmp_list.extend(line.split())
    base.close()
    return tmp_list

results = parse_base_data()
print(results)
