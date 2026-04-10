import sys
import tabix

# PRS calculation for multiple sclerosis.

# Parsing the base data (GWAS summary statistics) vcf file.

def parse_base_data():
    # First argument of the program is a base data vcf file.
    base = open(sys.argv[1],"r")
    lines = base.readlines()
    tmp_list = []
    base_SNPs = []
    for line in lines:
        # Ignoring metadata.
        if not line.startswith("#"):
            tmp_list = line.split()
            base_SNPs.append(tmp_list)
    base.close()
    return base_SNPs

# Filtering the genome positions in target data
# to have only the genome positions that are in the GWAS file.

def filter_target_data(base_SNPs):
    # The second argument of the program is a compressed vcf
    # target data file.
    target_pos = []
    target = tabix.open(sys.argv[2])
    for SNP in base_SNPs:
        chrom = SNP[0]
        pos = int(SNP[1])
        iterator = target.query(chrom, pos-1, pos)
        target_pos.append(next(iterator))
    return target_pos

def main():
    base = parse_base_data()
    target = filter_target_data(base)
    print(target)

if __name__=="__main__":
    main()
