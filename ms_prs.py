import sys
import tabix
import math
import numpy as np

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

# Recoding target data sample genotypes.

def recode_genotype(target_data):
    tmp = []
    for SNP in target_data:
        for i in range(9, len(SNP)):
            tmp = SNP[i].split(":")
            if tmp[0] == "0/0":
                SNP[i] = 1
            elif tmp[0] == "0/1":
                SNP[i] = 2
            elif tmp[0] == "1/1":
                SNP[i] = 3
            else:
                SNP[i] = 0

def calculate_pair_LD(SNP1, SNP2):
    r2 = np.corrcoef(SNP1, SNP2)[0, 1] ** 2
    return r2

# Checking for correlation between SNPs.
# Removing the less significant SNP if a pair of SNPs is correlated.

def LD_clump(base_data, target_data, r2_threshold, window):
    i = 0
    while i < len(base_data) - 1:
        # Checking if the SNPs aren't too far apart
        # and on the same chromosome.
        if int(base_data[i + 1][1]) - int(base_data[i][1]) < window and \
        base_data[i][0] == base_data[i+1][0]:
            r2 = calculate_pair_LD(target_data[i][9:], target_data[i + 1][9:])
            if r2 > r2_threshold:
                tmp1 = base_data[i][9].split(":")
                tmp2 = base_data[i + 1][9].split(":")
                LP1 = float(tmp1[2])
                LP2 = float(tmp2[2])
                if LP1 < LP2:
                    del base_data[i]
                    del target_data[i]
                else:
                    del base_data[i + 1]
                    del target_data[i + 1]
            else:
                i += 1
        else:
           i += 1

# Keeping SNPs that are below a certain p-value threshold.

def filter_by_pval_threshold(base_SNPs, threshold):
    tmp = []
    filtered_base_SNPs = []
    # Transforming the p-value threshold to match LP transformation.
    transformed_threshold = math.log10(threshold) * (-1)
    for SNP in base_SNPs:
        # Extracting the LP (-log10 p-value for effect estimate).
        tmp = SNP[9].split(":")
        LP = float(tmp[2])
        if LP > transformed_threshold:
            filtered_base_SNPs.append(SNP)
    return filtered_base_SNPs

# Shrinking all effect sizes according to their LD using
# the Bayesian regression formula.

def shrink_effect_sizes(base_data, target_data, window, shrinkage):
    i = 0
    k = 0
    starting_point = 0
    LD_matrix = []
    row = []
    ES_before = []
    while i < starting_point + window - 1:
        j = i + 1
        # Storing original effect sizes of SNPs into an array.
        tmp = base_data[i][9].split(":")
        ES = float(tmp[0])
        ES_before.append(ES)
        # Calculating LD for i SNP with other SNPs in the window.
        while j < starting_point + window:
            r = math.sqrt(calculate_pair_LD(target_data[i][9:], target_data[j][9:]))
            row.append(r)
            j += 1
        # Creating LD matrix.
        LD_matrix.append(row)
        row = []
        i += 1
    # Applying the Bayesian formula: (LD + shrinkage * I)^(-1)
    LD_matrix = np.array(LD_matrix)
    identity_matrix = np.identity(window)
    identity_matrix = shrinkage * identity_matrix
    LD_matrix = np.add(LD_matrix, identity_matrix)
    LD_matrix = np.linalg.inv(LD_matrix)
    # Recalculating the effect sizes using the LD matrix.
    # Adding them to base data.
    while k < starting_point + window - 1:
        ES_modified = np.dot(LD_matrix[k], ES_before)
        base_data[k].append(ES_modified)

def calculate_PRS_score(base_data, target_data):
    scores = []
    score = 0
    for i in range(9,len(target_data)):
        for j in range(0, len(base_data)):
            tmp = base_data[j][9].split(":")
            ES = float(tmp[0])
            if target_data[j][i] != 0:
                score += ES * (target_data[j][i] - 1)
        score.append(score)
        score = 0

def main():
    base = parse_base_data()
    target = filter_target_data(base)
    recode_genotype(target)
    print("--------------------------------")
    print(target)
    print("--------------------------------")
    LD_clump(base, target, 0.1, 250000)
    print(base)
    print("--------------------------------")
    filtered = filter_by_pval_threshold(base, 0.05)
    print(filtered)

if __name__=="__main__":
    main()
