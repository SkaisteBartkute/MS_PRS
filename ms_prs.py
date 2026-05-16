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
    i = 0
    target_pos = []
    target = tabix.open(sys.argv[2])
    while i < len(base_SNPs):
        chrom = base_SNPs[i][0]
        pos = int(base_SNPs[i][1])
        iterator = target.query(chrom, pos-1, pos)
        if not iterator:
            print("Target chromosome position not found.")
            del base_SNPs[i]
        else:
            target_pos.append(next(iterator))
            i += 1
    return target_pos

# Checking whether the fields of SNPs from base and target data match.
def match_SNPs(base, target):
    i = 0
    while i < len(base):
        if (base[i][0] == target[i][0] and
            base[i][1] == target[i][1] and
            base[i][2] == target[i][2] and
            base[i][3] == target[i][3]):
            i += 1
        else:
            #print(f"Removing {base[i][1]} from base data and  \
            #{target[i][1]} target data. SNPs didn't match.")
            #print(base[i], target[i])
            del base[i]
            del target[i]

# Counting the number of SNPs for each chromosome.
def count_chromosome_SNPs(base):
    i = 0
    count = 0
    SNP_count = []
    while i < len(base) - 1:
        if (i + 1 == len(base) - 1):
            tmp = [base[i][0], count + 2]
            SNP_count.append(tmp)
        elif (base[i][0] == base[i + 1][0]):
            count += 1
        else:
            tmp = [base[i][0], count + 1]
            SNP_count.append(tmp)
            count = 0
        i += 1
    return SNP_count

# Recoding target data sample genotypes.

def recode_genotype(target_data):
    tmp = []
    for SNP in target_data:
        for i in range(9, len(SNP)):
            tmp = SNP[i].split(":")
            if tmp[0] == "0/0":
                SNP[i] = 0
            elif tmp[0] == "0/1" or tmp[0] == "1/0":
                SNP[i] = 1
            elif tmp[0] == "1/1":
                SNP[i] = 2
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
                print(base_data[i], base_data[i + 1])
                LP1 = float(tmp1[2])
                LP2 = float(tmp2[2])
                if LP1 < LP2:
                    print(base_data[i])
                    del base_data[i]
                    del target_data[i]
                else:
                    print(base_data[i + 1])
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

def shrink_effect_sizes(base_data, target_data, window, shrinkage, starting_point):
    LD_matrix = []
    row = []
    ES_before = []
    for i in range(starting_point, starting_point + window):
        # Storing original effect sizes of SNPs into an array.
        tmp = base_data[i][9].split(":")
        ES = float(tmp[0])
        ES_before.append(ES)
        # Calculating LD for i SNP with other SNPs in the window.
        for j in range(starting_point, starting_point + window):
            print(target_data[i][1], target_data[j][1])
            r = math.sqrt(calculate_pair_LD(target_data[i][9:], target_data[j][9:]))
            print(r)
            row.append(r)
        # Creating LD matrix.
        LD_matrix.append(row)
        row = []
    # Applying the Bayesian formula: (LD + shrinkage * I)^(-1)
    LD_matrix = np.array(LD_matrix)
    print(LD_matrix)
    identity_matrix = np.identity(window)
    print(identity_matrix)
    identity_matrix = shrinkage * identity_matrix
    print(identity_matrix)
    LD_matrix = np.add(LD_matrix, identity_matrix)
    print(LD_matrix)
    LD_matrix = np.linalg.inv(LD_matrix)
    print(LD_matrix)
    # Recalculating the effect sizes using the LD matrix.
    # Adding them to base data.
    print(ES_before)
    idx = 0
    for k in range(starting_point, starting_point + window):
        ES_modified = np.dot(LD_matrix[idx], ES_before)
        print(ES_modified)
        base_data[k].append(ES_modified)
        idx += 1

# To account for the fact that SNP count in each chromosome isn't always
# dividible by the window size, the last window size of the chromosome
# needs to be recalculated (last_window).

# Performing effect size shrinkage for each chromosome separately.

def iterate_over_chromosomes(base_data, target_data, SNP_count, window, shrinkage):
    starting_point = 0
    for chrom in SNP_count:
        last_window = chrom[1] % window
        cycles = chrom[1] // window
        if cycles != 0:
           for i in range(cycles):
               shrink_effect_sizes(base_data, target_data, window, shrinkage, \
                                   starting_point)
               starting_point += window
        if last_window != 0:
            shrink_effect_sizes(base_data, target_data, last_window, shrinkage, \
                                starting_point)
            starting_point += last_window

def calculate_PRS_score(base_data, target_data):
    scores = []
    score = 0
    for i in range(9,len(target_data[0])):
        for j in range(0, len(base_data)):
            tmp = base_data[j][9].split(":")
            ES = float(tmp[0])
            score += ES * target_data[j][i]
        scores.append(score)
        score = 0
    return scores

def thresholding_by_pvalue_PRS(p_value):
    base = parse_base_data()
    filtered = filter_by_pval_threshold(base, p_value)
    target = filter_target_data(filtered)
    SNP_count = count_chromosome_SNPs(filtered)
    match_SNPs(filtered, target)
    recode_genotype(target)
    scores = calculate_PRS_score(filtered, target)
    print(*scores, sep = ', ')
    SNP_count = count_chromosome_SNPs(filtered)
    print(*SNP_count)
    # Printing used SNPs for PRS calculation.
    for SNP in filtered:
        print(SNP[0], SNP[1], SNP[2])
    print(filtered[len(filtered)-1],target[len(filtered)-1])

def choose_command_line_option(option):
    if option == "p_val_threshold":
        p_val = sys.argv[4]
        if p_val:
            thresholding_by_pvalue_PRS(float(p_val))
        else:
            print("Missing p-value threshold.")
    else:
        print("Missing options.")

def main():
    choose_command_line_option(sys.argv[3])

if __name__=="__main__":
    main()
