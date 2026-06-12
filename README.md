# MS_PRS

Implementation of different Polygenic Risk Score (PRS) calculation methods. The program takes effect sizes of SNPs from base data 
'.vcf' files and calculates PRS from target data '.vcf' file for each individual. The method of PRS calculation differs in
the way the SNPs are chosen for the calculation.

## Functions: 

**'parse_base_data'** - every SNP record is parsed into fields CHROM, POS, ID and etc.

**'filter_target_data'** - collecting the same position records from target data that are in base data.

**'match_SNPs'** - quality control: the base data and target data fields CHROM, POS, ID, REF and ALT should match
for record pairs from base and target data. If these fields do not match, the record pair is discarded.

**'recode_genotype'** - the target data genotypes are recoded into the number of alternative alleles.

**'count_chromosome_SNPs'** - counts the number of SNPs in each chromosome.

**'filter_by_pval_threshold'** - only SNPs that are below a certain p-value threshold are left for PRS calculation.

**'calculate_pair_LD'** - calculates Pearson's correlation between two SNPs from target data.

**'LD_clump'** - if correlation between two SNPs exceeds a specified threhold (r2), the less significant SNP is removed from 
base data. Two SNPs are compared if they are on the same chromosome and the distance between them doesn't exceed 250kb.

**'shrink_effect_sizes'** - calculates posterior mean effects using an LD matrix and 'shrinkage' parameter. Default LD matrix
window is 50 SNPs.

**'iterate_over_chromosomes'** - calculates the start and end points in the base data to calculate LD matrices in 
'shrink_effect_sizes.'
**'calculate_PRS_score'** - sums the SNP effect sizes that are multiplied by the recoded genotype (PRS is calculated for each
individual in target data).
**'calculate_shrunk_PRS_score'** - uses the posterior mean effects for calculation.
**'calculate_PRS_SE'** - calculates SE for each PRS score from individual effect size standard effors.
**'thresholding_by_pvalue_PRS'** - full pipeline to calculate PRS using p-value thresholding.
**'LD_clumping_PRS'** - full pipeline to calculate PRS using LD clumping.
**'thresholding_and_shrinking_PRS'** - full pipeline to calculate PRS using p-value thresholding and shrinking.
**'choose_command_line_option'** - implements switching between PRS calculation methods and passing parameters to these 
methods.
