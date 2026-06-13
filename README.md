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

**'calculate_PRS_SE'** - calculates SE for each PRS score from individual effect size standard errors.

**'thresholding_by_pvalue_PRS'** - full pipeline to calculate PRS using p-value thresholding.

**'LD_clumping_PRS'** - full pipeline to calculate PRS using LD clumping.

**'thresholding_and_clumping_PRS'** - full pipeline to calculate PRS using p-value thresholding and LD clumping.

**'thresholding_and_shrinking_PRS'** - full pipeline to calculate PRS using p-value thresholding and shrinking.

**'choose_command_line_option'** - implements switching between PRS calculation methods and passing parameters to these 
methods.

## Dependencies and versions:

Python - version 3.8.10

The program uses packages 'tabix' and 'numpy'.

## Usage of the program:

The first input file of the program - a sorted GWAS summary statistics vcf file. The second input file of the program
is a sorted merged annotated target data vcf.gz file with genotypes of individuals. The target data should have also have
an indexing 'tbi' file in the same directory it's called from.

The program has four different modes:

**p_val_threshold** - full pipeline to calculate PRS using p-value thresholding. The given parameter is p-value.

`python3 ms_prs.py base_data.sorted.vcf target_data.vcf.gz p_val_threshold 0.05`

**LD_clump** - full pipeline to calculate PRS using LD clumping. The given parameter is r2 threshold.

`python3 ms_prs.py base_data.sorted.vcf target_data.vcf.gz LD_clump 0.1`

**p_val_LD_clump** - full pipeline to calculate PRS using p-value thresholding and LD clumping. The first given parameter i
p-value, the second is r2 threshold.

`python3 ms_prs.py base_data.sorted.vcf target_data.vcf.gz p_val_LD_clump 0.05 0.1`

**p_val_shrinkage** - full pipeline to calculate PRS using p-value thresholding and shrinking. The first given parameter is
p-value, the second is shrinkage.

`python3 ms_prs.py base_data.sorted.vcf target_data.vcf.gz p_val_shrinkage 0.05 2.0`

## File format:

### Base data (GWAS summary statistics) file format:

`CHROM POS ID REF ALT QUAL FILTER FORMAT [data displayed like format]`

**FORMAT** is **'ES:SE:LP:ID'**.

**ES** - *'effect size estimate relative to the alternate allele'.*

**SE** - *'standard error of effect size estimate'.*

**LP** - *'-log10 p-value for effect estimate'.*

**ID** - *'study variant identifier'.*

### Target data (individual genotype) file format for n individuals:

`CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [genotype]*n`

**FORMAT** is **'GT:PL'**.

**GT** - *'genotype'.*

**PL** - *'List of Phred scaled genotype likelyhoods'.*

## R-script:

The R-script analyzes calculated PRS scores and normalizes them using 'z-score' normalization.
Outputs the summary statistics of the scores and creates graphs comparing the normalized scores between different PRS
calculation methods.

## Testing of the program:

**synthetic_vcf_generator.py** is an AI generated program that creates synthetic target data files according to the
supplied GWAS summary statistics base data file.

'synthetic.vcf.gz' and 'synthetic.vcf.gz.tbi' were created from OpenGWAS.io database GWAS summary statistics file with ID
'ieu-b-18' (this file was aligned to the GRCh38 genome, indexed and sorted, first 250 lines and last 10 lines of this file
were concatenated, after that the last record with X chromosome was removed because of ambiguity while calculating PRS).

'synthetic.vcf.gz' and 'synthetic.vcf.gz.tbi' are just examples.

### Creating synthetic files from given GWAS summary statistics:

`-i` is the input file.

`-o` is the output file.

`-n` is the number of individuals for PRS calculation. 

`python3 synthetic_vcf_generator.py -i base_data.sorted.vcf -o synthetic.vcf -n 5`

`bgzip synthetic.vcf`

`tabix -p vcf synthetic.vcf.gz`

The program is then tested with the GWAS summary statistics it was created from and the synthetic target data.
