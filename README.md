# snp_to_HardyWeinberg.R

Rim Gubaev, 2019

The script snp_to_HardyWeinberg.R allows one to visualise observed and theoretical genotype frequencies calculated according to Hardy-Weinberg equilibrium using the information on observed genotype frequencies. The input file for the script represents a tab-delimited table produced by vcftools using the following command:
```
vcftools --vcf sample.vcf --hardy --out sample.hwe
```
The input data for vcftools represent a vcf file obtained using GATK/samtools or other SNP-calling tools or vcf file obtained using DNA-array genotyping. Thr output of vcftools is a table that containes information on observed genotype frequencies for each genetic variant present in vcf file as well as p-values obtained from chi-square goodenes of fit test. The lower the p-value the lower the probability that the observed genetic variant is in the Hardy-Weinberg equilibrium.

The script is then 
