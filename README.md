# snp_to_HardyWeinberg.R

Rim Gubaev, 2019

The script snp_to_HardyWeinberg.R allows one to visualise observed and theoretical genotype frequencies calculated according to Hardy-Weinberg equilibrium (HWE) using the information on observed genotype frequencies. The input file for the script represents a tab-delimited table produced by vcftools using the following command:
```
vcftools --vcf sample.vcf --hardy --out sample.hwe
```
The input data for vcftools represent a vcf file obtained using GATK/samtools or other SNP-calling tools or vcf file obtained using DNA-array genotyping. The output of vcftools is a table that contains information on observed genotype frequencies for each genetic variant present in vcf file as well as p-values obtained from Chi-Square goodness of fit test. The lower the p-value the lower the probability that the observed genetic variant is in the Hardy-Weinberg equilibrium.

The using the table obtained by vcftools script produces the following picture:

![](https://raw.githubusercontent.com/RimGubaev/snp_to_HardyWeinberg/master/Hardy-Weinberg.png)

Panel **A** is a dot-plot representing observed and theoretical genotype frequencies. Lines correspond to theoretical genotype frequencies, dots correspond to observed genotype frequencies.
Panel **B** corresponds to a histogram of p-values obtained form Chi-Square goodness of fit test. Red threshold indicated the p-value equal to 0.05.

Email: rimgubaev@gmail.com
