# set the working directory
setwd("~/GitHub/Hardy-weinberg/ToGitHub/")

# upload required package
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)

# read the plink output file which was obtain by folowing command: /path/to/vcftools --vcf sample.vcf --hardy --out sample.hwe  
sample.hwe <- read.table("sample.hwe.hwe", header = T, sep = "\t")

# set the threshold for Chi-square test for goodnes of fit
alpha.hwe <- 0.05

# calculate theoretical genotype and allele frequencies according to Hardy-Weinberg equllibrium
aa <-  cbind(cbind(seq(0, 1, 0.01), seq(0, 1, 0.01)^2), "aa") # calculate "a" allele and "aa" genotype frequencies 
AA <-  cbind(cbind(seq(0, 1, 0.01), (1-seq(0, 1, 0.01))^2), "AA") # calculate "A" allele and "AA" genotype frequencies
Aa <-  cbind(cbind(seq(0, 1, 0.01), 2*(seq(0, 1, 0.01)) * (1-seq(0, 1, 0.01))), "Aa") # calculate "A" allele and "Aa" genotype frequencies

# merge theoretical genotype and allele frequencies
hwe.simulations <- data.frame(rbind(Aa, AA, aa), stringsAsFactors = F) # merge theoretical probabilities for all genotypes
hwe.simulations$X3 <- as.factor(hwe.simulations$X3) # genotype name
hwe.simulations$X1 <- as.numeric(hwe.simulations$X1) # allele frequency
hwe.simulations$X2 <- as.numeric(hwe.simulations$X2) # genotype frequency calculated according to Hardy-Weinberg equillibrium

# read observed genotype counts
hwe.plot <- as.data.frame(str_split_fixed(as.character(sample.hwe$OBS.HOM1.HET.HOM2.), "/", 3), stringsAsFactors = F)
hwe.plot$V1 <- as.numeric(hwe.plot$V1)
hwe.plot$V2 <- as.numeric(hwe.plot$V2)
hwe.plot$V3 <- as.numeric(hwe.plot$V3)
colnames(hwe.plot) <- c("p2", "pq2", "q2") 

# calculate observed genotype frequencies 
hwe.plot$p <- (2*hwe.plot$p2 + hwe.plot$pq2) / (2*(hwe.plot$p2 + hwe.plot$q2 + hwe.plot$pq2)) #obseved "A" allele frequency
hwe.plot$AA <- hwe.plot$p2/(hwe.plot$p2 + hwe.plot$q2 + hwe.plot$pq2) # observed "AA" genotype frequencies
hwe.plot$aa <- hwe.plot$q2/(hwe.plot$p2 + hwe.plot$q2 + hwe.plot$pq2) # observed "aa" genotype frequencies
hwe.plot$Aa <- hwe.plot$pq2/(hwe.plot$p2 + hwe.plot$q2 + hwe.plot$pq2) # observed "Aa" genotype frequencies 

# make a dataframe of observed genotype and allele frequencies similar to hwe.simulations dataframe
df <- data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors = F)
for (i in 1:nrow(hwe.plot)){
    df <- rbind(c(hwe.plot$p[i], hwe.plot$AA[i], 2), df)
    df <- rbind(c(hwe.plot$p[i], hwe.plot$aa[i], 1), df)
    df <- rbind(c(hwe.plot$p[i], hwe.plot$Aa[i], 3), df)
}
colnames(df) <- c("X1", "X2", "X3")
df$X3 <- ifelse(df$X3 == 1, "AA", ifelse(df$X3 == 2, "aa", "Aa"))

# make a dotplot of expected and observed genotype and allele frequences
HWEplot <- ggplot(NULL, aes(x = X1, y = X2, color = X3)) + geom_point(data = df, aes(), alpha = 0.5) +
    geom_line(data = hwe.simulations, size=1.2) +
    theme_bw(base_size = 16) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    labs(x="Allele frequencies",y="Genotype frequencies") + 
    theme(legend.position = "bottom") +
    scale_color_discrete(name = "Genotypes", labels = c("AA", "Aa", "aa"))

# make a histogram of p-values obtained using Chi-square test for goodnes of fit
# threshhold indicates the selected alpha level
# the plot on the right side of the threshhold corrsponds to SNPs within Hardy-Weinberg equillibrum  
CHISQplot <- ggplot(sample.hwe, aes(log10(P_HWE))) +
    geom_histogram(binwidth = 1, alpha=0.3, fill='white', colour='black') +
    theme_bw(base_size = 16) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    geom_vline(xintercept = log10(alpha.hwe), color = "red", size=0.8) +
    xlab("log10(P-value for deviation from HWE)") + ylab("Number of SNPs")

# at this step you can save the plots using RStudio interface or go further 

# this is optional step for joining two plots above together
CHISQplot <- ggarrange(HWEplot, CHISQplot,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# save joined plots as png file
png(filename="Hardy-Weinberg.png", width = 1200, height = 600, units = "px")
print(CHISQplot)
dev.off()

# optional qqplot for P-value for deviation from HWE distribution
ggplot(sample.hwe, aes(sample = P_HWE)) +
    stat_qq() + stat_qq_line() + theme_bw()