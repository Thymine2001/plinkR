# plinkR

R wrapper for PLINK/PLINK2 with automatic output capture, genotype reading (BED/BIM/FAM or PED/MAP), and report loading (missingness, PCA, etc.). Includes a built-in C++ BED reader (adapted from plink2R) so you can read binary outputs directly.

## Installation
```r
# Install from GitHub
install.packages("remotes")
remotes::install_github("Thymine2001/plinkR")
```

## Quick start
```r
library(plinkR)

# Check PLINK availability
check_plink()

# Missingness (reads imiss/lmiss into $report)
mis <- plink("--bfile data/mydata --missing")
mis$report$imiss  # sample missingness
mis$report$lmiss  # SNP missingness

# PCA (reads eigenvec/eigenval into $report)
pca <- plink("--bfile data/mydata --pca")


# Make-bed and read binary outputs
a <- plink("--bfile data/mydata --make-bed --out data/mydata_filtered")
dim(a$out_geno$bed)   # samples x SNPs
head(colnames(a$out_geno$bim)) # CHR SNP CM POS A1 A2
```
