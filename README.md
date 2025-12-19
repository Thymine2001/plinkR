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
pca <- plink("--bfile path/data/mydata --pca")


# Make-bed and read binary outputs
a <- plink("--bfile path/data/mydata --make-bed --out data/mydata_filtered")
dim(a$out_geno$bed)   # samples x SNPs
head(colnames(a$out_geno$bim)) # CHR SNP CM POS A1 A2

# And all other parameters same as plink
plink("--file path/data --maf --recode --out outfile")
path/data indicates the plink file path and format
or setwd(path/) then
plink("--file data --maf --recode --out outfile")

And so on

a <- plink("--file gen --pca --out out_geno")
saveRDS(a, file = "out_geno_plink_result.rds")

```

# Acknowledgements

The functionality for reading PLINK binary genotype files (BED/BIM/FAM) in plinkR is based on C++ code adapted from the plink2R project by Gabriel Abraham:

Gabriel Abraham plink2R: Read PLINK BED/BIM/FAM files into R.
GitHub repository: https://github.com/gabraham/plink2R
