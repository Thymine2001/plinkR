# plinkR

R wrapper for PLINK/PLINK2 with automatic output capture, genotype reading (BED/BIM/FAM or PED/MAP), and report loading (missingness, PCA, etc.). Includes bidirectional PLINK ↔ BLUPF90 format conversion with allele restoration. Features a built-in C++ BED reader (adapted from plink2R) so you can read binary outputs directly.

## Features

- **PLINK Wrapper**: Execute any PLINK command from R with automatic output parsing
- **Auto PLINK Detection**: Finds system PLINK or uses bundled version from linkbreedeR
- **Report Parsing**: Automatically reads missingness, PCA, frequency reports into R
- **Binary File Reading**: Fast C++ reader for BED/BIM/FAM files
- **BLUPF90 Conversion**: Bidirectional conversion between PLINK and BLUPF90 formats

## Installation

```r
# Install from GitHub
install.packages("remotes")
remotes::install_github("Thymine2001/plinkR")

# Optional: Install linkbreedeR for bundled PLINK binary
remotes::install_github("Thymine2001/linkbreedeR")
```

## Quick Start

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

# Other PLINK commands
plink("--file path/data --maf --recode --out outfile")
```

## BLUPF90 Format Conversion

plinkR provides bidirectional conversion between PLINK and BLUPF90 genotype formats, with automatic allele restoration using BIM files.

### PLINK → BLUPF90

```r
# Convert PLINK files to BLUPF90 format
# Input: pigs.ped + pigs.map (or pigs.bed + pigs.bim + pigs.fam)
# Output: geno.txt + geno.map + geno.bim

plink_to_blupf90("pigs", out_prefix = "geno")

# Default: output prefix matches input
plink_to_blupf90("pigs")
# Creates: pigs.txt, pigs.map, pigs.bim
```

**Output files:**
- `<prefix>.txt` - BLUPF90 continuous genotype format (ID + 0/1/2 string)
- `<prefix>.map` - SNP position file (SNP_ID CHR POS)
- `<prefix>.bim` - Allele information for reverse conversion

### BLUPF90 → PLINK

```r
# Convert BLUPF90 files back to PLINK PED/MAP
# Auto-detects: geno.txt, geno.map, geno.bim

blupf90_to_plink("geno")
# Creates: geno.ped, geno.map

# Custom output prefix
blupf90_to_plink("geno", out_prefix = "output")
# Creates: output.ped, output.map
```

**Allele Restoration:**
- If `.bim` file exists: Original alleles (A/C/G/T) are restored
- If no `.bim` file: Generic alleles (A/B) are used

### Round-trip Example

```r
# Complete round-trip conversion
plink_to_blupf90("original", out_prefix = "blupf90_geno")
blupf90_to_plink("blupf90_geno", out_prefix = "reconstructed")

# Verify: reconstructed.ped should match original.ped
```

### Function Parameters

#### plink_to_blupf90()

| Parameter | Description | Default |
|-----------|-------------|---------|
| `prefix` | Input PLINK file prefix | (required) |
| `out_prefix` | Output file prefix | same as input |
| `exe` | Path to PLINK executable | auto-detect |
| `plink_dir` | Directory containing PLINK | NULL |
| `verbose` | Print progress messages | TRUE |

#### blupf90_to_plink()

| Parameter | Description | Default |
|-----------|-------------|---------|
| `blupf90_prefix` | Input BLUPF90 file prefix | (required) |
| `geno_file` | Explicit genotype file path | `<prefix>.txt` |
| `map_file` | Explicit map file path | `<prefix>.map` |
| `bim_file` | BIM file for allele info | `<prefix>.bim` (if exists) |
| `out_prefix` | Output PLINK prefix | same as input |
| `allele_letters` | Alleles when no BIM | c("A", "B") |
| `verbose` | Print progress messages | TRUE |

## PLINK Detection Order

plinkR automatically finds PLINK in this order:

1. User-provided `exe` parameter
2. User-provided `plink_dir` directory
3. Cached option from previous session
4. System PATH (plink2, then plink)
5. Bundled PLINK from linkbreedeR package

## Acknowledgements

The functionality for reading PLINK binary genotype files (BED/BIM/FAM) in plinkR is based on C++ code adapted from the plink2R project by Gabriel Abraham:

- Gabriel Abraham plink2R: Read PLINK BED/BIM/FAM files into R.
- GitHub repository: https://github.com/gabraham/plink2R
