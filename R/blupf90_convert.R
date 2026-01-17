#' @title PLINK <-> BLUPF90 Format Conversion
#'
#' @description
#' Bidirectional conversion between PLINK (PED/BED) and BLUPF90 genotype formats.
#' Uses plinkR's automatic PLINK path detection (system PATH or bundled from linkbreedeR).
#'
#' @name blupf90-conversion
NULL

# ============================================================================
# Rcpp accelerated parsing (compiled on load if available)
# ============================================================================

.blupf90_rcpp_available <- FALSE

.compile_blupf90_rcpp <- function() {
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    return(FALSE)
  }
  
  tryCatch({
    Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace Rcpp;

// Fast parser for plink .raw file to BLUPF90 format
// [[Rcpp::export]]
int parse_and_write_plink_raw_cpp(
  const std::string& raw_file,
  const std::string& out_file,
  const std::string& format,
  bool replace_9_with_5 = true
) {
  std::ifstream infile(raw_file.c_str());
  if (!infile.is_open()) {
    Rcpp::stop("Cannot open input file: " + raw_file);
  }
  
  std::ofstream outfile(out_file.c_str());
  if (!outfile.is_open()) {
    infile.close();
    Rcpp::stop("Cannot open output file: " + out_file);
  }
  
  std::string line;
  int line_count = 0;
  bool skip_header = true;
  
  while (std::getline(infile, line)) {
    if (skip_header) {
      skip_header = false;
      continue;
    }
    if (line.empty()) continue;
    
    std::istringstream iss(line);
    std::vector<std::string> fields;
    std::string field;
    while (iss >> field) {
      fields.push_back(field);
    }
    if (fields.size() < 7) continue;
    
    std::string sample_id = fields[1];
    
    if (format == "continuous") {
      outfile << sample_id << " ";
      for (size_t i = 6; i < fields.size(); ++i) {
        std::string geno = fields[i];
        if (replace_9_with_5 && (geno == "9" || geno == "NA")) {
          outfile << "5";
        } else {
          outfile << geno;
        }
      }
      outfile << std::endl;
    } else {
      outfile << sample_id;
      for (size_t i = 6; i < fields.size(); ++i) {
        std::string geno = fields[i];
        if (replace_9_with_5 && (geno == "9" || geno == "NA")) {
          outfile << " 5";
        } else {
          outfile << " " << geno;
        }
      }
      outfile << std::endl;
    }
    line_count++;
  }
  
  infile.close();
  outfile.close();
  return line_count;
}

// Fast PED file writer for BLUPF90 to PLINK conversion
// [[Rcpp::export]]
int write_ped_fast_cpp(
  Rcpp::StringVector sample_ids,
  Rcpp::IntegerMatrix geno_matrix,
  Rcpp::StringVector allele1,
  Rcpp::StringVector allele2,
  std::string out_file
) {
  int n_samples = geno_matrix.nrow();
  int n_snps = geno_matrix.ncol();
  
  std::ofstream outfile(out_file.c_str());
  if (!outfile.is_open()) {
    Rcpp::stop("Cannot open output file: " + out_file);
  }
  
  for (int i = 0; i < n_samples; ++i) {
    std::string sid = Rcpp::as<std::string>(sample_ids[i]);
    outfile << sid << " " << sid << " 0 0 0 -9";
    
    for (int j = 0; j < n_snps; ++j) {
      int g = geno_matrix(i, j);
      std::string a1 = Rcpp::as<std::string>(allele1[j]);
      std::string a2 = Rcpp::as<std::string>(allele2[j]);
      
      std::string g1, g2;
      if (g == 0) {
        g1 = a2; g2 = a2;
      } else if (g == 1) {
        g1 = a1; g2 = a2;
      } else if (g == 2) {
        g1 = a1; g2 = a1;
      } else {
        g1 = "0"; g2 = "0";
      }
      outfile << " " << g1 << " " << g2;
    }
    outfile << std::endl;
  }
  
  outfile.close();
  return n_samples;
}
    ')
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

# ============================================================================
# plink_to_blupf90: Convert PLINK to BLUPF90 format
# ============================================================================

#' Convert PLINK format to BLUPF90 format
#'
#' Converts PLINK genotype data (PED/MAP or BED/BIM/FAM) to BLUPF90 format.
#' Uses plinkR's automatic PLINK path detection.
#'
#' @param prefix Character. Base name of PLINK files (without extension).
#' @param out_prefix Character. Output file prefix. Creates <prefix>.txt, <prefix>.map, <prefix>.bim.
#'        Default: same as input prefix.
#' @param exe Character. Path to plink executable (optional, uses plinkR auto-detection if NULL)
#' @param plink_dir Character. Directory containing plink (optional)
#' @param allow_extra_chr Logical. Allow non-standard chromosome IDs. Default: TRUE
#' @param allow_no_sex Logical. Allow unknown sex. Default: TRUE
#' @param nonfounders Logical. Include non-founders. Default: TRUE
#' @param tempdir Character. Temporary file directory. Default: current directory
#' @param verbose Logical. Print progress. Default: TRUE
#'
#' @return Invisibly returns a list with geno_file, map_file, bim_file, n_snps, n_individuals
#'
#' @examples
#' \dontrun{
#' # Basic usage - output files match input prefix
#' plink_to_blupf90("mydata")
#' # Creates: mydata.txt, mydata.map, mydata.bim
#'
#' # Custom output prefix
#' plink_to_blupf90("mydata", out_prefix = "geno")
#' # Creates: geno.txt, geno.map, geno.bim
#' }
#'
#' @export
plink_to_blupf90 <- function(
  prefix,
  out_prefix = NULL,
  exe = NULL,
  plink_dir = NULL,
  allow_extra_chr = TRUE,
  allow_no_sex = TRUE,
  nonfounders = TRUE,
  tempdir = ".",
  verbose = TRUE
) {
  
  # Default out_prefix to input prefix
  if (is.null(out_prefix)) {
    out_prefix <- prefix
  }
  
  # Initialize Rcpp if not done

if (!.blupf90_rcpp_available) {
    .blupf90_rcpp_available <<- .compile_blupf90_rcpp()
  }
  
  # Get PLINK path using plinkR's detection
  plink_exe <- .get_plink_path(verbose = FALSE, user_exe = exe, user_plink_dir = plink_dir)
  if (is.null(plink_exe)) {
    stop(
      "PLINK not found. Use check_plink() to diagnose.\n",
      "Options: exe='/path/to/plink' or install linkbreedeR for bundled PLINK."
    )
  }
  
  # Normalize output prefix (remove .txt extension if present)
  out_prefix <- sub("\\.txt$", "", out_prefix, ignore.case = TRUE)
  geno_out_file <- paste0(out_prefix, ".txt")
  map_out_file <- paste0(out_prefix, ".map")
  bim_out_file <- paste0(out_prefix, ".bim")
  
  if (verbose) cat("\n========== PLINK to BLUPF90 Conversion ==========\n")
  if (verbose) cat("Using PLINK:", plink_exe, "\n")
  
  # Detect input files
  has_bed <- file.exists(paste0(prefix, ".bed"))
  has_ped <- file.exists(paste0(prefix, ".ped"))
  
  if (!has_ped && !has_bed) {
    stop("Neither ", prefix, ".ped nor ", prefix, ".bed found!")
  }
  
  # Determine MAP/BIM file
  map_file <- if (has_bed) paste0(prefix, ".bim") else paste0(prefix, ".map")
  if (!file.exists(map_file)) stop("Map file not found: ", map_file)
  
  # Read MAP/BIM for chromosome detection
  if (verbose) cat("Step 1: Detecting chromosome range from", basename(map_file), "...\n")
  
  map_data <- read.table(map_file, colClasses = c("integer", "character", "numeric", "numeric"))
  max_chr <- max(as.numeric(map_data[[1]]), na.rm = TRUE)
  
  if (verbose) {
    cat("  Max chromosome:", max_chr, "\n")
    cat("  SNP count:", nrow(map_data), "\n")
  }
  
  # Build PLINK command
  if (verbose) cat("Step 2: Running plink --recode A ...\n")
  
  temp_prefix <- file.path(tempdir, ".plink_temp")
  
  plink_args <- c(
    "--file", prefix,
    "--allow-extra-chr",
    "--chr-set", max_chr,
    "--allow-no-sex",
    "--nonfounders",
    "--recode", "A",
    "--out", temp_prefix
  )
  
  if (!allow_extra_chr) plink_args <- plink_args[plink_args != "--allow-extra-chr"]
  if (!allow_no_sex) plink_args <- plink_args[plink_args != "--allow-no-sex"]
  if (!nonfounders) plink_args <- plink_args[plink_args != "--nonfounders"]
  
  # Run PLINK
  cmd <- paste(plink_exe, paste(plink_args, collapse = " "))
  result <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  if (result != 0) stop("PLINK command failed. Check input files.")
  
  raw_file <- paste0(temp_prefix, ".raw")
  if (!file.exists(raw_file)) stop("PLINK output file not found: ", raw_file)
  
  # Process .raw file
  if (verbose) cat("Step 3: Processing plink output...\n")
  
  if (.blupf90_rcpp_available && exists("parse_and_write_plink_raw_cpp")) {
    if (verbose) cat("  Using Rcpp acceleration...\n")
    n_processed <- parse_and_write_plink_raw_cpp(
      raw_file = raw_file,
      out_file = geno_out_file,
      format = "continuous",
      replace_9_with_5 = TRUE
    )
    if (verbose) cat("  Processed:", n_processed, "individuals\n")
  } else {
    # Pure R fallback
    if (verbose) cat("  Using pure R (slower)...\n")
    raw_lines <- readLines(raw_file)
    data_lines <- raw_lines[-1]
    
    geno_matrix <- do.call(rbind, strsplit(data_lines, " "))
    sample_ids <- geno_matrix[, 2]
    geno_data <- geno_matrix[, 7:ncol(geno_matrix), drop = FALSE]
    
    geno_data[geno_data == "9"] <- "5"
    geno_data[geno_data == "NA"] <- "5"
    
    geno_strings <- apply(geno_data, 1, paste, collapse = "")
    output_content <- paste(sample_ids, geno_strings)
    writeLines(output_content, geno_out_file)
  }
  
  if (!file.exists(geno_out_file)) stop("Failed to write output file")
  
  # Read .raw header for BIM generation
  raw_header_line <- readLines(raw_file, n = 1)
  
  # Create MAP file
  if (verbose) cat("Step 4: Creating BLUPF90 map and BIM files ...\n")
  
  map_output <- cbind(
    SNP_ID = map_data[[2]],
    CHR = map_data[[1]],
    POS = map_data[[4]]
  )
  
  map_text <- c(
    "SNP_ID CHR POS",
    apply(map_output, 1, function(x) paste(x, collapse = " "))
  )
  writeLines(map_text, map_out_file)
  
  # Create BIM file
  orig_bim_file <- paste0(prefix, ".bim")
  orig_ped_file <- paste0(prefix, ".ped")
  
  if (file.exists(orig_bim_file)) {
    file.copy(orig_bim_file, bim_out_file, overwrite = TRUE)
    if (verbose) cat("  Copied original BIM file\n")
  } else if (file.exists(orig_ped_file)) {
    # Generate BIM using PLINK --make-bed
    if (verbose) cat("  Generating BIM from PED file...\n")
    
    temp_bed_prefix <- file.path(tempdir, ".plink_temp_bed")
    bed_cmd <- paste(plink_exe, "--file", prefix,
                     "--allow-extra-chr --chr-set", max_chr,
                     "--allow-no-sex --nonfounders",
                     "--make-bed --out", temp_bed_prefix)
    bed_result <- system(bed_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    temp_bim <- paste0(temp_bed_prefix, ".bim")
    if (bed_result == 0 && file.exists(temp_bim)) {
      file.copy(temp_bim, bim_out_file, overwrite = TRUE)
      
      # Cleanup temp files
      for (ext in c(".bed", ".bim", ".fam", ".log", ".nosex")) {
        f <- paste0(temp_bed_prefix, ext)
        if (file.exists(f)) file.remove(f)
      }
      if (verbose) cat("  Generated BIM using PLINK --make-bed\n")
    } else {
      # Fallback: synthetic from .raw header
      header_fields <- strsplit(raw_header_line, "\\s+")[[1]]
      snp_fields <- header_fields[7:length(header_fields)]
      
      snp_ids <- sapply(snp_fields, function(s) {
        parts <- strsplit(s, "_")[[1]]
        if (length(parts) >= 2) paste(parts[1:(length(parts)-1)], collapse = "_") else s
      })
      allele_counted <- sapply(snp_fields, function(s) {
        parts <- strsplit(s, "_")[[1]]
        if (length(parts) >= 2) parts[length(parts)] else "B"
      })
      
      bim_df <- data.frame(
        CHR = map_data[[1]], SNP = snp_ids, CM = 0, POS = map_data[[4]],
        A1 = allele_counted, A2 = ifelse(allele_counted == "A", "B", "A"),
        stringsAsFactors = FALSE
      )
      write.table(bim_df, file = bim_out_file, quote = FALSE, sep = "\t",
                  row.names = FALSE, col.names = FALSE)
      if (verbose) cat("  Generated BIM from .raw header (synthetic)\n")
    }
  } else {
    # Fallback from .raw header
    header_fields <- strsplit(raw_header_line, "\\s+")[[1]]
    snp_fields <- header_fields[7:length(header_fields)]
    
    snp_ids <- sapply(snp_fields, function(s) {
      parts <- strsplit(s, "_")[[1]]
      if (length(parts) >= 2) paste(parts[1:(length(parts)-1)], collapse = "_") else s
    })
    allele_counted <- sapply(snp_fields, function(s) {
      parts <- strsplit(s, "_")[[1]]
      if (length(parts) >= 2) parts[length(parts)] else "B"
    })
    
    bim_df <- data.frame(
      CHR = map_data[[1]], SNP = snp_ids, CM = 0, POS = map_data[[4]],
      A1 = allele_counted, A2 = ifelse(allele_counted == "A", "B", "A"),
      stringsAsFactors = FALSE
    )
    write.table(bim_df, file = bim_out_file, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
    if (verbose) cat("  Generated BIM from .raw header (synthetic)\n")
  }
  
  if (verbose) {
    cat("  Map file:", map_out_file, "\n")
    cat("  BIM file:", bim_out_file, "\n")
  }
  
  # Cleanup
  if (verbose) cat("Step 5: Cleaning up ...\n")
  
  for (ext in c(".raw", ".log", ".nosex")) {
    f <- paste0(temp_prefix, ext)
    if (file.exists(f)) file.remove(f)
  }
  
  n_individuals <- as.integer(system(paste("wc -l <", geno_out_file), intern = TRUE))
  
  if (verbose) {
    cat("\n========== Conversion Complete ==========\n")
    cat("✓ Genotype:", geno_out_file, "\n")
    cat("✓ Map:     ", map_out_file, "\n")
    cat("✓ BIM:     ", bim_out_file, "\n")
    cat("Individuals:", n_individuals, "\n")
    cat("SNPs:       ", nrow(map_data), "\n")
    cat("==========================================\n\n")
  }
  
  invisible(list(
    geno_file = geno_out_file,
    map_file = map_out_file,
    bim_file = bim_out_file,
    n_snps = nrow(map_data),
    n_individuals = n_individuals,
    max_chr = max_chr
  ))
}

# ============================================================================
# blupf90_to_plink: Convert BLUPF90 back to PLINK format
# ============================================================================

#' Convert BLUPF90 format to PLINK PED/MAP
#'
#' Converts BLUPF90 genotype files back to PLINK format.
#' Uses BIM file (if available) to restore original allele letters.
#'
#' @param blupf90_prefix Character. Prefix for BLUPF90 files. Auto-detects <prefix>.txt, .map, .bim
#' @param geno_file Character. Explicit genotype file path (optional, overrides prefix)
#' @param map_file Character. Explicit map file path (optional, overrides prefix)
#' @param bim_file Character. BIM file for allele restoration (optional, auto-detected)
#' @param out_prefix Character. Output PLINK file prefix (default: same as blupf90_prefix)
#' @param allele_letters Character vector of length 2. Alleles to use when no BIM. Default: c("A","B")
#' @param ref_allele Character. Which allele is reference: "A1" or "A2". Default: "A1"
#' @param verbose Logical. Print progress. Default: TRUE
#'
#' @return Invisibly returns a list with ped_file, map_file, n_individuals, n_snps, used_bim
#'
#' @examples
#' \dontrun
#' # Simple usage with auto-detection
#' blupf90_to_plink("geno")
#' # Looks for: geno.txt, geno.map, geno.bim
#'
#' # Explicit files
#' blupf90_to_plink("geno", geno_file = "custom.txt", out_prefix = "output")
#' }
#'
#' @export
blupf90_to_plink <- function(
  blupf90_prefix,
  geno_file = NULL,
  map_file = NULL,
  bim_file = NULL,
  out_prefix = NULL,
  allele_letters = c("A", "B"),
  ref_allele = "A1",
  verbose = TRUE
) {
  
  # Initialize Rcpp if not done
  if (!.blupf90_rcpp_available) {
    .blupf90_rcpp_available <<- .compile_blupf90_rcpp()
  }
  
  # Auto-detect files from prefix
  if (is.null(geno_file)) geno_file <- paste0(blupf90_prefix, ".txt")
  if (is.null(map_file)) map_file <- paste0(blupf90_prefix, ".map")
  if (is.null(bim_file)) {
    bim_candidate <- paste0(blupf90_prefix, ".bim")
    if (file.exists(bim_candidate)) bim_file <- bim_candidate
  }
  if (is.null(out_prefix)) out_prefix <- blupf90_prefix
  
  # Validate inputs
  if (!file.exists(geno_file)) stop("Genotype file not found: ", geno_file)
  if (!file.exists(map_file)) stop("Map file not found: ", map_file)
  if (!ref_allele %in% c("A1", "A2")) stop("ref_allele must be 'A1' or 'A2'")
  if (length(allele_letters) != 2) stop("allele_letters must have length 2")
  
  a_ref <- allele_letters[1]
  a_alt <- allele_letters[2]
  
  if (verbose) cat("\n========== BLUPF90 → PLINK (PED/MAP) =========\n")
  
  # Read map file
  first_line <- readLines(map_file, n = 1)
  first_fields <- strsplit(trimws(first_line), "\\s+")[[1]]
  has_header <- !grepl("^[0-9]+$", first_fields[1]) ||
                tolower(first_fields[1]) %in% c("snp_id", "snp", "chr", "chrom")
  
  map <- read.table(map_file, header = has_header, stringsAsFactors = FALSE)
  
  # Normalize columns
  if (ncol(map) >= 4) {
    names(map)[1:4] <- c("CHR", "SNP_ID", "CM", "POS")
    map <- map[, c("SNP_ID", "CHR", "POS")]
  } else if (ncol(map) >= 3) {
    names(map)[1:3] <- c("SNP_ID", "CHR", "POS")
    map <- map[, 1:3]
  } else {
    stop("Map file must have at least 3 columns")
  }
  
  n_snps <- nrow(map)
  if (verbose) cat("SNPs:", n_snps, "\n")
  
  # Prepare output paths
  map_out <- file.path(dirname(out_prefix), paste0(basename(out_prefix), ".map"))
  ped_out <- file.path(dirname(out_prefix), paste0(basename(out_prefix), ".ped"))
  
  # Write MAP file
  snp_ids <- as.character(map$SNP_ID)
  map_df <- data.frame(
    CHR = map$CHR, SNP_ID = snp_ids, CM = 0, POS = map$POS,
    stringsAsFactors = FALSE
  )
  write.table(map_df, file = map_out, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  if (verbose) cat("MAP file:", map_out, "\n")
  
  # Determine alleles (from BIM or default)
  use_bim <- !is.null(bim_file) && file.exists(bim_file)
  
  if (use_bim) {
    if (verbose) cat("Using BIM to restore allele letters...\n")
    bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
    colnames(bim) <- c("CHR", "SNP", "CM", "POS", "A1", "A2")
    
    # Match SNPs
    bim_idx <- match(snp_ids, bim$SNP)
    
    if (any(is.na(bim_idx))) {
      warning("Some SNPs not found in BIM, using default alleles")
      missing <- is.na(bim_idx)
      allele_A1 <- ifelse(missing, a_ref, bim$A1[bim_idx])
      allele_A2 <- ifelse(missing, a_alt, bim$A2[bim_idx])
    } else {
      allele_A1 <- bim$A1[bim_idx]
      allele_A2 <- bim$A2[bim_idx]
    }
  } else {
    allele_A1 <- rep(a_ref, n_snps)
    allele_A2 <- rep(a_alt, n_snps)
  }
  
  # Read genotype file
  geno_lines <- readLines(geno_file)
  n_individuals <- length(geno_lines)
  
  # Parse genotypes
  sample_ids <- character(n_individuals)
  geno_matrix <- matrix(NA_integer_, nrow = n_individuals, ncol = n_snps)
  
  for (i in seq_len(n_individuals)) {
    parts <- strsplit(geno_lines[i], "\\s+")[[1]]
    sample_ids[i] <- parts[1]
    
    if (length(parts) == 2) {
      # Continuous format: ID + string
      geno_str <- parts[2]
      geno_vals <- as.integer(strsplit(geno_str, "")[[1]])
    } else {
      # Standard format: ID + separated values
      geno_vals <- as.integer(parts[-1])
    }
    
    geno_vals[geno_vals == 5] <- NA_integer_
    geno_matrix[i, ] <- geno_vals[1:n_snps]
  }
  
  # Write PED file
  if (.blupf90_rcpp_available && exists("write_ped_fast_cpp")) {
    if (verbose) cat("Using Rcpp fast writer for PED...\n")
    geno_matrix[is.na(geno_matrix)] <- -1
    write_ped_fast_cpp(sample_ids, geno_matrix, allele_A1, allele_A2, ped_out)
  } else {
    if (verbose) cat("Using pure R for PED...\n")
    
    ped_lines <- character(n_individuals)
    for (i in seq_len(n_individuals)) {
      sid <- sample_ids[i]
      geno_row <- geno_matrix[i, ]
      
      allele_pairs <- character(n_snps)
      for (j in seq_len(n_snps)) {
        g <- geno_row[j]
        a1 <- allele_A1[j]
        a2 <- allele_A2[j]
        
        if (is.na(g) || g < 0 || g > 2) {
          allele_pairs[j] <- "0 0"
        } else if (g == 0) {
          allele_pairs[j] <- paste(a2, a2)
        } else if (g == 1) {
          allele_pairs[j] <- paste(a1, a2)
        } else {
          allele_pairs[j] <- paste(a1, a1)
        }
      }
      
      ped_lines[i] <- paste(c(sid, sid, "0", "0", "0", "-9", allele_pairs), collapse = " ")
    }
    writeLines(ped_lines, ped_out)
  }
  
  if (verbose) {
    cat("PED file:", ped_out, "\n")
    cat("\n========== Conversion Complete =========\n\n")
  }
  
  invisible(list(
    ped_file = ped_out,
    map_file = map_out,
    n_individuals = n_individuals,
    n_snps = n_snps,
    used_bim = use_bim
  ))
}
