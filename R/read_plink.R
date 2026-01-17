#' Read PLINK outputs
#'
#' Internal function to read PLINK binary (BED/BIM/FAM) or text (PED/MAP) files.
#' Automatically detects format and assigns proper column names.
#'
#' @param root Character; file path prefix (without extension)
#' @param format Character; one of "auto", "binary", "text", or "all"
#' @param impute Character; imputation method ("none", "avg", or "random")
#' @param use_fread Logical; use data.table::fread (default TRUE)
#' @param verbose Logical; print progress (default FALSE)
#'
#' @return List with bed, bim, fam, ped, map components (NULL if not present)
#' @keywords internal
read_plink <- function(root, format = c("auto", "binary", "text", "all"), impute = c("none", "avg", "random"), use_fread = TRUE, verbose = FALSE) {
  proot <- path.expand(root)
  format <- match.arg(format)
  impute <- match.arg(impute)
  impute_int <- switch(impute, "none" = 0L, "avg" = 1L, "random" = 2L)

  bedfile <- paste0(proot, ".bed")
  famfile <- paste0(proot, ".fam")
  bimfile <- paste0(proot, ".bim")
  pedfile <- paste0(proot, ".ped")
  mapfile <- paste0(proot, ".map")

  result <- list(bed = NULL, fam = NULL, bim = NULL, ped = NULL, map = NULL)

  if (format == "auto") {
    if (file.exists(pedfile) || file.exists(mapfile)) {
      format <- "text"
    } else if (file.exists(bedfile) || file.exists(famfile) || file.exists(bimfile)) {
      format <- "binary"
    }
  }

  rd <- function(path) {
    if (use_fread) .fast_read_table(path) else .fast_read_table(path)
  }

  if (format %in% c("binary", "all")) {
    result$bim <- rd(bimfile)
    if (!is.null(result$bim) && ncol(result$bim) >= 6) {
      colnames(result$bim)[1:6] <- c("CHR", "SNP", "CM", "POS", "A1", "A2")
    }
    result$fam <- rd(famfile)
    if (!is.null(result$fam) && ncol(result$fam) >= 6) {
      colnames(result$fam)[1:6] <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")
    }
    if (file.exists(bedfile) && !is.null(result$fam)) {
      result$bed <- tryCatch(
        read_plink_cpp(bedfile, famfile, impute_int, verbose),
        error = function(e) NULL
      )
      if (!is.null(result$bed) && !is.null(result$bim) && !is.null(result$fam)) {
        colnames(result$bed) <- result$bim[, 2]
        rownames(result$bed) <- paste(result$fam[, 1], result$fam[, 2], sep = ":")
      }
    }
  }

  if (format %in% c("text", "all")) {
    result$map <- rd(mapfile)
    if (!is.null(result$map) && ncol(result$map) >= 4) {
      colnames(result$map)[1:4] <- c("CHR", "SNP", "CM", "POS")
    }
    result$ped <- rd(pedfile)
    if (!is.null(result$ped) && ncol(result$ped) >= 6) {
      # First 6 are FID IID PAT MAT SEX PHENO
      base_names <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")
      cn <- colnames(result$ped)
      cn[1:6] <- base_names
      # If MAP available and column count matches 6 + 2*nsnp, name genotype pairs
      if (!is.null(result$map)) {
        snp_ids <- result$map[, 2]
        expected_cols <- 6 + 2 * length(snp_ids)
        if (ncol(result$ped) == expected_cols) {
          geno_names <- as.vector(rbind(paste0(snp_ids, "_A1"), paste0(snp_ids, "_A2")))
          cn <- c(base_names, geno_names)
        }
      }
      colnames(result$ped) <- cn
    }
  }

  result
}
