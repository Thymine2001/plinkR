#' Check PLINK Installation Status
#'
#' Verifies whether PLINK is installed and accessible, displays version information.
#'
#' @param verbose Logical; print detailed information (default TRUE)
#'
#' @return Invisibly returns TRUE if PLINK found, FALSE otherwise
#' @export
#'
#' @examples
#' \dontrun{
#' check_plink()
#' }
check_plink <- function(verbose = TRUE) {
  exe <- getOption("plinkR_plink_exe", default = "")
  if (exe == "") {
    exe <- Sys.which("plink2")
    if (exe == "") exe <- Sys.which("plink")
  }
  
  if (exe == "") {
    if (verbose) {
      cat("\u2717 PLINK not found\n\n")
      cat("Installation instructions:\n")
      cat("  1. Download from: https://www.cog-genomics.org/plink/\n")
      cat("  2. Extract and add to PATH, or\n")
      cat("  3. Use exe= or plink_dir= parameter\n\n")
      cat("Usage examples:\n")
      cat("  plink('--bfile data --make-bed', exe='/path/to/plink')\n")
      cat("  plink('--bfile data --make-bed', plink_dir='/path/to/dir')\n")
    }
    invisible(FALSE)
  } else {
    if (verbose) {
      cat("\u2713 PLINK found at:", exe, "\n")
      ver <- tryCatch(
        system2(exe, "--version", stdout = TRUE, stderr = TRUE),
        error = function(e) "unknown"
      )
      cat("  Version:", paste(ver, collapse = " "), "\n")
    }
    invisible(TRUE)
  }
}
