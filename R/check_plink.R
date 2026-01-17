#' Check PLINK Installation Status
#'
#' Verifies whether PLINK is installed and accessible, displays version information.
#' Checks both system-wide PLINK and bundled PLINK from linkbreedeR package.
#'
#' @param verbose Logical; print detailed information (default TRUE)
#'
#' @return Invisibly returns a list with components:
#'   - found: Logical, whether PLINK is available
#'   - exe: Character, path to PLINK executable
#'   - source: Character, "system", "bundled", or "none"
#'   - version: Character, PLINK version
#'
#' @details
#' This function checks for PLINK availability in the following order:
#'   1. System PATH (plink2, then plink)
#'   2. Bundled PLINK from linkbreedeR package
#'
#' If PLINK is not found, provides installation instructions.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' check_plink()
#' }
check_plink <- function(verbose = TRUE) {
  result <- .check_plink_availability(verbose = verbose)
  invisible(result)
}
