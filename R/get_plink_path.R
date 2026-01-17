#' Get PLINK Executable Path with Fallback
#'
#' Attempts to locate PLINK executable. First searches the system PATH,
#' then falls back to bundled PLINK from linkbreedeR package if available.
#'
#' @param verbose Logical; print status messages (default TRUE)
#' @param user_exe Character; explicit path to PLINK executable (optional)
#' @param user_plink_dir Character; directory containing PLINK executable (optional)
#'
#' @return Character string with path to PLINK executable, or NULL if not found
#'
#' @details
#' Search order:
#'   1. User-provided exe parameter
#'   2. User-provided plink_dir parameter
#'   3. Stored option from .onLoad
#'   4. System PATH (plink2, then plink)
#'   5. linkbreedeR bundled plink
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'   path <- .get_plink_path()
#'   path <- .get_plink_path(user_exe = "/usr/bin/plink")
#' }
#'
.get_plink_path <- function(verbose = FALSE, user_exe = NULL, user_plink_dir = NULL) {
  
  # 1. Check user-provided explicit exe
  if (!is.null(user_exe)) {
    if (file.exists(user_exe)) {
      if (verbose) cat("Using user-provided PLINK: ", user_exe, "\n")
      return(user_exe)
    } else {
      stop("Specified exe '", user_exe, "' does not exist")
    }
  }
  
  # 2. Check user-provided plink_dir
  if (!is.null(user_plink_dir)) {
    user_plink_dir <- normalizePath(user_plink_dir, mustWork = FALSE)
    plink2_candidate <- file.path(user_plink_dir, "plink2")
    plink_candidate <- file.path(user_plink_dir, "plink")
    
    if (file.exists(plink2_candidate)) {
      if (verbose) cat("Using PLINK from plink_dir: ", plink2_candidate, "\n")
      return(plink2_candidate)
    } else if (file.exists(plink_candidate)) {
      if (verbose) cat("Using PLINK from plink_dir: ", plink_candidate, "\n")
      return(plink_candidate)
    } else {
      stop("No plink or plink2 executable found in '", user_plink_dir, "'")
    }
  }
  
  # 3. Check stored option
  exe <- getOption("plinkR_plink_exe", default = "")
  if (exe != "" && file.exists(exe)) {
    if (verbose) cat("Using cached PLINK from option: ", exe, "\n")
    return(exe)
  }
  
  # 4. Search system PATH
  exe <- Sys.which("plink2")
  if (exe != "") {
    if (verbose) cat("Found PLINK in PATH: ", exe, "\n")
    return(exe)
  }
  
  exe <- Sys.which("plink")
  if (exe != "") {
    if (verbose) cat("Found PLINK in PATH: ", exe, "\n")
    return(exe)
  }
  
  # 5. Try linkbreedeR bundled plink
  tryCatch({
    if (requireNamespace("linkbreedeR", quietly = TRUE)) {
      plink_path <- linkbreedeR::get_tool_path("plink")
      if (file.exists(plink_path)) {
        if (verbose) cat("Using bundled PLINK from linkbreedeR: ", plink_path, "\n")
        return(plink_path)
      }
    }
  }, error = function(e) {
    # linkbreedeR not available or tool not found
    if (verbose) cat("linkbreedeR PLINK not available:", conditionMessage(e), "\n")
  })
  
  # Not found anywhere
  NULL
}

#' Check PLINK Availability with Fallback
#'
#' Enhanced version that reports both system and bundled PLINK availability
#'
#' @param verbose Logical; print detailed information
#'
#' @return List with components:
#'   - found: Logical, whether PLINK is available
#'   - exe: Character, path to PLINK executable
#'   - source: Character, source of PLINK ("system", "bundled", "none")
#'   - version: Character, PLINK version if available
#'
#' @keywords internal
.check_plink_availability <- function(verbose = TRUE) {
  exe <- .get_plink_path(verbose = FALSE)
  
  if (!is.null(exe)) {
    version <- tryCatch(
      system2(exe, "--version", stdout = TRUE, stderr = TRUE),
      error = function(e) "unknown"
    )
    version <- paste(version, collapse = " ")
    
    # Determine source
    source <- if (grepl("extbin", exe, fixed = TRUE)) "bundled" else "system"
    
    if (verbose) {
      cat("\u2713 PLINK found\n")
      cat("  Path: ", exe, "\n")
      cat("  Source: ", source, "\n")
      cat("  Version: ", version, "\n")
    }
    
    return(list(
      found = TRUE,
      exe = exe,
      source = source,
      version = version
    ))
  } else {
    if (verbose) {
      cat("\u2717 PLINK not found\n\n")
      cat("Options:\n")
      cat("  1. Install PLINK from: https://www.cog-genomics.org/plink/\n")
      cat("  2. Install linkbreedeR package for bundled PLINK: devtools::install_github(\"Thymine2001/linkbreedeR\")\n")
      cat("  3. Provide path via: plink(..., exe='/path/to/plink')\n")
      cat("  4. Provide dir via:  plink(..., plink_dir='/path/to/dir')\n")
    }
    
    return(list(
      found = FALSE,
      exe = NULL,
      source = "none",
      version = NA_character_
    ))
  }
}
