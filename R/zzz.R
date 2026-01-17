# Package initialization and PLINK detection
.onLoad <- function(lib, pkg) {
  # Try to find PLINK executable in system PATH
  # Use direct Sys.which() instead of .get_plink_path() to avoid dependency issues
  exe <- Sys.which("plink2")
  if (exe == "") {
    exe <- Sys.which("plink")
  }
  
  if (exe != "") {
    # Store the exe path as default option
    options(plinkR_plink_exe = exe)
  } else {
    options(plinkR_plink_exe = "")
  }
}

.onAttach <- function(lib, pkg) {
  # Show message only when attaching (library()), not when loading as dependency
  plink_status <- .check_plink_availability(verbose = FALSE)
  
  if (!plink_status$found) {
    packageStartupMessage(
      "\u2717 PLINK not found\n\n",
      "Options:\n",
      "  1. Install PLINK from: https://www.cog-genomics.org/plink/\n",
      "  2. Install linkbreedeR package for bundled PLINK: devtools::install_github('Thymine2001/linkbreedeR')\n",
      "  3. Provide path via: plink(..., exe='/path/to/plink')\n",
      "  4. Provide dir via:  plink(..., plink_dir='/path/to/dir')\n\n",
      "Run check_plink() for details."
    )
  } else {
    packageStartupMessage(
      "plinkR: PLINK found at ", plink_status$exe, " (", plink_status$source, ")"
    )
  }
}

