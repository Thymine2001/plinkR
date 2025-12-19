# Package initialization and PLINK detection
.onLoad <- function(lib, pkg) {
  # Try to find PLINK executable
  plink2_exe <- Sys.which("plink2")
  plink_exe <- Sys.which("plink")
  
  if (plink2_exe != "" || plink_exe != "") {
    # PLINK found, store the exe path as default option
    found_exe <- if (plink2_exe != "") plink2_exe else plink_exe
    options(plinkR_plink_exe = found_exe)
  } else {
    options(plinkR_plink_exe = "")
  }
}

.onAttach <- function(lib, pkg) {
  # Show message only when attaching (library()), not when loading as dependency
  plink_exe <- getOption("plinkR_plink_exe", default = "")
  
  if (plink_exe == "") {
    packageStartupMessage(
      "plinkR: PLINK not found in PATH.\n",
      "  Install from https://www.cog-genomics.org/plink/ or use exe=/path/to/plink\n",
      "  Run check_plink() for details."
    )
  }
}
