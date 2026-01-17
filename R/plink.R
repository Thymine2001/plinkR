#' Run PLINK/PLINK2 from R
#'
#' Execute PLINK commands, capture outputs, and automatically read result files.
#' Without --out, infers output prefix from PLINK's log. Reads genotype files
#' (BED/BIM/FAM or PED/MAP) and report tables (--missing, --pca, etc.) into R.
#'
#' This function will automatically use bundled PLINK from linkbreedeR if system PLINK
#' is not found.
#'
#' @param cmd Character string with PLINK command (e.g., "--bfile input --make-bed --out output")
#' @param read_outputs What to read: "all" (read all available formats),
#'   "none" (skip reading), "binary" (BED/BIM/FAM), "text" (PED/MAP), or "auto" (infer from --recode/--make-bed)
#' @param use_fread Logical; use data.table::fread for faster reading (default TRUE)
#' @param exe Character; explicit path to PLINK executable. If NULL, searches PATH and bundled versions
#' @param plink_dir Character; directory containing plink executable
#'
#' @return A list of class plink_result with components:
#'   \itemize{
#'     \item cmd: Full command executed
#'     \item exe: Path to PLINK executable used
#'     \item args: Parsed argument vector
#'     \item workdir: Working directory
#'     \item exit_status: Exit code (0 = success)
#'     \item error: Error message if failed, NULL if succeeded
#'     \item start_time, end_time, elapsed: Timing information
#'     \item log: Contents of .log file
#'     \item log_path: Path to .log file
#'     \item outputs: Vector of all output files discovered
#'     \item report: Named list of report tables (imiss, lmiss, eigenvec, etc.)
#'     \item out_geno: List with bed, bim, fam, ped, map (genotype data)
#'     \item session_info: R session information
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' res <- plink("--bfile mydata --missing")
#' res$report$imiss  # sample missingness
#' res$report$lmiss  # SNP missingness
#'
#' # With explicit PLINK path
#' res <- plink("--bfile mydata --make-bed --out filtered",
#'              exe = "/usr/local/bin/plink")
#' }
#'
#' @export
plink <- function(cmd, read_outputs = c("all", "none", "binary", "text", "auto"), 
                  use_fread = TRUE, exe = NULL, plink_dir = NULL) {
  if (!is.character(cmd) || length(cmd) != 1) {
    stop("plink() requires a single command string")
  }

  read_outputs <- match.arg(read_outputs)

  # Determine PLINK executable path with fallback to linkbreedeR
  if (is.null(exe) && is.null(plink_dir)) {
    exe <- .get_plink_path(verbose = FALSE)
    if (is.null(exe)) {
      stop(
        "PLINK executable not found. Please:\n",
        "  1. Install PLINK from: https://www.cog-genomics.org/plink/\n",
        "  2. Add PLINK to your PATH, or\n",
        "  3. Install linkbreedeR package (provides bundled PLINK):\n",
        "     devtools::install_github('Thymine2001/linkbreedeR'), or\n",
        "  4. Pass plink_dir = '/path/to/plink/dir', or\n",
        "  5. Pass exe = '/full/path/to/plink'\n\n",
        "Check installation status with: check_plink()"
      )
    }
  } else {
    # User provided explicit exe or plink_dir, use validation logic
    if (!is.null(exe)) {
      if (!file.exists(exe)) {
        stop("Specified exe '", exe, "' does not exist")
      }
    } else if (!is.null(plink_dir)) {
      plink_dir <- normalizePath(plink_dir, mustWork = FALSE)
      plink2_candidate <- file.path(plink_dir, "plink2")
      plink_candidate <- file.path(plink_dir, "plink")
      if (file.exists(plink2_candidate)) {
        exe <- plink2_candidate
      } else if (file.exists(plink_candidate)) {
        exe <- plink_candidate
      } else {
        stop("No plink or plink2 executable found in '", plink_dir, "'")
      }
    }
  }

  args <- .split_shell_tokens(cmd)

  out_idx <- which(args == "--out")
  out_prefix <- if (length(out_idx) == 1 && out_idx < length(args)) args[out_idx + 1] else NULL

  # Check if this is a help/version command that should print directly to terminal
  is_interactive_cmd <- any(args %in% c("--help", "-h", "--version", "-v", "--help-full"))
  
  start_time <- Sys.time()
  if (is_interactive_cmd) {
    # For help/version commands, print directly to terminal
    exit_status <- suppressWarnings(system2(exe, args = args, stdout = "", stderr = ""))
    out <- character()  # No output to capture
  } else {
    # For normal commands, capture output for parsing
    out <- suppressWarnings(tryCatch(
      system2(exe, args = args, stdout = TRUE, stderr = TRUE),
      error = function(e) {
        structure(character(), status = 1L, msg = paste("system2 error:", conditionMessage(e)))
      }
    ))
    exit_status <- attr(out, "status"); if (is.null(exit_status)) exit_status <- 0L
  }
  end_time <- Sys.time()

  # If no explicit --out, infer prefix from stdout ('Logging to ...') or fallback
  if (is.null(out_prefix)) {
    combined <- paste(out, collapse = "\n")
    m <- regexec("Logging to ([^\\r\\n]+?)\\.log\\.", combined, perl = TRUE)
    mm <- regmatches(combined, m)
    if (length(mm) && length(mm[[1]]) == 2) {
      log_file <- mm[[1]][2]
      # Ensure simple filename (no path)
      log_name <- basename(log_file)
      out_prefix <- sub("\\.log$", "", log_name)
    } else {
      # Fallback to executable-based default
      default_base <- if (grepl("plink2", basename(exe), fixed = TRUE)) "plink2" else "plink"
      out_prefix <- default_base
    }
  }

  # Prepare error information
  error_message <- NULL
  
  # Check for errors and print detailed error messages
  if (exit_status != 0L) {
    cat("\n==== PLINK ERROR (exit status ", exit_status, ") ====\n", sep = "")
    cat("Command: ", paste(c(exe, args), collapse = " "), "\n\n", sep = "")
    
    # Print stdout/stderr combined output
    if (length(out) > 0) {
      cat("Output:\n")
      cat(paste(out, collapse = "\n"), "\n\n")
    }
    
    # Try to read log file if it exists
    log_path_temp <- if (!is.null(out_prefix)) paste0(out_prefix, ".log") else NULL
    if (!is.null(log_path_temp) && file.exists(log_path_temp)) {
      log_content <- readLines(log_path_temp, warn = FALSE)
      cat("Log file (", log_path_temp, "):\n", sep = "")
      cat(paste(log_content, collapse = "\n"), "\n")
    }
    
    cat("=============================================\n\n")
    error_message <- paste0("PLINK command failed with exit status ", exit_status, 
                           ". See error messages above.")
    # 中断脚本但不显示额外的错误框
    stop(call. = FALSE)
  }

  log_path <- if (!is.null(out_prefix)) paste0(out_prefix, ".log") else NULL
  log <- if (!is.null(log_path) && file.exists(log_path)) readLines(log_path, warn = FALSE) else character()

  outputs <- if (!is.null(out_prefix)) .discover_outputs(out_prefix) else character()

  out_format <- "all"
  if (read_outputs != "auto") {
    out_format <- read_outputs
  } else {
    # Auto-detect format from command arguments
    if ("--recode" %in% args) {
      out_format <- "text"
    } else if ("--make-bed" %in% args) {
      out_format <- "binary"
    } else {
      # Default to all if neither recode nor make-bed specified
      out_format <- "all"
    }
  }

  out_geno <- if (!is.null(out_prefix) && read_outputs != "none") {
    read_plink(out_prefix, format = out_format, use_fread = use_fread)
  } else {
    list(bed = NULL, fam = NULL, bim = NULL, ped = NULL, map = NULL)
  }

  # Collect reports based on args and discovered outputs
  report <- if (!is.null(out_prefix)) {
    .collect_reports(out_prefix, args, use_fread = use_fread,
                     start_time = start_time, end_time = end_time)
  } else {
    list()
  }

  structure(
    list(
      cmd = paste(c(exe, args), collapse = " "),
      exe = exe,
      args = args,
      workdir = getwd(),
      outputs = outputs,
      report = report,
      out_geno = out_geno,
      error = error_message,
      start_time = start_time,
      end_time = end_time,
      elapsed = as.numeric(difftime(end_time, start_time, units = "secs")),
      log = log,
      log_path = log_path,
      session_info = utils::sessionInfo()
    ),
    class = "plink_result"
  )
}
