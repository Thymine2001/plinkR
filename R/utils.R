# Fast table reader
.fast_read_table <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(tryCatch(as.data.frame(data.table::fread(path, showProgress = FALSE, ...)), error = function(e) NULL))
  }
  tryCatch(utils::read.table(path, header = FALSE, sep = "", stringsAsFactors = FALSE, ...), error = function(e) NULL)
}

.split_shell_tokens <- function(cmd) {
  toks <- tryCatch(
    scan(text = cmd, what = character(), quiet = TRUE, sep = " ", quote = "\"'"),
    error = function(e) character(0)
  )
  unname(toks)
}

.discover_outputs <- function(out_prefix) {
  dir <- dirname(out_prefix)
  base <- basename(out_prefix)
  files <- list.files(dir, full.names = TRUE)
  files <- files[startsWith(basename(files), base)]
  files[file.exists(files)]
}

# Heuristic report reader with header detection
.read_report_table <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
  if (is_gz) {
    # Use base reader for gz
    con <- gzfile(path, open = "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    first <- tryCatch(readLines(con, n = 1, warn = FALSE), error = function(e) "")
    has_header <- length(first) > 0 && grepl("[A-Za-z]", first)
    # Reopen to read all
    con <- gzfile(path, open = "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    return(tryCatch(utils::read.table(con, header = has_header, sep = "", stringsAsFactors = FALSE), error = function(e) NULL))
  } else {
    first <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) "")
    has_header <- length(first) > 0 && grepl("[A-Za-z]", first)
    return(.fast_read_table(path, header = has_header, ...))
  }
}

# Find candidate report files restricted by time window and pick newest per suffix
.find_report_files <- function(out_prefix, suffixes, start_time, end_time) {
  if (length(suffixes) == 0) return(character())
  cand <- unlist(lapply(suffixes, function(sfx) {
    c(paste0(out_prefix, ".", sfx), paste0(out_prefix, ".", sfx, ".gz"))
  }), use.names = FALSE)
  cand <- cand[file.exists(cand)]
  if (!length(cand)) return(character())

  if (!is.null(start_time) && !is.null(end_time)) {
    info <- file.info(cand)
    mt <- info$mtime
    in_win <- !is.na(mt) & mt >= (start_time - 5) & mt <= (end_time + 5)
    cand <- cand[in_win]
    if (!length(cand)) return(character())
  }

  # Group by suffix and choose newest
  keys <- sub(paste0("^", out_prefix, "\\."), "", cand)
  keys <- sub("\\.gz$", "", keys)
  split_idx <- split(seq_along(cand), keys)
  newest <- unlist(lapply(split_idx, function(idx) {
    group <- cand[idx]
    group[which.max(file.info(group)$mtime)]
  }), use.names = TRUE)
  unname(newest)
}

# Collect common PLINK report files based on args and discovered outputs
.collect_reports <- function(out_prefix, args, use_fread = TRUE, start_time = NULL, end_time = NULL) {
  # Known mapping from flags to suffixes
  flag_map <- list(
    "--missing" = c("imiss", "lmiss"),
    "--check-sex" = c("sexcheck"),
    "--hardy" = c("hwe"),
    "--het" = c("het"),
    "--freq" = c("frq", "frq.count", "frq.strat"),
    "--pca" = c("eigenvec", "eigenval"),
    "--ibc" = c("ibc"),
    "--genome" = c("genome"),
    "--indep" = c("prune.in", "prune.out"),
    "--indep-pairwise" = c("prune.in", "prune.out"),
    "--indep-pairphase" = c("prune.in", "prune.out"),
    "--mendel" = c("mendel", "imendel", "lmendel"),
    "--assoc" = c("assoc", "qassoc"),
    "--linear" = c("assoc.linear"),
    "--logistic" = c("assoc.logistic"),
    "--clump" = c("clumped")
  )

  requested_suffixes <- unique(unlist(flag_map[intersect(names(flag_map), args)], use.names = FALSE))

  # Also include any discovered outputs with known suffixes
  discovered <- .discover_outputs(out_prefix)
  known_suffixes <- unique(unlist(flag_map, use.names = FALSE))
  discovered_suffixes <- character()
  if (length(discovered)) {
    bases <- basename(discovered)
    m <- regexec(paste0("^", basename(out_prefix), "\\.(.+)$"), bases)
    parts <- regmatches(bases, m)
    for (p in parts) {
      if (length(p) == 2 && p[2] %in% known_suffixes) {
        discovered_suffixes <- c(discovered_suffixes, p[2])
      }
    }
  }

  suffixes <- unique(c(requested_suffixes, discovered_suffixes))
  files <- .find_report_files(out_prefix, suffixes, start_time, end_time)
  res <- list()
  for (fp in files) {
    key <- sub(paste0("^", out_prefix, "\\.|\\.gz$"), "", fp)
    tbl <- .read_report_table(fp)
    # Special handling for eigenvec: add FID/IID + PC names
    if (!is.null(tbl) && key == "eigenvec") {
      k <- ncol(tbl)
      if (k >= 2) {
        pcs <- if (k > 2) paste0("PC", seq_len(k - 2)) else character()
        colnames(tbl) <- c("FID", "IID", pcs)
      }
    }
    if (!is.null(tbl)) res[[key]] <- tbl
  }
  res
}
