# =============================================================================
# utils-nm.R
# Package-level imports, global variable declarations, and internal helpers
# shared by several nmtools functions.
# =============================================================================

#' @importFrom magrittr %>%
#' @importFrom rlang .data %||%
NULL


# Suppress R CMD check NOTEs for unquoted column names used in dplyr pipes
utils::globalVariables(c(
  # NONMEM reserved columns
  "ID", "TIME", "DV", "MDV", "EVID", "AMT", "CMT", "BLQ",
  # Internal computed columns
  ".DT", ".t0", ".dose_ref", ".blq", ".y", ".y_plot", ".ABS",
  # dplyr summarise outputs
  "lo", "hi", "med", "n",
  # Other
  "DOSE", "OCC", "STUDY"
))


# =============================================================================
# Shared internal helpers
# =============================================================================

# NULL-coalescing `%||%` comes from rlang (returns `b` only when `a` is NULL).
# Variant that also replaces a single NA (used by nm_gof):
`%|na|%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# Read a text file and strip Windows carriage returns.
.nm_read_lines <- function(path) sub("\r$", "", readLines(path, warn = FALSE))

# Remove a $RECORD block (header line to the next $ record or EOF). [verified]
.nm_strip_record <- function(ctl_lines, record) {
  pat  <- paste0("^\\s*\\$", record)
  drop <- rep(FALSE, length(ctl_lines))
  i <- 1L
  while (i <= length(ctl_lines)) {
    if (grepl(pat, ctl_lines[i], ignore.case = TRUE)) {
      drop[i] <- TRUE; j <- i + 1L
      while (j <= length(ctl_lines) && !grepl("^\\s*\\$", ctl_lines[j])) { drop[j] <- TRUE; j <- j + 1L }
      i <- j
    } else i <- i + 1L
  }
  ctl_lines[!drop]
}

# Absolute $DATA path (handles quoted paths, relative->ctl dir). [verified]
.nm_extract_data_path <- function(ctl) {
  lines <- .nm_read_lines(ctl)
  i <- grep("^\\s*\\$DATA", lines, ignore.case = TRUE)[1]
  if (is.na(i)) stop("No $DATA record found in control file: ", ctl)
  l <- sub("^\\s*\\$DATA\\s+", "", lines[i], ignore.case = TRUE)
  l <- strsplit(l, ";")[[1]][1]; l <- trimws(l)
  if (grepl('^["\']', l)) {
    q <- substr(l, 1, 1); f <- sub(paste0("^", q, "([^", q, "]*)", q, ".*$"), "\\1", l)
  } else f <- strsplit(l, "\\s+")[[1]][1]
  if (!grepl("^([A-Za-z]:|/|\\\\)", f)) f <- file.path(dirname(ctl), f)
  normalizePath(f, winslash = "/", mustWork = FALSE)
}

# Replace ONLY the data-file token, keeping every IGNORE/ACCEPT option. [verified]
.nm_patch_data_line <- function(line, new_file) {
  m <- regmatches(line, regexec('^(\\s*\\$DATA\\s+)("[^"]*"|\'[^\']*\'|\\S+)(.*)$',
                                line, ignore.case = TRUE))[[1]]
  if (length(m) == 0) return(paste0("$DATA ", new_file))
  paste0(m[2], new_file, m[4])
}


# Locate the directory holding the output files of a run, whatever the
# layout: <run>.reg (Wings / nm_run default), <run>.nm7, a folder named
# <run>, or the analysis directory itself (flat layout, outputs next to the
# control stream). Returns the first candidate containing all files
# <run_name><need>, or stops with the list of paths tried.
.nm_find_run_dir <- function(analysis_dir, run_name, need = ".ext",
                             subdirs = c(".reg", ".nm7")) {
  cand <- c(file.path(analysis_dir, paste0(run_name, subdirs)),
            file.path(analysis_dir, run_name),
            analysis_dir)
  for (d in cand) {
    if (all(file.exists(file.path(d, paste0(run_name, need))))) return(d)
  }
  stop("Files ", paste0(run_name, need, collapse = ", "),
       " not found for run '", run_name, "'. Folders tried:\n - ",
       paste(cand, collapse = "\n - "), call. = FALSE)
}
