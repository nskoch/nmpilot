# =============================================================================
# utils-nm.R
# Package-level imports, global variable declarations, and internal helpers
# shared by several nmpilot functions.
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

# -----------------------------------------------------------------------------
# Deprecated argument names (action 11, interface harmonised before 1.0)
# -----------------------------------------------------------------------------

# Old argument names given by name land in `...`. `map` gives, for each old
# name, the new one. Returns a named list (new name -> value) and warns once
# per old name; any other name in `...` is an error, as for an unused
# argument.
.nm_deprecated_args <- function(dots, map, fun) {
  if (!length(dots)) return(list())
  nm <- names(dots)
  if (is.null(nm) || any(!nzchar(nm)))
    stop("Unnamed extra argument(s) in ", fun, "().", call. = FALSE)
  bad <- setdiff(nm, names(map))
  if (length(bad))
    stop("unused argument(s) in ", fun, "(): ", paste(bad, collapse = ", "), call. = FALSE)
  out <- list()
  for (o in nm) {
    n <- map[[o]]
    warning("In ", fun, "(), `", o, "` is deprecated; use `", n, "` instead.", call. = FALSE)
    out[[n]] <- dots[[o]]
  }
  out
}

# TRUE when `x` looks like the nmfe launcher rather than a control stream,
# used to recognise the former argument order of nm_run() and nm_profile()
# (nmfe_path first).
.nm_looks_like_nmfe <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) &&
    grepl("^nmfe", basename(x), ignore.case = TRUE) &&
    !grepl("\\.(ctl|mod)$", x, ignore.case = TRUE)
}

# A finished run given as one path (`run`): its control stream
# (".../run2.ctl"), its run directory (".../run2.reg", ".nm7", ".../run2/")
# or its stem (".../run2"). Returns the analysis directory and the run name,
# as expected by the functions that look for <run>.reg, <run>.nm7, <run>/ or
# the flat layout.
.nm_run_arg <- function(run, extensions = c("nm7", "reg")) {
  if (!is.character(run) || length(run) != 1L || is.na(run) || !nzchar(trimws(run)))
    stop("`run` must be the path of a run: its control stream, its run directory or its stem.",
         call. = FALSE)
  p <- sub("/+$", "", gsub("\\\\", "/", trimws(run)))
  ext_pat <- paste0("\\.(", paste(extensions, collapse = "|"), ")$")
  if (grepl("\\.(ctl|mod|ext|lst)$", p, ignore.case = TRUE)) {
    name <- tools::file_path_sans_ext(basename(p))
    dir <- dirname(p)
    # control stream copied in its run directory (run2.reg/run2.ctl): the
    # analysis directory is one level up. In a folder named like the run
    # (run2/run2.ctl), the files are found there as a flat layout.
    if (grepl(ext_pat, basename(dir), ignore.case = TRUE))
      dir <- dirname(dir)
    return(list(analysis_dir = dir, run_name = name))
  }
  if (grepl(ext_pat, p, ignore.case = TRUE))
    return(list(analysis_dir = dirname(p),
                run_name = sub(ext_pat, "", basename(p), ignore.case = TRUE)))
  list(analysis_dir = dirname(p), run_name = basename(p))
}

# Former name of the argument that chooses the figure in a plot() method
# (`which` or `area`, now `type`): taken from `...` with a warning.
.nm_plot_type_arg <- function(type, dots, old, fun) {
  if (is.null(dots[[old]])) return(type)
  warning("In ", fun, "(), `", old, "` is deprecated; use `type` instead.", call. = FALSE)
  dots[[old]]
}
