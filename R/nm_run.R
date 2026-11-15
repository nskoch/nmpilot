# =============================================================================
# nm_run.R
# Launch one or more NONMEM runs with .reg directory management
#
# Assumptions:
#   - The user has written the control streams, including $DATA.
#   - The data file is NOT copied into .reg directories.
#   - nmfe76.bat (or equivalent) is accessible on the machine.
# =============================================================================


# =============================================================================
# Internal parsing functions (.lst)
# Prefixed with a dot to indicate they are internal (not exported).
# =============================================================================

.nm_read_lst_lines <- function(lst_file) {
  if (!file.exists(lst_file)) return(character(0))
  readLines(lst_file, warn = FALSE)
}

.nm_lst_ofv <- function(lst_file) {
  x <- .nm_read_lst_lines(lst_file)
  if (length(x) == 0) return(NA_real_)

  # Priority 1: "WITHOUT CONSTANT" (reported value in modern NONMEM)
  idx_woc <- grep("OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT", x, fixed = TRUE)
  if (length(idx_woc) > 0) {
    line <- x[idx_woc[length(idx_woc)]]
    num  <- regmatches(line, regexpr("-?[0-9]+\\.?[0-9]*([Ee][\\+\\-]?[0-9]+)?", line))
    if (length(num) > 0) return(round(as.numeric(num), 3))
  }

  # Priority 2: generic "OBJECTIVE FUNCTION VALUE" (older NONMEM)
  idx <- grep("OBJECTIVE FUNCTION VALUE", x, fixed = TRUE)
  if (length(idx) == 0) return(NA_real_)
  line <- x[idx[length(idx)]]
  num  <- regmatches(line, regexpr("-?[0-9]+\\.?[0-9]*([Ee][\\+\\-]?[0-9]+)?", line))
  if (length(num) == 0) NA_real_ else round(as.numeric(num), 3)
}

.nm_lst_status <- function(lst_file) {
  empty <- data.frame(
    minimization_status = NA_character_,
    min_ok              = NA,
    cov_ok              = NA,
    run_ok              = NA,
    cov_status          = NA_character_,
    stringsAsFactors    = FALSE
  )

  x <- .nm_read_lst_lines(lst_file)
  if (length(x) == 0) return(empty)

  # Detect NMTRAN/control stream errors (no minimization attempted)
  if (any(grepl("AN ERROR WAS FOUND IN THE CONTROL STATEMENTS", x, fixed = TRUE)) ||
      any(grepl("AN ERROR WAS FOUND ON LINE", x, fixed = TRUE))) {
    return(data.frame(
      minimization_status = "nmtran_error",
      min_ok              = FALSE,
      cov_ok              = NA,
      run_ok              = FALSE,
      cov_status          = NA_character_,
      stringsAsFactors    = FALSE
    ))
  }

  # Minimization status
  min_idx <- grep("MINIMIZATION", x, fixed = TRUE)
  if (length(min_idx) == 0) return(empty)
  min_line <- x[min_idx[1]]

  if (grepl("SUCCESSFUL", min_line, fixed = TRUE)) {
    min_status <- "successful"
    min_ok     <- TRUE
  } else if (grepl("TERMINATED", min_line, fixed = TRUE)) {
    min_status <- "terminated"
    min_ok     <- FALSE
  } else if (grepl("ROUNDING ERRORS", min_line, fixed = TRUE)) {
    min_status <- "rounding_errors"
    min_ok     <- FALSE
  } else {
    min_status <- "unknown"
    min_ok     <- FALSE
  }

  # Covariance step
  cov_idx <- grep("COVARIANCE STEP", x, fixed = TRUE)
  if (length(cov_idx) == 0) {
    cov_ok     <- NA
    cov_status <- "not_requested"
  } else {
    cov_line <- x[cov_idx[length(cov_idx)]]
    if (grepl("ABORTED", cov_line, fixed = TRUE)) {
      # Modern NONMEM: covariance step aborted
      cov_ok     <- FALSE
      cov_status <- "aborted"
    } else if (grepl("SUCCESSFUL", cov_line, fixed = TRUE)) {
      # Modern NONMEM: covariance step successful
      cov_ok     <- TRUE
      cov_status <- "successful"
    } else if (grepl("OMITTED:.*NO", cov_line)) {
      # Older NONMEM: "COVARIANCE STEP OMITTED: NO" = step was run
      # Check for specific covariance eigenvalue warning (not WRESCHOL)
      eigen_warn <- any(grepl("EIGENVALUE OF COV", x, fixed = TRUE))
      cov_ok     <- !eigen_warn
      cov_status <- if (eigen_warn) "eigenvalue_issue" else "successful"
    } else if (grepl("OMITTED:.*YES", cov_line)) {
      # Older NONMEM: "COVARIANCE STEP OMITTED: YES" = step was skipped
      cov_ok     <- NA
      cov_status <- "not_requested"
    } else {
      cov_ok     <- FALSE
      cov_status <- "unknown"
    }
  }

  if (isTRUE(min_ok) && identical(cov_ok, FALSE)) {
    min_status <- "min_success_cov_failed"
  }

  run_ok <- isTRUE(min_ok) && !identical(cov_ok, FALSE)

  data.frame(
    minimization_status = min_status,
    min_ok              = min_ok,
    cov_ok              = cov_ok,
    run_ok              = run_ok,
    cov_status          = cov_status,
    stringsAsFactors    = FALSE
  )
}

# Internal helper: check if a run is in a terminal state
.nm_is_run_terminal <- function(lst_file, accept_failed = TRUE) {
  if (!file.exists(lst_file)) return(FALSE)
  st <- .nm_lst_status(lst_file)
  if (is.na(st$min_ok)) return(FALSE)
  if (accept_failed) return(!is.na(st$minimization_status))
  isTRUE(st$run_ok)
}

# Internal helper: detect if a .ctl is in Wings extended format
.nm_is_wings_ctl <- function(ctl_file) {
  lines <- readLines(ctl_file, warn = FALSE)
  lines <- gsub("\r", "", lines, fixed = TRUE)

  # Build mapping from $THETA/$OMEGA/$SIGMA comments
  theta_map <- .nm_build_mapping(lines, "^\\s*\\$THETA", "THETA")
  eta_map   <- .nm_build_eta_mapping(lines)
  err_map   <- .nm_build_err_mapping(lines)
  all_names <- c(names(theta_map), names(eta_map), names(err_map))

  if (length(all_names) == 0L) return(FALSE)

  # Check if any name appears in model records ($PK, $PRED, $ERROR, $DES)
  model_pat <- "^\\s*\\$(PK|PRED|ERROR|DES|AES|MIX)\\b"
  in_model  <- FALSE
  for (ln in lines) {
    if (grepl(model_pat, ln, ignore.case = TRUE)) { in_model <- TRUE; next }
    if (in_model && grepl("^\\s*\\$[A-Za-z]", ln)) in_model <- FALSE
    if (!in_model) next
    code_part <- sub(";.*$", "", ln)
    for (nm in all_names) {
      if (grepl(paste0("\\b", nm, "\\b"), code_part)) return(TRUE)
    }
  }
  FALSE
}

# Internal helper: fix $DATA path to absolute in a .ctl file
.nm_fix_data_path <- function(ctl_file, analysis_dir, reg_dir) {
  lines <- readLines(ctl_file, warn = FALSE)
  lines <- gsub("\r", "", lines, fixed = TRUE)

  data_idx <- grep("^\\s*\\$DATA", lines, ignore.case = TRUE)
  if (length(data_idx) == 0L) return(invisible(NULL))

  for (idx in data_idx) {
    ln <- lines[idx]

    # Extract the data filename: second token on the $DATA line
    tokens    <- strsplit(trimws(ln), "\\s+")[[1]]
    if (length(tokens) < 2L) next
    data_file <- tokens[2]  # first token is $DATA, second is filename

    # Skip if already absolute (drive letter or UNC or Unix root)
    if (grepl("^[A-Za-z]:", data_file) ||
        grepl("^/",         data_file) ||
        grepl("^\\\\",  data_file)) next

    # Resolve relative path from reg_dir (Wings uses "..\" relative to .reg)
    abs_path <- normalizePath(
      file.path(reg_dir, data_file),
      mustWork = FALSE
    )
    # Use forward slashes (NONMEM accepts both on Windows)
    abs_path   <- gsub("\\", "/", abs_path, fixed = TRUE)
    lines[idx] <- sub(data_file, abs_path, ln, fixed = TRUE)

    message("  $DATA path resolved: ", data_file, " -> ", abs_path)
  }

  writeLines(lines, ctl_file)
  invisible(ctl_file)
}

# Internal helper: build one log row
.nm_run_log_row <- function(run_name, reg_dir, action, exit_status, lst_file) {
  ofv        <- .nm_lst_ofv(lst_file)
  status_tbl <- .nm_lst_status(lst_file)

  data.frame(
    run                 = run_name,
    reg_dir             = reg_dir,
    action              = action,
    exit_status         = exit_status,
    lst_exists          = file.exists(lst_file),
    ofv                 = ofv,
    minimization_status = status_tbl$minimization_status,
    min_ok              = status_tbl$min_ok,
    cov_ok              = status_tbl$cov_ok,
    run_ok              = status_tbl$run_ok,
    cov_status          = status_tbl$cov_status,
    stringsAsFactors    = FALSE
  )
}


# =============================================================================
# Main exported function
# =============================================================================

#' Run one or more NONMEM models
#'
#' Sequentially launches a list of NONMEM control streams via the `nmfe`
#' executable (direct call, no Wings for NONMEM, no PsN required). For each
#' run, a subdirectory `<run>.reg` is created inside `analysis_dir` to store
#' the control stream and all NONMEM output files. The data file is not copied:
#' the `$DATA` record in each control stream must point to the data file using
#' an absolute path.
#'
#' @param nmfe_path Full path to the nmfe executable
#'   (e.g. `"C:/nm76g64/run/nmfe76.bat"`).
#' @param ctl_files Character vector of paths to `.ctl` files. May be absolute
#'   paths or file names relative to `analysis_dir`.
#' @param analysis_dir Root directory of the analysis. One `.reg` subdirectory
#'   will be created here for each run.
#' @param skip_terminated Logical. If `TRUE` (default), runs whose `.lst` file
#'   already exists and reports a terminal state are skipped.
#' @param accept_failed Logical. Only used when `skip_terminated = TRUE`.
#'   If `TRUE` (default), runs that terminated with an error are also skipped
#'   (not re-launched). If `FALSE`, only runs with `run_ok = TRUE` are skipped.
#' @param overwrite_ctl Logical. If `TRUE` (default), the `.ctl` file is
#'   overwritten in the `.reg` directory if already present.
#' @param verbose Logical. If `TRUE` (default), progress messages are printed
#'   to the console.
#'
#' @return A `data.frame` (run log) with one row per run and the following
#'   columns:
#'   \describe{
#'     \item{run}{Run name (without extension).}
#'     \item{reg_dir}{Full path to the `.reg` directory.}
#'     \item{action}{One of `"launched"`, `"skipped_done"`, or
#'       `"skipped_no_ctl"`.}
#'     \item{exit_status}{Shell return code from `nmfe` (0 = success).}
#'     \item{lst_exists}{Logical: does the `.lst` file exist after the run?}
#'     \item{ofv}{Objective function value (OFV).}
#'     \item{minimization_status}{Minimization status string: `"successful"`,
#'       `"terminated"`, `"rounding_errors"`, `"min_success_cov_failed"`, or
#'       `"unknown"`.}
#'     \item{min_ok}{Logical: did minimization succeed?}
#'     \item{cov_ok}{Logical: did the covariance step succeed?}
#'     \item{run_ok}{Logical: is the run globally acceptable
#'       (minimization and covariance both successful)?}
#'     \item{cov_status}{Detailed covariance step status: `"successful"`,
#'       `"aborted"`, `"not_requested"`, or `"unknown"`.}
#'   }
#'
#' @details
#' The `.reg` directory convention mirrors the Wings for NONMEM output
#' structure: a run named `run001.ctl` produces a directory `run001.reg/`
#' containing `run001.ctl`, `run001.lst`, `run001.ext`, etc.
#'
#' Run status is parsed from the `.lst` file after completion.
#' The `skip_terminated` mechanism reads any existing `.lst` before launching,
#' allowing safe re-use of `nm_run()` on a partially completed batch.
#'
#' @examples
#' \dontrun{
#' log <- nm_run(
#'   nmfe_path       = "C:/nm76g64/run/nmfe76.bat",
#'   ctl_files       = c("run001.ctl", "run002.ctl", "run003.ctl"),
#'   analysis_dir    = "C:/analyses/my_pk_project",
#'   skip_terminated = TRUE,
#'   accept_failed   = TRUE,
#'   verbose         = TRUE
#' )
#'
#' # Show only failed runs
#' log[!log$run_ok, ]
#' }
#'
#' @seealso [nm_extract_runs()] for parsing NONMEM output files after a run.
#'
#' @export
nm_run <- function(
    nmfe_path,
    ctl_files,
    analysis_dir,
    skip_terminated = TRUE,
    accept_failed   = TRUE,
    overwrite_ctl   = TRUE,
    verbose         = TRUE
) {

  # --- Initial checks ---------------------------------------------------------

  if (!file.exists(nmfe_path)) {
    stop("nmfe executable not found: ", nmfe_path)
  }
  if (!dir.exists(analysis_dir)) {
    stop("Analysis directory does not exist: ", analysis_dir)
  }
  if (length(ctl_files) == 0) {
    stop("No .ctl files provided.")
  }

  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)

  # --- Loop over runs ---------------------------------------------------------

  results <- lapply(ctl_files, function(ctl) {

    # Resolve .ctl path
    if (!file.exists(ctl)) {
      ctl_full <- file.path(analysis_dir, ctl)
    } else {
      ctl_full <- normalizePath(ctl, mustWork = FALSE)
    }

    run_name <- tools::file_path_sans_ext(basename(ctl_full))
    reg_dir  <- file.path(analysis_dir, paste0(run_name, ".reg"))
    lst_file <- file.path(reg_dir, paste0(run_name, ".lst"))
    # Note: if Wings CTL is converted, run_name keeps original name

    if (verbose) message(">>> Run: ", run_name)

    # Check .ctl source exists
    if (!file.exists(ctl_full)) {
      if (verbose) message("    [SKIP] .ctl not found: ", ctl_full)
      return(.nm_run_log_row(run_name, reg_dir, "skipped_no_ctl", NA, lst_file))
    }

    # Skip if already terminated
    if (skip_terminated && .nm_is_run_terminal(lst_file, accept_failed)) {
      if (verbose) message("    [SKIP] Run already terminated.")
      return(.nm_run_log_row(run_name, reg_dir, "skipped_done", NA, lst_file))
    }

    # Create .reg directory
    if (!dir.exists(reg_dir)) dir.create(reg_dir, recursive = TRUE)

    # Auto-convert Wings extended CTL if detected
    ctl_to_use <- ctl_full
    if (.nm_is_wings_ctl(ctl_full)) {
      if (verbose) message("    [Wings CTL detected] Converting to standard NONMEM...")
      # Convert to a temp file, keeping original run name
      tmp_ctl <- tempfile(fileext = ".ctl")
      nm_convert_ctl(ctl_full, out_file = tmp_ctl,
                     overwrite = TRUE, verbose = FALSE)
      ctl_to_use <- tmp_ctl
      if (verbose) message("    Wings CTL converted successfully.")
    }

    # Copy .ctl into .reg (always using original run_name for the dest file)
    ctl_dest <- file.path(reg_dir, paste0(run_name, ".ctl"))
    if (overwrite_ctl || !file.exists(ctl_dest)) {
      file.copy(ctl_to_use, ctl_dest, overwrite = TRUE)
    }

    # Fix $DATA path to absolute (handles relative paths like ..\data.csv)
    .nm_fix_data_path(ctl_dest, analysis_dir, reg_dir)

    # Build and launch nmfe command
    lst_dest <- file.path(reg_dir, paste0(run_name, ".lst"))
    cmd <- paste(
      shQuote(nmfe_path),
      shQuote(ctl_dest),
      shQuote(lst_dest)
    )

    if (verbose) message("    Command: ", cmd)
    exit_status <- system(cmd, wait = TRUE)
    if (verbose) message("    Exit status: ", exit_status)

    .nm_run_log_row(run_name, reg_dir, "launched", exit_status, lst_file)
  })

  # Assemble log
  do.call(rbind, results)
}
