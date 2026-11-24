# =============================================================================
# nm_run.R
# Launch one or more NONMEM runs, each in its own run directory (<run>.reg by
# default) or directly in the folder of its control stream (run_subdir = "").
#
# Assumptions:
#   - The user has written the control streams, including $DATA.
#   - The data file is NOT copied into run directories.
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
  # several $EST / $PROB: OFV of the last estimation step
  x <- .nm_lst_step_lines(x)

  # Priority 1: "WITHOUT CONSTANT" (reported value in modern NONMEM)
  idx_woc <- grep("OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT", x, fixed = TRUE)
  if (length(idx_woc) > 0) {
    line <- x[idx_woc[length(idx_woc)]]
    num  <- regmatches(line, regexpr("-?[0-9]+\\.?[0-9]*([Ee][\\+\\-]?[0-9]+)?", line))
    if (length(num) > 0) return(round(as.numeric(num), 3))
  }

  # Priority 2: generic "OBJECTIVE FUNCTION VALUE" (older NONMEM)
  idx <- grep("OBJECTIVE FUNCTION VALUE", x, fixed = TRUE)
  if (length(idx) == 0) {
    # Priority 3: "#OBJV:" line (the only OFV line of an evaluation run,
    # MAXEVAL=0, which prints no "WITHOUT CONSTANT" line)
    idx_ov <- grep("^\\s*#OBJV:", x)
    if (length(idx_ov) == 0) return(NA_real_)
    line <- gsub("\\*", " ", sub("^\\s*#OBJV:", "", x[idx_ov[length(idx_ov)]]))
    num  <- regmatches(line, regexpr("-?[0-9]+\\.?[0-9]*([Ee][\\+\\-]?[0-9]+)?", line))
    return(if (length(num) == 0) NA_real_ else round(as.numeric(num), 3))
  }
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

  # Minimization status of the last estimation step (several $EST: the
  # earlier steps only prepare the last one; EM methods print
  # "OPTIMIZATION/PORTION ... COMPLETED" instead of "MINIMIZATION")
  stp <- .nm_lst_step_status(x)
  if (is.na(stp$code)) {
    # Simulation-only run ($SIMULATION ONLYSIM, e.g. the _sim run of
    # nm_vpc_ctl): no estimation step (no #METH line), no OFV; finished when
    # NONMEM printed its final lines (#CPUT / Stop Time)
    if (any(grepl("SIMULATION STEP PERFORMED", x, fixed = TRUE)) &&
        !any(grepl("^\\s*#METH:", x)) &&
        any(grepl("^\\s*(#CPUT:|Stop Time:)", x))) {
      return(data.frame(
        minimization_status = "simulation",
        min_ok              = NA,
        cov_ok              = NA,
        run_ok              = TRUE,
        cov_status          = "not_requested",
        stringsAsFactors    = FALSE
      ))
    }
    # Evaluation run (MAXEVAL=0): no minimization, an OFV is printed
    if (any(grepl("(EVALUATION)", x, fixed = TRUE)) && any(grepl("^\\s*#OBJV:", x))) {
      return(data.frame(
        minimization_status = "evaluation",
        min_ok              = NA,
        cov_ok              = NA,
        run_ok              = TRUE,
        cov_status          = "not_requested",
        stringsAsFactors    = FALSE
      ))
    }
    return(empty)
  }
  min_status <- stp$code
  min_ok <- switch(stp$code, successful = TRUE, not_tested = NA, FALSE)

  # Covariance step
  # only the status lines ("COVARIANCE STEP OMITTED: NO", "... ABORTED"...):
  # the sentence "AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE
  # OUTPUT" printed after "HOWEVER, PROBLEMS OCCURRED" is not one
  cov_idx <- grep("^\\s*0?\\s*COVARIANCE STEP\\s+(OMITTED|ABORTED|SUCCESSFUL|WAS)", x)
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

# Internal helper: warn when the analysis directory is inside a folder kept in
# sync by a cloud client (Dropbox, OneDrive, Google Drive, iCloud, Box). The
# client can lock the files NM-TRAN and the compiler are writing ("the process
# cannot access the file because it is being used by another process") and
# make a run fail. Silenced by options(nmpilot.warn_synced = FALSE).
.nm_warn_synced_dir <- function(dir) {
  if (isFALSE(getOption("nmpilot.warn_synced", TRUE))) return(invisible(FALSE))
  pat <- "(?i)(dropbox|onedrive|google ?drive|icloud ?drive|box ?sync|/box/|\\\\box\\\\|pcloud)"
  hit <- grepl(pat, dir, perl = TRUE)
  if (hit) warning(
    "The analysis directory is inside a synchronised folder (", dir, "). ",
    "The sync client can lock the files written by NM-TRAN and the compiler ",
    "and make a run fail: prefer a local folder (e.g. under C:/nm76g64/run/), ",
    "or pause the synchronisation. Silence with options(nmpilot.warn_synced = FALSE).",
    call. = FALSE)
  invisible(hit)
}

# Internal helper: check if a run is in a terminal state
.nm_is_run_terminal <- function(lst_file, accept_failed = TRUE) {
  if (!file.exists(lst_file)) return(FALSE)
  st <- .nm_lst_status(lst_file)
  if (st$minimization_status %in% c("evaluation", "simulation")) return(TRUE)
  if (is.na(st$min_ok)) return(FALSE)
  if (accept_failed) return(!is.na(st$minimization_status))
  isTRUE(st$run_ok)
}

# Internal helper: detect if a .ctl is in Wings extended format
.nm_is_wings_ctl <- function(ctl_file) {
  # Wings extended syntax: the model records use labels of the $THETA/$OMEGA/
  # $SIGMA comments (POPCL, eta_cl) that are never assigned in the code and
  # are not $INPUT items. A standard control stream whose comments repeat
  # parameter names (";CL" next to CL = THETA(1)*...) is not Wings.
  lines <- readLines(ctl_file, warn = FALSE)
  nrow(.nm_ip_expand_labels(lines)$substitutions) > 0L
}

# Internal helper: fix $DATA path to absolute in a .ctl file
.nm_fix_data_path <- function(ctl_file, src_dir, reg_dir, verbose = TRUE) {
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

    # A relative path is looked up first from the folder of the original
    # control stream (plain nmfe usage: $DATA data.csv), then from the run
    # directory (Wings usage: $DATA ../data.csv). If neither exists, the
    # control-stream folder is used and NM-TRAN will report the missing file.
    cand <- c(file.path(src_dir, data_file), file.path(reg_dir, data_file))
    hit <- cand[file.exists(cand)]
    abs_path <- normalizePath(if (length(hit)) hit[1] else cand[1], mustWork = FALSE)
    # Use forward slashes (NONMEM accepts both on Windows)
    abs_path   <- gsub("\\", "/", abs_path, fixed = TRUE)
    lines[idx] <- sub(data_file, abs_path, ln, fixed = TRUE)

    if (verbose) message("  $DATA path resolved: ", data_file, " -> ", abs_path)
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
#' executable (direct call, no Wings for NONMEM, no PsN required). By
#' default each run gets its own subdirectory `<run>.reg` inside
#' `analysis_dir`, holding a copy of the control stream and all NONMEM output
#' files; with `run_subdir = ""` the run is made directly in the folder of
#' its control stream, as when nmfe is called by hand. The data file is not
#' copied. A relative `$DATA` path is accepted: it is looked up from the
#' folder of the control stream (`$DATA data.csv`), then from the run
#' directory (`$DATA ../data.csv`, Wings style), and written as an absolute
#' path in the copied control stream.
#'
#' @param ctl Character vector of paths to `.ctl` files. May be absolute
#'   paths or file names relative to `analysis_dir`. (Formerly `ctl_files`,
#'   still accepted with a warning; the former order, `nmfe_path` first, is
#'   also recognised.)
#' @param nmfe_path Full path to the nmfe executable
#'   (e.g. `"C:/nm76g64/run/nmfe76.bat"` on Windows, `"/opt/nm760/run/nmfe76"`
#'   on Linux or macOS, where it is run directly and must be executable).
#' @param analysis_dir Root directory of the analysis. One run subdirectory
#'   is created here for each run (see `run_subdir`). A warning is issued
#'   when it lies in a folder synchronised by a cloud client (Dropbox,
#'   OneDrive, Google Drive, iCloud, Box): the client can lock the files
#'   written by NM-TRAN and make a run fail (silence it with
#'   `options(nmpilot.warn_synced = FALSE)`).
#' @param run_subdir Suffix of the run subdirectories: `".reg"` (default,
#'   Wings-compatible layout `run001.reg/`), any other suffix (e.g. `".run"`),
#'   or `""` to run in the folder of the control stream itself, without
#'   copying it (outputs next to the `.ctl`). Control streams in Wings label
#'   syntax need a subdirectory, because the converted copy is written there.
#' @param skip_terminated Logical. If `TRUE` (default), runs whose `.lst` file
#'   already exists and reports a terminal state are skipped.
#' @param accept_failed Logical. Only used when `skip_terminated = TRUE`.
#'   If `TRUE` (default), runs that terminated with an error are also skipped
#'   (not re-launched). If `FALSE`, only runs with `run_ok = TRUE` are skipped.
#' @param overwrite_ctl Logical. If `TRUE` (default), the `.ctl` file is
#'   overwritten in the `.reg` directory if already present.
#' @param n_cores Number of runs launched at the same time (default 1: one
#'   after the other, NONMEM output shown in the console). A whole number,
#'   `"auto"` (physical cores minus one, see [nm_cores()]) or `"ask"`
#'   (interactive sessions: the cores are shown and the number is asked once per R
#'   session, then reused; `options(nmpilot.n_cores = NULL)` asks again). In
#'   parallel, each run is started in the background from its run directory;
#'   the console output of nmfe goes to `<run>_nmpilot.log`. Runs written in
#'   the same folder (`run_subdir = ""`) are always run one after the other.
#'   If R is interrupted, the runs already started go on in the background.
#' @param verbose Logical. If `TRUE` (default), progress messages are printed
#'   to the console.
#' @param ... Former argument names, accepted with a warning until version
#'   1.0: `ctl_files` (now `ctl`).
#'
#' @return A `data.frame` (run log) with one row per run and the following
#'   columns:
#'   \describe{
#'     \item{run}{Run name (without extension).}
#'     \item{reg_dir}{Full path to the run directory (`<run>.reg`, or the folder of the control stream when `run_subdir = ""`).}
#'     \item{action}{One of `"launched"`, `"skipped_done"`, or
#'       `"skipped_no_ctl"`.}
#'     \item{exit_status}{Shell return code from `nmfe` (0 = success).}
#'     \item{lst_exists}{Logical: does the `.lst` file exist after the run?}
#'     \item{ofv}{Objective function value (OFV).}
#'     \item{minimization_status}{Minimization status string: `"successful"`,
#'       `"terminated"`, `"rounding_errors"`, `"min_success_cov_failed"`,
#'       `"evaluation"` (run without minimization, `MAXEVAL=0`; `min_ok` is
#'       then `NA` and `run_ok` `TRUE`), `"simulation"` (simulation-only
#'       run, `$SIMULATION ONLYSIM`, finished; `ofv` and `min_ok` are then
#'       `NA` and `run_ok` `TRUE`), or `"unknown"`.}
#'     \item{min_ok}{Logical: did minimization succeed?}
#'     \item{cov_ok}{Logical: did the covariance step succeed?}
#'     \item{run_ok}{Logical: is the run globally acceptable
#'       (minimization and covariance both successful)?}
#'     \item{cov_status}{Detailed covariance step status: `"successful"`,
#'       `"aborted"`, `"not_requested"`, or `"unknown"`.}
#'   }
#'
#' @details
#' Wings for NONMEM is not needed. The default `.reg` layout mirrors the
#' Wings output structure (a run named `run001.ctl` produces a directory
#' `run001.reg/` containing `run001.ctl`, `run001.lst`, `run001.ext`, etc.),
#' so that projects started with Wings can be continued; the other functions
#' of the package read runs in any of these layouts.
#'
#' Run status is parsed from the `.lst` file after completion.
#' The `skip_terminated` mechanism reads any existing `.lst` before launching,
#' allowing safe re-use of `nm_run()` on a partially completed batch.
#'
#' @examples
#' \dontrun{
#' log <- nm_run(
#'   ctl             = c("run001.ctl", "run002.ctl", "run003.ctl"),
#'   nmfe_path       = "C:/nm76g64/run/nmfe76.bat",
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
    ctl,
    nmfe_path,
    analysis_dir,
    skip_terminated = TRUE,
    accept_failed   = TRUE,
    overwrite_ctl   = TRUE,
    run_subdir      = ".reg",
    verbose         = TRUE,
    n_cores         = 1L,
    ...
) {
  dep <- .nm_deprecated_args(list(...), c(ctl_files = "ctl"), "nm_run")
  if (!is.null(dep$ctl)) {
    # nm_run(nmfe_path = ..., ctl_files = ..., ...): the nmfe path is then
    # in `ctl` when given by position
    if (missing(nmfe_path) && !missing(ctl)) nmfe_path <- ctl
    ctl <- dep$ctl
  }
  if (!missing(ctl) && !missing(nmfe_path) && .nm_looks_like_nmfe(ctl) &&
      !.nm_looks_like_nmfe(nmfe_path)) {
    warning("nm_run(): the control streams now come first, nm_run(ctl, nmfe_path, ...); ",
            "the former order was recognised.", call. = FALSE)
    tmp <- ctl; ctl <- nmfe_path; nmfe_path <- tmp
  }
  ctl_files <- ctl

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
  .nm_warn_synced_dir(analysis_dir)
  if (is.null(run_subdir)) run_subdir <- ""
  flat <- !nzchar(run_subdir)
  n_cores <- .nm_resolve_cores(n_cores)
  if (n_cores > 1L && flat && length(ctl_files) > 1L) {
    warning("run_subdir = \"\": the runs share a folder and cannot run at the same time; ",
            "they are run one after the other.", call. = FALSE)
    n_cores <- 1L
  }
  jobs <- list()   # runs prepared for the parallel queue

  # --- Loop over runs ---------------------------------------------------------

  results <- lapply(ctl_files, function(ctl) {

    # Resolve .ctl path
    if (!file.exists(ctl)) {
      ctl_full <- file.path(analysis_dir, ctl)
    } else {
      ctl_full <- normalizePath(ctl, mustWork = FALSE)
    }

    run_name <- tools::file_path_sans_ext(basename(ctl_full))
    reg_dir  <- if (flat) dirname(ctl_full) else file.path(analysis_dir, paste0(run_name, run_subdir))
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
    wings <- .nm_is_wings_ctl(ctl_full)
    if (wings && flat)
      stop("Control stream ", basename(ctl_full), " uses Wings label syntax: ",
           "it must be converted, which needs a run subdirectory. ",
           "Use run_subdir = \".reg\", or convert it first with nm_convert_ctl().",
           call. = FALSE)
    if (wings) {
      if (verbose) message("    [Wings CTL detected] Converting to standard NONMEM...")
      # Convert to a temp file, keeping original run name
      tmp_ctl <- tempfile(fileext = ".ctl")
      nm_convert_ctl(ctl_full, out_file = tmp_ctl,
                     overwrite = TRUE, verbose = FALSE)
      ctl_to_use <- tmp_ctl
      if (verbose) message("    Wings CTL converted successfully.")
    }

    if (flat) {
      # Run in place: the user's control stream is used as it is.
      ctl_dest <- ctl_full
    } else {
      # Copy .ctl into the run directory (always using the original run name)
      ctl_dest <- file.path(reg_dir, paste0(run_name, ".ctl"))
      if (overwrite_ctl || !file.exists(ctl_dest)) {
        file.copy(ctl_to_use, ctl_dest, overwrite = TRUE)
      }
      # Write the $DATA path as an absolute path in the copy
      .nm_fix_data_path(ctl_dest, dirname(ctl_full), reg_dir, verbose = verbose)
    }

    # Build and launch nmfe command
    lst_dest <- file.path(reg_dir, paste0(run_name, ".lst"))
    cmd <- paste0(shQuote(nmfe_path), " ", basename(ctl_dest), " ",
                  basename(lst_dest), "   (in ", reg_dir, ")")

    if (verbose) message("    Command: ", cmd)
    if (n_cores > 1L) {
      # launched later by the parallel queue
      jobs[[length(jobs) + 1L]] <<- list(run = run_name, dir = reg_dir)
      return(.nm_run_log_row(run_name, reg_dir, "launched", NA, lst_file))
    }
    # nmfe is run from inside the .reg directory: NM-TRAN writes its work
    # files in the current directory.
    exit_status <- .nm_launch(nmfe_path, reg_dir, run_name)
    if (verbose) message("    Exit status: ", exit_status)

    .nm_run_log_row(run_name, reg_dir, "launched", exit_status, lst_file)
  })

  # Parallel launch of the prepared runs, then log rows rebuilt from the .lst
  if (n_cores > 1L && length(jobs)) {
    if (verbose) message(sprintf(">>> Launching %d run(s), %d at a time", length(jobs), n_cores))
    st <- .nm_run_queue(jobs, nmfe_path, n_cores,
                        prepare = function(job, slot) list(dir = job$dir, name = job$run),
                        verbose = verbose)
    for (k in seq_along(jobs)) {
      i <- which(vapply(results, function(r) r$run == jobs[[k]]$run &&
                          r$action == "launched", logical(1)))[1]
      results[[i]] <- .nm_run_log_row(jobs[[k]]$run, jobs[[k]]$dir, "launched", st[k],
                                      file.path(jobs[[k]]$dir, paste0(jobs[[k]]$run, ".lst")))
    }
  }

  # Assemble log
  do.call(rbind, results)
}
