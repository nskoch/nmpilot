# =============================================================================
# Parallel NONMEM runs: background launcher and run queue
# =============================================================================
#
# NM-TRAN and the compiler write their work files in the current directory,
# so two simultaneous runs must use two different directories. Each run is
# started in the background by a small launcher script written in its
# directory (`<run>_nmpilot.bat` on Windows, `<run>_nmpilot.sh` elsewhere):
# it calls nmfe exactly as the sequential launcher does, sends the console
# output to `<run>_nmpilot.log` and writes the exit status of nmfe to the
# flag file `<run>.nmpilot_done`. R keeps at most `n_cores` runs going and
# polls the flag files. With n_cores = 1 the functions keep their sequential
# launcher (console output of NONMEM shown as before).
# =============================================================================

#' Number of processor cores available for parallel NONMEM runs
#'
#' Reports the logical and physical cores of the computer, and the number of
#' parallel runs that `n_cores = "auto"` would use (physical cores minus one,
#' at least 1). Several runs share the memory and the disk: on a laptop, or
#' with large data sets, fewer runs than cores can be faster.
#'
#' @param verbose Print the result.
#' @return Invisibly, a list with `logical`, `physical` and `auto`.
#' @seealso [nm_run()], [nm_bootstrap()], [nm_jackknife()], [nm_profile()]
#' @examples
#' nm_cores()
#' @export
nm_cores <- function(verbose = TRUE) {
  lg <- suppressWarnings(parallel::detectCores(logical = TRUE))
  ph <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(lg)) lg <- 1L
  if (is.na(ph)) ph <- lg
  auto <- max(1L, as.integer(ph) - 1L)
  if (verbose)
    message(sprintf("Cores: %d logical, %d physical. n_cores = \"auto\" would use %d parallel run%s.",
                    lg, ph, auto, if (auto > 1) "s" else ""))
  invisible(list(logical = as.integer(lg), physical = as.integer(ph), auto = auto))
}

# n_cores argument: a positive integer, "auto" or "ask" (interactive only).
.nm_resolve_cores <- function(n_cores) {
  if (is.null(n_cores)) return(1L)
  info <- nm_cores(verbose = FALSE)
  if (is.character(n_cores)) {
    n_cores <- match.arg(n_cores, c("auto", "ask"))
    if (n_cores == "auto") return(info$auto)
    # "ask": the answer is kept for the rest of the R session
    prev <- getOption("nmpilot.n_cores")
    if (!is.null(prev)) {
      message(sprintf("n_cores = %d (answer kept for this session; options(nmpilot.n_cores = NULL) to ask again)",
                      as.integer(prev)))
      return(as.integer(prev))
    }
    if (!interactive())
      stop("n_cores = \"ask\" needs an interactive session; give a number or \"auto\".",
           call. = FALSE)
    repeat {
      ans <- readline(sprintf(
        "This computer has %d logical / %d physical cores. Number of parallel NONMEM runs [%d]: ",
        info$logical, info$physical, info$auto))
      k <- if (!nzchar(trimws(ans))) info$auto else suppressWarnings(as.integer(ans))
      if (!is.na(k) && k >= 1L) { n_cores <- k; options(nmpilot.n_cores = k); break }
      message("Please type a whole number of at least 1.")
    }
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1L || is.na(n_cores) || n_cores < 1)
    stop("n_cores must be a whole number >= 1, \"auto\" or \"ask\".", call. = FALSE)
  n_cores <- as.integer(n_cores)
  if (n_cores > info$logical)
    warning(sprintf("n_cores = %d exceeds the %d logical cores of this computer.",
                    n_cores, info$logical), call. = FALSE)
  n_cores
}

# Start one NONMEM run in the background in `run_dir`; returns the flag path.
.nm_launch_bg <- function(nmfe_path, run_dir, run_name) {
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  done <- file.path(run_dir, paste0(run_name, ".nmpilot_done"))
  unlink(done)
  ctl <- paste0(run_name, ".ctl"); lst <- paste0(run_name, ".lst")
  log <- paste0(run_name, "_nmpilot.log")
  if (.Platform$OS.type == "windows") {
    nmfe_win <- normalizePath(nmfe_path, winslash = "\\", mustWork = FALSE)
    script <- file.path(run_dir, paste0(run_name, "_nmpilot.bat"))
    # No redirection inside the batch: "call nmfe ... 2>&1" makes cmd keep the
    # first "2>trash.tmp" of nmfe in force (later copies fail), and running nmfe
    # in a child cmd lost the last line. R itself opens the log files, so cmd
    # inherits plain handles, as when nmfe is run from the console.
    writeLines(c("@echo off",
                 "cd /d \"%~dp0\"",
                 sprintf("call \"%s\" %s %s", nmfe_win, ctl, lst),
                 # redirection first: "echo 1> f" would redirect handle 1
                 sprintf("> %s echo %%ERRORLEVEL%%", basename(done))),
               script)
    old <- setwd(run_dir); on.exit(setwd(old), add = TRUE)
    system2("cmd.exe", args = c("/c", basename(script)), wait = FALSE,
            stdout = log, stderr = paste0(run_name, "_nmpilot.err"))
  } else {
    script <- file.path(run_dir, paste0(run_name, "_nmpilot.sh"))
    writeLines(c("#!/bin/sh",
                 "cd \"$(dirname \"$0\")\"",
                 sprintf("\"%s\" %s %s > %s 2>&1", nmfe_path, ctl, lst, log),
                 sprintf("echo $? > %s", basename(done))),
               script)
    old <- setwd(run_dir); on.exit(setwd(old), add = TRUE)
    system2("sh", args = basename(script), wait = FALSE, stdout = FALSE, stderr = FALSE)
  }
  done
}

# Read the exit status written by the launcher (NA while not complete).
.nm_read_done <- function(done) {
  if (!file.exists(done)) return(NA_integer_)
  x <- tryCatch(readLines(done, warn = FALSE), error = function(e) character(0))
  x <- trimws(x[nzchar(trimws(x))])
  if (!length(x)) return(NA_integer_)
  v <- suppressWarnings(as.integer(x[1]))
  if (is.na(v)) NA_integer_ else v
}

# Run a list of jobs with at most n_cores at a time.
#   jobs     list; each element is passed to prepare()
#   prepare  function(job, slot) -> list(dir, name) of the run to launch, or
#            NULL to skip the job (e.g. already done)
#   finish   function(job, slot, status, seconds), called in R when a run ends
# Returns the exit statuses (NA for skipped jobs). If R is interrupted, the
# runs already started go on in the background.
.nm_run_queue <- function(jobs, nmfe_path, n_cores, prepare, finish = NULL,
                          poll = 1, verbose = TRUE, label = "runs") {
  n <- length(jobs)
  status <- rep(NA_integer_, n)
  if (!n) return(status)
  free <- seq_len(n_cores)
  running <- list()           # name: job index -> list(slot, done, t0, name, dir)
  nxt <- 1L; n_done <- 0L; n_skip <- 0L
  on.exit({
    if (length(running))
      message(sprintf("Interrupted: %d run(s) still going on in the background (%s); ",
                      length(running),
                      paste(vapply(running, `[[`, "", "name"), collapse = ", ")),
              "a later call with resume/skip of terminated runs picks up their results.")
  }, add = TRUE)
  say <- function(...) if (verbose) message(...)
  while (nxt <= n || length(running)) {
    while (nxt <= n && length(free)) {
      slot <- free[1]
      spec <- prepare(jobs[[nxt]], slot)
      if (is.null(spec)) { n_skip <- n_skip + 1L; nxt <- nxt + 1L; next }
      done <- .nm_launch_bg(nmfe_path, spec$dir, spec$name)
      running[[as.character(nxt)]] <- list(slot = slot, done = done, t0 = Sys.time(),
                                           name = spec$name, dir = spec$dir)
      free <- free[-1]
      say(sprintf("  started %s (slot %d)", spec$name, slot))
      nxt <- nxt + 1L
    }
    if (!length(running)) next
    Sys.sleep(poll)
    for (key in names(running)) {
      r <- running[[key]]
      st <- .nm_read_done(r$done)
      if (is.na(st)) next
      i <- as.integer(key)
      secs <- as.numeric(difftime(Sys.time(), r$t0, units = "secs"))
      status[i] <- st
      running[[key]] <- NULL
      free <- sort(c(free, r$slot))
      n_done <- n_done + 1L
      if (!is.null(finish)) finish(jobs[[i]], r$slot, st, secs)
      say(sprintf("  finished %s (exit %d, %.0fs) - %d/%d %s done, %d running",
                  r$name, st, secs, n_done + n_skip, n, label, length(running)))
    }
  }
  status
}
