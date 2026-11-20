# =============================================================================
# nm_vpc_ctl.R
# Generate _orig.ctl and _sim.ctl from an estimation run, then launch both
# via nm_run(). Used as the first step before nm_vpc().
#
# Differences between _orig and _sim:
#   _orig: $ESTIMATION MAXEVAL=0  /  NSIM = 1 in $PK
#   _sim:  $SIMULATION (seed) ONLYSIM SUBPROBLEMS=n_sim  /  NSIM=IREP in $PK
# =============================================================================

#' @importFrom utils read.table
NULL


# =============================================================================
# Internal helpers
# =============================================================================

# Read .ext and return the final estimates row (-1000000000)
.nm_vpc_read_ext_final <- function(ext_file) {
  raw <- tryCatch(
    utils::read.table(ext_file, skip = 1, header = TRUE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(raw)) stop("Cannot read .ext file: ", ext_file)
  iter_col <- as.character(raw[[1]])
  final    <- raw[iter_col == "-1000000000", , drop = FALSE]
  if (nrow(final) == 0) stop("No final estimates row (-1000000000) in: ", ext_file)
  final
}

# Extract THETA, OMEGA, SIGMA values from the final row
.nm_vpc_get_estimates <- function(final_row) {
  theta_cols <- grep("^THETA", names(final_row), value = TRUE)
  omega_cols <- grep("^OMEGA", names(final_row), value = TRUE)
  sigma_cols <- grep("^SIGMA", names(final_row), value = TRUE)

  list(
    theta = stats::setNames(as.numeric(final_row[1, theta_cols]), theta_cols),
    omega = stats::setNames(as.numeric(final_row[1, omega_cols]), omega_cols),
    sigma = stats::setNames(as.numeric(final_row[1, sigma_cols]), sigma_cols)
  )
}

# Format THETA block lines
.nm_vpc_write_theta <- function(theta_vals) {
  c("$THETA",
    vapply(theta_vals, function(v) sprintf("  %g FIX", v), character(1)),
    "")
}

# Format OMEGA or SIGMA block lines (diagonal only from .ext vector)
# .ext stores upper triangle: OMEGA(1,1) OMEGA(2,1) OMEGA(2,2) ...
.nm_vpc_write_matrix <- function(record, vals) {
  # Find dimension from number of elements: n*(n+1)/2 = length
  n_elem <- length(vals)
  dim    <- round((-1 + sqrt(1 + 8 * n_elem)) / 2)

  lines <- paste0("$", record, " BLOCK(", dim, ")")
  k <- 1
  for (i in seq_len(dim)) {
    row_vals <- vapply(seq_len(i), function(j) {
      sprintf("%g", vals[k + j - 1])
    }, character(1))
    lines <- c(lines, paste(" ", paste(row_vals, collapse = " "), "FIX"))
    k <- k + i
  }
  c(lines, "")
}

# Remove estimation-related records from a control stream
# Removes: $THETA $OMEGA $SIGMA $ESTIMATION $COVARIANCE $TABLE (all abbrevs)
.nm_vpc_remove_records <- function(ctl) {
  pattern <- "^\\s*\\$(THE|OME|SIG|EST|COV|TAB)[A-Za-z]*"
  starts  <- grep(pattern, ctl, ignore.case = TRUE)
  if (length(starts) == 0) return(ctl)

  all_rec <- grep("^\\s*\\$[A-Za-z]", ctl, ignore.case = TRUE)
  remove  <- unlist(lapply(starts, function(s) {
    next_rec <- all_rec[all_rec > s][1]
    end      <- if (is.na(next_rec)) length(ctl) else next_rec - 1
    s:end
  }))
  ctl[-unique(remove)]
}

# Set NSIM line in $PK block
.nm_vpc_set_nsim <- function(ctl, mode = c("orig", "sim")) {
  mode        <- match.arg(mode)
  replacement <- if (mode == "orig") "  NSIM = 1" else "  NSIM = IREP"

  pk_start <- grep("^\\s*\\$PK\\b", ctl, ignore.case = TRUE)[1]
  if (is.na(pk_start)) return(ctl)

  all_rec <- grep("^\\s*\\$[A-Za-z]", ctl, ignore.case = TRUE)
  pk_end  <- all_rec[all_rec > pk_start][1] - 1
  if (is.na(pk_end)) pk_end <- length(ctl)

  pk_lines <- ctl[pk_start:pk_end]
  nsim_idx <- grep("^\\s*NSIM\\s*=", pk_lines, ignore.case = TRUE)

  if (length(nsim_idx) > 0) {
    ctl[pk_start + nsim_idx[1] - 1] <- replacement
  } else {
    warning("NSIM line not found in $PK - not modified.")
  }
  ctl
}

# Build _orig.ctl: MAXEVAL=0 estimation, NSIM=1
.nm_vpc_build_orig_ctl <- function(ctl_base, estimates, run_base) {
  ctl <- .nm_vpc_remove_records(ctl_base)
  ctl <- .nm_vpc_set_nsim(ctl, "orig")

  # Append fixed parameter blocks + MAXEVAL=0 estimation + TABLE
  orig_file <- paste0(run_base, "_orig.fit")
  c(ctl,
    "; --- Fixed estimates for orig run ---",
    .nm_vpc_write_theta(estimates$theta),
    .nm_vpc_write_matrix("OMEGA", estimates$omega),
    .nm_vpc_write_matrix("SIGMA", estimates$sigma),
    "$ESTIMATION MAXEVAL=0",
    "",
    "$TABLE NSIM ID TIME AMT MDV DV PRED IPRED",
    paste0("  FILE=", orig_file),
    "  NOAPPEND NOPRINT NOHEADER",
    ""
  )
}

# Build _sim.ctl: $SIMULATION ONLYSIM, NSIM=IREP
.nm_vpc_build_sim_ctl <- function(ctl_base, estimates, run_base, seed, n_sim) {
  ctl <- .nm_vpc_remove_records(ctl_base)
  ctl <- .nm_vpc_set_nsim(ctl, "sim")

  sim_file <- paste0(run_base, "_sim.fit")
  c(ctl,
    "; --- Fixed estimates for sim run ---",
    .nm_vpc_write_theta(estimates$theta),
    .nm_vpc_write_matrix("OMEGA", estimates$omega),
    .nm_vpc_write_matrix("SIGMA", estimates$sigma),
    sprintf("$SIMULATION (%d) ONLYSIM SUBPROBLEMS=%d",
            as.integer(seed), as.integer(n_sim)),
    "",
    "$TABLE NSIM ID TIME AMT MDV DV PRED IPRED",
    paste0("  FILE=", sim_file),
    "  NOAPPEND NOPRINT NOHEADER",
    ""
  )
}


# =============================================================================
# Main exported function
# =============================================================================

#' Generate VPC control streams and launch NONMEM runs
#'
#' From an estimation run, generates
#' two NONMEM control streams required for a VPC:
#' \itemize{
#'   \item \code{<run>_orig.ctl}: re-runs the model at \code{MAXEVAL=0} with
#'     final estimates fixed, producing observed predictions (\code{_orig.fit}).
#'   \item \code{<run>_sim.ctl}: simulates \code{n_sim} replicates using
#'     \code{$SIMULATION ONLYSIM SUBPROBLEMS=n_sim}, producing
#'     \code{_sim.fit}.
#' }
#' Both control streams are then launched via \code{\link{nm_run}}.
#'
#' @param run_name Name of the estimation run, without extension
#'   (e.g. \code{"run001"}). Its \code{.ctl}, \code{.ext} and \code{.lst}
#'   files are looked for in \code{run001.reg/}, \code{run001.nm7/},
#'   \code{run001/}, then directly in \code{analysis_dir} (flat layout).
#' @param analysis_dir Root directory of the analysis, containing the
#'   estimation run (in a subdirectory or directly).
#' @param nmfe_path Full path to the nmfe executable
#'   (e.g. \code{"C:/nm76g64/run/nmfe76.bat"}).
#' @param n_sim Number of simulation replicates (default: 500).
#' @param seed Random seed for \code{$SIMULATION} (default: 12345).
#' @param skip_terminated Passed to \code{\link{nm_run}}. If \code{TRUE}
#'   (default), already completed VPC runs are skipped.
#' @param run_subdir Layout of the two VPC runs, passed to
#'   \code{\link{nm_run}}. By default (\code{NULL}) the layout of the
#'   estimation run is reused: \code{""} (flat) if its files were found
#'   directly in \code{analysis_dir}, \code{".reg"} otherwise.
#' @param verbose Passed to \code{\link{nm_run}} and used locally for messages.
#' @param n_cores Passed to \code{\link{nm_run}}: 2 runs the \code{_orig} and
#'   \code{_sim} control streams at the same time (default 1).
#'
#' @return A \code{data.frame} run log (from \code{\link{nm_run}}) with one
#'   row per VPC run (\code{_orig} and \code{_sim}).
#'
#' @details
#' The \code{$PK} block of the source control stream must contain a line of
#' the form \code{NSIM = <value>} which will be replaced by \code{NSIM = 1}
#' in \code{_orig.ctl} and \code{NSIM = IREP} in \code{_sim.ctl}.
#'
#' Final parameter estimates are read from the \code{.ext} file of the
#' estimation run (row \code{-1000000000}) and written as \code{FIX} values
#' in both generated control streams.
#'
#' The generated \code{.ctl} files are written to \code{analysis_dir} (not
#' inside a run subdirectory). They can be inspected or edited before running
#' if needed. Wings for NONMEM is not required.
#'
#' @examples
#' \dontrun{
#' log <- nm_vpc_ctl(
#'   run_name     = "run001",
#'   analysis_dir = "C:/analyses/my_pk_project",
#'   nmfe_path    = "C:/nm76g64/run/nmfe76.bat",
#'   n_sim        = 500,
#'   seed         = 12345
#' )
#' print(log)
#' # Then run the VPC:
#' # fig <- nm_vpc("C:/analyses/my_pk_project", run_name = "run001")
#' }
#'
#' @seealso \code{\link{nm_run}}
#'
#' @export
nm_vpc_ctl <- function(
    run_name,
    analysis_dir,
    nmfe_path,
    n_sim            = 500,
    seed             = 12345,
    skip_terminated  = TRUE,
    run_subdir       = NULL,
    verbose          = TRUE,
    n_cores          = 1L
) {

  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)

  # --- Locate estimation run files -------------------------------------------
  reg_dir  <- .nm_find_run_dir(analysis_dir, run_name,
                               need = c(".ctl", ".ext", ".lst"))
  if (is.null(run_subdir))
    run_subdir <- if (normalizePath(reg_dir) == analysis_dir) "" else ".reg"

  ctl_src  <- file.path(reg_dir, paste0(run_name, ".ctl"))
  ext_src  <- file.path(reg_dir, paste0(run_name, ".ext"))
  lst_src  <- file.path(reg_dir, paste0(run_name, ".lst"))

  # --- Read control stream and estimates -------------------------------------
  if (verbose) message("Reading estimation run: ", run_name)

  ctl_base  <- readLines(ctl_src, warn = FALSE)
  final_row <- .nm_vpc_read_ext_final(ext_src)
  estimates <- .nm_vpc_get_estimates(final_row)

  # --- Generate _orig.ctl ----------------------------------------------------
  orig_name <- paste0(run_name, "_orig")
  orig_ctl  <- file.path(analysis_dir, paste0(orig_name, ".ctl"))
  orig_lines <- .nm_vpc_build_orig_ctl(ctl_base, estimates, run_name)
  writeLines(orig_lines, orig_ctl)
  if (verbose) message("Written: ", orig_ctl)

  # --- Generate _sim.ctl -----------------------------------------------------
  sim_name  <- paste0(run_name, "_sim")
  sim_ctl   <- file.path(analysis_dir, paste0(sim_name, ".ctl"))
  sim_lines <- .nm_vpc_build_sim_ctl(ctl_base, estimates, run_name, seed, n_sim)
  writeLines(sim_lines, sim_ctl)
  if (verbose) message("Written: ", sim_ctl)

  # --- Launch both runs via nm_run -------------------------------------------
  if (verbose) message("Launching VPC runs...")

  nm_run(
    nmfe_path       = nmfe_path,
    ctl_files       = c(orig_ctl, sim_ctl),
    analysis_dir    = analysis_dir,
    skip_terminated = skip_terminated,
    run_subdir      = run_subdir,
    verbose         = verbose,
    n_cores         = n_cores
  )
}
