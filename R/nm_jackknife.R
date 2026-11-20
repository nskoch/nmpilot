#  ---------------------------------------------------------------------------
#  nm_jackknife.R  --  Leave-one-subject-out jackknife of a NONMEM model
#  Part of the nmtools package (PsN-/WFN-independent).
#
#  Exported :  nm_jackknife()
#  S3       :  print.nm_jackknife(), plot.nm_jackknife()
#
#  Method (Simon N, nonmem-jackknife-influence, doi:10.5281/zenodo.22257994):
#  remove one subject (ID) at a time, re-estimate the model with NONMEM and
#  measure that subject's INFLUENCE:
#    * primary  : |dOFV| = |OFV(-i) - OFV(ref)|      (likelihood influence);
#    * secondary: relative THETA shifts (th(-i)-th_ref)/th_ref, their L2 norm
#                 ||dtheta||, the max |shift|, and an SE-standardised L2 norm
#                 ||(th(-i)-th_ref)/SE_ref|| when the reference $COV is usable.
#  Subjects are ranked by |dOFV|. The metrics and the reading of OFV/THETA/SE
#  from the .lst are those of the published script; data reading, $DATA
#  patching and the NONMEM launch use the shared nmtools helpers (as
#  nm_bootstrap()): the data lines are copied verbatim, so number formatting
#  and conditional IGNORE rows are preserved.
#  ---------------------------------------------------------------------------

# Signed integer/decimal with optional scientific exponent
.NM_JK_NUM_PAT <- "-?[0-9]*\\.?[0-9]+(?:[Ee][+-]?[0-9]+)?"

# THETA estimates (se = FALSE) or their standard errors (se = TRUE) from a .lst.
.nm_jk_theta_lst <- function(lst_lines, se = FALSE) {
  # several $EST: estimates of the last estimation step
  lst_lines <- .nm_lst_step_lines(lst_lines)
  if (se) {
    h0 <- grep("STANDARD ERROR OF ESTIMATE", lst_lines)
    if (length(h0) == 0) return(numeric(0))
    th <- grep("THETA - VECTOR OF FIXED EFFECTS PARAMETERS", lst_lines)
    th <- th[th > h0[1]]
  } else {
    th <- grep("THETA - VECTOR OF FIXED EFFECTS PARAMETERS", lst_lines)
  }
  if (length(th) == 0) return(numeric(0))
  start <- th[1] + 1
  omega <- grep("OMEGA", lst_lines); omega <- omega[omega > start]
  stop_line <- if (length(omega)) omega[1] - 1 else length(lst_lines)
  vals <- numeric(0)
  for (ln in lst_lines[start:stop_line]) {
    if (!nzchar(trimws(ln))) next
    if (grepl("\\bTH\\s*[0-9]", ln)) next
    for (t in strsplit(trimws(ln), "\\s+")[[1]]) {
      vals <- c(vals,
                if (grepl(paste0("^", .NM_JK_NUM_PAT, "$"), t, perl = TRUE)) as.numeric(t)
                else NA_real_)
    }
  }
  vals
}

# Relative THETA shifts, their L2 norm, max |shift| and SE-standardised norm.
.nm_jk_theta_metrics <- function(th, ref_theta, ref_se = NULL) {
  k   <- length(ref_theta)
  out <- list(rel = rep(NA_real_, k), norm = NA_real_, max = NA_real_,
              z_norm = NA_real_)
  if (k == 0 || length(th) != k) return(out)
  rel <- ifelse(ref_theta != 0, (th - ref_theta) / ref_theta, NA_real_)
  out$rel <- rel
  fin <- is.finite(rel)
  if (any(fin)) {
    out$norm <- sqrt(sum(rel[fin]^2))
    out$max  <- max(abs(rel[fin]))
  }
  if (!is.null(ref_se) && length(ref_se) == k) {
    z  <- ifelse(is.finite(ref_se) & ref_se > 0, (th - ref_theta) / ref_se, NA_real_)
    zf <- is.finite(z)
    if (any(zf)) out$z_norm <- sqrt(sum(z[zf]^2))
  }
  out
}

# Safe file-name fragment for a subject ID.
.nm_jk_name <- function(id) paste0("jk_", gsub("[^A-Za-z0-9.-]", "_", id))


#' Leave-one-subject-out jackknife of a NONMEM model
#'
#' Removes one subject at a time from the data set, re-estimates the model with
#' NONMEM and measures how much each subject influences the fit. The primary
#' metric is the likelihood influence \eqn{|\Delta OFV| = |OFV_{(-i)} -
#' OFV_{ref}|}; the secondary metrics are the relative THETA shifts
#' \eqn{(\theta_{(-i)} - \theta_{ref}) / \theta_{ref}}, their L2 norm, their
#' maximum, and an L2 norm standardised by the reference standard errors when
#' the reference covariance step is usable. Subjects are ranked by
#' \eqn{|\Delta OFV|}. The analysis is descriptive (an influence diagnostic,
#' not a formal outlier test).
#'
#' The method is that of the standalone script nonmem-jackknife-influence
#' (Simon N, doi:10.5281/zenodo.22257994). The data lines are copied verbatim
#' (positional reading through `$INPUT`), so number formatting and rows
#' removed by a conditional `IGNORE=(...)` are preserved. Runs already
#' minimised are skipped on a second call (resume), unless
#' `overwrite_existing_runs = TRUE`.
#'
#' @param ctl Path to the reference control stream (`.ctl`), already run.
#' @param nmfe_path Path to the nmfe launcher (e.g.
#'   `"C:/nm76g64/run/nmfe76.bat"`). Needed when `run = TRUE`.
#' @param lst Reference output. Default: the `.lst` of the same name next to
#'   `ctl`.
#' @param reference_data Data file. `NULL` (default): read from `$DATA`.
#' @param id_col Label of the subject column in `$INPUT` (default `"ID"`).
#' @param comment_char Data comment/ignore character. `NULL` (default): the
#'   single-character `IGNORE` of `$DATA` (`"#"`, `"@"`, a letter...); `""`:
#'   none.
#' @param sep Field separator of the data file (default `","`).
#' @param output_dir Working directory (default `"<ctl_dir>/<run>_jackknife"`).
#'   Holds `jk_data/`, `jk_ctl/`, `jk_run/` (one subdirectory per run),
#'   `figures/` and the two CSV files.
#' @param overwrite_existing_runs `FALSE` (default): skip the runs already
#'   minimised (resume).
#' @param strip_tables Remove `$TABLE` from the jackknife control streams
#'   (default `TRUE`).
#' @param strip_cov Remove `$COVARIANCE` from the jackknife runs (default
#'   `FALSE`). `TRUE` is faster and does not change the influence metrics
#'   (the reference SE come from the reference run).
#' @param dtheta_threshold Relative THETA shift flagged as large (default
#'   0.01).
#' @param save_figures Write the three figures to `output_dir/figures`
#'   (default `TRUE`).
#' @param run `FALSE`: only write the data sets and control streams (no
#'   NONMEM).
#' @param verbose Progress messages (default `TRUE`).
#' @param n_cores Number of jackknife runs at the same time (default 1). A
#'   whole number,
#'   `"auto"` (physical cores minus one, see [nm_cores()]) or `"ask"`
#'   (interactive sessions: the cores are shown and the number is asked once per R
#'   session, then reused; `options(nmtools.n_cores = NULL)` asks again).
#'   Each run is then started in the background in its own directory; if R
#'   is interrupted, the runs already started go on in the background and a
#'   later call with the same arguments picks up their results.
#'
#' @return An object of class `nm_jackknife`: the reference OFV, THETA and SE,
#'   `all_runs` (every replicate, including the failed ones), `summary`
#'   (successful runs ranked by \eqn{|\Delta OFV|}, with `rank`,
#'   `rank_dtheta` and, when available, `rank_dtheta_z`) and the settings.
#'   `jackknife_all_runs.csv` and `jackknife_summary.csv` are written to
#'   `output_dir`.
#'
#' @examples
#' \dontrun{
#' jk <- nm_jackknife("C:/study/run2.reg/run2.ctl",
#'                    nmfe_path = "C:/nm76g64/run/nmfe76.bat")
#' jk
#' plot(jk, "bar")
#' plot(jk, "scatter")
#' }
#' @seealso [nm_bootstrap()], [nm_profile()]
#' @export
nm_jackknife <- function(ctl,
                         nmfe_path               = NULL,
                         lst                     = NULL,
                         reference_data          = NULL,
                         id_col                  = "ID",
                         comment_char            = NULL,
                         sep                     = ",",
                         output_dir              = NULL,
                         overwrite_existing_runs = FALSE,
                         strip_tables            = TRUE,
                         strip_cov               = FALSE,
                         dtheta_threshold        = 0.01,
                         save_figures            = TRUE,
                         run                     = TRUE,
                         verbose                 = TRUE,
                         n_cores                 = 1L) {

  say <- function(...) if (isTRUE(verbose)) message(...)

  if (!file.exists(ctl)) stop("Control stream not found: ", ctl, call. = FALSE)
  reference_ctl <- normalizePath(ctl, winslash = "/", mustWork = TRUE)
  ctl_dir <- dirname(reference_ctl)
  base    <- tools::file_path_sans_ext(basename(reference_ctl))
  if (is.null(lst)) lst <- file.path(ctl_dir, paste0(base, ".lst"))
  reference_lst <- normalizePath(lst, winslash = "/", mustWork = FALSE)
  if (isTRUE(run) && !file.exists(reference_lst))
    stop("Reference .lst not found: ", reference_lst, call. = FALSE)
  if (isTRUE(run) && (is.null(nmfe_path) || !file.exists(nmfe_path)))
    stop("A valid `nmfe_path` is required when run = TRUE.", call. = FALSE)
  if (is.null(output_dir)) output_dir <- file.path(ctl_dir, paste0(base, "_jackknife"))

  dirs <- file.path(output_dir, c("jk_data", "jk_ctl", "jk_run", "figures"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  jk_data_dir <- dirs[1]; jk_ctl_dir <- dirs[2]; jk_run_dir <- dirs[3]; fig_dir <- dirs[4]

  # --- Read the data (positional, formatting preserved) ---------------------
  ctl_lines <- .nm_read_lines(reference_ctl)
  if (is.null(reference_data)) reference_data <- .nm_extract_data_path(reference_ctl)
  input_nm <- .nm_extract_input_names(reference_ctl)
  id_idx   <- .nm_col_index(input_nm, id_col)
  if (is.na(id_idx))
    stop(sprintf("ID column '%s' not found in $INPUT: %s", id_col,
                 paste(input_nm, collapse = " ")), call. = FALSE)
  if (is.null(comment_char)) comment_char <- .nm_data_comment_char(ctl_lines)
  raw <- .nm_read_data(reference_data, comment_char = comment_char, sep = sep, id_idx = id_idx)
  ids <- unique(raw$id)
  N   <- length(ids)
  say(sprintf("Dataset: %d subjects, %d records (ID = col #%d '%s', comment = '%s').",
              N, length(raw$lines), id_idx, id_col, comment_char))
  if (N < 2L) stop("The jackknife needs at least two subjects.", call. = FALSE)

  # --- Leave-one-out data sets + control streams -----------------------------
  base_ctl <- ctl_lines
  if (strip_tables) base_ctl <- .nm_strip_record(base_ctl, "TABLE")
  if (isTRUE(strip_cov)) base_ctl <- .nm_strip_record(base_ctl, "COV")
  di <- grep("^\\s*\\$DATA", base_ctl, ignore.case = TRUE)[1]
  if (is.na(di)) stop("No $DATA record found in the reference control stream.", call. = FALSE)

  names_jk <- .nm_jk_name(ids)
  if (anyDuplicated(names_jk))
    stop("Subject IDs give duplicated file names after cleaning.", call. = FALSE)
  for (k in seq_len(N)) {
    keep <- raw$id != ids[k]
    writeLines(raw$lines[keep], file.path(jk_data_dir, paste0(names_jk[k], ".csv")))
    new_ctl <- base_ctl
    new_ctl[di] <- .nm_patch_data_line(new_ctl[di], paste0(names_jk[k], ".csv"))
    writeLines(new_ctl, file.path(jk_ctl_dir, paste0(names_jk[k], ".ctl")))
  }
  say(N, " jackknife data sets and control streams written to ", output_dir)

  obj <- list(call = match.call(), reference_ctl = reference_ctl,
              reference_lst = reference_lst, reference_data = reference_data,
              output_dir = output_dir, n_subjects = N, ids = ids,
              ref_ofv = NA_real_, ref_theta = numeric(0), ref_se = numeric(0),
              k_theta = 0L, all_runs = NULL, summary = NULL,
              dtheta_threshold = dtheta_threshold, run = FALSE)
  if (!isTRUE(run)) {
    say("run = FALSE: data sets and control streams only, NONMEM not launched.")
    return(structure(obj, class = "nm_jackknife"))
  }

  # --- Reference run ---------------------------------------------------------
  ref_lines  <- readLines(reference_lst, warn = FALSE)
  ref_status <- .nm_parse_status(reference_lst)
  ref_ofv    <- .nm_extract_ofv(reference_lst)
  ref_theta  <- .nm_jk_theta_lst(ref_lines, se = FALSE)
  ref_se     <- .nm_jk_theta_lst(ref_lines, se = TRUE)
  k_theta    <- length(ref_theta)
  if (is.na(ref_ofv))
    stop("Could not read the reference OFV (#OBJV) from ", reference_lst, call. = FALSE)
  say("Reference OFV: ", ref_ofv)
  if (!isTRUE(ref_status$min_ok))
    warning("The reference run did not minimise successfully: dOFV are relative ",
            "to a non-converged reference.", call. = FALSE)
  if (k_theta == 0)
    warning("No THETA vector read from the reference .lst: dtheta metrics are NA.",
            call. = FALSE)
  if (length(ref_se) != k_theta || all(is.na(ref_se))) {
    say("No usable THETA SE in the reference run: the SE-standardised distance is NA.")
    ref_se <- rep(NA_real_, k_theta)
  }

  # --- NONMEM runs (resume-aware), each in its own subdirectory of jk_run/ ---
  n_cores <- .nm_resolve_cores(n_cores)
  prep_run <- function(k) {
    name <- names_jk[k]
    rdir <- file.path(jk_run_dir, name)
    lst_file <- file.path(rdir, paste0(name, ".lst"))
    if (!overwrite_existing_runs && file.exists(lst_file) &&
        isTRUE(.nm_parse_status(lst_file)$min_ok)) {
      say("Skipping (already minimised): ", name)
      return(NULL)
    }
    dir.create(rdir, showWarnings = FALSE)
    file.copy(file.path(jk_ctl_dir,  paste0(name, ".ctl")), rdir, overwrite = TRUE)
    file.copy(file.path(jk_data_dir, paste0(name, ".csv")), rdir, overwrite = TRUE)
    list(dir = rdir, name = name)
  }
  warn_status <- function(status, name)
    if (!identical(as.integer(status), 0L))
      warning("NONMEM returned a non-zero exit status (", status, ") for ", name, ".",
              call. = FALSE)
  if (n_cores > 1L) {
    say(sprintf("Running %d jackknife runs, %d at a time", N, n_cores))
    .nm_run_queue(as.list(seq_len(N)), nmfe_path, n_cores, verbose = verbose,
                  label = "jackknife runs",
                  prepare = function(k, slot) prep_run(k),
                  finish = function(k, slot, st, secs) warn_status(st, names_jk[k]))
  } else {
    for (k in seq_len(N)) {
      spec <- prep_run(k)
      if (is.null(spec)) next
      say(sprintf("Running %s (%d/%d)", spec$name, k, N))
      warn_status(.nm_launch(nmfe_path, spec$dir, spec$name), spec$name)
    }
  }

  # --- Results ---------------------------------------------------------------
  rows <- lapply(seq_len(N), function(k) {
    f <- file.path(jk_run_dir, names_jk[k], paste0(names_jk[k], ".lst"))
    if (!file.exists(f)) {
      st <- data.frame(min_ok = FALSE, cov_ok = FALSE, cov_warning = FALSE, nmtran_error = NA)
      ofv <- NA_real_; th <- numeric(0)
    } else {
      st  <- .nm_parse_status(f)
      ofv <- .nm_extract_ofv(f)
      th  <- .nm_jk_theta_lst(readLines(f, warn = FALSE), se = FALSE)
    }
    tm <- .nm_jk_theta_metrics(th, ref_theta, ref_se)
    row <- data.frame(id_removed = ids[k], file = basename(f), ofv = ofv,
                      dtheta_norm = tm$norm, dtheta_max = tm$max, dtheta_znorm = tm$z_norm,
                      st, run_ok = st$min_ok, stringsAsFactors = FALSE)
    if (k_theta > 0)
      row <- cbind(row, stats::setNames(as.data.frame(as.list(tm$rel)),
                                        paste0("dtheta_rel_", seq_len(k_theta))))
    row
  })
  res_all <- do.call(rbind, rows)
  res_all$ref_ofv  <- ref_ofv
  res_all$dOFV     <- res_all$ofv - ref_ofv
  res_all$abs_dOFV <- abs(res_all$dOFV)
  res_all$large_dtheta <- !is.na(res_all$dtheta_max) & res_all$dtheta_max > dtheta_threshold
  utils::write.csv(res_all, file.path(output_dir, "jackknife_all_runs.csv"), row.names = FALSE)

  n_failed <- sum(!res_all$run_ok); n_no_ofv <- sum(is.na(res_all$dOFV))
  if (n_failed > 0 || n_no_ofv > 0)
    say("Excluded from the ranking: ", n_failed, " non-minimised and ",
        n_no_ofv, " missing-OFV run(s).")

  res <- res_all[res_all$run_ok & is.finite(res_all$dOFV), , drop = FALSE]
  if (!nrow(res))
    stop("No usable run (all failed or without OFV). See ", jk_run_dir, call. = FALSE)
  res <- res[order(-res$abs_dOFV), , drop = FALSE]
  res$rank <- seq_len(nrow(res))
  res$rank_dtheta <- rank(-res$dtheta_norm, ties.method = "first", na.last = "keep")
  if (any(is.finite(res$dtheta_znorm)))
    res$rank_dtheta_z <- rank(-res$dtheta_znorm, ties.method = "first", na.last = "keep")
  rownames(res) <- NULL
  utils::write.csv(res, file.path(output_dir, "jackknife_summary.csv"), row.names = FALSE)
  say("Summary written to: ", file.path(output_dir, "jackknife_summary.csv"))

  obj[c("ref_ofv", "ref_theta", "ref_se", "k_theta", "all_runs", "summary", "run")] <-
    list(ref_ofv, ref_theta, ref_se, k_theta, res_all, res, TRUE)
  obj <- structure(obj, class = "nm_jackknife")

  if (isTRUE(save_figures)) {
    nm_export_plot(plot(obj, "bar"),  file.path(fig_dir, "abs_dOFV_bar.png"),
                   width = 6, height = 5, verbose = FALSE)
    nm_export_plot(plot(obj, "rank"), file.path(fig_dir, "rank_dOFV.png"),
                   width = 7, height = 4, verbose = FALSE)
    if (k_theta > 0 && any(is.finite(res$dtheta_norm)))
      nm_export_plot(plot(obj, "scatter"), file.path(fig_dir, "dOFV_vs_dtheta.png"),
                     width = 6, height = 5, verbose = FALSE)
    say("Figures saved to: ", fig_dir)
  }
  obj
}


#' @param x An `nm_jackknife` object.
#' @param digits Significant digits in the printed table (default 4).
#' @param ... Ignored.
#' @rdname nm_jackknife
#' @export
print.nm_jackknife <- function(x, digits = 4, ...) {
  cat("<nm_jackknife>\n")
  cat(sprintf("  reference : %s\n", basename(x$reference_ctl)))
  cat(sprintf("  run dir   : %s\n", x$output_dir))
  cat(sprintf("  subjects  : %d\n", x$n_subjects))
  if (!isTRUE(x$run) || is.null(x$summary)) {
    cat("  (dry run: data sets and control streams only)\n")
    return(invisible(x))
  }
  # OFV: three decimals, as NONMEM prints it (signif() gave -1224)
  cat(sprintf("  ref OFV   : %.3f\n", x$ref_ofv))
  cat(sprintf("  runs used : %d (of %d)\n", nrow(x$summary), nrow(x$all_runs)))
  top <- utils::head(x$summary[, c("rank", "id_removed", "dOFV", "abs_dOFV",
                                   "dtheta_norm", "dtheta_max")], 10)
  num <- vapply(top, is.numeric, logical(1)); top[num] <- lapply(top[num], signif, digits)
  cat("  most influential subjects (by |dOFV|):\n")
  print(top, row.names = FALSE)
  invisible(x)
}

#' Plots for `nm_jackknife` objects
#'
#' @param x An `nm_jackknife` object.
#' @param type `"bar"` (default): \eqn{|\Delta OFV|} per subject; `"rank"`:
#'   ranked \eqn{|\Delta OFV|}; `"scatter"`: \eqn{|\Delta OFV|} against the
#'   relative THETA shift norm.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.nm_jackknife <- function(x, type = c("bar", "rank", "scatter"), ...) {
  type <- match.arg(type)
  res <- x$summary
  if (is.null(res) || !nrow(res)) stop("No summary to plot (dry run?).", call. = FALSE)
  res$id_removed <- as.character(res$id_removed)

  if (type == "bar") {
    return(ggplot2::ggplot(res, ggplot2::aes(
             x = stats::reorder(factor(.data$id_removed), .data$abs_dOFV),
             y = .data$abs_dOFV)) +
             ggplot2::geom_col(fill = "steelblue") + ggplot2::coord_flip() +
             ggplot2::theme_bw() +
             ggplot2::labs(x = "Subject removed", y = "|\u0394OFV|",
                           title = "Jackknife: influence by subject"))
  }
  if (type == "rank") {
    return(ggplot2::ggplot(res, ggplot2::aes(x = .data$rank, y = .data$abs_dOFV,
                                             label = .data$id_removed)) +
             ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
             ggplot2::geom_text(vjust = -0.8, size = 3) + ggplot2::theme_bw() +
             ggplot2::labs(x = "Rank", y = "|\u0394OFV|",
                           title = "Jackknife: influence ranking"))
  }
  if (!(x$k_theta > 0 && any(is.finite(res$dtheta_norm))))
    stop("No THETA shift to plot (the scatter needs the reference THETA vector).",
         call. = FALSE)
  ggplot2::ggplot(res, ggplot2::aes(x = .data$dtheta_norm, y = .data$abs_dOFV,
                                    label = .data$id_removed)) +
    ggplot2::geom_point(size = 2) + ggplot2::geom_text(vjust = -0.8, size = 3) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "\u2016\u0394\u03b8\u2016 (relative parameter shift)", y = "|\u0394OFV|",
                  title = "Jackknife: likelihood influence vs parameter sensitivity")
}
