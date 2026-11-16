# =============================================================================
# nm_compare.R
# Paired comparison of two NONMEM runs (CWRES, NPDE, iOFV)
# Part of the nmtools package
# =============================================================================

#' Compare two NONMEM runs
#'
#' Paired comparison of two NONMEM runs at the observation level (CWRES, NPDE)
#' and at the subject level (individual OFV from `.phi`). Provides a summary
#' of residual diagnostics, per-subject iOFV deltas, and identifies which
#' subjects drive the OFV difference.
#'
#' Observations are merged on `join_keys` (default `c("ID", "TIME")`). If the
#' two datasets differ (e.g. different BLQ handling), only the intersection
#' is retained and a warning is issued.
#'
#' @param run1,run2 Path to a run directory or run prefix (without extension).
#'   Each location must contain at least a `sdtab*` table. `.phi` and `.lst`
#'   files are used if available for iOFV and global OFV respectively.
#' @param labels Character vector of length 2. Display labels for the two runs.
#' @param join_keys Character vector. Columns used to pair observations between
#'   runs. Default `c("ID", "TIME")`. Add `"OCC"` for BOV models, `"CMT"` for
#'   multi-DV models.
#'
#' @param n_par Numbers of estimated parameters of the two runs (length 2).
#'   Default: read from the `.ext` files (parameters not fixed).
#' @param nested `TRUE` if one model is nested in the other: a likelihood
#'   ratio test is added to the Vuong and Clarke tests.
#' @param covariates Covariates (columns of the data set) against which the
#'   per-subject delta iOFV is examined.
#' @param data Data set (path or data.frame) holding the covariates. Default:
#'   the `$DATA` file of run 1.
#' @param n_boot Number of subject resamples for the bootstrap of the delta
#'   iOFV.
#' @param seed Seed of the bootstrap.
#'
#' @details
#' **Tests on the individual contributions.** With `d_i = (iOFV_1,i -
#' iOFV_2,i)/2`, the difference of the individual log-likelihoods, the Vuong
#' (1989) statistic is `z = (sum(d_i) - c) / (sqrt(N) * sd(d_i))` and the
#' Clarke (2007) test is a sign test on `d_i - c/N`, where `N` is the number
#' of subjects and `c` the correction for the number of parameters: none,
#' `k_2 - k_1` (AIC-like) or `(k_2 - k_1) log(N) / 2` (BIC-like). Both apply
#' to non-nested models fitted to the same data; the subjects are the
#' independent units. The individual contributions are those of the
#' estimation method (FOCE: the iOFV of the `.phi` file).
#'
#' **Bootstrap.** The subjects are resampled with their delta iOFV, without
#' re-estimation: the distribution measures how much the difference depends
#' on a few subjects, conditionally on the estimates (optimistic for the
#' larger model).
#'
#' **Influence.** `n_reverse` is the number of subjects with the largest
#' contributions whose removal reverses the sign of dOFV (as `dOFV.vs.id` in
#' xpose4).
#'
#' @references
#' Vuong QH. Likelihood ratio tests for model selection and non-nested
#' hypotheses. Econometrica 1989;57:307-33.
#'
#' Clarke KA. A simple distribution-free test for nonnested model selection.
#' Political Analysis 2007;15:347-63.
#'
#' @return An object of class `nm_compare`, a list with:
#'   \describe{
#'     \item{`data_obs`}{Observation-level paired data.frame with `_1`/`_2`
#'       suffixed columns and delta columns (`dCWRES`, `dAbsCWRES`, `dNPDE`,
#'       `dAbsNPDE`).}
#'     \item{`data_subj`}{Per-subject summary including `iOFV_1`, `iOFV_2`,
#'       `diOFV`, and residual aggregates.}
#'     \item{`metrics`}{Global metrics table (OFV, RMSE(CWRES), Mean |CWRES|,
#'       proportions of outliers, NPDE moments). Wilcoxon paired p-value
#'       stored as attribute `wilcox_absCWRES_p`.}
#'     \item{`meta`}{Labels, global OFV, numbers of parameters, sample sizes,
#'       dataset compatibility check.}
#'     \item{`tests`}{Vuong, Clarke (and LRT) tests, subject bootstrap and
#'       influence (`NULL` without `.phi` files).}
#'     \item{`parameters`}{Estimates, SE, RSE and shrinkage of the two runs
#'       matched by parameter name, with their ratio.}
#'     \item{`etas`}{Paired EBEs of the ETAs present in both runs.}
#'     \item{`covariate_tests`}{Delta iOFV versus each covariate.}
#'     \item{`run_info`}{Run-level criteria (RSE, condition number,
#'       shrinkage, status).}
#'   }
#'
#' @examples
#' \dontrun{
#'   cmp <- nm_compare("runs/run42", "runs/run43",
#'                     labels = c("Base", "+ CLcr"))
#'   print(cmp)
#'   plot(cmp)
#'   plot(cmp, which = "iofv_caterpillar")
#'   summary(cmp)
#'   plot(cmp, which = c("gain_map", "overlay_gof"))
#' }
#'
#' @seealso [nm_extract_runs()], [plot.nm_compare()], [print.nm_compare()]
#' @importFrom rlang warn .data
#' @importFrom utils read.table
#' @importFrom stats sd wilcox.test
#' @importFrom tools file_path_sans_ext
#' @export
nm_compare <- function(run1, run2,
                       labels = c("Run 1", "Run 2"),
                       join_keys = c("ID", "TIME"),
                       n_par = NULL, nested = FALSE,
                       covariates = NULL, data = NULL,
                       n_boot = 2000L, seed = 2026L) {

  stopifnot(length(labels) == 2L, is.character(join_keys))

  d1 <- .nm_cmp_read_run(run1)
  d2 <- .nm_cmp_read_run(run2)

  check <- .nm_cmp_check_compat(d1$obs, d2$obs, keys = join_keys)
  if (!check$identical) {
    rlang::warn(sprintf(
      paste0("Datasets differ: %d common obs, %d only in run1, %d only in ",
             "run2. Comparison restricted to intersection."),
      check$n_common, check$n_only1, check$n_only2
    ))
  }

  obs_cols <- c(join_keys, "DV", "PRED", "IPRED", "CWRES", "NPDE")
  obs_cols <- intersect(obs_cols, intersect(names(d1$obs), names(d2$obs)))

  data_obs <- merge(
    d1$obs[, obs_cols, drop = FALSE],
    d2$obs[, obs_cols, drop = FALSE],
    by = join_keys, suffixes = c("_1", "_2")
  )

  if (all(c("CWRES_1", "CWRES_2") %in% names(data_obs))) {
    data_obs$dCWRES    <- data_obs$CWRES_2 - data_obs$CWRES_1
    data_obs$dAbsCWRES <- abs(data_obs$CWRES_2) - abs(data_obs$CWRES_1)
  }
  if (all(c("NPDE_1", "NPDE_2") %in% names(data_obs))) {
    data_obs$dNPDE    <- data_obs$NPDE_2 - data_obs$NPDE_1
    data_obs$dAbsNPDE <- abs(data_obs$NPDE_2) - abs(data_obs$NPDE_1)
  }

  data_subj <- .nm_cmp_summarize_subj(data_obs, d1$phi, d2$phi)
  metrics   <- .nm_cmp_global_metrics(data_obs, d1$ofv, d2$ofv)

  # Run-level information, parameters, ETAs
  i1 <- .nm_cmp_runinfo(run1, d1$paths)
  i2 <- .nm_cmp_runinfo(run2, d2$paths)
  info_list <- lapply(list(i1, i2), function(i) {
    if (is.null(i)) return(NULL)
    pr <- i$params[i$params$type %in% c("THETA", "OMEGA", "SIGMA") & !i$params$fixed, , drop = FALSE]
    list(n_par = i$run$n_par, cond = i$cond, max_eta_shrink = i$run$max_eta_shrink,
         rse_max = if (nrow(pr) && any(is.finite(pr$rse_pct))) max(pr$rse_pct, na.rm = TRUE) else NA_real_,
         min_status = i$run$min_status, cov_status = i$run$cov_status)
  })
  if (is.null(n_par)) {
    n_par <- vapply(info_list, function(i) if (is.null(i) || is.null(i$n_par)) NA_real_ else as.numeric(i$n_par), 0)
  }
  stopifnot(length(n_par) == 2L)
  ofv <- c(d1$ofv, d2$ofv)
  if (any(!is.finite(ofv)) && !is.null(d1$phi) && !is.null(d2$phi) &&
      all(c("OBJ") %in% intersect(names(d1$phi), names(d2$phi)))) {
    ofv[!is.finite(ofv)] <- c(sum(d1$phi$OBJ), sum(d2$phi$OBJ))[!is.finite(ofv)]
  }
  tests <- if (all(c("iOFV_1", "iOFV_2") %in% names(data_subj)))
    .nm_cmp_tests(data_subj, n_par, ofv, nrow(data_obs), nested, n_boot, seed, labels) else NULL
  parameters <- .nm_cmp_params(i1$params, i2$params)
  etas <- .nm_cmp_etas(d1$phi, d2$phi, i1$params, i2$params)
  cov_tests <- cov_obj <- NULL
  if (length(covariates)) {
    cov <- .nm_w_covariates(data, if (is.null(i1)) NA else i1$ctl)
    miss <- setdiff(covariates, names(cov))
    if (length(miss)) stop("Covariate(s) not in the data: ", paste(miss, collapse = ", "), call. = FALSE)
    cov_tests <- .nm_cmp_cov_tests(data_subj, cov, covariates)
    cov_obj <- list(data = cov[, c("ID", covariates), drop = FALSE], names = covariates)
  }

  out <- list(
    data_obs  = data_obs,
    data_subj = data_subj,
    metrics   = metrics,
    tests     = tests,
    parameters = parameters,
    etas      = etas,
    covariate_tests = cov_tests,
    covariates = cov_obj,
    run_info  = info_list,
    meta = list(
      labels        = labels,
      ofv           = ofv,
      n_par         = n_par,
      n_obs         = nrow(data_obs),
      n_subj        = nrow(data_subj),
      compatibility = check,
      join_keys     = join_keys
    )
  )
  class(out) <- "nm_compare"
  out
}


# -----------------------------------------------------------------------------
# Print method
# -----------------------------------------------------------------------------

#' @rdname nm_compare
#' @param x An `nm_compare` object.
#' @param ... Unused.
#' @export
print.nm_compare <- function(x, ...) {
  cat("=== NONMEM run comparison ===\n")
  cat(sprintf("Run 1 : %s  (OFV = %s)\n",
              x$meta$labels[1], .nm_cmp_fmt(x$meta$ofv[1])))
  cat(sprintf("Run 2 : %s  (OFV = %s)\n",
              x$meta$labels[2], .nm_cmp_fmt(x$meta$ofv[2])))
  if (all(is.finite(x$meta$ofv))) {
    cat(sprintf("Delta OFV : %+.3f\n", x$meta$ofv[2] - x$meta$ofv[1]))
  }
  cat(sprintf("N subjects : %d | N paired obs : %d\n",
              x$meta$n_subj, x$meta$n_obs))
  cat("\n--- Global metrics ---\n")
  print(x$metrics, row.names = FALSE)
  p <- attr(x$metrics, "wilcox_absCWRES_p")
  if (!is.null(p)) {
    cat(sprintf("\nPaired Wilcoxon on |CWRES| : p = %.4g\n", p))
  }
  if (!is.null(x$tests)) {
    cat("\n--- Tests on individual contributions (iOFV) ---\n")
    tt <- x$tests$table
    tt$p <- format.pval(tt$p, digits = 3)
    print(tt[, c("test", "correction", "statistic", "p", "favours")], row.names = FALSE)
    b <- x$tests$bootstrap
    cat(sprintf("\nSubject bootstrap (%d, no refit): dOFV 95%% CI [%.1f; %.1f], P(%s better) = %.3f\n",
                b$n_boot, b$ci[1], b$ci[3], x$meta$labels[2], b$p_run2_better))
    cat(sprintf("Subjects to remove to reverse dOFV: %s of %d\n",
                format(x$tests$influence$n_reverse), x$tests$influence$n_subjects))
  }
  if (!is.null(x$covariate_tests)) {
    cat("\n--- Delta iOFV versus covariates ---\n")
    ct <- x$covariate_tests; ct$p <- signif(ct$p, 3); ct$p_adj <- signif(ct$p_adj, 3)
    print(ct, row.names = FALSE)
  }
  if (!is.null(x$parameters)) {
    cat("\n--- Parameters (matched by name) ---\n")
    pp <- x$parameters[x$parameters$type %in% c("THETA", "OMEGA", "SIGMA"),
                       c("type", "name", "value_1", "rse_pct_1", "value_2", "rse_pct_2", "ratio"), drop = FALSE]
    num <- vapply(pp, is.numeric, logical(1)); pp[num] <- lapply(pp[num], signif, 4)
    names(pp)[3:6] <- c(paste("est", x$meta$labels[1]), "RSE%", paste("est", x$meta$labels[2]), "RSE% ")
    print(pp, row.names = FALSE)
  }
  cat("\nsummary() gives the decision criteria; plot() the paired diagnostics.\n")
  invisible(x)
}


# -----------------------------------------------------------------------------
# Plot method
# -----------------------------------------------------------------------------

#' Diagnostic plots for `nm_compare` objects
#'
#' @param x An `nm_compare` object.
#' @param which Character. One or several of `"cwres_scatter"`,
#'   `"npde_scatter"`, `"iofv_caterpillar"`, `"delta_density"`, `"qq_npde"`,
#'   `"gain_map"` (Delta|CWRES| versus TIME and PRED, loess with a 95 %
#'   subject-bootstrap band), `"overlay_gof"` (diagnostics of both runs
#'   superimposed), `"bootstrap"` (subject bootstrap of dOFV),
#'   `"iofv_covariate"` (delta iOFV versus the covariates),
#'   `"parameters"` (estimates with 95 % Wald intervals), `"eta_scatter"`
#'   (EBEs of the two runs), or `"all"` (default).
#' @param n_band Number of subject resamples for the bands of `"gain_map"`
#'   and `"overlay_gof"`.
#' @param ... Unused.
#'
#' @return A `ggplot` object (single plot) or a list of `ggplot` objects
#'   (multiple plots). Requires the `ggplot2` package.
#'
#' @export
plot.nm_compare <- function(x, which = "all", n_band = 200L, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot.nm_compare(). ",
         "Install it with install.packages('ggplot2').", call. = FALSE)
  }

  all_plots <- c("cwres_scatter", "npde_scatter", "iofv_caterpillar",
                 "delta_density", "qq_npde", "gain_map", "overlay_gof",
                 "bootstrap", "iofv_covariate", "parameters", "eta_scatter")
  if (identical(which, "all")) which <- all_plots
  which <- match.arg(which, all_plots, several.ok = TRUE)

  d   <- x$data_obs
  s   <- x$data_subj
  lab <- x$meta$labels
  plots <- list()

  if ("cwres_scatter" %in% which &&
      all(c("CWRES_1", "CWRES_2") %in% names(d))) {
    plots$cwres_scatter <- ggplot2::ggplot(
      d, ggplot2::aes(x = .data$CWRES_1, y = .data$CWRES_2)
    ) +
      ggplot2::geom_hline(yintercept = c(-2, 2),
                          linetype = "dotted", colour = "grey50") +
      ggplot2::geom_vline(xintercept = c(-2, 2),
                          linetype = "dotted", colour = "grey50") +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red") +
      ggplot2::geom_point(alpha = 0.4, size = 1) +
      ggplot2::labs(
        title = "CWRES: paired comparison",
        subtitle = "Points below diagonal: run 2 improves fit",
        x = paste0("CWRES (", lab[1], ")"),
        y = paste0("CWRES (", lab[2], ")")
      ) +
      ggplot2::coord_equal() +
      ggplot2::theme_bw()
  }

  if ("npde_scatter" %in% which &&
      all(c("NPDE_1", "NPDE_2") %in% names(d))) {
    plots$npde_scatter <- ggplot2::ggplot(
      d, ggplot2::aes(x = .data$NPDE_1, y = .data$NPDE_2)
    ) +
      ggplot2::geom_hline(yintercept = c(-1.96, 1.96),
                          linetype = "dotted", colour = "grey50") +
      ggplot2::geom_vline(xintercept = c(-1.96, 1.96),
                          linetype = "dotted", colour = "grey50") +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red") +
      ggplot2::geom_point(alpha = 0.4, size = 1) +
      ggplot2::labs(
        title = "NPDE: paired comparison",
        x = paste0("NPDE (", lab[1], ")"),
        y = paste0("NPDE (", lab[2], ")")
      ) +
      ggplot2::coord_equal() +
      ggplot2::theme_bw()
  }

  if ("iofv_caterpillar" %in% which && "diOFV" %in% names(s)) {
    s2 <- s[order(s$diOFV), , drop = FALSE]
    s2$rank <- seq_len(nrow(s2))
    s2$improved <- s2$diOFV < 0
    plots$iofv_caterpillar <- ggplot2::ggplot(
      s2, ggplot2::aes(x = .data$rank, y = .data$diOFV,
                       fill = .data$improved)
    ) +
      ggplot2::geom_col() +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::scale_fill_manual(
        values = c(`TRUE` = "#2c7fb8", `FALSE` = "#d95f0e"),
        labels = c(`TRUE`  = paste(lab[2], "better"),
                   `FALSE` = paste(lab[1], "better"))
      ) +
      ggplot2::labs(
        title    = "Per-subject delta iOFV (run 2 - run 1)",
        subtitle = "Negative = run 2 improves individual fit",
        x = "Subjects (sorted)", y = "Delta iOFV", fill = NULL
      ) +
      ggplot2::theme_bw()
  }

  if ("delta_density" %in% which && "dAbsCWRES" %in% names(d)) {
    plots$delta_density <- ggplot2::ggplot(
      d, ggplot2::aes(x = .data$dAbsCWRES)
    ) +
      ggplot2::geom_density(fill = "#2c7fb8", alpha = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::labs(
        title    = "Distribution of Delta|CWRES|",
        subtitle = "< 0: run 2 reduces residual magnitude",
        x = "|CWRES run2| - |CWRES run1|", y = "Density"
      ) +
      ggplot2::theme_bw()
  }

  if ("qq_npde" %in% which &&
      all(c("NPDE_1", "NPDE_2") %in% names(d))) {
    qq_df <- rbind(
      data.frame(run = lab[1], npde = d$NPDE_1),
      data.frame(run = lab[2], npde = d$NPDE_2)
    )
    plots$qq_npde <- ggplot2::ggplot(
      qq_df, ggplot2::aes(sample = .data$npde, colour = .data$run)
    ) +
      ggplot2::stat_qq(alpha = 0.5) +
      ggplot2::stat_qq_line() +
      ggplot2::labs(
        title = "QQ-plot NPDE vs N(0,1)",
        x = "Theoretical quantiles", y = "Sample quantiles", colour = NULL
      ) +
      ggplot2::theme_bw()
  }

  if ("gain_map" %in% which) plots$gain_map <- .nm_cmp_plot_gain(x, n_band)
  if ("overlay_gof" %in% which) plots$overlay_gof <- .nm_cmp_plot_overlay(x, n_band)
  if ("bootstrap" %in% which) plots$bootstrap <- .nm_cmp_plot_bootstrap(x)
  if ("iofv_covariate" %in% which) {
    pc <- .nm_cmp_plot_cov(x)
    if (inherits(pc, "ggplot")) plots$iofv_covariate <- pc
    else if (length(pc)) plots[paste0("iofv_covariate_", names(pc))] <- pc
  }
  if ("parameters" %in% which) plots$parameters <- .nm_cmp_plot_params(x)
  if ("eta_scatter" %in% which) plots$eta_scatter <- .nm_cmp_plot_etas(x)
  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 1L) plots[[1]] else plots
}


# =============================================================================
# Internal helpers (prefix `.nm_cmp_`, not exported)
# =============================================================================

#' Read a NONMEM run (sdtab + .phi + OFV from .lst)
#' @keywords internal
#' @noRd
.nm_cmp_read_run <- function(path) {
  info <- .nm_cmp_resolve_paths(path)

  if (is.null(info$sdtab) || !file.exists(info$sdtab)) {
    stop("No sdtab found for run: ", path, call. = FALSE)
  }

  obs <- .nm_cmp_read_table(info$sdtab)
  phi <- if (!is.null(info$phi) && file.exists(info$phi)) {
    .nm_cmp_read_table(info$phi)
  } else NULL
  ofv <- if (!is.null(info$lst) && file.exists(info$lst)) {
    .nm_cmp_extract_ofv(info$lst)
  } else NA_real_

  if ("EVID" %in% names(obs)) {
    obs <- obs[obs$EVID == 0, , drop = FALSE]
  } else if ("MDV" %in% names(obs)) {
    obs <- obs[obs$MDV == 0, , drop = FALSE]
  }

  list(obs = obs, phi = phi, ofv = ofv, paths = info)
}

#' Resolve sdtab / phi / lst paths from a directory or prefix
#' @keywords internal
#' @noRd
.nm_cmp_resolve_paths <- function(path) {
  if (dir.exists(path)) {
    sdtab_files <- list.files(path, pattern = "^sdtab",      full.names = TRUE)
    phi_files   <- list.files(path, pattern = "\\.phi$",     full.names = TRUE)
    lst_files   <- list.files(path, pattern = "\\.(lst|res)$", full.names = TRUE)
    list(
      sdtab = if (length(sdtab_files)) sdtab_files[1] else NULL,
      phi   = if (length(phi_files))   phi_files[1]   else NULL,
      lst   = if (length(lst_files))   lst_files[1]   else NULL
    )
  } else {
    dir  <- dirname(path)
    base <- tools::file_path_sans_ext(basename(path))
    list(
      sdtab = if (file.exists(path)) path
              else file.path(dir, paste0("sdtab", base)),
      phi   = file.path(dir, paste0(base, ".phi")),
      lst   = file.path(dir, paste0(base, ".lst"))
    )
  }
}

#' Read a NONMEM-formatted table (skip `TABLE NO.` header if present)
#' @keywords internal
#' @noRd
.nm_cmp_read_table <- function(file) {
  first <- readLines(file, n = 1L, warn = FALSE)
  skip  <- if (grepl("^TABLE NO", first)) 1L else 0L
  utils::read.table(file, skip = skip, header = TRUE,
                    stringsAsFactors = FALSE, check.names = FALSE)
}

#' Extract final OFV from a `.lst` file
#' @keywords internal
#' @noRd
.nm_cmp_extract_ofv <- function(lst_file) {
  lines <- readLines(lst_file, warn = FALSE)
  idx   <- grep("#OBJV", lines)
  if (!length(idx)) return(NA_real_)
  last  <- lines[idx[length(idx)]]
  val   <- suppressWarnings(as.numeric(
    regmatches(last, regexpr("-?\\d+\\.\\d+", last))
  ))
  if (!length(val)) NA_real_ else val
}

#' Check observation-level compatibility between two runs
#' @keywords internal
#' @noRd
.nm_cmp_check_compat <- function(obs1, obs2, keys) {
  k1 <- do.call(paste, c(obs1[, keys, drop = FALSE], sep = "_"))
  k2 <- do.call(paste, c(obs2[, keys, drop = FALSE], sep = "_"))
  list(
    identical = setequal(k1, k2),
    n_common  = length(intersect(k1, k2)),
    n_only1   = length(setdiff(k1, k2)),
    n_only2   = length(setdiff(k2, k1))
  )
}

#' Per-subject aggregation of residuals and iOFV merge
#' @keywords internal
#' @noRd
.nm_cmp_summarize_subj <- function(data_obs, phi1, phi2) {

  agg_one <- function(d) {
    data.frame(
      ID             = d$ID[1],
      n_obs          = nrow(d),
      rmse_CWRES_1   = if ("CWRES_1" %in% names(d)) sqrt(mean(d$CWRES_1^2)) else NA_real_,
      rmse_CWRES_2   = if ("CWRES_2" %in% names(d)) sqrt(mean(d$CWRES_2^2)) else NA_real_,
      mean_absNPDE_1 = if ("NPDE_1"  %in% names(d)) mean(abs(d$NPDE_1))     else NA_real_,
      mean_absNPDE_2 = if ("NPDE_2"  %in% names(d)) mean(abs(d$NPDE_2))     else NA_real_
    )
  }

  subj <- do.call(rbind, lapply(split(data_obs, data_obs$ID), agg_one))
  rownames(subj) <- NULL

  if (!is.null(phi1) && "OBJ" %in% names(phi1)) {
    subj <- merge(subj, phi1[, c("ID", "OBJ")], by = "ID", all.x = TRUE)
    names(subj)[names(subj) == "OBJ"] <- "iOFV_1"
  }
  if (!is.null(phi2) && "OBJ" %in% names(phi2)) {
    subj <- merge(subj, phi2[, c("ID", "OBJ")], by = "ID", all.x = TRUE)
    names(subj)[names(subj) == "OBJ"] <- "iOFV_2"
  }
  if (all(c("iOFV_1", "iOFV_2") %in% names(subj))) {
    subj$diOFV <- subj$iOFV_2 - subj$iOFV_1
  }
  subj
}

#' Global metrics across all paired observations
#' @keywords internal
#' @noRd
.nm_cmp_global_metrics <- function(data_obs, ofv1, ofv2) {

  m <- data.frame(metric = character(), run1 = numeric(),
                  run2 = numeric(), delta = numeric(),
                  stringsAsFactors = FALSE)
  add <- function(name, v1, v2) {
    m[nrow(m) + 1L, ] <<- list(name, v1, v2, v2 - v1)
  }

  add("OFV", ofv1, ofv2)

  if (all(c("CWRES_1", "CWRES_2") %in% names(data_obs))) {
    add("RMSE(CWRES)",
        sqrt(mean(data_obs$CWRES_1^2)),
        sqrt(mean(data_obs$CWRES_2^2)))
    add("Mean |CWRES|",
        mean(abs(data_obs$CWRES_1)),
        mean(abs(data_obs$CWRES_2)))
    add("Prop |CWRES|>2",
        mean(abs(data_obs$CWRES_1) > 2),
        mean(abs(data_obs$CWRES_2) > 2))
  }
  if (all(c("NPDE_1", "NPDE_2") %in% names(data_obs))) {
    add("Mean NPDE",       mean(data_obs$NPDE_1),       mean(data_obs$NPDE_2))
    add("SD NPDE",         stats::sd(data_obs$NPDE_1),  stats::sd(data_obs$NPDE_2))
    add("Prop |NPDE|>1.96",
        mean(abs(data_obs$NPDE_1) > 1.96),
        mean(abs(data_obs$NPDE_2) > 1.96))
  }

  if (all(c("CWRES_1", "CWRES_2") %in% names(data_obs))) {
    wt <- suppressWarnings(stats::wilcox.test(
      abs(data_obs$CWRES_1), abs(data_obs$CWRES_2), paired = TRUE
    ))
    attr(m, "wilcox_absCWRES_p") <- wt$p.value
  }
  m
}

#' Pretty-format a numeric (used in print)
#' @keywords internal
#' @noRd
.nm_cmp_fmt <- function(x) {
  if (!is.finite(x)) "NA" else sprintf("%.3f", x)
}
