# =============================================================================
# nm_vpc.R
# VPC and pc-VPC from NONMEM _orig.fit and _sim.fit files.
#
# Methodology: Bergstrand et al. 2011 (CPT:PSP)
#   - MEDpred_ref computed once from pooled simulated PRED across all
#     replicates per bin, then applied identically to observed and simulated
#     data (Bergstrand 2011, Eq. 1, lb = 0 assumed).
#   - Bins defined on observed TIME; same breaks applied to simulated data.
#   - No mean line plotted (non-standard per Bergstrand 2011).
#
# Four binning methods:
#   1 = Manual breaks (user-supplied vector)
#   2 = Equiprobable quantiles
#   3 = Nominal time points (midpoint-derived breaks)
#   4 = Automatic rules via hist() (Sturges / Scott / FD)
# =============================================================================

#' @importFrom stats quantile median
#' @importFrom rlang .data
#' @importFrom utils read.table
NULL


# =============================================================================
# Internal helpers
# =============================================================================

# Read a NONMEM .fit table. Header lines ("TABLE NO. n", column names,
# repeated for each subproblem unless NOHEADER/ONEHEADER) are recognised and
# left out; the tables written by nm_vpc_ctl() have none (NOHEADER), so no
# line is skipped blindly: the first record is kept even when it is an
# observation. A data line starts with a number.
.nm_vpc_read_fit <- function(fit_file, col_names) {
  if (!file.exists(fit_file)) stop(".fit file not found: ", fit_file)
  L <- readLines(fit_file, warn = FALSE)
  hdr <- grepl("^\\s*TABLE NO", L) | !grepl("^\\s*[-+.0-9]", L)
  read <- function(...) utils::read.table(..., header = FALSE,
                                          col.names = col_names, na.strings = ".",
                                          stringsAsFactors = FALSE)
  if (!any(hdr)) read(fit_file) else read(text = L[!hdr])
}

# Compute bin breaks from observed TIME using the chosen method
.nm_vpc_make_breaks <- function(time_obs, method, manual = NULL,
                         n_quantile = 8, rule = "FD") {
  switch(as.character(method),
    "1" = {
      if (is.null(manual)) stop("bin_method=1 requires TTIME_seq_manual.")
      manual
    },
    "2" = {
      unique(quantile(time_obs, probs = seq(0, 1, length.out = n_quantile + 1),
                      na.rm = TRUE))
    },
    "3" = {
      nom <- sort(unique(time_obs))
      if (length(nom) < 2) stop("Not enough unique TIME values for nominal binning.")
      mids <- c(nom[1] - diff(nom[1:2]) / 2,
                (nom[-length(nom)] + nom[-1]) / 2,
                nom[length(nom)] + diff(nom[c(length(nom)-1, length(nom))]) / 2)
      mids
    },
    "4" = {
      h <- graphics::hist(time_obs, breaks = rule, plot = FALSE)
      h$breaks
    },
    stop("bin_method must be 1, 2, 3, or 4.")
  )
}

# Assign bin labels (median TIME per bin) and bin midpoints
.nm_vpc_assign_bins <- function(df, breaks, time_col = "TIME") {
  df$TTIME <- cut(df[[time_col]], breaks = breaks,
                  include.lowest = TRUE, labels = FALSE)
  # Bin midpoint = median observed TIME per bin. tapply() returns only the
  # non-empty bins: look the midpoint up by bin name, not by position (with
  # an empty bin, positions shift and the last bins get wrong times or NA).
  bin_mids <- tapply(df[[time_col]], df$TTIME, median, na.rm = TRUE)
  df$MED_T <- as.vector(bin_mids[as.character(df$TTIME)])
  df
}


# =============================================================================
# Main exported function
# =============================================================================

#' Visual Predictive Check (VPC or pc-VPC) from NONMEM simulation output
#'
#' Reads the \code{_orig.fit} and \code{_sim.fit} files produced by
#' \code{\link{nm_vpc_ctl}}, computes observed and simulated percentiles in
#' time bins, and returns a \code{ggplot2} VPC or pc-VPC plot.
#'
#' @param run Path of the estimation run, as given to
#'   \code{\link{nm_vpc_ctl}}: its control stream
#'   (\code{"C:/project/run001.ctl"}), its run directory or its stem
#'   (\code{"C:/project/run001"}). The function looks for
#'   \code{run001_orig.fit} and \code{run001_sim.fit}, looked for in
#'   \code{run001_orig.reg/} (resp. \code{_sim.reg/}), \code{.nm7/}, a folder
#'   of the same name, then directly in the folder of the run.
#' @param analysis_dir Deprecated (former form \code{run_name} +
#'   \code{analysis_dir}): when given, \code{run} is read as the run name in
#'   this directory, with a warning.
#' @param col_names Character vector of column names matching the \code{$TABLE}
#'   output columns in the \code{.fit} files
#'   (e.g. \code{c("NSIM","ID","TIME","AMT","MDV","DV","PRED","IPRED")}).
#' @param type Integer. \code{1} = standard VPC (default), \code{2} = pc-VPC
#'   (prediction-corrected, Bergstrand 2011).
#' @param pi_low Lower bound of the prediction interval (default: 0.05).
#' @param pi_up Upper bound of the prediction interval (default: 0.95).
#' @param ci_low Lower bound of the confidence interval around simulated
#'   percentiles (default: 0.025).
#' @param ci_up Upper bound of the confidence interval around simulated
#'   percentiles (default: 0.975).
#' @param bin_method Integer. Binning method: \code{1} = manual,
#'   \code{2} = equiprobable quantiles, \code{3} = nominal time points,
#'   \code{4} = automatic (default: \code{4}).
#' @param TTIME_seq_manual Numeric vector of bin breaks. Required when
#'   \code{bin_method = 1}.
#' @param n_bins_quantile Number of bins for quantile method
#'   (\code{bin_method = 2}). Default: 8.
#' @param bin_rule Rule for automatic binning (\code{bin_method = 4}):
#'   \code{"Sturges"}, \code{"Scott"}, or \code{"FD"} (default).
#' @param plot_title Character. Title of the VPC plot. If \code{NULL}
#'   (default), a title is generated from the run name and \code{type}.
#' @param plot_file Character. If provided, the plot is saved to this path
#'   (e.g. \code{"vpc_run001.jpeg"}) through \code{\link{nm_export_plot}}
#'   (\code{.jpg}, \code{.png}, \code{.pdf} or \code{.tiff}).
#' @param show_hist Logical. If \code{TRUE} (default), prints a diagnostic
#'   histogram of observed TIME with bin boundaries before the VPC plot.
#' @param ... Former argument name, accepted with a warning until version
#'   1.0: \code{run_name} (now \code{run}).
#'
#' @return A \code{ggplot2} object, displayed when not assigned.
#'
#' @details
#' \strong{pc-VPC correction (Bergstrand 2011, Eq. 1, lb = 0):}
#' \deqn{DV_{pc} = DV \times \tilde{PRED}_{bin} / PRED}
#' where \eqn{\tilde{PRED}_{bin}} (\code{MEDpred_ref}) is the median of all
#' simulated PRED values pooled across replicates within each bin. The same
#' reference is applied to both observed and simulated data.
#'
#' \strong{Required packages:} \code{ggplot2} and \code{dplyr} must be
#' installed (listed in \code{Imports}).
#'
#' @examples
#' \dontrun{
#' col_names <- c("NSIM", "ID", "TIME", "AMT", "MDV", "DV", "PRED", "IPRED")
#'
#' # Standard VPC
#' p <- nm_vpc(
#'   run          = "C:/analyses/my_pk_project/run001.ctl",
#'   col_names    = col_names,
#'   type         = 1,
#'   bin_method   = 4
#' )
#'
#' # pc-VPC with manual bins, saved to file
#' p <- nm_vpc(
#'   run              = "C:/analyses/my_pk_project/run001.ctl",
#'   col_names        = col_names,
#'   type             = 2,
#'   bin_method       = 1,
#'   TTIME_seq_manual = c(0, 1, 2, 4, 8, 12, 24),
#'   plot_file        = "pcvpc_run001.jpeg"
#' )
#' }
#'
#' @seealso \code{\link{nm_vpc_ctl}} for generating the simulation runs.
#'
#' @export
nm_vpc <- function(
    run,
    analysis_dir     = NULL,
    col_names,
    type             = 1,
    pi_low           = 0.05,
    pi_up            = 0.95,
    ci_low           = 0.025,
    ci_up            = 0.975,
    bin_method       = 4,
    TTIME_seq_manual = NULL,
    n_bins_quantile  = 8,
    bin_rule         = "FD",
    plot_title       = NULL,
    plot_file        = NULL,
    show_hist        = TRUE,
    ...
) {
  dep <- .nm_deprecated_args(list(...), c(run_name = "run"), "nm_vpc")
  if (!is.null(dep$run)) run <- dep$run
  if (!is.null(analysis_dir)) {
    warning("In nm_vpc(), `analysis_dir` is deprecated: give the path of the run, ",
            "e.g. run = \"", file.path(analysis_dir, run), ".ctl\".", call. = FALSE)
    run_name <- run
  } else {
    ra <- .nm_run_arg(run)
    run_name <- ra$run_name; analysis_dir <- ra$analysis_dir
  }

  # --- Checks -----------------------------------------------------------------
  if (!type %in% c(1, 2)) stop("type must be 1 (VPC) or 2 (pc-VPC).")
  if (!bin_method %in% 1:4) stop("bin_method must be 1, 2, 3, or 4.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Package 'dplyr' is required. Install it with install.packages('dplyr').")

  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)

  # --- Locate .fit files ------------------------------------------------------
  orig_name <- paste0(run_name, "_orig")
  sim_name  <- paste0(run_name, "_sim")
  orig_fit <- file.path(.nm_find_run_dir(analysis_dir, orig_name, ".fit"),
                        paste0(orig_name, ".fit"))
  sim_fit  <- file.path(.nm_find_run_dir(analysis_dir, sim_name, ".fit"),
                        paste0(sim_name, ".fit"))

  if (!file.exists(orig_fit)) stop("_orig.fit not found: ", orig_fit)
  if (!file.exists(sim_fit))  stop("_sim.fit not found: ", sim_fit)

  # --- Read data --------------------------------------------------------------
  message("Reading _orig.fit...")
  orig <- .nm_vpc_read_fit(orig_fit, col_names)
  message("Reading _sim.fit...")
  sim  <- .nm_vpc_read_fit(sim_fit,  col_names)

  # Filter observed rows (MDV == 0, NSIM == 1 for orig)
  obs <- orig[orig$MDV == 0 & orig$NSIM == 1, ]
  dd  <- sim[sim$MDV == 0, ]

  # --- Binning ----------------------------------------------------------------
  breaks <- .nm_vpc_make_breaks(obs$TIME, bin_method, TTIME_seq_manual,
                         n_bins_quantile, bin_rule)

  # Diagnostic histogram
  if (show_hist) {
    graphics::hist(obs$TIME, breaks = 30,
                   main = "Observed TIME distribution with bin boundaries",
                   xlab = "TIME", col = "lightgrey", border = "white")
    graphics::abline(v = breaks, col = "red", lty = 2)
  }

  # Assign bins
  obs_t <- .nm_vpc_assign_bins(obs, breaks)
  dd_2  <- .nm_vpc_assign_bins(dd,  breaks)

  # Bin midpoints for x-axis
  obs_tt <- unique(obs_t[, c("TTIME", "MED_T")])

  # --- Observed statistics ----------------------------------------------------
  if (type == 1) {
    # Standard VPC
    obs2 <- do.call(rbind, lapply(split(obs_t, obs_t$TTIME), function(b) {
      data.frame(
        TTIME = b$TTIME[1],
        Q025  = quantile(b$DV, pi_low,  na.rm = TRUE),
        MED   = median(b$DV,             na.rm = TRUE),
        Q975  = quantile(b$DV, pi_up,   na.rm = TRUE)
      )
    }))

  } else {
    # pc-VPC: MEDpred_ref from pooled simulated PRED (Bergstrand 2011, Eq. 1)
    ref_MEDpred <- do.call(rbind, lapply(split(dd_2, dd_2$TTIME), function(b) {
      data.frame(TTIME = b$TTIME[1], MEDpred_ref = median(b$PRED, na.rm = TRUE))
    }))

    obs_m <- merge(obs_t, ref_MEDpred, by = "TTIME", all.x = TRUE)
    obs_m$pcDV <- obs_m$DV * obs_m$MEDpred_ref / obs_m$PRED

    obs2 <- do.call(rbind, lapply(split(obs_m, obs_m$TTIME), function(b) {
      data.frame(
        TTIME = b$TTIME[1],
        Q025  = quantile(b$pcDV, pi_low, na.rm = TRUE),
        MED   = median(b$pcDV,            na.rm = TRUE),
        Q975  = quantile(b$pcDV, pi_up,  na.rm = TRUE)
      )
    }))
  }

  obs3 <- merge(obs2, obs_tt, by = "TTIME", all.x = TRUE)

  # --- Simulated statistics ---------------------------------------------------
  if (type == 1) {
    dd_4 <- do.call(rbind, lapply(
      split(dd_2, list(dd_2$NSIM, dd_2$TTIME)),
      function(b) {
        if (nrow(b) == 0) return(NULL)
        data.frame(
          NSIM   = b$NSIM[1],
          TTIME  = b$TTIME[1],
          Q_025  = quantile(b$DV, pi_low, na.rm = TRUE),
          MED_DV = median(b$DV,            na.rm = TRUE),
          Q_975  = quantile(b$DV, pi_up,  na.rm = TRUE)
        )
      }
    ))

  } else {
    dd_2m <- merge(dd_2, ref_MEDpred, by = "TTIME", all.x = TRUE)
    dd_2m$pcDV <- dd_2m$DV * dd_2m$MEDpred_ref / dd_2m$PRED

    dd_4 <- do.call(rbind, lapply(
      split(dd_2m, list(dd_2m$NSIM, dd_2m$TTIME)),
      function(b) {
        if (nrow(b) == 0) return(NULL)
        data.frame(
          NSIM   = b$NSIM[1],
          TTIME  = b$TTIME[1],
          Q_025  = quantile(b$pcDV, pi_low, na.rm = TRUE),
          MED_DV = median(b$pcDV,            na.rm = TRUE),
          Q_975  = quantile(b$pcDV, pi_up,  na.rm = TRUE)
        )
      }
    ))
  }

  # CI around simulated percentiles
  dd_5 <- merge(dd_4, obs_tt, by = "TTIME", all.x = TRUE)
  dd_6 <- do.call(rbind, lapply(split(dd_5, dd_5$TTIME), function(b) {
    data.frame(
      TTIME        = b$TTIME[1],
      MED_T        = b$MED_T[1],
      low_ciQ025   = quantile(b$Q_025,  ci_low, na.rm = TRUE),
      up_ciQ025    = quantile(b$Q_025,  ci_up,  na.rm = TRUE),
      low_ciMED_DV = quantile(b$MED_DV, ci_low, na.rm = TRUE),
      up_ciMED_DV  = quantile(b$MED_DV, ci_up,  na.rm = TRUE),
      low_ciQ975   = quantile(b$Q_975,  ci_low, na.rm = TRUE),
      up_ciQ975    = quantile(b$Q_975,  ci_up,  na.rm = TRUE)
    )
  }))

  # --- Plot -------------------------------------------------------------------
  y_label <- if (type == 1) "Concentration" else "Prediction-corrected concentration"
  if (is.null(plot_title)) {
    plot_title <- if (type == 1) {
      paste("VPC -", run_name)
    } else {
      paste("pc-VPC -", run_name)
    }
  }

  # Use DV or pcDV for observed points
  obs_pts <- if (type == 1) obs_t else obs_m

  grf <- ggplot2::ggplot() +
    # Observed points
    ggplot2::geom_point(
      data = obs_pts,
      mapping = ggplot2::aes(
        x = .data[["TIME"]],
        y = if (type == 1) .data[["DV"]] else .data[["pcDV"]]
      ),
      alpha = 0.2, size = 0.8
    ) +
    # CI ribbons (simulated)
    ggplot2::geom_ribbon(
      data = dd_6,
      mapping = ggplot2::aes(x = .data[["MED_T"]],
                             ymin = .data[["low_ciQ025"]],
                             ymax = .data[["up_ciQ025"]]),
      fill = "blue", alpha = 0.2
    ) +
    ggplot2::geom_ribbon(
      data = dd_6,
      mapping = ggplot2::aes(x = .data[["MED_T"]],
                             ymin = .data[["low_ciMED_DV"]],
                             ymax = .data[["up_ciMED_DV"]]),
      fill = "red", alpha = 0.2
    ) +
    ggplot2::geom_ribbon(
      data = dd_6,
      mapping = ggplot2::aes(x = .data[["MED_T"]],
                             ymin = .data[["low_ciQ975"]],
                             ymax = .data[["up_ciQ975"]]),
      fill = "blue", alpha = 0.2
    ) +
    # Observed percentile lines
    ggplot2::geom_line(
      data = obs3,
      mapping = ggplot2::aes(x = .data[["MED_T"]], y = .data[["Q025"]]),
      color = "blue", linewidth = 1.0
    ) +
    ggplot2::geom_line(
      data = obs3,
      mapping = ggplot2::aes(x = .data[["MED_T"]], y = .data[["MED"]]),
      color = "red", linewidth = 1.0
    ) +
    ggplot2::geom_line(
      data = obs3,
      mapping = ggplot2::aes(x = .data[["MED_T"]], y = .data[["Q975"]]),
      color = "blue", linewidth = 1.0
    ) +
    ggplot2::labs(
      title = plot_title,
      x     = "Time",
      y     = y_label,
      caption = paste0(
        "PI: [", pi_low * 100, "%, ", pi_up * 100, "%]  ",
        "CI: [", ci_low * 100, "%, ", ci_up * 100, "%]"
      )
    ) +
    ggplot2::theme_bw()

  if (!is.null(plot_file)) {
    nm_export_plot(grf, plot_file, width = 8, height = 6, dpi = 300)
  }

  grf
}
