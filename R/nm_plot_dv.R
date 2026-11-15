# =============================================================================
# nm_plot_dv.R
# -----------------------------------------------------------------------------
# Build a ggplot of DV vs TIME (or DV/dose vs TIME) from a data.frame already
# parsed by nm_read_csv().
#
# Expected columns: TIME, DV (or DV_norm), ID, and optionally MDV, EVID, CMT,
# BLQ. The `time_unit` and `lloq` attributes are used when present.
# =============================================================================

#' Plot DV vs TIME with median and interval band
#'
#' @param data         A data.frame from [nm_read_csv()].
#' @param y            Y variable: `"DV"` (default) or `"DV_norm"`.
#' @param time_var     X variable: `"TIME"` (default) or `"TAD"` (time after
#'                     dose, computed by [nm_read_csv()] when EVID is present).
#' @param strat        Name of a stratification column (e.g. `"DOSE"`).
#' @param ci_probs     Quantiles of the interval band (default 90% interval,
#'                     `c(0.05, 0.95)`).
#' @param cmt          CMT value(s) to keep (`NULL` = all).
#' @param xlim         X-axis limits `c(min, max)` (`NULL` = full range).
#' @param scales       `"linear"`, `"log"` or `"both"` (side-by-side panels).
#' @param blq_handling BLQ handling: `"show"` (grey crosses, default),
#'                     `"exclude"` or `"keep"`.
#' @param verbose      Print diagnostics (default `TRUE`).
#'
#' @return A ggplot (or a patchwork object if `scales = "both"`).
#'
#' @examples
#' \dontrun{
#' parsed <- nm_read_csv("data.csv", lloq = 0.5)
#' p <- nm_plot_dv(parsed$data, strat = "DOSE", cmt = 2, scales = "both")
#' print(p)
#' }
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate group_by summarise across all_of
#' @importFrom rlang sym .data
#' @importFrom stats median quantile
#' @export
nm_plot_dv <- function(data,
                       y            = c("DV", "DV_norm"),
                       time_var     = c("TIME", "TAD"),
                       strat        = NULL,
                       ci_probs     = c(0.05, 0.95),
                       cmt          = NULL,
                       xlim         = NULL,
                       scales       = c("linear", "log", "both"),
                       blq_handling = c("show", "exclude", "keep"),
                       verbose      = TRUE) {

  y            <- match.arg(y)
  time_var     <- match.arg(time_var)
  scales       <- match.arg(scales)
  blq_handling <- match.arg(blq_handling)
  stopifnot(length(ci_probs) == 2,
            all(ci_probs > 0 & ci_probs < 1),
            ci_probs[1] < ci_probs[2])
  if (!is.null(xlim))
    stopifnot(length(xlim) == 2, is.numeric(xlim), xlim[1] < xlim[2])

  # ---------------------------------------------------------------------------
  # Attributes
  # ---------------------------------------------------------------------------
  time_unit <- attr(data, "time_unit"); if (is.null(time_unit)) time_unit <- "h"
  lloq      <- attr(data, "lloq")

  if (!y %in% names(data))
    stop(sprintf("Column '%s' not found in the data.frame.", y))
  if (!time_var %in% names(data))
    stop(sprintf("Column '%s' not found in the data.frame.", time_var))
  if (!"ID" %in% names(data))
    stop("Column ID not found.")

  # ---------------------------------------------------------------------------
  # Keep observation records
  # ---------------------------------------------------------------------------
  obs <- data
  if ("MDV" %in% names(obs))
    obs <- obs %>% filter(MDV == 0)
  else if ("EVID" %in% names(obs))
    obs <- obs %>% filter(EVID == 0)

  # CMT filter
  cmt_label <- NULL
  if ("CMT" %in% names(obs)) {
    obs$CMT <- suppressWarnings(as.numeric(obs$CMT))
    avail <- sort(unique(obs$CMT[!is.na(obs$CMT) & obs$CMT != 0]))
    if (verbose)
      message(sprintf("Available CMT: %s", paste(avail, collapse = ", ")))
    if (!is.null(cmt)) {
      bad <- setdiff(cmt, avail)
      if (length(bad) > 0L)
        warning(sprintf("CMT not found: %s", paste(bad, collapse = ", ")))
      obs <- obs %>% filter(CMT %in% cmt)
      cmt_label <- paste0("CMT = ", paste(cmt, collapse = ","))
    }
  } else if (!is.null(cmt)) {
    warning("Argument `cmt` supplied but no CMT column: ignored.")
  }

  # Stratification
  if (!is.null(strat)) {
    if (!strat %in% names(obs)) {
      warning(sprintf("Column '%s' not found: stratification ignored.", strat))
      strat <- NULL
    } else {
      obs[[strat]] <- as.factor(obs[[strat]])
    }
  }

  # Apply blq_handling
  if (!"BLQ" %in% names(obs)) obs$BLQ <- FALSE
  if (blq_handling == "exclude") obs <- obs %>% filter(!BLQ)

  # ---------------------------------------------------------------------------
  # Plot builder
  # ---------------------------------------------------------------------------
  ylab <- if (y == "DV") "DV" else "DV / dose"

  build <- function(log_y) {
    pdf <- obs %>%
      mutate(.y = suppressWarnings(as.numeric(.data[[y]]))) %>%
      filter(!is.na(.y), !is.na(.data[[time_var]]))
    if (log_y) pdf <- pdf %>% filter(.y > 0 | BLQ)
    if (nrow(pdf) == 0L) return(NULL)

    show_blq <- blq_handling == "show" && any(pdf$BLQ)
    if (show_blq) {
      stat_df <- pdf %>% filter(!BLQ)
      blq_df  <- pdf %>% filter(BLQ)
    } else {
      stat_df <- pdf
      blq_df  <- NULL
    }

    grp <- c(time_var, if (!is.null(strat)) strat)
    summ <- stat_df %>%
      group_by(across(all_of(grp))) %>%
      summarise(med = median(.y, na.rm = TRUE),
                lo  = quantile(.y, probs = ci_probs[1], na.rm = TRUE),
                hi  = quantile(.y, probs = ci_probs[2], na.rm = TRUE),
                n   = sum(!is.na(.y)),
                .groups = "drop") %>%
      filter(n >= 2)

    if (is.null(strat)) {
      g <- ggplot2::ggplot(stat_df, ggplot2::aes(x = .data[[time_var]], y = .y)) +
        ggplot2::geom_line(ggplot2::aes(group = ID), alpha = 0.20, colour = "grey50") +
        ggplot2::geom_point(alpha = 0.35, colour = "grey30", size = 1) +
        ggplot2::geom_ribbon(data = summ, ggplot2::aes(x = .data[[time_var]], ymin = lo, ymax = hi),
                    inherit.aes = FALSE, fill = "firebrick", alpha = 0.20) +
        ggplot2::geom_line(data = summ, ggplot2::aes(x = .data[[time_var]], y = med),
                  inherit.aes = FALSE, colour = "firebrick", linewidth = 1.1)
    } else {
      ssym <- sym(strat)
      g <- ggplot2::ggplot(stat_df, ggplot2::aes(x = .data[[time_var]], y = .y,
                               colour = !!ssym, fill = !!ssym)) +
        ggplot2::geom_line(ggplot2::aes(group = ID), alpha = 0.20) +
        ggplot2::geom_point(alpha = 0.35, size = 1) +
        ggplot2::geom_ribbon(data = summ,
                    ggplot2::aes(x = .data[[time_var]], ymin = lo, ymax = hi, fill = !!ssym),
                    inherit.aes = FALSE, alpha = 0.20, colour = NA) +
        ggplot2::geom_line(data = summ,
                  ggplot2::aes(x = .data[[time_var]], y = med, colour = !!ssym),
                  inherit.aes = FALSE, linewidth = 1.1)
    }

    # BLQ crosses
    if (show_blq && nrow(blq_df) > 0L) {
      blq_df$.y_plot <- if (log_y && !is.null(lloq)) lloq / 2 else blq_df$.y
      g <- g + ggplot2::geom_point(data = blq_df,
                          ggplot2::aes(x = TIME, y = .y_plot),
                          inherit.aes = FALSE,
                          shape = 4, colour = "grey40",
                          size = 1.5, alpha = 0.6)
    }

    # LLOQ line
    if (!is.null(lloq) && y == "DV")
      g <- g + ggplot2::geom_hline(yintercept = lloq, linetype = "dashed",
                          colour = "grey40", linewidth = 0.4)

    ci_lbl <- sprintf("%d-%d%% interval",
                      round(ci_probs[1] * 100), round(ci_probs[2] * 100))
    sub_parts <- c(paste("Median +", ci_lbl),
                   if (!is.null(strat)) paste("stratified by", strat),
                   if (!is.null(cmt_label)) cmt_label,
                   if (show_blq) sprintf("BLQ: %d crosses", nrow(blq_df)))

    g <- g +
      ggplot2::labs(x = paste0("TIME (", time_unit, ")"),
           y = ylab,
           subtitle = paste(sub_parts, collapse = " | ")) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(legend.position = "bottom")

    if (!is.null(xlim)) g <- g + ggplot2::coord_cartesian(xlim = xlim)
    if (log_y) g <- g + ggplot2::scale_y_log10()
    g
  }

  # ---------------------------------------------------------------------------
  # Final assembly according to `scales`
  # ---------------------------------------------------------------------------
  lin <- build(FALSE)
  log <- build(TRUE)
  if (!is.null(lin)) lin <- lin + ggplot2::ggtitle(paste(ylab, "vs", time_var, "(linear)"))
  if (!is.null(log)) log <- log + ggplot2::ggtitle(paste(ylab, "vs", time_var, "(log10)"))

  grf <- switch(scales,
         "linear" = lin,
         "log"    = log,
         "both"   = {
           if (requireNamespace("patchwork", quietly = TRUE)) {
             patchwork::wrap_plots(lin, log, ncol = 2) +
               patchwork::plot_layout(guides = "collect") &
               ggplot2::theme(legend.position = "bottom")
           } else {
             warning("Package 'patchwork' not installed: returning the linear plot only.")
             lin
           }
         })

  grf
}
