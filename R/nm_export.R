# =============================================================================
# nm_export.R
# -----------------------------------------------------------------------------
# Export and diagnostic helpers for NONMEM-like datasets and plots.
#
# Contains:
#   - nm_export_plot() : export a ggplot/patchwork object to an image file.
#   - nm_export_data() : export a data.frame to CSV using the NONMEM
#                        convention for missing values (na = ".").
#   - nm_summary()     : one-row diagnostic data.frame (N subjects, % BLQ,
#                        time range...).
# =============================================================================

#' Export a ggplot to an image file
#'
#' Single export entry point for all nmtools plots: plot methods return the
#' ggplot object and leave saving to this function.
#'
#' @param plot    A ggplot (or patchwork) object.
#' @param file    Output path. The file type is taken from the extension:
#'   `.jpg`/`.jpeg`, `.png`, `.pdf` or `.tif`/`.tiff`. The directory is
#'   created if needed.
#' @param width   Width in inches (default 9).
#' @param height  Height in inches (default 5).
#' @param dpi     Resolution in dots per inch (default 300; ignored for PDF).
#' @param verbose Print the path of the written file (default `TRUE`).
#'
#' @return The path of the written file (invisibly), or `NULL` if `plot` is
#'   `NULL`.
#' @export
#' @examples
#' \dontrun{
#' p <- nm_plot_dv(nm_read_csv("data.csv")$data)
#' nm_export_plot(p, "figures/dv_vs_time.png")
#' }
nm_export_plot <- function(plot, file,
                           width = 9, height = 5, dpi = 300,
                           verbose = TRUE) {
  if (is.null(plot)) {
    warning("`plot` is NULL: nothing to export.")
    return(invisible(NULL))
  }
  if (!inherits(plot, c("gg", "ggplot", "patchwork")))
    stop("`plot` must be a ggplot or patchwork object.", call. = FALSE)

  ext <- tolower(tools::file_ext(file))
  device <- switch(ext,
                   jpg  = , jpeg = "jpeg",
                   png  = "png",
                   pdf  = "pdf",
                   tif  = , tiff = "tiff",
                   stop("Unsupported file extension '", ext, "' in `file`. ",
                        "Use .jpg, .png, .pdf or .tiff.", call. = FALSE))

  out_dir <- dirname(file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  ggplot2::ggsave(filename = file, plot = plot,
                  width = width, height = height,
                  dpi = dpi, device = device)
  if (verbose) message("  Plot written: ", file)
  invisible(file)
}

#' Export a data.frame to CSV (NONMEM convention: na = ".")
#'
#' @param data A data.frame, typically from [nm_read_csv()].
#' @param file Output path (.csv). The directory is created if needed.
#'
#' @return The path of the written file (invisibly).
#' @export
nm_export_data <- function(data, file) {
  out_dir <- dirname(file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  utils::write.csv(data, file, row.names = FALSE, na = ".")
  message("  CSV written: ", file)
  invisible(file)
}

#' Diagnostic summary of a NONMEM-like dataset
#'
#' @param data    A data.frame from [nm_read_csv()].
#' @param verbose Print the summary to the console (default `TRUE`).
#'
#' @return A one-row data.frame with: `n_subjects`, `n_records`,
#'   `n_observations`, `n_blq`, `pct_blq`, `time_min`, `time_max`,
#'   `time_unit`, `doses_uniques` (distinct positive doses), `n_cmt_obs`.
#' @export
nm_summary <- function(data, verbose = TRUE) {

  mapping <- attr(data, "mapping"); if (is.null(mapping)) mapping <- list()
  time_unit <- attr(data, "time_unit"); if (is.null(time_unit)) time_unit <- "h"

  dv_num <- suppressWarnings(as.numeric(data$DV))
  amt_num <- if ("AMT" %in% names(data)) suppressWarnings(as.numeric(data$AMT)) else NA
  doses <- if (!all(is.na(amt_num)))
    sort(unique(amt_num[!is.na(amt_num) & amt_num > 0])) else NA

  # Number of observations
  n_obs <- if ("MDV" %in% names(data)) sum(data$MDV == 0, na.rm = TRUE)
           else if ("EVID" %in% names(data)) sum(data$EVID == 0, na.rm = TRUE)
           else sum(!is.na(dv_num))

  # BLQ: counted on observation records only
  blq_flag <- if ("BLQ" %in% names(data)) data$BLQ else rep(FALSE, nrow(data))
  obs_mask <- if ("MDV" %in% names(data)) data$MDV == 0
              else if ("EVID" %in% names(data)) data$EVID == 0
              else rep(TRUE, nrow(data))
  n_blq <- sum(blq_flag & obs_mask, na.rm = TRUE)
  pct_blq <- if (n_obs > 0) round(100 * n_blq / n_obs, 1) else 0

  # Observation compartments
  n_cmt <- NA_integer_
  if ("CMT" %in% names(data)) {
    cmt_obs <- suppressWarnings(as.numeric(data$CMT[obs_mask]))
    cmt_obs <- cmt_obs[!is.na(cmt_obs) & cmt_obs != 0]
    n_cmt <- length(unique(cmt_obs))
  }

  res <- data.frame(
    n_subjects     = length(unique(data$ID)),
    n_records      = nrow(data),
    n_observations = n_obs,
    n_blq          = n_blq,
    pct_blq        = pct_blq,
    time_min       = min(data$TIME, na.rm = TRUE),
    time_max       = max(data$TIME, na.rm = TRUE),
    time_unit      = time_unit,
    doses_uniques  = if (all(is.na(doses))) NA_character_
                     else paste(doses, collapse = ", "),
    n_cmt_obs      = n_cmt,
    stringsAsFactors = FALSE)

  if (verbose) {
    message("=== Diagnostic summary ===")
    message(sprintf("  Subjects             : %d", res$n_subjects))
    message(sprintf("  Records              : %d", res$n_records))
    message(sprintf("  Observations         : %d", res$n_observations))
    if (res$n_blq > 0)
      message(sprintf("  BLQ                  : %d (%.1f%%)",
                      res$n_blq, res$pct_blq))
    message(sprintf("  TIME range           : %.2f - %.2f %s",
                    res$time_min, res$time_max, time_unit))
    if (!is.na(res$doses_uniques))
      message(sprintf("  Distinct doses       : %s", res$doses_uniques))
    if (!is.na(res$n_cmt_obs))
      message(sprintf("  Observation CMTs     : %d", res$n_cmt_obs))
  }

  invisible(res)
}
