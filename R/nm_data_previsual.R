# =============================================================================
# nm_data_previsual.R
# -----------------------------------------------------------------------------
# This file contains two exported functions:
#   - nm_data_previsual() : read and pre-visualise a NONMEM dataset.
#   - nm_add_rate()       : add a RATE column (estimated infusion duration or
#                           rate) without altering the rest of the dataset.
# -----------------------------------------------------------------------------
# Read and pre-visualise a NONMEM-like database in CSV format.
#
# Features (nm_data_previsual):
#   - Identify NONMEM reserved columns (ID, DV, TIME, DATE, MDV, EVID, AMT,
#     RATE, II, SS, ADDL, CMT, DAT1, DAT2, DAT3...).
#   - By default, the TIME column is used AS IS (numeric, already relative).
#     A relative time is rebuilt from a clock-time column (hh:mm[:ss]) and/or
#     a DATE column ONLY on explicit request (reconstruct_time = TRUE).
#     Detecting a DATE/DAT1/DAT2/DAT3 column NEVER triggers the conversion.
#   - ggplot of DV vs TIME: individual profiles, median, interval band
#     (configurable quantiles), optional stratification.
#
# NONMEM conventions (Guide IV - III.B.2 and Guide VIII - $INPUT):
#   DATE : month-day-year   (mm-dd-yy)
#   DAT1 : day-month-year   (dd-mm-yy)
#   DAT2 : year-month-day   (yy-mm-dd)
#   DAT3 : year-day-month   (yy-dd-mm)
#   A DATE column may also hold a plain integer day counter (1, 2, 3...):
#   in that case use date_type = "continuous".
# =============================================================================

# -----------------------------------------------------------------------------
# nm_data_previsual()
# -----------------------------------------------------------------------------
#' Pre-visualise a NONMEM-like dataset
#'
#' @param file        Path to the CSV file (or an already loaded data.frame).
#' @param sep         Column separator (default `","`).
#' @param na          Strings to be treated as missing values.
#' @param time_unit   Output unit for TIME: `"h"` (default) or `"d"`.
#' @param reconstruct_time If `TRUE`, rebuild a relative TIME (t0 = first dose
#'                    per subject) from a clock-time column and, optionally,
#'                    a DATE column. If `FALSE` (default), the TIME column is
#'                    used as is (`as.numeric`). Detecting a
#'                    DATE/DAT1/DAT2/DAT3 column NEVER triggers the
#'                    reconstruction automatically.
#' @param clock_col   Name of the column holding the clock time (`hh:mm:ss`,
#'                    `hh:mm` or decimal hours), used when
#'                    `reconstruct_time = TRUE`. `NULL` (default) = column
#'                    `"TIME"`.
#' @param date_type   How the DATE column is interpreted during the
#'                    reconstruction (only if `reconstruct_time = TRUE`):
#'                    - `NULL` (default): clock time only, anchored on a single day
#'                    - `"continuous"`: DATE = integer day counter (1, 2, 3...)
#'                    - `"DATE"`: calendar date month-day-year (mm-dd-yy)
#'                    - `"DAT1"`: calendar date day-month-year (dd-mm-yy)
#'                    - `"DAT2"`: calendar date year-month-day (yy-mm-dd)
#'                    - `"DAT3"`: calendar date year-day-month (yy-dd-mm)
#' @param dv_col      Name of the DV column to use (manual override).
#' @param strat       Name of a stratification column for the plot
#'                    (e.g. `"DOSE"`, `"OCC"`, `"STUDY"`). `NULL` = none.
#' @param ci_probs    Two quantiles defining the interval band around the
#'                    median (default `c(0.05, 0.95)`, i.e. a 90% interval).
#' @param dose_norm   If `TRUE`, also compute DV_norm = DV / AMT (dose given
#'                    to the subject) and produce the matching plot.
#'                    Requires an identified AMT column.
#' @param scales      `"linear"`, `"log"` or `"both"`: Y axes produced.
#'                    `"both"` gives side-by-side linear and log panels
#'                    (patchwork).
#' @param cmt         CMT value(s) to keep for the plot. Requires an
#'                    identified CMT column. `NULL` (default) = all observed
#'                    compartments (CMT != 0 when EVID is available).
#'                    E.g. `cmt = 2` for the central compartment.
#' @param xlim        X-axis limits `c(min, max)`. `NULL` (default) = full
#'                    data range.
#' @param blq_handling How BLQ values are handled:
#'                    - `"show"` (default): BLQ drawn as grey crosses below a
#'                      dashed line at LLOQ; excluded from the median and
#'                      quantiles
#'                    - `"exclude"`: BLQ removed from the plot and statistics
#'                    - `"half_lloq"`: DV replaced by LLOQ/2 (requires `lloq`)
#'                    - `"keep"`: nothing done, plotted as is
#' @param lloq        Numeric LLOQ. Needed to identify BLQ records when no
#'                    flag column is found, and for `"half_lloq"`.
#'                    `NULL` = no threshold.
#' @param blq_col     Name of the censoring flag column (1 = BLQ,
#'                    0 = quantified). If `NULL`, auto-detected among `"BLQ"`,
#'                    `"BQL"`, `"LLQ"`, `"CENS"`.
#' @param export_jpg  If `TRUE`, export the plots as JPG files in `output_dir`.
#' @param export_csv  If `TRUE`, export the final data.frame (relative TIME,
#'                    DV_norm) as CSV in `output_dir`.
#' @param output_dir  Output directory (default: directory of the source file).
#'                    Formerly `out_dir`, still accepted with a warning.
#' @param jpg_width   Width in inches of the JPG export (default 9).
#' @param jpg_height  Height in inches of the JPG export (default 5).
#' @param jpg_dpi     Resolution in dpi (default 300).
#' @param verbose     Print the identification report (default `TRUE`).
#' @param ...         Former argument name, accepted with a warning until
#'                    version 1.0: `out_dir` (now `output_dir`).
#'
#' @return An object of class `nm_data_previsual` (a list, printed as a short
#'   summary):
#'   - `data`: data.frame with TIME (relative if rebuilt) and DV_norm
#'   - `mapping`: reserved labels -> original column names
#'   - `plot_dv`: ggplot of DV vs TIME (according to `scales`)
#'   - `plot_norm`: ggplot of DV/dose vs TIME (`NULL` if `dose_norm = FALSE`)
#'   - `summary`: diagnostic data.frame (N subjects, N obs, time range, % BLQ...)
#'   - `files`: paths of the written files (`NULL` if nothing was exported)
#'
#' @examples
#' \dontrun{
#' # TIME already numeric/relative: used as is (default behaviour)
#' res <- nm_data_previsual("data.csv")
#'
#' # Rebuild from a clock-time column TIME_R + DATE = day counter
#' res <- nm_data_previsual("data.csv", reconstruct_time = TRUE,
#'                          clock_col = "TIME_R", date_type = "continuous")
#'
#' # Rebuild from a calendar DATE in day-month-year order (DAT1)
#' res <- nm_data_previsual("data.csv", reconstruct_time = TRUE,
#'                          clock_col = "TIME", date_type = "DAT1")
#'
#' res <- nm_data_previsual("data.csv", strat = "DOSE", dose_norm = TRUE,
#'                          scales = "both", lloq = 0.5, export_jpg = TRUE)
#' print(res$plot_dv); print(res$summary)
#' }
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate group_by summarise across all_of left_join bind_rows select first n_distinct
#' @importFrom lubridate parse_date_time
#' @importFrom rlang sym .data
#' @importFrom utils read.csv write.csv
#' @importFrom stats median quantile
#' @importFrom tools file_path_sans_ext
#' @export
nm_data_previsual <- function(file,
                              sep              = ",",
                              na               = c("", "NA", "."),
                              time_unit        = c("h", "d"),
                              reconstruct_time = FALSE,
                              clock_col        = NULL,
                              date_type        = NULL,
                              dv_col           = NULL,
                              strat            = NULL,
                              ci_probs         = c(0.05, 0.95),
                              dose_norm        = FALSE,
                              scales           = c("linear", "log", "both"),
                              cmt              = NULL,
                              xlim             = NULL,
                              blq_handling     = c("show", "exclude", "half_lloq", "keep"),
                              lloq             = NULL,
                              blq_col          = NULL,
                              export_jpg       = FALSE,
                              export_csv       = FALSE,
                              output_dir       = NULL,
                              jpg_width        = 9,
                              jpg_height       = 5,
                              jpg_dpi          = 300,
                              verbose          = TRUE,
                              ...) {
  dep <- .nm_deprecated_args(list(...), c(out_dir = "output_dir"), "nm_data_previsual")
  if (!is.null(dep$output_dir)) output_dir <- dep$output_dir
  out_dir <- output_dir

  time_unit    <- match.arg(time_unit)
  scales       <- match.arg(scales)
  blq_handling <- match.arg(blq_handling)
  stopifnot(length(ci_probs) == 2,
            all(ci_probs > 0 & ci_probs < 1),
            ci_probs[1] < ci_probs[2])
  if (!is.null(xlim))
    stopifnot(length(xlim) == 2, is.numeric(xlim), xlim[1] < xlim[2])
  if (blq_handling == "half_lloq" && is.null(lloq))
    stop("blq_handling='half_lloq' requires the `lloq` argument.")
  if (!is.null(date_type))
    date_type <- match.arg(date_type,
                           c("continuous", "DATE", "DAT1", "DAT2", "DAT3"))

  # ---------------------------------------------------------------------------
  # 1. Raw read (or data.frame supplied directly)
  # ---------------------------------------------------------------------------
  if (is.data.frame(file)) {
    df <- file
  } else {
    df <- utils::read.csv(file, sep = sep, na.strings = na,
                          stringsAsFactors = FALSE, check.names = FALSE)
    # Clean column names: remove leading '#' (NONMEM comment marker on first col)
    names(df) <- sub("^#+", "", names(df))
  }

  # ---------------------------------------------------------------------------
  # 2. Identify NONMEM reserved words
  # ---------------------------------------------------------------------------
  reserved <- c("ID", "L1", "L2", "DV", "MDV", "TIME", "DATE",
                "DAT1", "DAT2", "DAT3", "EVID", "AMT", "RATE",
                "SS", "II", "ADDL", "CMT", "PCMT")

  cols_upper <- toupper(names(df))
  mapping <- list()
  for (r in reserved) {
    hit <- which(cols_upper == r)
    if (length(hit) == 1L) {
      mapping[[r]] <- names(df)[hit]
      names(df)[hit] <- r  # normalise
    } else if (length(hit) > 1L) {
      warning(sprintf("Several columns match the label '%s': %s",
                      r, paste(names(df)[hit], collapse = ", ")))
    }
  }

  # Manual override for DV column
  if (!is.null(dv_col)) {
    if (!dv_col %in% names(df)) {
      warning(sprintf("dv_col '%s' not found in dataset.", dv_col))
    } else {
      mapping[["DV"]] <- dv_col
      names(df)[names(df) == dv_col] <- "DV"
    }
  }

  # Auto-detect DV if still not found: look for column whose name starts with "DV"
  if (!"DV" %in% names(mapping)) {
    dv_candidates <- names(df)[grepl("^DV", toupper(names(df))) &
                                !grepl("^MDV", toupper(names(df)))]
    if (length(dv_candidates) == 1L) {
      if (verbose) message(sprintf("  DV auto-detected: '%s'", dv_candidates[1]))
      mapping[["DV"]] <- dv_candidates[1]
      names(df)[names(df) == dv_candidates[1]] <- "DV"
    } else if (length(dv_candidates) > 1L) {
      # Interactive selection when multiple DV columns found
      message("Multiple DV columns detected:")
      for (i in seq_along(dv_candidates))
        message(sprintf("  %d: %s", i, dv_candidates[i]))
      choice <- readline(prompt = "Select DV column (enter number): ")
      idx <- suppressWarnings(as.integer(choice))
      if (!is.na(idx) && idx >= 1L && idx <= length(dv_candidates)) {
        chosen <- dv_candidates[idx]
        message(sprintf("  DV set to: '%s'", chosen))
        mapping[["DV"]] <- chosen
        names(df)[names(df) == chosen] <- "DV"
      } else {
        warning("Invalid selection - no DV column assigned. Use dv_col= to specify.")
      }
    }
  }

  if (verbose) {
    message("=== Identified NONMEM columns ===")
    if (length(mapping) == 0L) {
      message("  (none)")
    } else {
      for (k in names(mapping))
        message(sprintf("  %-6s -> %s", k, mapping[[k]]))
    }
  }

  # ---------------------------------------------------------------------------
  # 3. Minimal checks
  # ---------------------------------------------------------------------------
  if (!"TIME" %in% names(mapping))
    stop("No TIME column identified.")
  if (!"ID" %in% names(mapping)) {
    warning("No ID column identified: a single subject is assumed.")
    df$ID <- 1L
    mapping$ID <- "ID"
  }
  if (!"DV" %in% names(mapping))
    warning("No DV column identified: no plot will be produced.")

  # ---------------------------------------------------------------------------
  # 4. TIME handling
  # ---------------------------------------------------------------------------
  # Default: numeric TIME used as is (already relative).
  # Reconstruction ONLY if reconstruct_time = TRUE.
  # DATE/DAT1/DAT2/DAT3 columns are detected but never used automatically.
  clock_source <- if (!is.null(clock_col)) clock_col else "TIME"
  date_label   <- intersect(c("DATE", "DAT1", "DAT2", "DAT3"), names(mapping))
  has_date     <- length(date_label) >= 1L
  if (has_date) date_label <- date_label[1L]

  if (verbose && has_date && !reconstruct_time)
    message(sprintf(
      "Column '%s' detected but IGNORED (reconstruct_time=FALSE).",
      mapping[[date_label]]))

  if (!reconstruct_time) {
    # -- Default case: numeric TIME used as is -------------------------------
    if (!clock_source %in% names(df))
      stop(sprintf("Time column '%s' not found.", clock_source))
    time_chr <- as.character(df[[clock_source]])
    if (any(grepl(":", time_chr)))
      warning(sprintf(
        paste0("Column '%s' seems to contain clock times (hh:mm). ",
               "Set reconstruct_time=TRUE (with clock_col=/date_type=) ",
               "to rebuild a relative time."),
        clock_source))
    df$TIME <- suppressWarnings(as.numeric(time_chr))
    if (verbose)
      message("TIME format: numeric -> used as is (relative time).")

  } else {
    # -- Rebuild a relative time ----------------------------------------------
    if (!clock_source %in% names(df))
      stop(sprintf("reconstruct_time=TRUE: clock-time column '%s' not found.",
                   clock_source))
    if (verbose)
      message(sprintf(
        "Rebuilding relative time from '%s'%s.",
        clock_source,
        if (!is.null(date_type)) sprintf(" + DATE (type '%s')", date_type)
        else " (clock time only)"))

    clock_raw <- df[[clock_source]]

    # Helper: clock time -> decimal hours
    clock_to_hours <- function(x) {
      x <- as.character(x)
      out <- vapply(x, function(v) {
        if (is.na(v) || v == "") return(NA_real_)
        if (grepl(":", v)) {
          p <- suppressWarnings(as.numeric(strsplit(v, ":", fixed = TRUE)[[1]]))
          p <- c(p, rep(0, max(0L, 3L - length(p))))
          return(p[1] + p[2] / 60 + p[3] / 3600)
        }
        suppressWarnings(as.numeric(v))  # already decimal
      }, numeric(1))
      unname(out)
    }

    # Helper: clock time -> "hh:mm:ss" string (for calendar parsing)
    clean_time <- function(x) {
      x <- as.character(x)
      out <- vapply(x, function(v) {
        if (is.na(v) || v == "") return(NA_character_)
        if (grepl(":", v)) {
          n_colon <- lengths(regmatches(v, gregexpr(":", v)))
          if (n_colon == 1L) return(paste0(v, ":00"))
          return(v)
        }
        num <- suppressWarnings(as.numeric(v))
        if (is.na(num)) return(NA_character_)
        hh <- floor(num)
        mm <- floor((num - hh) * 60)
        ss <- round(((((num - hh) * 60) - mm)) * 60)
        sprintf("%02d:%02d:%02d", hh, mm, ss)
      }, character(1))
      unname(out)
    }

    if (is.null(date_type) || date_type == "continuous") {
      # ---- Clock time only, or DATE = integer day counter -------------------
      hrs <- clock_to_hours(clock_raw)
      if (!is.null(date_type) && date_type == "continuous") {
        if (!has_date)
          stop("date_type='continuous' but no DATE column detected.")
        day <- suppressWarnings(as.numeric(df[[date_label]]))
        if (all(is.na(day)))
          stop(sprintf("DATE column '%s' is not numeric: ",
                       mapping[[date_label]]),
               "date_type='continuous' is not possible.")
        df$.ABS <- day * 24 + hrs
      } else {
        if (verbose)
          message("  No DATE used: anchored on a single day ",
                  "(beware of records crossing midnight).")
        df$.ABS <- hrs
      }
    } else {
      # ---- Calendar DATE (DATE / DAT1 / DAT2 / DAT3) ------------------------
      if (!has_date)
        stop(sprintf("date_type='%s' but no DATE column detected.",
                     date_type))
      date_orders <- switch(date_type,
             "DATE" = c("mdy HMS", "mdY HMS"),
             "DAT1" = c("dmy HMS", "dmY HMS"),
             "DAT2" = c("ymd HMS", "Ymd HMS"),
             "DAT3" = c("ydm HMS", "Ydm HMS"))
      dt_str <- paste(df[[date_label]], clean_time(clock_raw))
      DT <- lubridate::parse_date_time(dt_str, orders = date_orders,
                                       quiet = TRUE)
      n_fail <- sum(is.na(DT))
      if (n_fail > 0L && verbose)
        message(sprintf("  Warning: %d DATE/time value(s) could not be parsed.",
                        n_fail))
      origin  <- suppressWarnings(min(DT, na.rm = TRUE))
      df$.ABS <- as.numeric(difftime(DT, origin, units = "hours"))
    }

    # t0 per subject = first dose (EVID 1/4 or AMT > 0), otherwise first record
    has_evid <- "EVID" %in% names(mapping)
    has_amt  <- "AMT"  %in% names(mapping)
    if (has_evid || has_amt) {
      ref <- df %>%
        group_by(ID) %>%
        filter(
          (if (has_evid) EVID %in% c(1, 4) else FALSE) |
          (if (has_amt) {
             .a <- suppressWarnings(as.numeric(AMT)); !is.na(.a) & .a > 0
           } else FALSE)
        ) %>%
        summarise(.t0 = suppressWarnings(min(.ABS, na.rm = TRUE)),
                  .groups = "drop")

      no_dose <- setdiff(unique(df$ID), ref$ID)
      if (length(no_dose) > 0L) {
        if (verbose)
          message(sprintf("  %d subject(s) without dose: t0 = first record.",
                          length(no_dose)))
        fb <- df %>% filter(ID %in% no_dose) %>%
          group_by(ID) %>%
          summarise(.t0 = suppressWarnings(min(.ABS, na.rm = TRUE)),
                    .groups = "drop")
        ref <- bind_rows(ref, fb)
      }
    } else {
      ref <- df %>%
        group_by(ID) %>%
        summarise(.t0 = suppressWarnings(min(.ABS, na.rm = TRUE)),
                  .groups = "drop")
      if (verbose)
        message("  No EVID/AMT: t0 = first record per subject.")
    }

    df <- df %>%
      left_join(ref, by = "ID") %>%
      mutate(TIME = .ABS - .t0) %>%
      select(-.ABS, -.t0)
  }

  # Unit conversion
  if (time_unit == "d") df$TIME <- df$TIME / 24

  # ---------------------------------------------------------------------------
  # 5. DV_norm = DV / administered dose (if requested)
  # ---------------------------------------------------------------------------
  if (dose_norm) {
    if (!"AMT" %in% names(mapping)) {
      warning("dose_norm=TRUE but no AMT column: normalisation skipped.")
      dose_norm <- FALSE
    } else {
      # Reference dose per subject = first dose > 0
      # (simple case: single dose or identical doses. With heterogeneous
      # multiple doses, the first one is used and the user is told.)
      dose_ref <- df %>%
        mutate(AMT = suppressWarnings(as.numeric(AMT))) %>%
        filter(!is.na(AMT), AMT > 0) %>%
        group_by(ID) %>%
        summarise(.dose_ref = first(AMT),
                  .n_doses  = n_distinct(AMT),
                  .groups   = "drop")

      n_multi <- sum(dose_ref$.n_doses > 1, na.rm = TRUE)
      if (n_multi > 0L && verbose)
        message(sprintf(
          "  Note: %d subject(s) with heterogeneous doses; DV_norm uses the first dose.",
          n_multi))

      df <- df %>%
        left_join(dose_ref %>% select(ID, .dose_ref), by = "ID") %>%
        mutate(DV_norm = suppressWarnings(as.numeric(DV)) / .dose_ref) %>%
        select(-.dose_ref)
    }
  }

  # ---------------------------------------------------------------------------
  # 6. Build the DV (and DV_norm) vs TIME ggplots
  # ---------------------------------------------------------------------------
  build_plot <- function(data, yvar, ylab, log_y) {

    plot_df <- data %>%
      mutate(.y = suppressWarnings(as.numeric(.data[[yvar]])),
             .blq = if (".blq" %in% names(data)) .blq else FALSE) %>%
      filter(!is.na(.y), !is.na(TIME))

    # Log scale: drop .y <= 0
    if (log_y) plot_df <- plot_df %>% filter(.y > 0)

    if (nrow(plot_df) == 0L) return(NULL)

    # In "show" mode, BLQ are excluded from the statistics but drawn as crosses
    show_blq <- blq_handling == "show" && any(plot_df$.blq)
    if (show_blq) {
      stat_df  <- plot_df %>% filter(!.blq)
      blq_df   <- plot_df %>% filter(.blq)
    } else {
      stat_df  <- plot_df
      blq_df   <- NULL
    }

    grp_vars <- c("TIME", if (!is.null(strat)) strat)
    summ <- stat_df %>%
      group_by(across(all_of(grp_vars))) %>%
      summarise(
        med = median(.y, na.rm = TRUE),
        lo  = quantile(.y, probs = ci_probs[1], na.rm = TRUE),
        hi  = quantile(.y, probs = ci_probs[2], na.rm = TRUE),
        n   = sum(!is.na(.y)),
        .groups = "drop") %>%
      filter(n >= 2)

    if (is.null(strat)) {
      g <- ggplot2::ggplot(stat_df, ggplot2::aes(x = TIME, y = .y)) +
        ggplot2::geom_line(ggplot2::aes(group = ID), alpha = 0.20, colour = "grey50") +
        ggplot2::geom_point(alpha = 0.35, colour = "grey30", size = 1) +
        ggplot2::geom_ribbon(data = summ,
                    ggplot2::aes(x = TIME, ymin = lo, ymax = hi),
                    inherit.aes = FALSE,
                    fill = "firebrick", alpha = 0.20) +
        ggplot2::geom_line(data = summ,
                  ggplot2::aes(x = TIME, y = med),
                  inherit.aes = FALSE,
                  colour = "firebrick", linewidth = 1.1)
    } else {
      strat_sym <- sym(strat)
      g <- ggplot2::ggplot(stat_df, ggplot2::aes(x = TIME, y = .y, colour = !!strat_sym,
                               fill = !!strat_sym)) +
        ggplot2::geom_line(ggplot2::aes(group = ID), alpha = 0.20) +
        ggplot2::geom_point(alpha = 0.35, size = 1) +
        ggplot2::geom_ribbon(data = summ,
                    ggplot2::aes(x = TIME, ymin = lo, ymax = hi,
                        fill = !!strat_sym),
                    inherit.aes = FALSE, alpha = 0.20, colour = NA) +
        ggplot2::geom_line(data = summ,
                  ggplot2::aes(x = TIME, y = med, colour = !!strat_sym),
                  inherit.aes = FALSE, linewidth = 1.1)
    }

    # Draw BLQ as grey crosses
    if (show_blq && !is.null(blq_df) && nrow(blq_df) > 0L) {
      # On the log scale, BLQ with DV <= 0 were already dropped; they are
      # placed at lloq/2 to keep them visible when lloq is supplied
      if (log_y && !is.null(lloq)) {
        blq_df$.y_plot <- lloq / 2
      } else {
        blq_df$.y_plot <- blq_df$.y
      }
      g <- g + ggplot2::geom_point(data = blq_df,
                          ggplot2::aes(x = TIME, y = .y_plot),
                          inherit.aes = FALSE,
                          shape = 4, colour = "grey40",
                          size = 1.5, alpha = 0.6)
    }

    # LLOQ line (dashed horizontal) if lloq is supplied and yvar == "DV"
    if (!is.null(lloq) && yvar == "DV") {
      g <- g + ggplot2::geom_hline(yintercept = lloq, linetype = "dashed",
                          colour = "grey40", linewidth = 0.4)
    }

    ci_lbl <- sprintf("%d-%d%% interval",
                      round(ci_probs[1] * 100),
                      round(ci_probs[2] * 100))

    sub_parts <- c(paste("Median +", ci_lbl),
                   if (!is.null(strat)) paste("stratified by", strat),
                   if (!is.null(cmt_label)) cmt_label,
                   if (show_blq) sprintf("BLQ: %d grey crosses", nrow(blq_df)))

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

  # Keep observation records first
  obs_df <- df
  if ("MDV" %in% names(mapping)) {
    obs_df <- obs_df %>% filter(MDV == 0)
  } else if ("EVID" %in% names(mapping)) {
    obs_df <- obs_df %>% filter(EVID == 0)
  }

  # CMT filter
  cmt_label <- NULL
  if ("CMT" %in% names(mapping)) {
    obs_df <- obs_df %>%
      mutate(CMT = suppressWarnings(as.numeric(CMT)))
    avail_cmt <- sort(unique(obs_df$CMT[!is.na(obs_df$CMT) & obs_df$CMT != 0]))
    if (verbose)
      message(sprintf("Available observation compartments: %s",
                      paste(avail_cmt, collapse = ", ")))

    if (!is.null(cmt)) {
      bad <- setdiff(cmt, avail_cmt)
      if (length(bad) > 0L)
        warning(sprintf("CMT not present in the observations: %s",
                        paste(bad, collapse = ", ")))
      obs_df <- obs_df %>% filter(CMT %in% cmt)
      cmt_label <- paste0("CMT = ", paste(cmt, collapse = ","))
    } else if (length(avail_cmt) > 1L) {
      if (verbose)
        message("  Several CMT present: all plotted. Use `cmt = ...` to filter.")
    }
  } else if (!is.null(cmt)) {
    warning("Argument `cmt` supplied but no CMT column identified: ignored.")
  }

  # BLQ detection and flag
  # A .blq column (TRUE/FALSE) is added to obs_df
  obs_df$.blq <- FALSE

  # 1. Auto-detect a censoring flag column
  #    Known names: BLQ, BQL, LLQ, CENS (1 = censored/BLQ, 0 = quantified)
  if (is.null(blq_col)) {
    blq_candidates <- intersect(c("BLQ", "BQL", "LLQ", "CENS"),
                                toupper(names(obs_df)))
    if (length(blq_candidates) >= 1L) {
      blq_col <- names(obs_df)[toupper(names(obs_df)) == blq_candidates[1]][1]
      if (verbose)
        message(sprintf("Censoring column auto-detected: '%s'.", blq_col))
    }
  }

  if (!is.null(blq_col) && blq_col %in% names(obs_df)) {
    obs_df$.blq <- suppressWarnings(as.numeric(obs_df[[blq_col]])) == 1
    obs_df$.blq[is.na(obs_df$.blq)] <- FALSE
    if (verbose)
      message(sprintf("BLQ flag from column '%s': %d BLQ record(s).",
                      blq_col, sum(obs_df$.blq)))
  } else if (!is.null(lloq)) {
    dv_num <- suppressWarnings(as.numeric(obs_df$DV))
    obs_df$.blq <- !is.na(dv_num) & dv_num < lloq
    if (verbose)
      message(sprintf("BLQ identified from lloq=%g: %d record(s).",
                      lloq, sum(obs_df$.blq)))
  }

  # Apply blq_handling
  pct_blq <- if (nrow(obs_df) > 0L) mean(obs_df$.blq) * 100 else 0
  if (any(obs_df$.blq)) {
    if (blq_handling == "exclude") {
      obs_df <- obs_df %>% filter(!.blq)
    } else if (blq_handling == "half_lloq") {
      # Keep the reference dose per subject to recompute DV_norm
      if ("DV_norm" %in% names(obs_df)) {
        dose_ref_blq <- obs_df$DV[obs_df$.blq] / obs_df$DV_norm[obs_df$.blq]
        obs_df$DV[obs_df$.blq] <- lloq / 2
        obs_df$DV_norm[obs_df$.blq] <- (lloq / 2) / dose_ref_blq
      } else {
        obs_df$DV[obs_df$.blq] <- lloq / 2
      }
    }
    # "show" and "keep": .blq is kept for use in build_plot
  }

  if (!is.null(strat)) {
    if (!strat %in% names(obs_df)) {
      warning(sprintf("Column '%s' not found: stratification ignored.", strat))
      strat <- NULL
    } else {
      obs_df[[strat]] <- as.factor(obs_df[[strat]])
    }
  }

  # Helper: assemble linear/log panels according to `scales`
  assemble <- function(yvar, ylab) {
    if (!yvar %in% names(obs_df)) return(NULL)
    lin <- build_plot(obs_df, yvar, ylab, log_y = FALSE)
    log <- build_plot(obs_df, yvar, ylab, log_y = TRUE)
    if (!is.null(lin)) lin <- lin + ggplot2::ggtitle(paste(ylab, "vs TIME (linear)"))
    if (!is.null(log)) log <- log + ggplot2::ggtitle(paste(ylab, "vs TIME (log10)"))
    switch(scales,
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
  }

  plot_dv   <- if ("DV" %in% names(mapping)) assemble("DV", "DV") else NULL
  plot_norm <- if (dose_norm) assemble("DV_norm", "DV / dose") else NULL

  # ---------------------------------------------------------------------------
  # 7. Diagnostic summary
  # ---------------------------------------------------------------------------
  dv_num <- suppressWarnings(as.numeric(df$DV))
  amt_num <- if ("AMT" %in% names(mapping)) suppressWarnings(as.numeric(df$AMT)) else NA
  doses_uniq <- if (!all(is.na(amt_num))) sort(unique(amt_num[!is.na(amt_num) & amt_num > 0])) else NA

  summary_df <- data.frame(
    n_subjects     = length(unique(df$ID)),
    n_records      = nrow(df),
    n_observations = if ("MDV" %in% names(mapping)) sum(df$MDV == 0, na.rm = TRUE)
                     else if ("EVID" %in% names(mapping)) sum(df$EVID == 0, na.rm = TRUE)
                     else sum(!is.na(dv_num)),
    n_blq          = sum(obs_df$.blq),
    pct_blq        = round(pct_blq, 1),
    time_min       = min(df$TIME, na.rm = TRUE),
    time_max       = max(df$TIME, na.rm = TRUE),
    time_unit      = time_unit,
    unique_doses  = if (all(is.na(doses_uniq))) NA_character_
                     else paste(doses_uniq, collapse = ", "),
    n_cmt_obs      = if ("CMT" %in% names(mapping))
                       length(unique(obs_df$CMT[!is.na(obs_df$CMT)]))
                     else NA_integer_,
    stringsAsFactors = FALSE)

  if (verbose) {
    message("=== Diagnostic summary ===")
    message(sprintf("  Subjects             : %d", summary_df$n_subjects))
    message(sprintf("  Records              : %d", summary_df$n_records))
    message(sprintf("  Observations         : %d", summary_df$n_observations))
    if (summary_df$n_blq > 0)
      message(sprintf("  BLQ                  : %d (%.1f%%)",
                      summary_df$n_blq, summary_df$pct_blq))
    message(sprintf("  TIME range           : %.2f - %.2f %s",
                    summary_df$time_min, summary_df$time_max, time_unit))
    if (!is.na(summary_df$unique_doses))
      message(sprintf("  Distinct doses       : %s", summary_df$unique_doses))
  }

  # ---------------------------------------------------------------------------
  # 8. JPG / CSV export (optional)
  # ---------------------------------------------------------------------------
  files <- NULL
  if (export_jpg || export_csv) {
    if (is.null(out_dir)) {
      out_dir <- if (is.data.frame(file)) getwd()
                 else dirname(normalizePath(file, mustWork = FALSE))
    }
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    base <- if (is.data.frame(file)) "dataset"
            else tools::file_path_sans_ext(basename(file))
  }

  if (export_jpg) {
    save_jpg <- function(plot_obj, suffix) {
      if (is.null(plot_obj)) return(NULL)
      fp <- file.path(out_dir, sprintf("%s_%s.jpg", base, suffix))
      nm_export_plot(plot_obj, fp, width = jpg_width, height = jpg_height,
                     dpi = jpg_dpi, verbose = verbose)
      fp
    }
    f1 <- save_jpg(plot_dv,   "DV_vs_TIME")
    f2 <- save_jpg(plot_norm, "DVnorm_vs_TIME")
    files <- c(files, f1, f2)
  }

  if (export_csv) {
    fp_csv <- file.path(out_dir, sprintf("%s_processed.csv", base))
    utils::write.csv(df, fp_csv, row.names = FALSE, na = ".")
    if (verbose) message("  CSV written: ", fp_csv)
    files <- c(files, fp_csv)
  }

  structure(list(data      = df,
                 mapping   = mapping,
                 plot_dv   = plot_dv,
                 plot_norm = plot_norm,
                 summary   = summary_df,
                 files     = files),
            class = c("nm_data_previsual", "list"))
}

#' @export
print.nm_data_previsual <- function(x, ...) {
  cat("Data set preview:", nrow(x$data), "rows\n")
  if (length(x$mapping)) {
    m <- unlist(x$mapping)
    m <- m[!is.na(m)]
    if (length(m)) cat("  columns:", paste(sprintf("%s = %s", names(m), m), collapse = ", "), "\n")
  }
  cat("  plots  :", paste(c(if (!is.null(x$plot_dv)) "plot_dv", if (!is.null(x$plot_norm)) "plot_norm"),
                          collapse = ", "), "\n")
  if (length(x$files)) cat("  files  :", paste(x$files, collapse = ", "), "\n")
  if (!is.null(x$summary)) { cat("\n"); print(as.data.frame(x$summary), row.names = FALSE) }
  invisible(x)
}


# =============================================================================
# nm_add_rate()
# -----------------------------------------------------------------------------
# Adds a RATE column to an existing NONMEM dataset WITHOUT modifying the
# original columns, then writes the result to a NEW CSV file.
#
# Rule: RATE = base_value (0) everywhere, except on dosing records
# (EVID 1/4, or AMT > 0 when there is no EVID) where RATE = dose_value.
#   dose_value = -2 (default): zero-order input DURATION estimated through
#                              a Dn parameter (e.g. D1) defined in $PK.
#   dose_value = -1          : RATE estimated through an Rn parameter.
#
# Multiple-route datasets (e.g. a subject receives an IV AND an oral dose):
# `dose_cmt` restricts dose_value to doses given in one (or several)
# compartment(s). Other doses (e.g. IV bolus) keep base_value.
#
# No-alteration principle: the file is read ENTIRELY AS CHARACTER (no numeric
# conversion, no NA substitution), so existing fields are written back
# unchanged. Only the RATE column is added.
# =============================================================================

#' Add a RATE column to a NONMEM dataset
#'
#' Reads the CSV verbatim (all columns as character) so that nothing is
#' reformatted, adds a RATE column based on the dosing records, and writes a
#' new file. The original columns (values, format, '#ID' header) are kept
#' unchanged.
#'
#' @param file        Path to the original CSV.
#' @param output_file    Name (or path) of the output CSV. `NULL` (default) ->
#'                    `"<base>_rate.csv"` in the same directory as `file`.
#'                    If only a file name is given (no directory), it is
#'                    placed next to the source file. (Formerly `out_file`, still accepted with a warning.)
#' @param ... Former argument name, accepted with a warning until version 1.0.
#' @param sep         Column separator (default `","`).
#' @param dose_value  RATE value on the targeted dosing records: -2 (default,
#'                    estimated duration) or -1 (estimated rate).
#' @param base_value  RATE value on all other records (default 0).
#' @param after       Name of the column after which RATE is inserted
#'                    (default `"AMT"`). `NULL` = added as the last column.
#' @param dose_cmt    CMT value(s) identifying the route(s) that receive
#'                    `dose_value`. `NULL` (default): all dosing records receive
#'                    `dose_value` (single-route dataset). E.g. `dose_cmt = 1`
#'                    to affect only the oral depot, leaving IV boluses (other
#'                    CMT) at `base_value`.
#' @param cmt_col     Name of the CMT column (default: auto-detected). Required
#'                    only when `dose_cmt` is supplied.
#' @param evid_col    Name of the EVID column (default: auto-detected).
#' @param amt_col     Name of the AMT column (default: auto-detected; used only
#'                    when no EVID column is found).
#' @param overwrite   If `TRUE`, overwrite an existing RATE column (default
#'                    `FALSE`: the function stops if RATE already exists).
#' @param verbose     Print a diagnostic report (default `TRUE`).
#'
#' @return The path of the written file (invisibly).
#'
#' @examples
#' \dontrun{
#' # Single route: RATE=-2 on all doses -> data_nmol_rate.csv
#' nm_add_rate("C:/nm74g64/run/Bienvenu/data_nmol.csv")
#'
#' # Estimated rate rather than duration
#' nm_add_rate("data.csv", dose_value = -1)
#'
#' # Multiple routes: -2 only on the oral depot (CMT=1),
#' # IV boluses (CMT=2) keep RATE=0
#' nm_add_rate("data.csv", dose_cmt = 1)
#'
#' # Chosen output name (placed next to the source file)
#' nm_add_rate("data.csv", output_file = "data_dur.csv")
#' }
#' @importFrom utils read.csv write.table
#' @importFrom tools file_path_sans_ext
#' @export
nm_add_rate <- function(file,
                        output_file = NULL,
                        sep        = ",",
                        dose_value = -2,
                        base_value = 0,
                        after      = "AMT",
                        dose_cmt   = NULL,
                        cmt_col    = NULL,
                        evid_col   = NULL,
                        amt_col    = NULL,
                        overwrite  = FALSE,
                        verbose    = TRUE,
                        ...) {
  dep <- .nm_deprecated_args(list(...), c(out_file = "output_file"), "nm_add_rate")
  if (!is.null(dep$output_file)) output_file <- dep$output_file
  out_file <- output_file

  if (!is.character(file) || !file.exists(file))
    stop("File not found: ", file)

  # ---------------------------------------------------------------------------
  # 1. Verbatim read: all character, no NA / quote conversion
  # ---------------------------------------------------------------------------
  df <- utils::read.csv(file, sep = sep, colClasses = "character",
                        check.names = FALSE, na.strings = character(0),
                        quote = "", comment.char = "")

  if (ncol(df) == 0L || nrow(df) == 0L)
    stop("The file contains no usable data.")

  nms_up <- toupper(names(df))

  # ---------------------------------------------------------------------------
  # 2. Existing RATE column?
  # ---------------------------------------------------------------------------
  if ("RATE" %in% nms_up) {
    if (!overwrite)
      stop("A RATE column already exists. Use overwrite=TRUE to replace it.")
    df[[ names(df)[nms_up == "RATE"][1] ]] <- NULL
    nms_up <- toupper(names(df))
    if (verbose) message("Existing RATE column replaced.")
  }

  # ---------------------------------------------------------------------------
  # 3. Locate dosing records + CMT column (multiple routes)
  # ---------------------------------------------------------------------------
  find_col <- function(user, label) {
    if (!is.null(user)) {
      if (!user %in% names(df)) stop(sprintf("Column '%s' not found.", user))
      return(user)
    }
    hit <- names(df)[nms_up == label]
    if (length(hit) >= 1L) hit[1] else NA_character_
  }
  evid_name <- find_col(evid_col, "EVID")
  amt_name  <- find_col(amt_col,  "AMT")
  cmt_name  <- find_col(cmt_col,  "CMT")

  if (!is.na(evid_name)) {
    ev      <- suppressWarnings(as.numeric(df[[evid_name]]))
    is_dose <- ev %in% c(1, 4)
    if (verbose)
      message(sprintf("Dosing records identified from EVID ('%s' = 1/4).",
                      evid_name))
  } else if (!is.na(amt_name)) {
    am      <- suppressWarnings(as.numeric(df[[amt_name]]))
    is_dose <- !is.na(am) & am > 0
    if (verbose)
      message(sprintf("No EVID: dosing identified from AMT ('%s' > 0).",
                      amt_name))
  } else {
    stop("Cannot identify dosing records: ",
         "no EVID or AMT column found (supply evid_col= or amt_col=).")
  }

  # ---------------------------------------------------------------------------
  # 4. Restrict by route (CMT) for multiple-route datasets
  # ---------------------------------------------------------------------------
  apply_dose <- is_dose
  if (!is.null(dose_cmt)) {
    if (is.na(cmt_name))
      stop("dose_cmt supplied but no CMT column found (supply cmt_col=).")
    cmt_num    <- suppressWarnings(as.numeric(df[[cmt_name]]))
    in_cmt     <- cmt_num %in% dose_cmt
    apply_dose <- is_dose & in_cmt
    n_target   <- sum(apply_dose, na.rm = TRUE)
    n_other    <- sum(is_dose & !in_cmt, na.rm = TRUE)
    if (n_target == 0L)
      warning(sprintf("No dosing record with CMT %s: RATE=%s will not be applied anywhere.",
                      paste(dose_cmt, collapse = ","), dose_value))
    if (verbose)
      message(sprintf(
        "CMT restriction %s: %d targeted dose(s) -> RATE=%s; %d other dose(s) -> RATE=%s.",
        paste(dose_cmt, collapse = ","), n_target, dose_value, n_other, base_value))
  }

  # ---------------------------------------------------------------------------
  # 5. Build the RATE column (as character, like the rest)
  # ---------------------------------------------------------------------------
  rate_num <- ifelse(apply_dose, dose_value, base_value)
  df$RATE  <- format(rate_num, trim = TRUE, scientific = FALSE)  # "-2" / "0"

  # ---------------------------------------------------------------------------
  # 6. Position of the RATE column
  # ---------------------------------------------------------------------------
  if (!is.null(after)) {
    if (!toupper(after) %in% toupper(names(df)) || toupper(after) == "RATE") {
      if (verbose)
        message(sprintf("Column '%s' not found: RATE added as the last column.",
                        after))
    } else {
      other     <- setdiff(names(df), "RATE")            # original order
      apos      <- match(toupper(after), toupper(other))
      new_order <- append(other, "RATE", after = apos)   # RATE right after
      df        <- df[, new_order, drop = FALSE]
    }
  }

  # ---------------------------------------------------------------------------
  # 7. Output file name
  # ---------------------------------------------------------------------------
  src_dir <- dirname(file)
  base    <- tools::file_path_sans_ext(basename(file))
  if (is.null(out_file)) {
    out_file <- file.path(src_dir, paste0(base, "_rate.csv"))
  } else if (dirname(out_file) == ".") {
    # only a file name supplied -> next to the source file
    out_file <- file.path(src_dir, out_file)
  }
  if (normalizePath(out_file, mustWork = FALSE) ==
      normalizePath(file,     mustWork = FALSE))
    stop("The output file is the same as the source file. ",
         "Choose another name with output_file=.")

  # ---------------------------------------------------------------------------
  # 8. Verbatim write (no quotes, '#ID' header kept)
  # ---------------------------------------------------------------------------
  utils::write.table(df, out_file, sep = sep, quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "")

  if (verbose) {
    message(sprintf("Total records     : %d", nrow(df)))
    message(sprintf("Records RATE=%-3s : %d", dose_value, sum(apply_dose, na.rm = TRUE)))
    message(sprintf("Records RATE=%-3s : %d", base_value, sum(!apply_dose, na.rm = TRUE)))
    message("File written: ", out_file)
  }

  invisible(out_file)
}
