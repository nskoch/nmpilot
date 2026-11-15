# =============================================================================
# nm_read_csv.R
# -----------------------------------------------------------------------------
# Read and normalise a NONMEM-like dataset in CSV format.
#
# Roles:
#   - Identify NONMEM reserved columns (ID, DV, TIME, DATE, MDV, EVID,
#     AMT, RATE, II, SS, ADDL, CMT, DAT1-3...).
#   - Detect the TIME format (clock time hh:mm[:ss] vs continuous) and convert
#     to relative time (t0 = first dose per subject), using DATE/DAT1/DAT2/DAT3
#     when present.
#   - Compute DV_norm = DV / first dose when AMT is available.
#   - Flag BLQ records (BLQ/CENS column auto-detected, or via lloq).
#
# NONMEM conventions (Guide IV - III.B.2):
#   DATE = mdy, DAT1 = dmy, DAT2 = ymd, DAT3 = ydm.
#
# Returns a clean data.frame and a mapping only. No plot.
# =============================================================================

#' Read and normalise a NONMEM-like dataset
#'
#' @param file       Path to the CSV file.
#' @param sep        Column separator (default `","`).
#' @param na         Strings to be treated as missing values.
#' @param time_unit  Output unit for TIME: `"h"` (default) or `"d"`.
#' @param lloq       LLOQ value. If supplied and no BLQ/CENS column exists,
#'                   records with DV < lloq are flagged as BLQ.
#' @param blq_col    Name of the BLQ flag column (1 = BLQ, 0 = quantified).
#'                   By default, `"BLQ"` or `"CENS"` is auto-detected.
#' @param date_format Date format forced in `lubridate::parse_date_time()`
#'                   (e.g. `"dmy"` for DD/MM/YYYY). By default, chosen from the
#'                   name of the DATE/DAT1-3 column.
#' @param dv_col     Name of the DV column if non-standard (e.g. `"DV_nmol/L"`).
#'                   If `NULL` (default), a column whose name starts with "DV"
#'                   is auto-detected. Required when several DV columns exist.
#' @param verbose    Print the identification report (default `TRUE`).
#'
#' @return A list:
#'   - `data`: data.frame with relative TIME, BLQ (logical), and DV_norm when
#'     AMT is available. Reserved columns are renamed in upper case.
#'   - `mapping`: reserved labels -> original column names.
#'   - `lloq`: LLOQ value used (`NULL` if not supplied).
#'
#' @importFrom utils read.csv write.csv
#' @importFrom magrittr %>%
#' @importFrom lubridate parse_date_time
#' @importFrom dplyr group_by filter summarise mutate select left_join bind_rows first n_distinct
#' @export
#' @examples
#' \dontrun{
#' res <- nm_read_csv("data.csv", lloq = 0.5)
#' head(res$data)
#' }
nm_read_csv <- function(file,
                        sep         = ",",
                        na          = c("", "NA", "."),
                        time_unit   = c("h", "d"),
                        lloq        = NULL,
                        blq_col     = NULL,
                        date_format = NULL,
                        dv_col      = NULL,
                        verbose     = TRUE) {

  time_unit <- match.arg(time_unit)

  # ---------------------------------------------------------------------------
  # 1. Raw read
  # ---------------------------------------------------------------------------
  df <- utils::read.csv(file, sep = sep, na.strings = na,
                        stringsAsFactors = FALSE, check.names = FALSE)

  # Clean column names: remove leading '#' (NONMEM comment marker on first col)
  names(df) <- sub("^#+", "", names(df))

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
      names(df)[hit] <- r
    } else if (length(hit) > 1L) {
      warning(sprintf("Several columns match the label '%s': %s",
                      r, paste(names(df)[hit], collapse = ", ")))
    }
  }

  if (verbose) {
    message("=== Identified NONMEM columns ===")
    if (length(mapping) == 0L) message("  (none)")
    else for (k in names(mapping))
      message(sprintf("  %-6s -> %s", k, mapping[[k]]))
  }

  # ---------------------------------------------------------------------------
  # 3. Minimal checks
  # ---------------------------------------------------------------------------
  # Manual override for DV column
  if (!is.null(dv_col)) {
    if (!dv_col %in% names(df)) {
      warning(sprintf("dv_col '%s' not found in dataset.", dv_col))
    } else {
      mapping[["DV"]] <- dv_col
      names(df)[names(df) == dv_col] <- "DV"
    }
  }

  # Auto-detect DV if still not found: look for column starting with "DV"
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

  if (!"TIME" %in% names(mapping))
    stop("No TIME column identified.")
  if (!"ID" %in% names(mapping)) {
    warning("No ID column identified: a single subject is assumed.")
    df$ID <- 1L
    mapping$ID <- "ID"
  }

  # ---------------------------------------------------------------------------
  # 4. Detect the TIME format
  # ---------------------------------------------------------------------------
  time_raw <- df$TIME
  is_clock <- is.character(time_raw) && any(grepl(":", time_raw))
  date_label <- intersect(c("DATE", "DAT1", "DAT2", "DAT3"), names(mapping))
  has_date <- length(date_label) >= 1L
  if (has_date) date_label <- date_label[1L]

  if (verbose) {
    message(sprintf("TIME format: %s",
                    if (is_clock) "clock time (hh:mm[:ss])" else "numeric"))
    if (has_date)
      message(sprintf("DATE column detected: %s", date_label))
  }

  # ---------------------------------------------------------------------------
  # 5. Convert TIME to relative time
  # ---------------------------------------------------------------------------
  if (is_clock || has_date) {

    date_orders <- if (!is.null(date_format)) {
      # User-supplied format takes priority
      date_format
    } else if (has_date) {
      switch(date_label,
             "DATE" = c("mdy", "mdY", "dmy", "dmY"),  # try both mdy and dmy
             "DAT1" = c("dmy", "dmY"),
             "DAT2" = c("ymd", "Ymd"),
             "DAT3" = c("ydm", "Ydm"))
    } else NULL

    # Normalise clock times
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
        ss <- round((((num - hh) * 60) - mm) * 60)
        sprintf("%02d:%02d:%02d", hh, mm, ss)
      }, character(1))
      unname(out)
    }

    t_clean <- clean_time(time_raw)

    if (has_date) {
      df$.DT <- lubridate::parse_date_time(
        paste(df[[date_label]], t_clean),
        orders = paste(date_orders, "HMS"),
        quiet  = TRUE)
    } else {
      df$.DT <- lubridate::parse_date_time(
        paste("1970-01-01", t_clean),
        orders = "Ymd HMS", quiet = TRUE)
    }

    n_fail <- sum(is.na(df$.DT))
    if (n_fail > 0L && verbose)
      message(sprintf("  %d TIME/DATE value(s) could not be parsed.", n_fail))

    # t0 per subject = first dose, otherwise first record
    if ("EVID" %in% names(mapping) || "AMT" %in% names(mapping)) {
      ref <- df %>%
        group_by(ID) %>%
        filter(
          (if ("EVID" %in% names(mapping)) EVID %in% c(1, 4) else FALSE) |
          (if ("AMT"  %in% names(mapping)) !is.na(AMT) & AMT > 0 else FALSE)
        ) %>%
        summarise(.t0 = suppressWarnings(min(.DT, na.rm = TRUE)),
                  .groups = "drop")

      no_dose <- setdiff(unique(df$ID), ref$ID)
      if (length(no_dose) > 0L) {
        if (verbose)
          message(sprintf("  %d subject(s) without dose: t0 = first record.",
                          length(no_dose)))
        fb <- df %>% filter(ID %in% no_dose) %>%
          group_by(ID) %>%
          summarise(.t0 = suppressWarnings(min(.DT, na.rm = TRUE)),
                    .groups = "drop")
        ref <- bind_rows(ref, fb)
      }
    } else {
      ref <- df %>%
        group_by(ID) %>%
        summarise(.t0 = suppressWarnings(min(.DT, na.rm = TRUE)),
                  .groups = "drop")
      if (verbose)
        message("  No EVID/AMT: t0 = first record per subject.")
    }

    df <- df %>%
      left_join(ref, by = "ID") %>%
      mutate(TIME = as.numeric(difftime(.DT, .t0, units = "hours"))) %>%
      select(-.DT, -.t0)
  } else {
    df$TIME <- suppressWarnings(as.numeric(df$TIME))
  }

  if (time_unit == "d") df$TIME <- df$TIME / 24

  # ---------------------------------------------------------------------------
  # 6. DV_norm = DV / dose_ref (when AMT is available)
  # ---------------------------------------------------------------------------
  if ("AMT" %in% names(mapping) && "DV" %in% names(mapping)) {
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
        "  %d subject(s) with heterogeneous doses; DV_norm uses the first one.",
        n_multi))

    df <- df %>%
      left_join(dose_ref %>% select(ID, .dose_ref), by = "ID") %>%
      mutate(DV_norm = suppressWarnings(as.numeric(DV)) / .dose_ref) %>%
      select(-.dose_ref)
  }

  # ---------------------------------------------------------------------------
  # 7. BLQ flag (logical column)
  # ---------------------------------------------------------------------------
  df$BLQ <- FALSE

  if (is.null(blq_col)) {
    cand <- intersect(c("BLQ", "CENS"), toupper(names(df)))
    if (length(cand) >= 1L) {
      # Do not pick up the BLQ column created just above
      cand_real <- setdiff(cand, "BLQ")
      if (length(cand_real) == 0L && "BLQ" %in% toupper(names(df))) {
        # A BLQ column already existed in the original file
        orig_blq <- names(df)[toupper(names(df)) == "BLQ"]
        orig_blq <- setdiff(orig_blq, "BLQ")
        if (length(orig_blq) >= 1L) blq_col <- orig_blq[1]
      } else if (length(cand_real) >= 1L) {
        blq_col <- names(df)[toupper(names(df)) == cand_real[1]][1]
      }
    }
  }

  if (!is.null(blq_col) && blq_col %in% names(df)) {
    flag <- suppressWarnings(as.numeric(df[[blq_col]])) == 1
    flag[is.na(flag)] <- FALSE
    df$BLQ <- flag
    if (verbose)
      message(sprintf("BLQ flag from column '%s': %d record(s).",
                      blq_col, sum(df$BLQ)))
  } else if (!is.null(lloq) && "DV" %in% names(mapping)) {
    dv_num <- suppressWarnings(as.numeric(df$DV))
    df$BLQ <- !is.na(dv_num) & dv_num < lloq
    if (verbose)
      message(sprintf("BLQ identified from lloq=%g: %d record(s).",
                      lloq, sum(df$BLQ)))
  }

  # ---------------------------------------------------------------------------
  # 8. TAD (time after last dose) per subject
  # ---------------------------------------------------------------------------
  if ("EVID" %in% names(mapping) && "TIME" %in% names(df)) {
    dose_times <- df[df$EVID %in% c(1, 4) & !is.na(df$TIME),
                     c("ID", "TIME")]
    names(dose_times)[2] <- ".dose_TIME"

    df$TAD <- NA_real_
    for (sid in unique(df$ID)) {
      idx_subj  <- which(df$ID == sid)
      dt_subj   <- dose_times$.dose_TIME[dose_times$ID == sid]
      times_subj <- df$TIME[idx_subj]
      tad_subj  <- vapply(times_subj, function(t) {
        if (is.na(t)) return(NA_real_)
        eligible <- dt_subj[!is.na(dt_subj) & dt_subj <= t]
        if (length(eligible) == 0L) return(NA_real_)
        t - max(eligible)
      }, numeric(1))
      df$TAD[idx_subj] <- tad_subj
    }
    if (verbose) message(sprintf("  TAD computed (%d values).", sum(!is.na(df$TAD))))
  }

  # ---------------------------------------------------------------------------
  # 9. Return attributes
  # ---------------------------------------------------------------------------
  attr(df, "time_unit") <- time_unit
  attr(df, "lloq")      <- lloq
  attr(df, "mapping")   <- mapping

  invisible(list(data    = df,
                 mapping = mapping,
                 lloq    = lloq))
}
