# =============================================================================
# nm_check_dataset.R
# -----------------------------------------------------------------------------
# Checks that a dataset can be read by NONMEM (point B5). Report only: the
# data are never modified.
# =============================================================================

#' Check that a dataset is compatible with NONMEM
#'
#' Reads a dataset as NONMEM would (every value as text, `.` or an empty cell
#' as missing) and reports the problems that make NM-TRAN refuse the data or
#' NONMEM misread them. The data are never modified.
#'
#' Checks:
#' * **columns**: ID, TIME and DV present (error); AMT absent (warning: no
#'   dose record);
#' * **missing ID or TIME** values (error);
#' * **non-numeric values** in a column (error; NONMEM reads numbers only,
#'   unless the column is dropped in `$INPUT`);
#' * **times**: TIME decreasing within a subject (error), except at a reset
#'   record (EVID 3 or 4);
#' * **EVID / MDV / AMT**: EVID outside 0 to 4 (error); dose record (EVID 1
#'   or 4) without a positive AMT (error); observation (EVID 0) with a
#'   positive AMT (error); MDV 0 with a missing DV (error); MDV 0 on a record
#'   that is not an observation (error); observation with MDV 1 but a DV
#'   (warning: the value is ignored); without EVID, a record with AMT > 0 is
#'   a dose and MDV must be 1;
#' * **duplicates**: identical records (warning); two observations of the
#'   same subject, time and compartment with different DV (warning);
#' * **outliers** in the observed DV: negative values (warning) and values
#'   far from the observations taken at the same time after the last dose
#'   (log scale, more than `outlier_k` median absolute deviations from their
#'   median; per compartment when CMT exists; groups of at least 10
#'   observations) (warning).
#'
#' The header line may start with `#` or `@` (`#ID,TIME,...`), as read with
#' `IGNORE=#` or `IGNORE=@` in `$DATA`.
#'
#' @param data Path of a CSV file, or a data frame (for example the `data`
#'   element of [nm_read_csv()]).
#' @param cols Optional named character vector giving the column of a
#'   standard item when its name differs, e.g. `c(ID = "SUBJ", DV = "CONC")`.
#'   Names are otherwise matched without regard to case.
#' @param sep Column separator of the CSV file (default `","`).
#' @param outlier_k Number of median absolute deviations (log DV) beyond which
#'   an observation is reported as an outlier (default 4).
#' @param verbose Print the report (default `TRUE`).
#'
#' @return An object of class `nm_dataset_check`: a list with `findings`
#'   (one row per problem found: `check`, `severity` "error" or "warning",
#'   `column`, `n_rows`, `rows` (row numbers in the data, header excluded),
#'   `ids` and `message`), `n_rows`, `n_subjects`, `columns` (standard item
#'   -> column used) and `ok` (`TRUE` when no error was found). Printed
#'   (invisibly returned) when `verbose = TRUE`.
#' @export
#' @examples
#' ex <- system.file("extdata", "example", package = "nmpilot")
#' chk <- nm_check_dataset(file.path(ex, "example_data.csv"))
#' chk$ok
#' chk$findings
nm_check_dataset <- function(data, cols = NULL, sep = ",", outlier_k = 4,
                             verbose = TRUE) {
  # ---- read as text, as NONMEM does ------------------------------------------
  if (is.character(data) && length(data) == 1L) {
    if (!file.exists(data)) stop("File not found: ", data, call. = FALSE)
    raw <- utils::read.csv(data, sep = sep, colClasses = "character", check.names = FALSE,
                           na.strings = character(0), strip.white = TRUE)
    source <- data
  } else if (is.data.frame(data)) {
    raw <- data
    # columns added by nm_read_csv() are not in the file
    if (!is.null(attr(data, "mapping")))
      raw <- raw[, setdiff(names(raw), c("DV_norm", "BLQ", "TAD")), drop = FALSE]
    raw[] <- lapply(raw, function(v) if (is.factor(v)) as.character(v) else v)
    source <- NA_character_
  } else {
    stop("`data` must be the path of a CSV file or a data frame.", call. = FALSE)
  }
  # header written "#ID" or "@ID" (skipped by IGNORE=# / IGNORE=@ in $DATA)
  if (ncol(raw)) names(raw)[1] <- sub("^[#@]\\s*", "", names(raw)[1])
  n <- nrow(raw)
  findings <- list()
  add <- function(check, severity, column, rows, message, ids = NULL) {
    rows <- sort(unique(rows))
    if (!length(rows) && !identical(check, "columns")) return(invisible())
    findings[[length(findings) + 1L]] <<- data.frame(
      check = check, severity = severity, column = column, n_rows = length(rows),
      rows = I(list(rows)),
      ids = if (is.null(ids)) "" else paste(utils::head(unique(ids), 20), collapse = ", "),
      message = message, stringsAsFactors = FALSE)
  }

  # ---- columns -----------------------------------------------------------------
  items <- c("ID", "TIME", "DV", "AMT", "EVID", "MDV", "CMT")
  col_of <- stats::setNames(rep(NA_character_, length(items)), items)
  up <- toupper(names(raw))
  for (it in items) {
    k <- which(up == it)
    if (length(k)) col_of[it] <- names(raw)[k[1]]
  }
  if (!is.null(cols)) {
    if (is.null(names(cols)) || any(!names(cols) %in% items))
      stop("`cols` must be named with standard items: ", paste(items, collapse = ", "), call. = FALSE)
    miss <- setdiff(cols, names(raw))
    if (length(miss)) stop("Column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
    col_of[names(cols)] <- cols
  }
  for (it in c("ID", "TIME", "DV"))
    if (is.na(col_of[it]))
      add("columns", "error", it, integer(0), sprintf("No %s column: NONMEM needs one.", it))
  if (is.na(col_of["AMT"]))
    add("columns", "warning", "AMT", integer(0), "No AMT column: the dataset has no dose record.")

  # ---- missing values and non-numeric values -----------------------------------
  is_miss <- function(v) is.na(v) | (is.character(v) & trimws(v) %in% c("", "."))
  num <- list()
  for (nm in names(raw)) {
    v <- raw[[nm]]
    if (is.numeric(v)) { num[[nm]] <- as.numeric(v); next }
    if (is.logical(v)) { num[[nm]] <- as.numeric(v); next }
    m <- is_miss(v)
    x <- suppressWarnings(as.numeric(trimws(as.character(v))))
    bad <- which(!m & is.na(x))
    if (length(bad))
      add("non_numeric", "error", nm, bad,
          sprintf("Non-numeric values in %s (e.g. \"%s\"): NONMEM reads numbers only; recode the column or drop it in $INPUT.",
                  nm, as.character(v[bad[1]])),
          ids = if (!is.na(col_of["ID"])) raw[[col_of["ID"]]][bad])
    x[m] <- NA_real_
    num[[nm]] <- x
  }
  get <- function(it) if (is.na(col_of[it])) NULL else num[[col_of[it]]]
  id <- get("ID"); time <- get("TIME"); dv <- get("DV"); amt <- get("AMT")
  evid <- get("EVID"); mdv <- get("MDV"); cmt <- get("CMT")
  idv <- if (!is.null(id)) id else rep(NA_real_, n)
  for (it in c("ID", "TIME")) {
    v <- get(it)
    if (!is.null(v)) {
      r <- which(is.na(v))
      add("missing_values", "error", col_of[it], r,
          sprintf("Missing %s: every record needs one.", it), ids = idv[r])
    }
  }

  # ---- times -------------------------------------------------------------------
  if (!is.null(id) && !is.null(time) && n > 1L) {
    same <- c(FALSE, id[-1] == id[-n])
    dec <- c(FALSE, time[-1] < time[-n])
    reset <- if (!is.null(evid)) evid %in% c(3, 4) else rep(FALSE, n)
    r <- which(same & dec & !reset & !is.na(same) & !is.na(dec))
    add("time_order", "error", col_of["TIME"], r,
        "TIME decreases within a subject (not at a reset record EVID 3 or 4): sort by ID and TIME.",
        ids = idv[r])
  }

  # ---- EVID / MDV / AMT --------------------------------------------------------
  amt_pos <- if (!is.null(amt)) !is.na(amt) & amt > 0 else rep(FALSE, n)
  not_mdv1 <- if (!is.null(mdv)) !mdv %in% 1 else rep(TRUE, n)
  if (!is.null(evid)) {
    r <- which(!is.na(evid) & !evid %in% 0:4)
    add("evid_mdv_amt", "error", col_of["EVID"], r, "EVID must be 0, 1, 2, 3 or 4.", ids = idv[r])
    r <- which(evid %in% c(1, 4) & !amt_pos)
    add("evid_mdv_amt", "error", col_of["AMT"] %|na|% "AMT", r,
        "Dose record (EVID 1 or 4) without a positive AMT.", ids = idv[r])
    r <- which(evid %in% 0 & amt_pos)
    add("evid_mdv_amt", "error", col_of["AMT"] %|na|% "AMT", r,
        "Observation record (EVID 0) with a positive AMT.", ids = idv[r])
    is_obs <- evid %in% 0
  } else {
    is_obs <- !amt_pos
  }
  if (!is.null(mdv)) {
    if (!is.null(dv)) {
      r <- which(mdv %in% 0 & is.na(dv))
      add("evid_mdv_amt", "error", col_of["MDV"], r, "MDV 0 with a missing DV.", ids = idv[r])
      r <- which(is_obs & mdv %in% 1 & !is.na(dv))
      add("evid_mdv_amt", "warning", col_of["MDV"], r,
          "Observation with MDV 1 but a DV: the value is ignored by NONMEM.", ids = idv[r])
    }
    r <- which(mdv %in% 0 & !is_obs)
    add("evid_mdv_amt", "error", col_of["MDV"], r,
        if (!is.null(evid)) "MDV 0 on a record that is not an observation (EVID 1 to 4): MDV must be 1."
        else "MDV 0 on a dose record (AMT > 0): MDV must be 1.", ids = idv[r])
  }

  # ---- duplicates --------------------------------------------------------------
  key <- do.call(paste, c(lapply(raw, as.character), sep = "\r"))
  r <- which(duplicated(key) | duplicated(key, fromLast = TRUE))
  add("duplicates", "warning", "", r, "Identical records.", ids = idv[r])
  if (!is.null(id) && !is.null(time) && !is.null(dv)) {
    ob <- which(is_obs & !is.na(dv) & not_mdv1)
    k <- paste(id[ob], time[ob], if (!is.null(cmt)) cmt[ob] else "", sep = "|")
    grp <- split(ob, k)
    r <- unlist(lapply(grp, function(g) if (length(g) > 1L && length(unique(dv[g])) > 1L) g))
    add("duplicates", "warning", col_of["DV"], r,
        "Several observations of the same subject, time and compartment with different DV.",
        ids = idv[r])
  }

  # ---- outliers ----------------------------------------------------------------
  if (!is.null(dv)) {
    ob <- which(is_obs & !is.na(dv) & not_mdv1)
    r <- ob[dv[ob] < 0]
    add("outliers", "warning", col_of["DV"], r, "Negative observed DV.", ids = idv[r])
    pos <- ob[dv[ob] > 0]
    # observations are compared with those taken at the same time after the
    # last dose (concentrations fall with time: a late sample is not an
    # outlier), per compartment when CMT exists
    tad <- if (!is.null(time)) time else rep(0, n)
    if (!is.null(id) && !is.null(time)) {
      last <- rep(NA_real_, n); cur <- NA_real_
      for (i in seq_len(n)) {
        if (i > 1L && !identical(id[i], id[i - 1L])) cur <- NA_real_
        if (amt_pos[i] && !is.na(time[i])) cur <- time[i]
        last[i] <- cur
      }
      tad <- ifelse(is.na(last), time, time - last)
    }
    tp <- tad[pos]
    ut <- unique(round(tp[!is.na(tp)], 6))
    bin <- if (length(ut) <= 30L) round(tp, 6) else {
      br <- unique(stats::quantile(tp, seq(0, 1, length.out = 11), na.rm = TRUE))
      if (length(br) < 2L) rep(1, length(tp)) else cut(tp, br, include.lowest = TRUE, labels = FALSE)
    }
    g <- paste(if (!is.null(cmt)) cmt[pos] else 1, bin)
    out <- unlist(lapply(split(pos, g), function(p) {
      if (length(p) < 10L) return(integer(0))
      l <- log(dv[p]); md <- stats::median(l); s <- stats::mad(l)
      if (!is.finite(s) || s <= 0) return(integer(0))
      p[abs(l - md) > outlier_k * s]
    }))
    add("outliers", "warning", col_of["DV"], out,
        sprintf("Observed DV far from the others taken at the same time after dose (more than %g median absolute deviations from their median, log scale%s).",
                outlier_k, if (!is.null(cmt)) ", per CMT" else ""), ids = idv[out])
  }

  fnd <- if (length(findings)) do.call(rbind, findings) else
    data.frame(check = character(0), severity = character(0), column = character(0),
               n_rows = integer(0), rows = I(list()), ids = character(0),
               message = character(0), stringsAsFactors = FALSE)
  if (nrow(fnd)) {
    fnd <- fnd[order(match(fnd$severity, c("error", "warning")),
                     match(fnd$check, c("columns", "missing_values", "non_numeric", "time_order",
                                        "evid_mdv_amt", "duplicates", "outliers"))), ]
    rownames(fnd) <- NULL
  }
  res <- structure(list(findings = fnd, n_rows = n,
                        n_subjects = if (!is.null(id)) length(unique(id[!is.na(id)])) else NA_integer_,
                        columns = col_of, ok = !any(fnd$severity == "error"), source = source),
                   class = "nm_dataset_check")
  if (verbose) print(res)
  invisible(res)
}

# compact row list: "3, 5, 8-10, ..."
.nm_cd_rows <- function(r, max = 12L) {
  if (!length(r)) return("")
  r <- sort(r); br <- c(0, which(diff(r) != 1), length(r))
  s <- vapply(seq_len(length(br) - 1), function(i) {
    a <- r[br[i] + 1]; b <- r[br[i + 1]]
    if (a == b) as.character(a) else paste0(a, "-", b)
  }, character(1))
  if (length(s) > max) paste0(paste(s[seq_len(max)], collapse = ", "), ", ...") else paste(s, collapse = ", ")
}

#' @export
print.nm_dataset_check <- function(x, ...) {
  f <- x$findings
  ne <- sum(f$severity == "error"); nw <- sum(f$severity == "warning")
  cat("<nm_dataset_check>", if (!is.na(x$source)) basename(x$source) else "", "\n")
  cat(sprintf("  %d records, %s subjects: %d error(s), %d warning(s)%s\n", x$n_rows,
              format(x$n_subjects), ne, nw,
              if (ne == 0L) " -> readable by NONMEM" else " -> to fix before NONMEM"))
  if (nrow(f)) {
    for (i in seq_len(nrow(f))) {
      cat(sprintf("  [%s] %s%s\n", toupper(f$severity[i]), f$message[i],
                  if (f$n_rows[i]) sprintf(" (%d record(s): rows %s%s)", f$n_rows[i],
                                           .nm_cd_rows(f$rows[[i]]),
                                           if (nzchar(f$ids[i])) paste0("; ID ", f$ids[i]) else "")
                  else ""))
    }
  }
  invisible(x)
}

#' @export
summary.nm_dataset_check <- function(object, ...) {
  f <- object$findings
  if (!nrow(f))
    return(data.frame(check = character(0), severity = character(0), n_findings = integer(0),
                      n_rows = integer(0)))
  k <- paste(f$severity, f$check, sep = "\r")
  s <- data.frame(check = f$check[!duplicated(k)], severity = f$severity[!duplicated(k)],
                  n_findings = as.integer(tapply(f$n_rows, k, length)[k[!duplicated(k)]]),
                  n_rows = as.integer(tapply(f$n_rows, k, sum)[k[!duplicated(k)]]),
                  stringsAsFactors = FALSE)
  rownames(s) <- NULL
  s
}
