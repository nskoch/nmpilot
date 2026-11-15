# =============================================================================
# nm_extract_runs.R
# -----------------------------------------------------------------------------
# Extract and compare several NONMEM runs.
#
# Design
#   - Parameter names: when a .ctl is present, nm_identify_parameters() uses
#     the $PK equations to name THETA/ETA (CL, V, KA...) and the $THETA
#     comment labels to name the covariate-effect and residual-error THETAs
#     (CL_WT, PROP, ADD...); the
#     .smr file is a fallback. nm_compare_runs() aligns on these names by
#     default; align_by = "index" gives a strict NONMEM-index comparison.
#   - Run filtering (run_filter / nm_filter_runs()):
#       all                   : all runs (default)
#       successful            : min_ok == TRUE
#       successful_no_warning : successful and no warning
#       successful_no_comment : successful and no status_comment
#       clean                 : successful, no warning and no status_comment
#     Filtering recomputes dOFV, rank_ofv, dAIC and dBIC on the kept runs.
#   - Covariance step status:
#       STANDARD ERROR OF ESTIMATE      -> COMPLETED
#       COVARIANCE STEP OMITTED: NO, no SE -> REQUESTED
#       COVARIANCE STEP OMITTED: YES    -> OMITTED
#       no mention of covariance        -> NOT_REQUESTED
#   .ext : primary NONMEM numeric source
#          THETA, OMEGA, SIGMA, OFV, SE, RSE, FIX
#   .lst : run diagnostics
#          minimisation, covariance, warnings, shrinkage, control OFV
#   .smr : optional Wings enrichment
#          THETA/ETA/ERR names, ETASD/EPSSD when available, n_obs
#   .shk : optional fallback for shrinkage
#
# The run directory name/extension does not matter:
# .nm7, .reg or a directory without extension are accepted.
#
# Main functions:
#   nm_extract_runs()
#   nm_filter_runs()
#   nm_compare_runs()
#   print.nm_runs()
#
# Parameter types in $params:
#   THETA      : structural parameters
#   OMEGA      : diagonal ETA variances
#   OMEGA_COV  : ETA covariances
#   ETA_SD     : sqrt(diagonal OMEGA), derived
#   ETA_COR    : correlation derived from an OMEGA block
#   SIGMA      : diagonal EPS variances
#   SIGMA_COV  : EPS covariances
#   EPS_SD     : sqrt(diagonal SIGMA), derived
#   EPS_COR    : correlation derived from a SIGMA block
# =============================================================================


# =============================================================================
# 1. MAIN FUNCTION
# =============================================================================

#' Extract the results of several NONMEM runs
#'
#' @param path Analysis directory, or the directory of a single run.
#' @param runs Optional names of the runs to keep.
#' @param verbose Print progress.
#'
#' @return An object of class `nm_runs` with:
#'   - `$runs`: one row per run
#'   - `$params`: one row per parameter
#'   - `$path`: the analysis directory
#'   - `$created`: time of extraction
#'
#' @details
#' When the control stream of a run is found, parameters are named from the
#' model code with [nm_identify_parameters()]: typical-value THETAs and ETAs
#' after the PK parameter they act on (`CL`, `V`...), covariate-effect and
#' residual-error THETAs after their `$THETA` comment label when it is a short
#' name (`CL_WT`, `PROP`, `ADD`...). Other THETAs keep their NONMEM name
#' (`THETA6`...).
#'
#' @export
nm_extract_runs <- function(path, runs = NULL, verbose = TRUE) {

  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  rt <- .nm_find_runs(path, runs)

  if (nrow(rt) == 0L) {
    stop("No NONMEM run (.ext/.lst/.smr) found under: ", path, call. = FALSE)
  }

  parsed <- purrr::pmap(
    rt,
    function(run, dir, ctl, smr, lst, ext, shk) {
      if (isTRUE(verbose)) message("  * run : ", run)
      .nm_extract_one(run, dir, ctl, smr, lst, ext, shk)
    }
  )

  runs_df <- purrr::map_dfr(parsed, "run")
  params_df <- purrr::map_dfr(parsed, "params")

  safe_min <- function(x) {
    if (!length(x) || all(is.na(x))) return(NA_real_)
    min(x, na.rm = TRUE)
  }

  min_ofv <- safe_min(runs_df$ofv)
  min_aic <- safe_min(runs_df$aic)
  min_bic <- safe_min(runs_df$bic)

  runs_df <- runs_df %>%
    dplyr::arrange(.data$ofv) %>%
    dplyr::mutate(
      rank_ofv = dplyr::if_else(
        is.na(.data$ofv), NA_integer_,
        as.integer(rank(.data$ofv, ties.method = "first", na.last = "keep"))
      ),
      dOFV = dplyr::if_else(
        is.na(.data$ofv) | is.na(min_ofv), NA_real_, .data$ofv - min_ofv
      ),
      dAIC = dplyr::if_else(
        is.na(.data$aic) | is.na(min_aic), NA_real_, .data$aic - min_aic
      ),
      dBIC = dplyr::if_else(
        is.na(.data$bic) | is.na(min_bic), NA_real_, .data$bic - min_bic
      )
    ) %>%
    dplyr::relocate(dplyr::any_of(c(
      "run", "ofv", "dOFV", "rank_ofv",
      "min_status", "min_ok", "cov_status",
      "ext_available", "smr_available", "lst_available",
      "estimation_available"
    )))

  structure(
    list(
      runs = runs_df,
      params = params_df,
      path = path,
      created = Sys.time()
    ),
    class = "nm_runs"
  )
}


# =============================================================================
# 2. FILTERING AND COMPARISON
# =============================================================================

#' Filter NONMEM runs and recompute relative criteria
#'
#' @param x An object of class `nm_runs`.
#' @param run_filter Filter to apply: `"all"` (default), `"successful"`,
#'   `"successful_no_warning"`, `"successful_no_comment"` or `"clean"`.
#'
#' @return A filtered `nm_runs` object. The `$runs` and `$params` tables are
#'   filtered together; dOFV, rank_ofv, dAIC and dBIC are recomputed on the
#'   remaining runs.
#'
#' @export
nm_filter_runs <- function(
    x,
    run_filter = c(
      "all",
      "successful",
      "successful_no_warning",
      "successful_no_comment",
      "clean"
    )) {

  if (!inherits(x, "nm_runs")) {
    stop("`x` must be an object of class `nm_runs`.", call. = FALSE)
  }

  run_filter <- match.arg(run_filter)
  runs_df <- x$runs

  required <- c("run", "min_ok", "warnings", "status_comment")
  missing_cols <- setdiff(required, names(runs_df))
  if (length(missing_cols)) {
    stop(
      "Missing column(s) in x$runs: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  is_blank <- function(z) is.na(z) | trimws(as.character(z)) == ""

  keep <- switch(
    run_filter,
    all = rep(TRUE, nrow(runs_df)),
    successful = runs_df$min_ok %in% TRUE,
    successful_no_warning =
      runs_df$min_ok %in% TRUE & is_blank(runs_df$warnings),
    successful_no_comment =
      runs_df$min_ok %in% TRUE & is_blank(runs_df$status_comment),
    clean =
      runs_df$min_ok %in% TRUE &
      is_blank(runs_df$warnings) &
      is_blank(runs_df$status_comment)
  )

  runs_df <- runs_df[keep, , drop = FALSE]

  safe_min <- function(z) {
    if (!length(z) || all(is.na(z))) return(NA_real_)
    min(z, na.rm = TRUE)
  }

  min_ofv <- if ("ofv" %in% names(runs_df)) safe_min(runs_df$ofv) else NA_real_
  min_aic <- if ("aic" %in% names(runs_df)) safe_min(runs_df$aic) else NA_real_
  min_bic <- if ("bic" %in% names(runs_df)) safe_min(runs_df$bic) else NA_real_

  if ("ofv" %in% names(runs_df)) {
    runs_df$rank_ofv <- ifelse(
      is.na(runs_df$ofv),
      NA_integer_,
      as.integer(rank(runs_df$ofv, ties.method = "first", na.last = "keep"))
    )
    runs_df$dOFV <- ifelse(
      is.na(runs_df$ofv) | is.na(min_ofv),
      NA_real_,
      runs_df$ofv - min_ofv
    )
  }

  if ("aic" %in% names(runs_df)) {
    runs_df$dAIC <- ifelse(
      is.na(runs_df$aic) | is.na(min_aic),
      NA_real_,
      runs_df$aic - min_aic
    )
  }

  if ("bic" %in% names(runs_df)) {
    runs_df$dBIC <- ifelse(
      is.na(runs_df$bic) | is.na(min_bic),
      NA_real_,
      runs_df$bic - min_bic
    )
  }

  if ("ofv" %in% names(runs_df)) {
    runs_df <- runs_df[order(runs_df$ofv, na.last = TRUE), , drop = FALSE]
  }

  kept_runs <- runs_df$run
  params_df <- x$params
  if ("run" %in% names(params_df)) {
    params_df <- params_df[params_df$run %in% kept_runs, , drop = FALSE]
  }

  out <- x
  out$runs <- tibble::as_tibble(runs_df)
  out$params <- tibble::as_tibble(params_df)
  attr(out, "run_filter") <- run_filter
  out
}


#' Compare the parameters of several runs side by side
#'
#' @param x An `nm_runs` object or a params tibble.
#' @param what `"value"`, `"rse_pct"`, `"se"` or `"shrinkage"`.
#' @param types Parameter types to keep. `NULL` = all.
#' @param format If `TRUE`, cells are formatted as "value (RSE%)".
#' @param digits Significant digits when `format = TRUE`.
#' @param align_by `"name"` (default) for semantic alignment; `"index"` for a
#'   strict comparison by NONMEM index.
#' @param run_filter Run filter: `"all"` (default), `"successful"`,
#'   `"successful_no_warning"`, `"successful_no_comment"` or `"clean"`.
#'
#' @return A wide tibble with one row per parameter and one column per run.
#'
#' @examples
#' \dontrun{
#' nm_compare_runs(res, types = "THETA", format = TRUE)
#' nm_compare_runs(res, types = "THETA", format = TRUE, run_filter = "successful")
#' nm_compare_runs(res, types = "THETA", format = TRUE, run_filter = "clean")
#' res_clean <- nm_filter_runs(res, "clean")
#' res_clean$runs
#' }
#'
#' @export
nm_compare_runs <- function(
    x,
    what = c("value", "rse_pct", "se", "shrinkage"),
    types = NULL,
    format = FALSE,
    digits = 3,
    align_by = c("name", "index"),
    run_filter = c(
      "all",
      "successful",
      "successful_no_warning",
      "successful_no_comment",
      "clean"
    )) {

  what <- match.arg(what)
  align_by <- match.arg(align_by)
  run_filter <- match.arg(run_filter)

  if (inherits(x, "nm_runs")) {
    x_use <- nm_filter_runs(x, run_filter = run_filter)
    params <- x_use$params
  } else {
    if (run_filter != "all") {
      stop(
        "A `run_filter` other than \"all\" requires an object of class `nm_runs`.",
        call. = FALSE
      )
    }
    params <- x
  }

  needed <- c("run", "type", "index", "name")
  if (!all(needed %in% names(params))) {
    stop("`x` must contain the columns run/type/index/name.", call. = FALSE)
  }

  if (!is.null(types)) {
    unknown <- setdiff(types, unique(params$type))
    if (length(unknown)) {
      warning(
        "Type(s) not found in the table: ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    params <- dplyr::filter(params, .data$type %in% types)
  }

  if (nrow(params) == 0L) {
    if (align_by == "name") {
      return(tibble::tibble(type = character(), name = character()))
    }
    return(tibble::tibble(
      type = character(), index = integer(), name = character()
    ))
  }

  type_order <- c(
    "THETA",
    "OMEGA", "ETA_SD", "OMEGA_COV", "ETA_COR",
    "SIGMA", "EPS_SD", "SIGMA_COV", "EPS_COR"
  )

  if (isTRUE(format)) {
    params <- params %>%
      dplyr::mutate(
        .cell = dplyr::case_when(
          is.na(.data$value) ~ NA_character_,
          "fixed" %in% names(.) & .data$fixed & !is.na(.data$value) ~
            paste0(formatC(.data$value, format = "fg", digits = digits), " FIX"),
          "rse_pct" %in% names(.) & !is.na(.data$rse_pct) ~
            sprintf(
              "%s (%.0f%%)",
              formatC(.data$value, format = "fg", digits = digits),
              .data$rse_pct
            ),
          TRUE ~ formatC(.data$value, format = "fg", digits = digits)
        )
      )
    value_col <- ".cell"
  } else {
    if (!what %in% names(params)) {
      stop("Column `", what, "` does not exist in params.", call. = FALSE)
    }
    value_col <- what
  }

  params <- params %>%
    dplyr::mutate(.val = .data[[value_col]])

  if (align_by == "name") {
    dup <- params %>%
      dplyr::count(.data$run, .data$type, .data$name, name = ".n") %>%
      dplyr::filter(.data$.n > 1L)

    if (nrow(dup)) {
      stop(
        "Several parameters share the same name within a run. ",
        "Use `align_by = \"index\"` for a strict comparison by NONMEM index.",
        call. = FALSE
      )
    }

    out <- params %>%
      dplyr::select("type", "name", "run", ".val") %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(
        names_from = "run",
        values_from = ".val"
      ) %>%
      dplyr::arrange(
        match(.data$type, type_order),
        .data$name
      )
  } else {
    out <- params %>%
      dplyr::select("type", "index", "name", "run", ".val") %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(
        names_from = "run",
        values_from = ".val"
      ) %>%
      dplyr::arrange(
        match(.data$type, type_order),
        .data$index,
        .data$name
      )
  }

  out
}


# =============================================================================
# 3. PRINTING
# =============================================================================

#' @export
print.nm_runs <- function(x, ...) {
  cat("=== nm_runs : ", nrow(x$runs), " run(s) ===\n", sep = "")
  cat("Analysis: ", x$path, "\n\n", sep = "")

  cols <- intersect(
    c(
      "run", "ofv", "dOFV", "rank_ofv",
      "min_status", "cov_status",
      "estimation_available", "ext_available",
      "n_par", "aic", "bic",
      "max_eta_shrink", "warnings"
    ),
    names(x$runs)
  )

  compact <- x$runs[, cols, drop = FALSE]
  num <- vapply(compact, is.numeric, logical(1))
  compact[num] <- lapply(compact[num], function(v) round(v, 3))
  print(as.data.frame(compact), row.names = FALSE)
  invisible(x)
}


# =============================================================================
# 4. GENERAL TOOLS
# =============================================================================

.nm_read_lines_latin1 <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)

  # NONMEM/Wings on Windows: latin1/CP1252 is common.
  con <- file(path, encoding = "latin1")
  on.exit(close(con))
  lines <- tryCatch(readLines(con, warn = FALSE), error = function(e) NULL)
  if (is.null(lines)) return(NULL)

  suppressWarnings(iconv(lines, from = "latin1", to = "UTF-8", sub = "?"))
}


.nm_empty_params <- function() {
  tibble::tibble(
    type = character(),
    index = integer(),
    name = character(),
    raw_name = character(),
    value = numeric(),
    fixed = logical(),
    rse_pct = numeric(),
    se = numeric(),
    shrinkage = numeric(),
    source = character()
  )
}


# =============================================================================
# 5. RUN DISCOVERY
# =============================================================================

.nm_find_runs <- function(path, runs = NULL) {

  # The directory name (.nm7/.reg/other) does not matter.
  cand <- c(path, list.dirs(path, recursive = FALSE, full.names = TRUE))

  rows <- purrr::map(cand, function(d) {

    ctl <- list.files(d, pattern = "(?i)\\.ctl$", full.names = TRUE)
    lst <- list.files(d, pattern = "(?i)\\.lst$", full.names = TRUE)
    smr <- list.files(d, pattern = "(?i)\\.smr$", full.names = TRUE)
    ext <- list.files(d, pattern = "(?i)\\.ext$", full.names = TRUE)
    shk <- list.files(d, pattern = "(?i)\\.shk$", full.names = TRUE)

    if (!length(lst) && !length(smr) && !length(ext)) return(NULL)

    # Prefer the stem of the .ext, then .smr, then .lst.
    primary <- if (length(ext)) ext[1] else if (length(smr)) smr[1] else lst[1]
    stem <- tools::file_path_sans_ext(basename(primary))

    # If several files exist, prefer the one with the same stem.
    pick <- function(files) {
      if (!length(files)) return(NA_character_)
      stems <- tools::file_path_sans_ext(basename(files))
      hit <- which(stems == stem)
      if (length(hit)) files[hit[1]] else files[1]
    }

    tibble::tibble(
      run = stem,
      dir = d,
      ctl = pick(ctl),
      smr = pick(smr),
      lst = pick(lst),
      ext = pick(ext),
      shk = pick(shk)
    )
  })

  out <- purrr::compact(rows)

  if (!length(out)) {
    return(tibble::tibble(
      run = character(), dir = character(),
      ctl = character(), smr = character(), lst = character(),
      ext = character(), shk = character()
    ))
  }

  out <- dplyr::bind_rows(out) %>%
    dplyr::distinct(.data$run, .keep_all = TRUE)

  if (!is.null(runs)) {
    out <- dplyr::filter(out, .data$run %in% runs)
  }

  out
}


# =============================================================================
# 6. .EXT PARSER
# =============================================================================

# Read the last estimation block of the .ext containing ITERATION.
.nm_parse_ext <- function(path) {

  L <- .nm_read_lines_latin1(path)
  if (is.null(L) || !length(L)) return(NULL)

  header_idx <- stringr::str_which(L, "^\\s*ITERATION\\s+")
  if (!length(header_idx)) return(NULL)

  parsed_blocks <- lapply(header_idx, function(h) {

    # Block ends at the next TABLE NO. or the next ITERATION header.
    later_table <- which(
      seq_along(L) > h &
        stringr::str_detect(L, "^\\s*TABLE NO\\.")
    )
    later_header <- header_idx[header_idx > h]

    ends <- c(
      if (length(later_table)) later_table[1] - 1L else length(L),
      if (length(later_header)) later_header[1] - 1L else length(L)
    )
    e <- min(ends)

    hdr <- stringr::str_split(stringr::str_trim(L[h]), "\\s+")[[1]]

    dat_lines <- L[(h + 1L):e]
    dat_lines <- dat_lines[nzchar(stringr::str_trim(dat_lines))]
    if (!length(dat_lines)) return(NULL)

    toks <- stringr::str_split(stringr::str_trim(dat_lines), "\\s+")
    toks <- toks[lengths(toks) == length(hdr)]
    if (!length(toks)) return(NULL)

    mat <- do.call(rbind, toks)
    df <- as.data.frame(mat, stringsAsFactors = FALSE)
    names(df) <- hdr

    for (nm in names(df)) {
      df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
    }

    if (!"ITERATION" %in% names(df)) return(NULL)
    df
  })

  parsed_blocks <- purrr::compact(parsed_blocks)
  if (!length(parsed_blocks)) return(NULL)

  # Prefer the last block that has the NONMEM final row.
  has_final <- vapply(
    parsed_blocks,
    function(z) any(z$ITERATION == -1000000000),
    logical(1)
  )
  df <- if (any(has_final)) {
    parsed_blocks[[utils::tail(which(has_final), 1)]]
  } else {
    parsed_blocks[[length(parsed_blocks)]]
  }

  par_cols <- setdiff(names(df), c("ITERATION", "OBJ"))

  row_final <- df[df$ITERATION == -1000000000, , drop = FALSE]
  if (!nrow(row_final)) {
    # Fallback: last non-negative iteration.
    pos <- df[df$ITERATION >= 0, , drop = FALSE]
    if (nrow(pos)) row_final <- pos[nrow(pos), , drop = FALSE]
  }

  if (!nrow(row_final)) return(NULL)

  row_se <- df[df$ITERATION == -1000000001, , drop = FALSE]
  row_fix <- df[df$ITERATION == -1000000006, , drop = FALSE]

  values <- vapply(par_cols, function(nm) row_final[[nm]][1], numeric(1))

  ses <- stats::setNames(rep(NA_real_, length(par_cols)), par_cols)
  if (nrow(row_se)) {
    ses[] <- vapply(par_cols, function(nm) row_se[[nm]][1], numeric(1))
    # NONMEM uses 1e10 as a sentinel value for unavailable/FIX SE.
    ses[!is.finite(ses) | abs(ses) >= 1e9] <- NA_real_
  }

  fixed <- stats::setNames(rep(FALSE, length(par_cols)), par_cols)
  if (nrow(row_fix)) {
    fixed[] <- vapply(par_cols, function(nm) {
      z <- row_fix[[nm]][1]
      is.finite(z) && abs(z) > 0.5
    }, logical(1))
  } else {
    # Fallback when the FIX row is missing: sentinel SE in row -1000000001.
    if (nrow(row_se)) {
      fixed[] <- vapply(par_cols, function(nm) {
        z <- row_se[[nm]][1]
        is.finite(z) && abs(z) >= 1e9
      }, logical(1))
    }
  }

  rse <- stats::setNames(rep(NA_real_, length(par_cols)), par_cols)
  ok_rse <- !fixed & is.finite(values) & values != 0 & is.finite(ses)
  rse[ok_rse] <- abs(ses[ok_rse] / values[ok_rse]) * 100

  raw <- tibble::tibble(
    raw_name = par_cols,
    value = as.numeric(values),
    fixed = as.logical(fixed),
    se = as.numeric(ses),
    rse_pct = as.numeric(rse)
  )

  ofv <- if ("OBJ" %in% names(row_final)) row_final$OBJ[1] else NA_real_

  # Last real iteration.
  real_iter <- df$ITERATION[df$ITERATION >= 0]
  last_iteration <- if (length(real_iter)) max(real_iter, na.rm = TRUE) else NA_real_

  list(
    raw = raw,
    ofv = ofv,
    last_iteration = last_iteration,
    final_row_available = any(df$ITERATION == -1000000000),
    table = df
  )
}


# Convert raw EXT parameters to a tidy table and add derived quantities.
.nm_ext_params <- function(ext_obj, names_smr = NULL, semantic_map = NULL,
                           residual_map = NULL) {

  if (is.null(ext_obj) || is.null(ext_obj$raw) || !nrow(ext_obj$raw)) {
    return(.nm_empty_params())
  }

  raw <- ext_obj$raw

  # Optional names supplied by Wings.
  theta_names <- names_smr$theta %||% character(0)
  eta_names   <- names_smr$eta   %||% character(0)

  # Semantic names from the .ctl take priority over Wings names.
  theta_sem <- character(0); eta_sem <- character(0)
  if (!is.null(semantic_map) && nrow(semantic_map)) {
    mt <- semantic_map[semantic_map$element_type == "THETA" &
                         semantic_map$role == "typical_value" &
                         semantic_map$likely_pk_parameter, , drop = FALSE]
    # A THETA/ETA linked to several parameters (CL = THETA(2)*..., KE = CL/V)
    # is named after the closest one (smallest dependency depth); on ties,
    # the one matching the $THETA/$OMEGA comment label, then file order.
    pick <- function(m) {
      if (!nrow(m)) return(character(0))
      d <- if ("depth" %in% names(m)) m$depth else rep(0, nrow(m))
      lab <- if ("comment_matches_pk" %in% names(m)) !(m$comment_matches_pk %in% TRUE) else rep(TRUE, nrow(m))
      m <- m[order(m$index, d, lab, seq_len(nrow(m))), , drop = FALSE]
      m <- m[!duplicated(m$index), , drop = FALSE]
      stats::setNames(m$pk_parameter, as.character(m$index))
    }
    theta_sem <- pick(mt)
    me <- semantic_map[semantic_map$element_type == "ETA" &
                         semantic_map$likely_pk_parameter, , drop = FALSE]
    eta_sem <- pick(me)
  }
  # Covariate-effect THETAs (CL_WT...) and residual-error THETAs (PROP, ADD...)
  # take their $THETA comment label when it is a short name ("PROP",
  # "4 Prop (sigma_prop)" -> "Prop", "SD1 = proportional SD" -> "SD1"); a
  # free-text comment ("erreur proportionnelle"), a label already used by a
  # PK parameter or a duplicated label keeps the NONMEM name.
  theta_res <- character(0)
  lab_tab <- data.frame(theta = integer(0), label = character(0))
  if (!is.null(semantic_map) && nrow(semantic_map) && "label" %in% names(semantic_map)) {
    mc <- semantic_map[semantic_map$element_type == "THETA" &
                         semantic_map$role == "covariate_effect" &
                         !(as.character(semantic_map$index) %in% names(theta_sem)), , drop = FALSE]
    mc <- mc[!duplicated(mc$index), , drop = FALSE]
    lab_tab <- rbind(lab_tab, data.frame(theta = mc$index, label = as.character(mc$label)))
  }
  if (!is.null(residual_map) && nrow(residual_map))
    lab_tab <- rbind(lab_tab, data.frame(theta = residual_map$theta,
                                         label = as.character(residual_map$label)))
  lab_tab <- lab_tab[!duplicated(lab_tab$theta), , drop = FALSE]
  if (nrow(lab_tab)) {
    lab <- trimws(lab_tab$label)
    lab <- sub("^[0-9]+[.)]?\\s+", "", lab)          # leading index: "4 ", "6. "
    lab <- trimws(sub("[(=:;].*$", "", lab))           # text before "(", "=", ":"
    ok  <- !is.na(lab) & grepl("^[A-Za-z][A-Za-z0-9_.]{0,19}$", lab) &
      !(lab %in% theta_sem)
    ok  <- ok & !(lab %in% lab[ok][duplicated(lab[ok])])
    theta_res <- stats::setNames(lab[ok], as.character(lab_tab$theta[ok]))
  }
  eps_names   <- names_smr$err   %||% character(0)

  theta_label <- function(i) {
    si <- theta_sem[as.character(i)]
    ri <- theta_res[as.character(i)]
    if (length(si) && !is.na(si) && nzchar(si)) si
    else if (length(ri) && !is.na(ri) && nzchar(ri)) ri
    else if (i <= length(theta_names) && nzchar(theta_names[i])) theta_names[i]
    else paste0("THETA", i)
  }
  eta_label <- function(i) {
    si <- eta_sem[as.character(i)]
    if (length(si) && !is.na(si) && nzchar(si)) si
    else if (i <= length(eta_names) && nzchar(eta_names[i])) eta_names[i]
    else paste0("ETA", i)
  }
  eps_label <- function(i) {
    if (i <= length(eps_names) && nzchar(eps_names[i])) eps_names[i]
    else paste0("EPS", i)
  }

  out <- list()
  k <- 0L

  # THETA ---------------------------------------------------------------------
  th <- raw %>% dplyr::filter(stringr::str_detect(.data$raw_name, "^THETA\\d+$"))
  if (nrow(th)) {
    idx <- as.integer(stringr::str_extract(th$raw_name, "\\d+"))
    k <- k + 1L
    out[[k]] <- tibble::tibble(
      type = "THETA",
      index = idx,
      name = vapply(idx, theta_label, character(1)),
      raw_name = th$raw_name,
      value = th$value,
      fixed = th$fixed,
      rse_pct = th$rse_pct,
      se = th$se,
      shrinkage = NA_real_,
      source = "ext"
    )
  }

  # Matrices OMEGA / SIGMA -----------------------------------------------------
  parse_matrix <- function(prefix, diag_type, cov_type, sd_type, cor_type, labeller) {

    z <- raw %>%
      dplyr::filter(stringr::str_detect(
        .data$raw_name,
        paste0("^", prefix, "\\(\\d+,\\d+\\)$")
      ))

    if (!nrow(z)) return(NULL)

    ij <- stringr::str_match(z$raw_name, "\\((\\d+),(\\d+)\\)")
    ii <- as.integer(ij[, 2])
    jj <- as.integer(ij[, 3])

    diag_rows <- ii == jj
    cov_rows <- ii != jj

    res <- list()

    if (any(diag_rows)) {
      zd <- z[diag_rows, , drop = FALSE]
      id <- ii[diag_rows]

      raw_diag <- tibble::tibble(
        type = diag_type,
        index = id,
        name = vapply(id, labeller, character(1)),
        raw_name = zd$raw_name,
        value = zd$value,
        fixed = zd$fixed,
        rse_pct = zd$rse_pct,
        se = zd$se,
        shrinkage = NA_real_,
        source = "ext"
      )

      # Derived SD only when the variance is >= 0.
      sd_val <- ifelse(zd$value >= 0, sqrt(zd$value), NA_real_)

      derived_sd <- tibble::tibble(
        type = sd_type,
        index = id,
        name = vapply(id, labeller, character(1)),
        raw_name = zd$raw_name,
        value = sd_val,
        fixed = zd$fixed,
        # Derived transformation: no RSE is reported, so that it is not
        # mistaken for a direct NONMEM estimate.
        rse_pct = NA_real_,
        se = NA_real_,
        shrinkage = NA_real_,
        source = "derived_ext"
      )

      res <- c(res, list(raw_diag, derived_sd))
    }

    if (any(cov_rows)) {
      zc <- z[cov_rows, , drop = FALSE]
      ic <- ii[cov_rows]
      jc <- jj[cov_rows]

      pair_name <- paste0(
        vapply(ic, labeller, character(1)),
        " ~ ",
        vapply(jc, labeller, character(1))
      )

      cov_tbl <- tibble::tibble(
        type = cov_type,
        index = seq_len(nrow(zc)),
        name = pair_name,
        raw_name = zc$raw_name,
        value = zc$value,
        fixed = zc$fixed,
        rse_pct = zc$rse_pct,
        se = zc$se,
        shrinkage = NA_real_,
        source = "ext"
      )

      # Derived correlation.
      diag_map <- stats::setNames(
        z$value[diag_rows],
        as.character(ii[diag_rows])
      )

      cor_val <- vapply(seq_len(nrow(zc)), function(q) {
        vi <- diag_map[[as.character(ic[q])]]
        vj <- diag_map[[as.character(jc[q])]]
        cij <- zc$value[q]
        if (
          is.null(vi) || is.null(vj) ||
          !is.finite(vi) || !is.finite(vj) ||
          vi <= 0 || vj <= 0 || !is.finite(cij)
        ) return(NA_real_)
        cij / sqrt(vi * vj)
      }, numeric(1))

      cor_tbl <- tibble::tibble(
        type = cor_type,
        index = seq_len(nrow(zc)),
        name = pair_name,
        raw_name = zc$raw_name,
        value = cor_val,
        fixed = zc$fixed,
        rse_pct = NA_real_,
        se = NA_real_,
        shrinkage = NA_real_,
        source = "derived_ext"
      )

      res <- c(res, list(cov_tbl, cor_tbl))
    }

    dplyr::bind_rows(res)
  }

  om <- parse_matrix(
    "OMEGA", "OMEGA", "OMEGA_COV", "ETA_SD", "ETA_COR", eta_label
  )
  sg <- parse_matrix(
    "SIGMA", "SIGMA", "SIGMA_COV", "EPS_SD", "EPS_COR", eps_label
  )

  dplyr::bind_rows(out, om, sg)
}


# =============================================================================
# 7. .SMR PARSER (optional Wings file)
# =============================================================================

.nm_parse_smr <- function(path) {

  L <- .nm_read_lines_latin1(path)
  if (is.null(L) || !length(L)) return(NULL)

  vals_after <- function(pat) {
    i <- stringr::str_which(L, pat)
    if (!length(i)) return(character(0))
    rhs <- stringr::str_trim(stringr::str_replace(L[i[1]], pat, ""))
    toks <- stringr::str_split(rhs, "\\s+")[[1]]
    toks[nzchar(toks)]
  }

  pnum <- function(x) suppressWarnings(as.numeric(
    stringr::str_remove(x, "(?i)[a-z%]+$")
  ))
  pfix <- function(x) stringr::str_detect(x, "(?i)c$|fix")

  theta_nm <- vals_after("(?i)^\\s*THETA:(?!se%)\\s*")
  eta_nm   <- vals_after("(?i)^\\s*ETA:\\s*")
  err_nm   <- vals_after("(?i)^\\s*ERR:\\s*")

  theta_v  <- vals_after("(?i)^\\s*THETA\\s*=\\s*")
  etasd_v  <- vals_after("(?i)^\\s*ETASD\\s*=\\s*")
  epssd_v  <- vals_after("(?i)^\\s*EPSSD\\s*=\\s*")

  theta_se <- vals_after("(?i)^\\s*THETA:se%\\s*=\\s*")
  etasd_se <- vals_after("(?i)^\\s*ETASD:se%\\s*=\\s*")
  epssd_se <- vals_after("(?i)^\\s*EPSSD:se%\\s*=\\s*")

  mk <- function(type, nms, vv, ss) {
    n <- length(vv)
    if (!n) return(NULL)
    if (length(nms) != n) nms <- paste0(type, seq_len(n))
    rse <- if (length(ss) == n) pnum(ss) else rep(NA_real_, n)

    tibble::tibble(
      type = type,
      index = seq_len(n),
      name = nms,
      value = pnum(vv),
      fixed = pfix(vv),
      rse_pct = rse,
      se = rse * abs(pnum(vv)) / 100
    )
  }

  params <- dplyr::bind_rows(
    mk("THETA", theta_nm, theta_v, theta_se),
    mk("ETA_SD", eta_nm, etasd_v, etasd_se),
    mk("EPS_SD", err_nm, epssd_v, epssd_se)
  )

  joined <- paste(L, collapse = " ")
  n_obs <- suppressWarnings(as.integer(
    stringr::str_match(joined, "(?i)obs\\s*=\\s*(\\d+)")[, 2]
  ))

  # Fallback OFV from the SMR.
  ofv <- NA_real_
  ie <- stringr::str_which(L, "(?i)^\\s*ERR:\\s*")
  if (length(ie) && ie[1] < length(L)) {
    num <- stringr::str_extract(
      L[ie[1] + 1],
      "-?\\d+\\.?\\d*([eE][+-]?\\d+)?"
    )
    if (!is.na(num)) ofv <- as.numeric(num)
  }

  status <- stringr::str_subset(L, "(?i)MINIMIZATION|OPTIMIZATION")

  list(
    params = params,
    names = list(theta = theta_nm, eta = eta_nm, err = err_nm),
    n_obs = n_obs,
    ofv = ofv,
    status = if (length(status)) stringr::str_trim(status[1]) else NA_character_
  )
}


# =============================================================================
# 8. .LST PARSER
# =============================================================================

.nm_parse_lst <- function(path) {

  L <- .nm_read_lines_latin1(path)
  if (is.null(L)) return(NULL)

  # OFV
  ofv <- NA_real_

  i <- stringr::str_which(
    L,
    "(?i)OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT"
  )
  if (length(i)) {
    for (k in i[1]:min(i[1] + 2L, length(L))) {
      num <- stringr::str_extract(L[k], "-?\\d+\\.\\d+")
      if (!is.na(num)) {
        ofv <- as.numeric(num)
        break
      }
    }
  }

  if (is.na(ofv)) {
    i <- stringr::str_which(L, "#OBJV")
    if (length(i)) {
      ofv <- suppressWarnings(as.numeric(
        stringr::str_extract(L[i[length(i)]], "-?\\d+\\.\\d+")
      ))
    }
  }

  if (is.na(ofv)) {
    i <- stringr::str_which(L, "(?i)MINIMUM VALUE OF OBJECTIVE FUNCTION")
    if (length(i)) {
      for (k in i[1]:min(i[1] + 3L, length(L))) {
        num <- stringr::str_extract(L[k], "-?\\d+\\.\\d+")
        if (!is.na(num)) {
          ofv <- as.numeric(num)
          break
        }
      }
    }
  }

  # Minimisation status
  spat <- paste0(
    "(?i)MINIMIZATION (SUCCESSFUL|TERMINATED|ABORTED)",
    "|OPTIMIZATION WAS (COMPLETED|NOT COMPLETED)",
    "|REDUCED STOCHASTIC"
  )

  si <- stringr::str_which(L, spat)

  min_status <- NA_character_
  min_ok <- NA
  status_comment <- NA_character_

  if (length(si)) {
    ln <- L[si[1]]
    min_status <- stringr::str_trim(stringr::str_remove(ln, "^[01+]"))
    min_ok <- stringr::str_detect(ln, "(?i)SUCCESSFUL|COMPLETED") &&
      !stringr::str_detect(ln, "(?i)NOT COMPLETED")

    stop_pat <- "(?i)^\\s*(NO\\. OF|REFERENCE|ITERATION|MONITORING|OBJECTIVE|#|\\$)"
    j <- si[1] + 1L
    buf <- character(0)

    while (j <= length(L)) {
      cur <- L[j]
      if (!nzchar(stringr::str_trim(cur))) break
      if (stringr::str_detect(cur, "^[01]")) break
      if (stringr::str_detect(cur, stop_pat)) break
      buf <- c(buf, stringr::str_trim(cur))
      j <- j + 1L
      if (length(buf) > 8L) break
    }

    if (length(buf)) status_comment <- paste(buf, collapse = " ")
  }

  # Warnings
  flags <- c(
    "ROUNDING ERRORS"       = "(?i)DUE TO ROUNDING ERRORS|ROUNDING ERRORS \\(ERROR",
    "NEAR BOUNDARY"         = "(?i)NEAR ITS BOUNDARY",
    "R MATRIX SINGULAR"     = "(?i)R MATRIX.*SINGULAR",
    "S MATRIX SINGULAR"     = "(?i)S MATRIX.*SINGULAR",
    "COV EIGENVALUE"        = "(?i)EIGENVALUE OF COV",
    "NON-POSITIVE HESSIAN"  = "(?i)NON-?POSITIVE",
    "LARGE STANDARD ERRORS" = "(?i)LARGE STANDARD ERRORS",
    "ZERO GRADIENT"         = "(?i)ZERO GRADIENT",
    "PROGRAM TERMINATED"    = "(?i)PROGRAM TERMINATED"
  )

  w <- names(flags)[vapply(
    flags,
    function(p) any(stringr::str_detect(L, p)),
    logical(1)
  )]
  warnings <- if (length(w)) paste(w, collapse = " | ") else ""

  # Covariance step status.
  # On the NONMEM/Wings outputs tested:
  #   - if $COV is not requested, covariance is not mentioned at all;
  #   - "COVARIANCE STEP OMITTED: NO" means the step was requested;
  #   - "STANDARD ERROR OF ESTIMATE" confirms that it completed.
  cov_status <- dplyr::case_when(
    any(stringr::str_detect(
      L,
      "(?i)COVARIANCE STEP ABORTED"
    )) ~ "ABORTED",

    any(stringr::str_detect(
      L,
      "(?i)STANDARD ERROR OF ESTIMATE"
    )) ~ "COMPLETED",

    any(stringr::str_detect(
      L,
      "(?i)COVARIANCE STEP OMITTED:\\s*NO"
    )) ~ "REQUESTED",

    any(stringr::str_detect(
      L,
      "(?i)COVARIANCE STEP OMITTED:\\s*YES"
    )) ~ "OMITTED",

    TRUE ~ "NOT_REQUESTED"
  )

  list(
    ofv = ofv,
    min_status = min_status,
    min_ok = min_ok,
    status_comment = status_comment,
    warnings = warnings,
    cov_status = cov_status,
    eta_shrink = .nm_lst_shrink(L, "ETA.?SHRINK"),
    eps_shrink = .nm_lst_shrink(L, "EPS.?SHRINK")
  )
}


.nm_lst_shrink <- function(L, key) {
  idx <- stringr::str_which(L, paste0("(?i)", key))
  if (!length(idx)) return(numeric(0))

  nums <- unlist(lapply(idx, function(i) {
    suppressWarnings(as.numeric(
      stringr::str_extract_all(
        L[i],
        "-?\\d+\\.?\\d+([eE][+-]?\\d+)?"
      )[[1]]
    ))
  }))

  nums[!is.na(nums)]
}


# =============================================================================
# 9. .SHK PARSER
# =============================================================================

.nm_parse_shk <- function(path) {

  L <- .nm_read_lines_latin1(path)
  if (is.null(L)) return(list(eta = numeric(0), eps = numeric(0)))

  tb <- stringr::str_which(L, "(?i)^TABLE NO")
  if (length(tb)) L <- L[tb[length(tb)]:length(L)]

  row_of_type <- function(t) {
    rx <- sprintf("^\\s*%d\\s+", t)
    i <- stringr::str_which(L, rx)
    if (!length(i)) return(numeric(0))

    v <- suppressWarnings(as.numeric(
      stringr::str_split(stringr::str_trim(L[i[1]]), "\\s+")[[1]]
    ))
    v[-1]
  }

  eta <- row_of_type(1)
  eps <- row_of_type(5)
  if (!length(eps)) eps <- row_of_type(3)

  list(eta = eta, eps = eps)
}


# =============================================================================
# 10. SINGLE-RUN EXTRACTION
# =============================================================================

.nm_extract_one <- function(run, dir, ctl, smr, lst, ext, shk) {

  ctl_available <- !is.na(ctl) && file.exists(ctl)

  ext_available <- !is.na(ext) && file.exists(ext)
  smr_available <- !is.na(smr) && file.exists(smr)
  lst_available <- !is.na(lst) && file.exists(lst)
  shk_available <- !is.na(shk) && file.exists(shk)

  e <- if (ext_available) .nm_parse_ext(ext) else NULL
  s <- if (smr_available) .nm_parse_smr(smr) else NULL
  l <- if (lst_available) .nm_parse_lst(lst) else NULL

  # Semantic identification from the .ctl (nm_identify_parameters()).
  # If parsing fails, names fall back to the .smr / NONMEM defaults.
  semantic <- NULL
  if (ctl_available) {
    semantic <- tryCatch(nm_identify_parameters(ctl, likely_pk_only = FALSE),
                         error = function(e) NULL)
  }
  semantic_map <- if (!is.null(semantic)) semantic$mapping_detail_all else NULL
  residual_map <- if (!is.null(semantic)) semantic$residual_error else NULL

  # Numeric source: EXT first.
  if (!is.null(e)) {
    params <- .nm_ext_params(e, names_smr = s$names %||% NULL,
                             semantic_map = semantic_map,
                             residual_map = residual_map)
  } else if (!is.null(s) && nrow(s$params)) {
    # Wings fallback when EXT is missing or unusable.
    params <- s$params %>%
      dplyr::mutate(
        raw_name = paste0(.data$type, .data$index),
        shrinkage = NA_real_,
        source = "smr"
      ) %>%
      dplyr::select(dplyr::all_of(c(
        "type", "index", "name", "raw_name",
        "value", "fixed", "rse_pct", "se",
        "shrinkage", "source"
      )))
  } else {
    params <- .nm_empty_params()
  }

  # ETA shrinkage: LST, then SHK.
  eta_shr <- if (!is.null(l)) l$eta_shrink else numeric(0)
  eps_shr <- if (!is.null(l)) l$eps_shrink else numeric(0)

  if (!length(eta_shr) && shk_available) {
    sh <- .nm_parse_shk(shk)
    eta_shr <- sh$eta
    eps_shr <- sh$eps
  }

  if (nrow(params)) {
    eta_rows <- which(params$type == "ETA_SD")
    if (length(eta_rows) && length(eta_shr)) {
      n <- min(length(eta_rows), length(eta_shr))
      params$shrinkage[eta_rows[seq_len(n)]] <- eta_shr[seq_len(n)]
    }

    eps_rows <- which(params$type == "EPS_SD")
    if (length(eps_rows) && length(eps_shr)) {
      n <- min(length(eps_rows), length(eps_shr))
      params$shrinkage[eps_rows[seq_len(n)]] <- eps_shr[seq_len(n)]
    }
  }

  params <- dplyr::mutate(params, run = run, .before = 1)

  # OFV: EXT first, LST as check/fallback, SMR as last fallback.
  ext_ofv <- e$ofv %||% NA_real_
  lst_ofv <- l$ofv %||% NA_real_
  smr_ofv <- s$ofv %||% NA_real_

  ofv <- if (!is.na(ext_ofv)) ext_ofv
  else if (!is.na(lst_ofv)) lst_ofv
  else smr_ofv

  ofv_diff_ext_lst <- if (
    is.finite(ext_ofv) && is.finite(lst_ofv)
  ) ext_ofv - lst_ofv else NA_real_

  min_status <- l$min_status %||%
    (if (!is.null(s)) s$status else NA_character_)
  min_ok <- l$min_ok %||% NA
  status_comment <- l$status_comment %||% NA_character_
  warnings <- l$warnings %||% NA_character_
  cov_status <- l$cov_status %||% NA_character_
  cov_ok <- !is.na(cov_status) && cov_status == "COMPLETED"

  estimation_available <- nrow(params) > 0L

  # Count of parameters actually estimated by NONMEM:
  # THETA + OMEGA/OMEGA_COV + SIGMA/SIGMA_COV only, never derived ones.
  raw_types <- c("THETA", "OMEGA", "OMEGA_COV", "SIGMA", "SIGMA_COV")

  n_theta <- sum(params$type == "THETA")
  n_omega <- sum(params$type %in% c("OMEGA", "OMEGA_COV"))
  n_sigma <- sum(params$type %in% c("SIGMA", "SIGMA_COV"))
  n_eta <- sum(params$type == "ETA_SD")
  n_eps <- sum(params$type == "EPS_SD")

  par_mask <- params$type %in% raw_types
  n_par <- if (any(par_mask)) {
    sum(!params$fixed[par_mask], na.rm = TRUE)
  } else {
    # Fallback for the old SMR format.
    sum(!params$fixed, na.rm = TRUE)
  }

  n_obs <- if (!is.null(s)) s$n_obs else NA_integer_

  aic <- if (is.finite(ofv) && n_par > 0) {
    ofv + 2 * n_par
  } else NA_real_

  bic <- if (
    is.finite(ofv) && n_par > 0 &&
    !is.na(n_obs) && n_obs > 0
  ) {
    ofv + n_par * log(n_obs)
  } else NA_real_

  max_shr <- if (any(!is.na(params$shrinkage))) {
    max(params$shrinkage, na.rm = TRUE)
  } else NA_real_

  # RSE > 50% on estimated raw NONMEM parameters only.
  rse_mask <- par_mask & !params$fixed & !is.na(params$rse_pct)
  pct_rse <- if (any(rse_mask)) {
    round(mean(params$rse_pct[rse_mask] > 50) * 100, 1)
  } else NA_real_

  run_row <- tibble::tibble(
    run = run,
    ofv = ofv,
    ext_ofv = ext_ofv,
    lst_ofv = lst_ofv,
    ofv_diff_ext_lst = ofv_diff_ext_lst,

    min_status = min_status,
    min_ok = min_ok,
    status_comment = status_comment,
    warnings = warnings,
    cov_status = cov_status,
    cov_ok = cov_ok,

    ctl_available = ctl_available,
    semantic_identification = !is.null(semantic),
    ext_available = ext_available,
    smr_available = smr_available,
    lst_available = lst_available,
    shk_available = shk_available,
    estimation_available = estimation_available,
    ext_final_row = e$final_row_available %||% FALSE,
    ext_last_iteration = e$last_iteration %||% NA_real_,

    n_theta = n_theta,
    n_omega = n_omega,
    n_sigma = n_sigma,
    n_eta = n_eta,
    n_eps = n_eps,
    n_par = n_par,
    n_obs = n_obs,

    aic = round(aic, 3),
    bic = round(bic, 3),

    max_eta_shrink = max_shr,
    pct_rse_gt50 = pct_rse,

    dir = dir
  )

  list(run = run_row, params = params)
}
