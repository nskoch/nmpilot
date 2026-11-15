# =============================================================================
# nm_extract_run_table.R
# -----------------------------------------------------------------------------
# Read a NONMEM .ext file directly and export a parameter table to Word.
#
# - Numeric source: .ext only
# - Automatic detection of the run directory (.nm7, .reg, <run>/ or flat)
# - THETA: direct estimate
# - Diagonal OMEGA: also shown as ETA standard deviations = sqrt(OMEGA)
# - Diagonal SIGMA: also shown as EPS standard deviations = sqrt(SIGMA)
# - SE / RSE / 95% CI:
#     * THETA: directly from the .ext
#     * ETA/EPS SD: delta method from the .ext variances
# - FIX parameters detected from the NONMEM row -1000000006
#
# The .ext parser follows the one of nm_extract_runs.R.
# Word export requires the suggested packages flextable and officer.
# =============================================================================

# =============================================================================
# 1. FIND THE RUN DIRECTORY (.nm7 OR .reg)
# =============================================================================

.nm_rt_find_run_dir <- function(run_rep,
                                run_name,
                                extensions = c("nm7", "reg")) {

  # Run subdirectories (<run>.nm7, <run>.reg), a folder named <run>, then the
  # parent directory itself (flat layout: outputs next to the control stream)
  candidate_dirs <- c(
    file.path(run_rep, paste0(run_name, ".", extensions)),
    file.path(run_rep, run_name),
    run_rep
  )
  extensions <- c(extensions, "dir", "flat")

  candidate_ext <- file.path(
    candidate_dirs,
    paste0(run_name, ".ext")
  )

  exists_ext <- file.exists(candidate_ext)

  if (!any(exists_ext)) {
    stop(
      paste0(
        "No .ext file found for run '", run_name, "'.\n\n",
        "Paths tried:\n - ",
        paste(candidate_ext, collapse = "\n - ")
      ),
      call. = FALSE
    )
  }

  if (sum(exists_ext) > 1) {
    warning(
      paste0(
        "Several .ext files found for '", run_name, "'. ",
        "The first one in the order ",
        paste(extensions, collapse = " then "),
        " will be used."
      ),
      call. = FALSE
    )
  }

  idx <- which(exists_ext)[1]

  list(
    run_dir = candidate_dirs[idx],
    ext_file = candidate_ext[idx],
    extension = extensions[idx]
  )
}


# =============================================================================
# 2. ROBUST FILE READING
# =============================================================================

.nm_rt_read_lines <- function(path) {

  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  # NONMEM / Wings on Windows: latin1 / CP1252 is common
  con <- file(path, encoding = "latin1")
  on.exit(close(con))

  lines <- tryCatch(
    readLines(con, warn = FALSE),
    error = function(e) NULL
  )

  if (is.null(lines) || !length(lines)) {
    stop("Cannot read the .ext file: ", path, call. = FALSE)
  }

  suppressWarnings(
    iconv(lines, from = "latin1", to = "UTF-8", sub = "?")
  )
}


# =============================================================================
# 3. PARSE THE .EXT
# =============================================================================
#
# NONMEM uses in particular:
#   -1000000000 : final estimates
#   -1000000001 : standard errors
#   -1000000006 : FIX indicators
#
# =============================================================================

.nm_rt_parse_ext <- function(file_path) {

  L <- .nm_rt_read_lines(file_path)

  header_idx <- stringr::str_which(
    L,
    "^\\s*ITERATION\\s+"
  )

  if (!length(header_idx)) {
    stop(
      "No block starting with ITERATION found in the .ext file.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Read all ITERATION blocks: a .ext can contain several tables
  # ---------------------------------------------------------------------------

  parsed_blocks <- lapply(header_idx, function(h) {

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

    hdr <- stringr::str_split(
      stringr::str_trim(L[h]),
      "\\s+"
    )[[1]]

    if (h >= e) return(NULL)

    dat_lines <- L[(h + 1L):e]
    dat_lines <- dat_lines[
      nzchar(stringr::str_trim(dat_lines))
    ]

    if (!length(dat_lines)) return(NULL)

    toks <- stringr::str_split(
      stringr::str_trim(dat_lines),
      "\\s+"
    )

    # Keep only the rows with exactly the number of columns
    # given by the header.
    toks <- toks[lengths(toks) == length(hdr)]

    if (!length(toks)) return(NULL)

    mat <- do.call(rbind, toks)

    df <- as.data.frame(
      mat,
      stringsAsFactors = FALSE
    )

    names(df) <- hdr

    for (nm in names(df)) {
      df[[nm]] <- suppressWarnings(
        as.numeric(df[[nm]])
      )
    }

    if (!"ITERATION" %in% names(df)) {
      return(NULL)
    }

    df
  })

  parsed_blocks <- Filter(Negate(is.null), parsed_blocks)

  if (!length(parsed_blocks)) {
    stop(
      "ITERATION blocks were found but none could be read.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Prefer the last block that has the NONMEM final row
  # ---------------------------------------------------------------------------

  has_final <- vapply(
    parsed_blocks,
    function(z) any(z$ITERATION == -1000000000),
    logical(1)
  )

  if (any(has_final)) {
    df <- parsed_blocks[[utils::tail(which(has_final), 1)]]
  } else {
    df <- parsed_blocks[[length(parsed_blocks)]]
  }

  par_cols <- setdiff(
    names(df),
    c("ITERATION", "OBJ")
  )

  if (!length(par_cols)) {
    stop(
      "No parameter column found in the .ext.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Final estimates
  # ---------------------------------------------------------------------------

  row_final <- df[
    df$ITERATION == -1000000000,
    ,
    drop = FALSE
  ]

  final_row_available <- nrow(row_final) > 0

  if (!nrow(row_final)) {

    # Fallback: last positive/real NONMEM iteration
    pos <- df[
      df$ITERATION >= 0,
      ,
      drop = FALSE
    ]

    if (nrow(pos)) {
      row_final <- pos[nrow(pos), , drop = FALSE]
    }
  }

  if (!nrow(row_final)) {
    stop(
      "Cannot find the final estimates in the .ext.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # SE and FIX
  # ---------------------------------------------------------------------------

  row_se <- df[
    df$ITERATION == -1000000001,
    ,
    drop = FALSE
  ]

  row_fix <- df[
    df$ITERATION == -1000000006,
    ,
    drop = FALSE
  ]

  values <- vapply(
    par_cols,
    function(nm) row_final[[nm]][1],
    numeric(1)
  )

  # SE
  ses <- stats::setNames(
    rep(NA_real_, length(par_cols)),
    par_cols
  )

  if (nrow(row_se)) {

    ses[] <- vapply(
      par_cols,
      function(nm) row_se[[nm]][1],
      numeric(1)
    )

    # NONMEM sentinel value for unavailable / FIX SE
    ses[
      !is.finite(ses) | abs(ses) >= 1e9
    ] <- NA_real_
  }

  # FIX
  fixed <- stats::setNames(
    rep(FALSE, length(par_cols)),
    par_cols
  )

  if (nrow(row_fix)) {

    fixed[] <- vapply(
      par_cols,
      function(nm) {
        z <- row_fix[[nm]][1]
        is.finite(z) && abs(z) > 0.5
      },
      logical(1)
    )

  } else if (nrow(row_se)) {

    # Fallback when the FIX row does not exist:
    # sentinel SE >= 1e9
    fixed[] <- vapply(
      par_cols,
      function(nm) {
        z <- row_se[[nm]][1]
        is.finite(z) && abs(z) >= 1e9
      },
      logical(1)
    )
  }

  # RSE on the native .ext scale
  rse <- stats::setNames(
    rep(NA_real_, length(par_cols)),
    par_cols
  )

  ok_rse <- (
    !fixed &
      is.finite(values) &
      values != 0 &
      is.finite(ses)
  )

  rse[ok_rse] <- abs(
    ses[ok_rse] / values[ok_rse]
  ) * 100

  raw <- tibble::tibble(
    raw_name = par_cols,
    value = as.numeric(values),
    fixed = as.logical(fixed),
    se = as.numeric(ses),
    rse_pct = as.numeric(rse)
  )

  # ---------------------------------------------------------------------------
  # OFV
  # ---------------------------------------------------------------------------

  ofv <- if ("OBJ" %in% names(row_final)) {
    row_final$OBJ[1]
  } else {
    NA_real_
  }

  # Last real iteration
  real_iter <- df$ITERATION[
    is.finite(df$ITERATION) &
      df$ITERATION >= 0
  ]

  last_iteration <- if (length(real_iter)) {
    max(real_iter, na.rm = TRUE)
  } else {
    NA_real_
  }

  list(
    raw = raw,
    ofv = ofv,
    last_iteration = last_iteration,
    final_row_available = final_row_available,
    table = df
  )
}


# =============================================================================
# 4. CONVERT PARAMETERS FOR THE TABLE
# =============================================================================
#
# The table follows the layout of the former .smr-based script:
#
# THETAi
# ETAi SD = sqrt(OMEGA(i,i))
# EPSi SD = sqrt(SIGMA(i,i))
#
# For a variance V with SD = sqrt(V), the delta method gives:
#   SE(SD) ~ SE(V) / [2 sqrt(V)]
#
# =============================================================================

.nm_rt_param_table <- function(ext_obj) {

  raw <- ext_obj$raw

  if (is.null(raw) || !nrow(raw)) {
    stop("No parameter available in the .ext.", call. = FALSE)
  }

  out <- list()

  # THETA ---------------------------------------------------------------------
  th <- raw %>%
    dplyr::filter(stringr::str_detect(.data$raw_name, "^THETA\\d+$"))

  if (nrow(th)) {
    idx <- as.integer(stringr::str_extract(th$raw_name, "\\d+"))

    out[[length(out) + 1L]] <- tibble::tibble(
      Type = "THETA",
      Index = idx,
      Parameter = paste0("THETA", idx),
      Estimate = th$value,
      SE = th$se,
      RSE_pct = th$rse_pct,
      Fixed = th$fixed,
      Raw_name = th$raw_name,
      Derived = FALSE
    )
  }

  # OMEGA: native NONMEM values ----------------------------------------------
  zom <- raw %>%
    dplyr::filter(
      stringr::str_detect(.data$raw_name, "^OMEGA\\(\\d+,\\d+\\)$")
    )

  if (nrow(zom)) {

    ij <- stringr::str_match(zom$raw_name, "\\((\\d+),(\\d+)\\)")
    ii <- as.integer(ij[, 2])
    jj <- as.integer(ij[, 3])

    diag_rows <- ii == jj
    off_rows  <- ii != jj

    # Diagonal OMEGA: variances estimated directly by NONMEM
    if (any(diag_rows)) {
      zd <- zom[diag_rows, , drop = FALSE]
      id <- ii[diag_rows]

      out[[length(out) + 1L]] <- tibble::tibble(
        Type = "OMEGA",
        Index = id,
        Parameter = paste0("OMEGA(", id, ",", id, ")"),
        Estimate = zd$value,
        SE = zd$se,
        RSE_pct = zd$rse_pct,
        Fixed = zd$fixed,
        Raw_name = zd$raw_name,
        Derived = FALSE
      )
    }

    # Off-diagonal OMEGA covariances: estimated directly by NONMEM
    if (any(off_rows)) {
      zc <- zom[off_rows, , drop = FALSE]
      ic <- ii[off_rows]
      jc <- jj[off_rows]

      out[[length(out) + 1L]] <- tibble::tibble(
        Type = "OMEGA_COV",
        Index = seq_len(nrow(zc)),
        Parameter = paste0("Cov(ETA", ic, ", ETA", jc, ")"),
        Estimate = zc$value,
        SE = zc$se,
        RSE_pct = zc$rse_pct,
        Fixed = zc$fixed,
        Raw_name = zc$raw_name,
        Derived = FALSE
      )
    }

    # ETA SD: equivalent of Wings ETASD = sqrt(OMEGA(i,i))
    if (any(diag_rows)) {
      zd <- zom[diag_rows, , drop = FALSE]
      id <- ii[diag_rows]

      eta_sd <- ifelse(
        is.finite(zd$value) & zd$value >= 0,
        sqrt(zd$value),
        NA_real_
      )

      # Delta method for the SE of sqrt(OMEGA)
      eta_sd_se <- ifelse(
        !zd$fixed & is.finite(zd$se) & is.finite(eta_sd) & eta_sd > 0,
        abs(zd$se) / (2 * eta_sd),
        NA_real_
      )

      eta_sd_rse <- ifelse(
        !zd$fixed & is.finite(eta_sd) & eta_sd > 0 & is.finite(eta_sd_se),
        abs(eta_sd_se / eta_sd) * 100,
        NA_real_
      )

      out[[length(out) + 1L]] <- tibble::tibble(
        Type = "ETA_SD",
        Index = id,
        Parameter = paste0("ETA", id, " SD (Wings)"),
        Estimate = eta_sd,
        SE = eta_sd_se,
        RSE_pct = eta_sd_rse,
        Fixed = zd$fixed,
        Raw_name = zd$raw_name,
        Derived = TRUE
      )
    }

    # ETA correlations derived from the OMEGA matrix
    if (any(off_rows)) {
      zc <- zom[off_rows, , drop = FALSE]
      ic <- ii[off_rows]
      jc <- jj[off_rows]

      diag_map <- stats::setNames(
        zom$value[diag_rows],
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

      out[[length(out) + 1L]] <- tibble::tibble(
        Type = "ETA_COR",
        Index = seq_len(nrow(zc)),
        Parameter = paste0("Cor(ETA", ic, ", ETA", jc, ")"),
        Estimate = cor_val,
        SE = NA_real_,
        RSE_pct = NA_real_,
        Fixed = FALSE,
        Raw_name = zc$raw_name,
        Derived = TRUE
      )
    }
  }

  # SIGMA: native NONMEM values + derived EPS SD ----------------------------
  zsg <- raw %>%
    dplyr::filter(
      stringr::str_detect(.data$raw_name, "^SIGMA\\(\\d+,\\d+\\)$")
    )

  if (nrow(zsg)) {

    ij <- stringr::str_match(zsg$raw_name, "\\((\\d+),(\\d+)\\)")
    ii <- as.integer(ij[, 2])
    jj <- as.integer(ij[, 3])
    diag_rows <- ii == jj

    if (any(diag_rows)) {
      zd <- zsg[diag_rows, , drop = FALSE]
      id <- ii[diag_rows]

      out[[length(out) + 1L]] <- tibble::tibble(
        Type = "SIGMA",
        Index = id,
        Parameter = paste0("SIGMA(", id, ",", id, ")"),
        Estimate = zd$value,
        SE = zd$se,
        RSE_pct = zd$rse_pct,
        Fixed = zd$fixed,
        Raw_name = zd$raw_name,
        Derived = FALSE
      )

      eps_sd <- ifelse(
        is.finite(zd$value) & zd$value >= 0,
        sqrt(zd$value),
        NA_real_
      )

      eps_sd_se <- ifelse(
        !zd$fixed & is.finite(zd$se) & is.finite(eps_sd) & eps_sd > 0,
        abs(zd$se) / (2 * eps_sd),
        NA_real_
      )

      eps_sd_rse <- ifelse(
        !zd$fixed & is.finite(eps_sd) & eps_sd > 0 & is.finite(eps_sd_se),
        abs(eps_sd_se / eps_sd) * 100,
        NA_real_
      )

      out[[length(out) + 1L]] <- tibble::tibble(
        Type = "EPS_SD",
        Index = id,
        Parameter = paste0("EPS", id, " SD (Wings)"),
        Estimate = eps_sd,
        SE = eps_sd_se,
        RSE_pct = eps_sd_rse,
        Fixed = zd$fixed,
        Raw_name = zd$raw_name,
        Derived = TRUE
      )
    }
  }

  params <- dplyr::bind_rows(out)

  if (!nrow(params)) {
    stop("No THETA, OMEGA or SIGMA recognised.", call. = FALSE)
  }

  params <- params %>%
    dplyr::mutate(
      CI95_lower = dplyr::if_else(
        !.data$Fixed & is.finite(.data$SE),
        .data$Estimate - 1.96 * .data$SE,
        NA_real_
      ),
      CI95_upper = dplyr::if_else(
        !.data$Fixed & is.finite(.data$SE),
        .data$Estimate + 1.96 * .data$SE,
        NA_real_
      ),
      CI95 = dplyr::case_when(
        .data$Fixed ~ "FIX",
        is.finite(.data$CI95_lower) & is.finite(.data$CI95_upper) ~
          paste0(
            "[",
            sprintf("%.4f", .data$CI95_lower),
            " ; ",
            sprintf("%.4f", .data$CI95_upper),
            "]"
          ),
        TRUE ~ NA_character_
      ),
      RSE = dplyr::case_when(
        .data$Fixed ~ "FIX",
        is.finite(.data$RSE_pct) ~
          paste0(round(.data$RSE_pct, 1), "%"),
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::arrange(
      match(
        .data$Type,
        c("THETA", "OMEGA", "OMEGA_COV", "ETA_SD", "ETA_COR", "SIGMA", "EPS_SD")
      ),
      .data$Index
    )

  params
}


# =============================================================================
# 5. BUILD THE FLEXTABLE
# =============================================================================

.nm_rt_flextable <- function(df,
                                       ofv,
                                       last_iteration = NA_real_,
                                       final_row_available = FALSE) {

  df_display <- df %>%
    dplyr::select(
      "Parameter",
      "Estimate",
      "CI95",
      "RSE"
    )

  ft <- flextable::flextable(df_display)

  ft <- flextable::colformat_double(
    ft,
    j = "Estimate",
    digits = 4,
    na_str = ""
  )

  ft <- flextable::set_header_labels(
    ft,
    Parameter = "Parameter",
    Estimate = "Estimate",
    CI95 = "95% CI",
    RSE = "RSE"
  )

  ft <- flextable::width(
    ft,
    j = "Parameter",
    width = 1.8
  )

  ft <- flextable::width(
    ft,
    j = "Estimate",
    width = 1.2
  )

  ft <- flextable::width(
    ft,
    j = "CI95",
    width = 2.2
  )

  ft <- flextable::width(
    ft,
    j = "RSE",
    width = 1.0
  )

  # Header
  ft <- flextable::bold(
    ft,
    part = "header"
  )

  ft <- flextable::bg(
    ft,
    part = "header",
    bg = "#4472C4"
  )

  ft <- flextable::color(
    ft,
    part = "header",
    color = "white"
  )

  ft <- flextable::align(
    ft,
    part = "header",
    align = "center"
  )

  # Body
  if (nrow(df_display) > 0) {

    ft <- flextable::bg(
      ft,
      i = seq(1, nrow(df_display), 2),
      bg = "#F2F2F2",
      part = "body"
    )
  }

  ft <- flextable::align(
    ft,
    j = c("Estimate", "CI95", "RSE"),
    align = "center",
    part = "body"
  )

  ft <- flextable::align(
    ft,
    j = "Parameter",
    align = "left",
    part = "body"
  )

  # Borders
  ft <- flextable::border_outer(
    ft,
    border = officer::fp_border(
      width = 2,
      color = "#4472C4"
    )
  )

  ft <- flextable::border_inner_h(
    ft,
    border = officer::fp_border(
      width = 1,
      color = "#CCCCCC"
    )
  )

  ft <- flextable::border_inner_v(
    ft,
    border = officer::fp_border(
      width = 1,
      color = "#CCCCCC"
    )
  )

  # Table footer
  ofv_txt <- if (is.finite(ofv)) {
    paste0(
      "Objective Function Value (OFV): ",
      sprintf("%.3f", ofv)
    )
  } else {
    "Objective Function Value (OFV): not available"
  }

  iter_txt <- if (is.finite(last_iteration)) {
    paste0(
      "Last estimation iteration: ",
      last_iteration
    )
  } else {
    "Last estimation iteration: not available"
  }

  final_txt <- if (isTRUE(final_row_available)) {
    "Final NONMEM estimation row (-1000000000): available"
  } else {
    "Final NONMEM estimation row (-1000000000): absent; last positive iteration used"
  }

  ft <- flextable::add_footer_lines(
    ft,
    values = c(
      ofv_txt,
      iter_txt,
      final_txt
    )
  )

  ft <- flextable::align(
    ft,
    part = "footer",
    align = "left"
  )

  ft <- flextable::italic(
    ft,
    part = "footer"
  )

  ft <- flextable::fontsize(
    ft,
    part = "footer",
    size = 9
  )

  ft
}


# =============================================================================
# 6. EXPORT A .EXT FILE TO WORD
# =============================================================================

.nm_rt_export_word <- function(file_path,
                               output_file = "nonmem_results.docx",
                               run_name = NULL) {

  message("Reading .ext file: ", file_path)

  ext_obj <- .nm_rt_parse_ext(file_path)

  message("Extracting parameters...")

  params <- .nm_rt_param_table(ext_obj)

  message("Building the flextable...")

  ft <- .nm_rt_flextable(
    df = params,
    ofv = ext_obj$ofv,
    last_iteration = ext_obj$last_iteration,
    final_row_available = ext_obj$final_row_available
  )

  message("Building the Word document...")

  doc <- officer::read_docx()

  title <- if (!is.null(run_name)) {
    paste0("popPK Results - ", run_name)
  } else {
    "popPK Results"
  }

  doc <- officer::body_add_par(
    doc,
    title,
    style = "heading 1"
  )

  doc <- officer::body_add_par(
    doc,
    paste(
      "Generated on:",
      format(Sys.Date(), "%d/%m/%Y")
    ),
    style = "Normal"
  )

  doc <- officer::body_add_par(
    doc,
    ""
  )

  doc <- flextable::body_add_flextable(
    doc,
    value = ft
  )

  # ---------------------------------------------------------------------------
  # Legend / definitions below the table
  # ---------------------------------------------------------------------------

  doc <- officer::body_add_par(doc, "")
  doc <- officer::body_add_par(
    doc,
    "Definitions and formulas",
    style = "heading 2"
  )

  doc <- officer::body_add_par(
    doc,
    "OMEGA(i,i) (NONMEM): variance of ETAi, estimated directly by NONMEM."
  )

  doc <- officer::body_add_par(
    doc,
    "ETAi SD (Wings): standard deviation of ETAi, corresponding to ETASD in Wings."
  )

  doc <- officer::body_add_par(
    doc,
    "ETAi SD (Wings) = sqrt[OMEGA(i,i)]"
  )

  doc <- officer::body_add_par(
    doc,
    "Cov(ETAi, ETAj): covariance between ETAi and ETAj. For an OMEGA BLOCK, this is the off-diagonal OMEGA(i,j) element estimated by NONMEM."
  )

  doc <- officer::body_add_par(
    doc,
    "Cov(ETAi, ETAj) = OMEGA(i,j)"
  )

  doc <- officer::body_add_par(
    doc,
    "Cor(ETAi, ETAj): correlation derived from the OMEGA covariance matrix."
  )

  doc <- officer::body_add_par(
    doc,
    "Cor(ETAi, ETAj) = OMEGA(i,j) / sqrt[OMEGA(i,i) x OMEGA(j,j)]"
  )

  doc <- officer::body_add_par(
    doc,
    "Note: with a diagonal OMEGA structure, off-diagonal covariances are constrained to zero and no ETA correlation is estimated. With BLOCK(n), NONMEM estimates the off-diagonal covariance terms; correlations shown here are calculated from the estimated OMEGA matrix."
  )

  print(
    doc,
    target = output_file
  )

  message("Table exported to: ", output_file)
  message("=== Summary ===")
  message("Parameters shown: ", nrow(params))
  if (is.finite(ext_obj$ofv)) {
    message("OFV: ", sprintf("%.3f", ext_obj$ofv))
  } else {
    message("OFV: not available")
  }
  message("Last iteration: ", ext_obj$last_iteration)
  message("Final row -1000000000: ",
          ifelse(ext_obj$final_row_available, "yes", "no"))

  invisible(
    list(
      flextable = ft,
      params = params,
      ext = ext_obj,
      output_file = output_file
    )
  )
}


# =============================================================================
# 7. CONVENIENCE: EXPORT ONE RUN DIRECTLY
# =============================================================================

.nm_rt_export_run <- function(run_rep,
                              run_name,
                              output_name = "nonmem_results.docx",
                              extensions = c("nm7", "reg")) {

  run_info <- .nm_rt_find_run_dir(
    run_rep = run_rep,
    run_name = run_name,
    extensions = extensions
  )

  message("Run directory detected: ", run_info$run_dir)

  output_file <- file.path(
    run_info$run_dir,
    output_name
  )

  res <- .nm_rt_export_word(
    file_path = run_info$ext_file,
    output_file = output_file,
    run_name = run_name
  )

  invisible(
    c(
      res,
      list(
        ext_file = run_info$ext_file,
        run_dir = run_info$run_dir,
        extension = run_info$extension
      )
    )
  )
}


# =============================================================================
# 8. MAIN FUNCTION
# =============================================================================

#' Parameter table of one NONMEM run, exported to Word
#'
#' Reads the `.ext` file of a run and writes a Word document with the
#' parameter estimates, 95% confidence intervals and RSE. ETA and EPS
#' standard deviations are derived from the variances (delta method for their
#' SE); ETA correlations are derived from the OMEGA matrix.
#'
#' Two call forms are accepted:
#' - `run_rep` (parent directory) + `run_name`;
#' - `run_path`, the full path to the run directory (`.nm7`, `.reg`) or to
#'   the run stem (e.g. `C:/project/run1` for `C:/project/run1.ext`).
#'
#' The `.ext` file is looked for in `<run>.nm7/`, `<run>.reg/`, a folder named
#' `<run>`, then directly in `run_rep` (flat layout, as produced by a plain
#' nmfe call). Wings for NONMEM is not required.
#'
#' @param run_rep     Parent directory containing the run directory, or the
#'   run files themselves (flat layout).
#' @param run_name    Run name (stem of the `.ext` file). A full path given
#'   here by mistake is converted automatically.
#' @param run_path    Full path to the run (alternative to `run_rep` +
#'   `run_name`).
#' @param output_name Name of the Word file, written in the run directory
#'   (default `"nonmem_results.docx"`).
#' @param extensions  Run directory extensions to try first, in order
#'   (default `c("nm7", "reg")`), before the folder `<run>` and the flat
#'   layout.
#'
#' @return Invisibly, a list with `flextable`, `params` (parameter table),
#'   `ext` (parsed .ext), `output_file`, `ext_file`, `run_dir` and
#'   `extension`.
#'
#' @examples
#' \dontrun{
#' res <- nm_extract_run_table(run_rep = "C:/nm76g64/run/para",
#'                             run_name = "run_1_prop_kaBbloc")
#' res <- nm_extract_run_table(run_path = "C:/nm76g64/run/para/run_1_prop_kaBbloc.reg")
#' res$params
#' }
#' @export

nm_extract_run_table <- function(run_rep = NULL,
                                 run_name = NULL,
                                 run_path = NULL,
                                 output_name = "nonmem_results.docx",
                                 extensions = c("nm7", "reg")) {

  for (pkg in c("flextable", "officer")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for the Word export. ",
           "Install it with install.packages(\"", pkg, "\").", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # CASE 1: run_path supplied
  # ---------------------------------------------------------------------------

  if (!is.null(run_path)) {

    run_path <- trimws(run_path)

    if (!nzchar(run_path)) {
      stop("`run_path` is empty.", call. = FALSE)
    }

    # Normalise slashes
    run_path <- gsub("\\\\", "/", run_path)

    # Remove a trailing slash
    run_path <- sub("/+$", "", run_path)

    # A path ending with .nm7 or .reg is taken as the Wings run directory
    # itself.
    ext_pat <- paste0("\\.(", paste(extensions, collapse = "|"), ")$")

    if (grepl(ext_pat, run_path, ignore.case = TRUE)) {

      run_dir <- run_path
      run_base <- basename(run_dir)

      run_name2 <- sub(
        ext_pat,
        "",
        run_base,
        ignore.case = TRUE
      )

      run_rep2 <- dirname(run_dir)

    } else {

      # Path to the run stem, without .nm7/.reg
      run_name2 <- basename(run_path)
      run_rep2 <- dirname(run_path)
    }

    return(
      .nm_rt_export_run(
        run_rep = run_rep2,
        run_name = run_name2,
        output_name = output_name,
        extensions = extensions
      )
    )
  }

  # ---------------------------------------------------------------------------
  # CASE 2: run_rep + run_name
  # ---------------------------------------------------------------------------

  if (is.null(run_rep) || !nzchar(trimws(run_rep))) {
    stop(
      "Supply either `run_path`, or `run_rep` + `run_name`.",
      call. = FALSE
    )
  }

  if (is.null(run_name) || !nzchar(trimws(run_name))) {
    stop(
      "Supply either `run_path`, or `run_rep` + `run_name`.",
      call. = FALSE
    )
  }

  run_rep <- trimws(run_rep)
  run_name <- trimws(run_name)

  # Tolerance: if run_name accidentally contains a full path,
  # it is converted automatically.
  run_name_norm <- gsub("\\\\", "/", run_name)

  if (grepl("/", run_name_norm, fixed = TRUE)) {

    run_name_norm <- sub("/+$", "", run_name_norm)

    ext_pat <- paste0("\\.(", paste(extensions, collapse = "|"), ")$")

    if (grepl(ext_pat, run_name_norm, ignore.case = TRUE)) {
      run_name2 <- sub(
        ext_pat,
        "",
        basename(run_name_norm),
        ignore.case = TRUE
      )
      run_rep2 <- dirname(run_name_norm)
    } else {
      run_name2 <- basename(run_name_norm)
      run_rep2 <- dirname(run_name_norm)
    }

    message(
      "Note: `run_name` contained a full path. ",
      "Converted automatically to run_rep/run_name."
    )

    return(
      .nm_rt_export_run(
        run_rep = run_rep2,
        run_name = run_name2,
        output_name = output_name,
        extensions = extensions
      )
    )
  }

  # Standard case
  .nm_rt_export_run(
    run_rep = run_rep,
    run_name = run_name,
    output_name = output_name,
    extensions = extensions
  )
}
