#' Goodness-of-fit diagnostic plots from a NONMEM $TABLE file
#'
#' Reads a NONMEM output table (space/tab-delimited, possibly with repeated
#' headers and multiple "TABLE NO." blocks), keeps the observation rows, and
#' builds the standard goodness-of-fit diagnostics: observed vs predictions and
#' residual plots (CWRES / NPDE / IWRES). Plots are only built when the required
#' columns are present. Optionally writes all plots to a Word (.docx) file in the
#' analysis directory.
#'
#' @param file Path to a NONMEM table file, OR a run directory (.reg / .nm7 /
#'   any folder) inside which a NONMEM table is auto-detected. Since the table is
#'   not necessarily named "sdtab", detection is by content ("TABLE NO." header),
#'   optionally filtered by \code{stem}.
#' @param pred Preference for the population prediction column used on x-axes.
#'   The first available name is used. Default \code{c("PRED","CPRED")}.
#'   Pass e.g. \code{pred = "CPRED"} to force one.
#' @param time Preference for the time column used on residual-vs-time plots.
#'   Default \code{c("TIME","TAD")}.
#' @param by Optional character vector of column names for faceting
#'   (\code{facet_wrap}), e.g. \code{c("SEX","CMT")}.
#' @param loess Add a loess trend (default TRUE).
#' @param reg Add a linear regression (default TRUE).
#' @param se Show confidence band around the linear regression (default FALSE).
#' @param alpha Point transparency (default 0.5).
#' @param log_pred Log10 axes for the DV-vs-prediction plots (default FALSE).
#' @param stem Optional filename stem to disambiguate auto-detection in a run
#'   directory (e.g. "sdtab", "patab"). Not required to be "sdtab".
#' @param save_docx If TRUE, write all plots to a .docx in the analysis
#'   directory (default FALSE).
#' @param docx_file Output .docx path. If NULL and \code{save_docx=TRUE}, a name
#'   is derived from the table file, written next to it.
#' @param width,height,dpi Image size (inches) and resolution for the docx
#'   export.
#'
#' @return An object of class \code{"nm_gof"}: a list with \code{plots} (named
#'   list of ggplot objects), \code{data} (the filtered observation rows),
#'   \code{file}, \code{columns} (resolved column mapping), \code{skipped}
#'   (plots that could not be built and why), and \code{by}.
#'
#' @examples
#' \dontrun{
#' g <- nm_gof("C:/nm74g64/run/lopi_04.reg/lopi_04_03_04.fit")
#' g                       # summary
#' plot(g)                 # all panels
#' plot(g, which = "dv_pred")
#' nm_gof(".../run.reg", pred = "CPRED", by = "SEX", save_docx = TRUE)
#' }
#' @export
nm_gof <- function(file,
                   pred = c("PRED", "CPRED"),
                   time = c("TIME", "TAD"),
                   by = NULL,
                   loess = TRUE,
                   reg = TRUE,
                   se = FALSE,
                   alpha = 0.5,
                   log_pred = FALSE,
                   stem = NULL,
                   save_docx = FALSE,
                   docx_file = NULL,
                   width = 6.5, height = 5, dpi = 150) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("nm_gof() requires 'ggplot2'. Install it or add to Suggests.", call. = FALSE)
  }

  ## ---- resolve the table file (explicit file or run/project directory) ----
  tab_file <- .nm_locate_table(file, stem = stem)

  ## ---- read + filter observation rows ----
  df  <- .nm_read_table(tab_file)
  obs <- .nm_filter_obs(df)
  if (nrow(obs) == 0L) {
    stop("No observation rows found (checked MDV==0, then EVID==0, then non-missing DV).",
         call. = FALSE)
  }

  ## ---- resolve the flexible columns ----
  pred_col <- .nm_resolve_col(obs, pred)
  time_col <- .nm_resolve_col(obs, time)
  has <- function(x) !is.na(x) && x %in% names(obs)

  ## validate 'by'
  if (!is.null(by)) {
    miss_by <- setdiff(by, names(obs))
    if (length(miss_by)) {
      warning("Faceting column(s) not found and ignored: ",
              paste(miss_by, collapse = ", "), call. = FALSE)
      by <- setdiff(by, miss_by)
    }
    if (length(by) == 0L) by <- NULL
  }

  ## ---- plot catalogue: id, y, x, reference type, title ----
  ipred_col <- .nm_resolve_col(obs, "IPRED")
  specs <- list(
    dv_pred   = list(y = "DV",    x = pred_col,  ref = "identity",
                     title = sprintf("DV vs %s", pred_col %|na|% "PRED")),
    dv_ipred  = list(y = "DV",    x = ipred_col, ref = "identity",
                     title = "DV vs IPRED"),
    cwres_time= list(y = "CWRES", x = time_col,  ref = "zero",
                     title = sprintf("CWRES vs %s", time_col %|na|% "TIME")),
    cwres_pred= list(y = "CWRES", x = pred_col,  ref = "zero",
                     title = sprintf("CWRES vs %s", pred_col %|na|% "PRED")),
    npde_time = list(y = "NPDE",  x = time_col,  ref = "zero",
                     title = sprintf("NPDE vs %s", time_col %|na|% "TIME")),
    npde_pred = list(y = "NPDE",  x = pred_col,  ref = "zero",
                     title = sprintf("NPDE vs %s", pred_col %|na|% "PRED")),
    iwres_pred= list(y = "IWRES", x = pred_col,  ref = "zero",
                     title = sprintf("IWRES vs %s", pred_col %|na|% "PRED"))
  )

  plots   <- list()
  skipped <- character(0)
  for (id in names(specs)) {
    sp <- specs[[id]]
    if (has(sp$y) && has(sp$x)) {
      is_dv <- identical(sp$ref, "identity")
      plots[[id]] <- .nm_gof_one(
        obs, xcol = sp$x, ycol = sp$y, ref = sp$ref, title = sp$title,
        by = by, loess = loess, reg = reg, se = se, alpha = alpha,
        log_xy = is_dv && log_pred
      )
    } else {
      why <- character(0)
      if (!has(sp$y)) why <- c(why, paste0("no ", sp$y))
      if (!has(sp$x)) why <- c(why, paste0("no ", sp$x %|na|% "x column"))
      skipped[id] <- paste(why, collapse = ", ")
    }
  }

  out <- structure(
    list(
      plots   = plots,
      data    = obs,
      file    = tab_file,
      columns = list(pred = pred_col, time = time_col, ipred = ipred_col,
                     available = names(obs)),
      skipped = skipped,
      by      = by
    ),
    class = "nm_gof"
  )

  if (isTRUE(save_docx)) {
    out$docx <- .nm_gof_docx(out, docx_file = docx_file,
                             width = width, height = height, dpi = dpi)
  }
  out
}

# ----------------------------------------------------------------------------
# internal: locate a NONMEM table given a file or directory
# ----------------------------------------------------------------------------
.nm_locate_table <- function(file, stem = NULL) {
  if (!file.exists(file)) stop("Path does not exist: ", file, call. = FALSE)
  if (!dir.exists(file)) return(normalizePath(file, mustWork = TRUE))

  cand <- list.files(file, full.names = TRUE, recursive = FALSE)
  cand <- cand[!dir.exists(cand)]
  if (!is.null(stem)) {
    cand <- cand[grepl(stem, basename(cand), fixed = TRUE)]
  }
  is_tab <- vapply(cand, function(f) {
    ln <- tryCatch(readLines(f, n = 1L, warn = FALSE), error = function(e) "")
    length(ln) && grepl("^\\s*TABLE NO", ln)
  }, logical(1))
  cand <- cand[is_tab]

  if (length(cand) == 0L)
    stop("No NONMEM table (\"TABLE NO.\" header) found in directory: ", file,
         if (!is.null(stem)) paste0(" (stem = '", stem, "')") else "", call. = FALSE)
  if (length(cand) > 1L)
    stop("Several NONMEM tables found; disambiguate with 'stem' or pass the file directly:\n  ",
         paste(basename(cand), collapse = "\n  "), call. = FALSE)
  normalizePath(cand, mustWork = TRUE)
}

# ----------------------------------------------------------------------------
# internal: robust NONMEM $TABLE reader
#   - drops "TABLE NO." lines and repeated header lines
#   - tolerates multiple concatenated tables (uses first header for names)
#   - fast numeric parse via scan()
# ----------------------------------------------------------------------------
.nm_read_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- sub("\r$", "", lines)
  lines <- lines[nzchar(trimws(lines))]

  ## a "meta" line contains a letter other than the e/E of scientific notation
  is_meta <- grepl("[A-DF-Za-df-z]", lines)
  is_tno  <- grepl("^\\s*TABLE NO", lines)

  hdr_idx <- which(is_meta & !is_tno)
  if (length(hdr_idx) == 0L)
    stop("No column header line found in table: ", path, call. = FALSE)

  col_names <- strsplit(trimws(lines[hdr_idx[1L]]), "\\s+")[[1]]
  ncol <- length(col_names)

  ## warn if a later header disagrees (different table structure concatenated)
  if (length(hdr_idx) > 1L) {
    other <- strsplit(trimws(lines[hdr_idx[2L]]), "\\s+")[[1]]
    if (!identical(other, col_names))
      warning("Multiple differing headers found; using the first one.", call. = FALSE)
  }

  data_lines <- lines[!is_meta]
  if (length(data_lines) == 0L)
    stop("No numeric data rows found in table: ", path, call. = FALSE)

  vals <- scan(text = data_lines, what = numeric(), quiet = TRUE)
  if (length(vals) %% ncol != 0L) {
    ## ragged rows: fall back to a per-row split keeping full-width rows only
    sp  <- strsplit(trimws(data_lines), "\\s+")
    keep <- lengths(sp) == ncol
    if (!all(keep))
      warning(sum(!keep), " malformed data row(s) dropped.", call. = FALSE)
    vals <- as.numeric(unlist(sp[keep], use.names = FALSE))
  }

  m <- matrix(vals, ncol = ncol, byrow = TRUE)
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  names(df) <- col_names
  df
}

# ----------------------------------------------------------------------------
# internal: keep observation rows (MDV==0 > EVID==0 > non-missing DV)
# ----------------------------------------------------------------------------
.nm_filter_obs <- function(df) {
  nm <- toupper(names(df))
  if ("MDV" %in% nm) {
    return(df[df[[which(nm == "MDV")[1]]] == 0, , drop = FALSE])
  }
  if ("EVID" %in% nm) {
    return(df[df[[which(nm == "EVID")[1]]] == 0, , drop = FALSE])
  }
  if ("DV" %in% nm) {
    dv <- df[[which(nm == "DV")[1]]]
    return(df[!is.na(dv), , drop = FALSE])
  }
  df
}

# ----------------------------------------------------------------------------
# internal: resolve first available column (case-insensitive), return real name
# ----------------------------------------------------------------------------
.nm_resolve_col <- function(df, candidates) {
  if (is.null(candidates)) return(NA_character_)
  nm <- names(df)
  up <- toupper(nm)
  for (c in candidates) {
    hit <- which(up == toupper(c))
    if (length(hit)) return(nm[hit[1]])
  }
  NA_character_
}

# ----------------------------------------------------------------------------
# internal: build one GOF ggplot
# ----------------------------------------------------------------------------
.nm_gof_one <- function(df, xcol, ycol, ref, title, by = NULL,
                        loess = TRUE, reg = TRUE, se = FALSE, alpha = 0.5,
                        log_xy = FALSE) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]])) +
    ggplot2::geom_point(alpha = alpha, size = 1.4, na.rm = TRUE)

  if (identical(ref, "identity")) {
    p <- p + ggplot2::geom_abline(intercept = 0, slope = 1,
                            linetype = "dashed", colour = "grey40")
  } else {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40")
  }

  if (isTRUE(reg)) {
    p <- p + ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = se,
                            colour = "#c0392b", linewidth = 0.7, na.rm = TRUE)
  }
  if (isTRUE(loess)) {
    p <- p + ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                            colour = "#2471a3", linewidth = 0.7, na.rm = TRUE)
  }

  if (isTRUE(log_xy)) {
    p <- p + ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
  }

  if (!is.null(by) && length(by)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", paste(by, collapse = "+"))))
  }

  p + ggplot2::labs(title = title, x = xcol, y = ycol) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 11,
        margin = ggplot2::margin(b = 6)
      ),
      axis.title = ggplot2::element_text(
        size = 10,
        face = "bold"
      ),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 8)
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 8)
      ),
      axis.text = ggplot2::element_text(
        size = 9
      ),
      plot.margin = ggplot2::margin(
        t = 8,
        r = 10,
        b = 18,
        l = 10,
        unit = "pt"
      )
    )
}

# ----------------------------------------------------------------------------
# internal: write all plots to a .docx in the analysis directory
# ----------------------------------------------------------------------------
.nm_gof_docx <- function(x, docx_file = NULL, width = 6.5, height = 5, dpi = 150) {
  if (!requireNamespace("officer", quietly = TRUE))
    stop("save_docx = TRUE requires the 'officer' package.", call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("save_docx = TRUE requires 'ggplot2'.", call. = FALSE)

  if (is.null(docx_file)) {
    base <- tools::file_path_sans_ext(basename(x$file))
    docx_file <- file.path(dirname(x$file), paste0(base, "_GOF.docx"))
  }

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Goodness-of-fit diagnostics", style = "heading 1")
  doc <- officer::body_add_par(doc, paste("Table:", basename(x$file)), style = "Normal")
  doc <- officer::body_add_par(
    doc, sprintf("Observation rows: %d%s", nrow(x$data),
                 if (!is.null(x$by)) paste0("  |  faceted by: ",
                                            paste(x$by, collapse = ", ")) else ""),
    style = "Normal")

  for (id in names(x$plots)) {
    p <- x$plots[[id]]
    png <- tempfile(fileext = ".png")
    ggplot2::ggsave(png, plot = p, width = width, height = height,
                    dpi = dpi, units = "in")
    doc <- officer::body_add_par(doc, p$labels$title %|na|% id, style = "heading 2")
    doc <- officer::body_add_img(doc, src = png, width = width, height = height)
    unlink(png)
  }

  print(doc, target = docx_file)
  normalizePath(docx_file, mustWork = TRUE)
}

# ----------------------------------------------------------------------------
# S3 methods
# ----------------------------------------------------------------------------
#' @export
print.nm_gof <- function(x, ...) {
  cat("<nm_gof>\n")
  cat("  table   :", x$file, "\n")
  cat("  obs rows:", nrow(x$data), "\n")
  cat("  pred    :", x$columns$pred %|na|% "-",
      " | time:", x$columns$time %|na|% "-",
      " | ipred:", x$columns$ipred %|na|% "-", "\n")
  if (!is.null(x$by)) cat("  facet by:", paste(x$by, collapse = ", "), "\n")
  cat("  plots   :", if (length(x$plots)) paste(names(x$plots), collapse = ", ") else "none", "\n")
  if (length(x$skipped))
    cat("  skipped :",
        paste(sprintf("%s (%s)", names(x$skipped), x$skipped), collapse = "; "), "\n")
  if (!is.null(x$docx)) cat("  docx    :", x$docx, "\n")
  invisible(x)
}

#' @export
plot.nm_gof <- function(x, which = NULL, ncol = 2, ...) {
  if (length(x$plots) == 0L) {
    message("No plots to display."); return(invisible(x))
  }
  sel <- if (is.null(which)) names(x$plots) else intersect(which, names(x$plots))
  if (length(sel) == 0L) {
    message("None of the requested plots exist. Available: ",
            paste(names(x$plots), collapse = ", "))
    return(invisible(x))
  }
  ps <- x$plots[sel]
  if (length(ps) == 1L) {
    print(ps[[1]]); return(invisible(x))
  }
  if (requireNamespace("patchwork", quietly = TRUE)) {
    print(patchwork::wrap_plots(ps, ncol = ncol))
  } else {
    for (p in ps) print(p)
  }
  invisible(x)
}

