# Summary table of all the runs of a project, exported to Word (action 8):
# one row per run with its parent, OFV, dOFV, number of parameters, AIC, BIC,
# minimisation and covariance status, maximal shrinkage and share of RSE > 50%.

#' Summary table of all the runs of a project
#'
#' Builds one table with one row per run of a project and, when
#' `output_file` is given, writes it to a Word document. Each run is shown
#' with its description (text of `$PROBLEM`), the run it was derived from
#' (its parent, read in the control stream as in [nm_model_diagram()]), its
#' OFV, the change of OFV from its parent and from the best run, its number
#' of estimated parameters, its AIC and BIC, its minimisation and covariance
#' status, its largest ETA shrinkage and the share of parameters estimated
#' with a relative standard error above 50%.
#'
#' @param x Result of [nm_extract_runs()], or the path of an analysis
#'   directory, which is then read with `nm_extract_runs(x)`.
#' @param parents Optional links that complete or replace those read in the
#'   control streams, as in [nm_model_diagram()]: a named character vector
#'   (`c(run2 = "run1")`) or a data frame with columns `run` and `parent`.
#' @param sort Order of the rows: `"tree"` (default; each run right after its
#'   parent, as in [nm_model_diagram()]), `"ofv"` (best OFV first) or `"run"`
#'   (run names, numbers sorted as numbers).
#' @param output_file Path of the Word file (.docx). `NULL` (default): nothing
#'   is written.
#' @param overwrite Logical. Replace an existing `output_file` (default
#'   `FALSE`).
#' @param title Title written above the table in the Word file.
#' @param verbose Print a message when the file is written.
#'
#' @return An object of class `nm_runs_table`: a data frame with one row per
#'   run and the columns `run`, `description`, `parent`, `ofv`, `dofv_parent`
#'   (OFV of the run minus OFV of its parent), `dofv_best` (minus the lowest
#'   OFV of the project), `n_par`, `aic`, `bic`, `minimization`, `covariance`,
#'   `max_eta_shrink` (%), `pct_rse_gt50` (% of the parameters with an RSE
#'   above 50%) and `status` (as in [nm_model_diagram()]). The Word file, when
#'   written, has the same rows in landscape format; its path is kept in the
#'   attribute `output_file`.
#'
#' @examples
#' ex <- system.file("extdata", "example", package = "nmpilot")
#' tab <- nm_runs_table(ex, parents = c(run2 = "run1"))
#' tab
#' \donttest{
#' if (requireNamespace("flextable", quietly = TRUE) &&
#'     requireNamespace("officer", quietly = TRUE)) {
#'   f <- file.path(tempdir(), "runs_summary.docx")
#'   nm_runs_table(ex, parents = c(run2 = "run1"), output_file = f, overwrite = TRUE)
#' }
#' }
#' @seealso [nm_extract_runs()], [nm_model_diagram()], [nm_compare_runs()],
#'   [nm_extract_run_table()] (table of the estimates of one run)
#' @export
nm_runs_table <- function(x, parents = NULL, sort = c("tree", "ofv", "run"),
                          output_file = NULL, overwrite = FALSE,
                          title = "Summary of the runs", verbose = TRUE) {
  sort <- match.arg(sort)
  if (is.character(x) && length(x) == 1L) {
    if (!dir.exists(x)) stop("Directory not found: ", x, call. = FALSE)
    x <- nm_extract_runs(x, verbose = FALSE)
  }
  r <- if (is.list(x) && !is.null(x$runs)) x$runs else x
  if (!is.data.frame(r) || !all(c("run", "ofv") %in% names(r)))
    stop("`x` must be the result of nm_extract_runs() or an analysis directory.", call. = FALSE)
  r <- as.data.frame(r, stringsAsFactors = FALSE)
  if (!is.null(output_file)) {
    if (!is.character(output_file) || length(output_file) != 1L || !nzchar(output_file))
      stop("`output_file` must be a single file path.", call. = FALSE)
    if (!grepl("\\.docx$", output_file, ignore.case = TRUE))
      stop("`output_file` must end with .docx.", call. = FALSE)
    if (file.exists(output_file) && !isTRUE(overwrite))
      stop("Output file already exists (use overwrite = TRUE): ", output_file, call. = FALSE)
  }

  # Parents, dOFV to the parent, description and status: same rules as the
  # diagram of the models.
  d <- nm_model_diagram(x, parents = parents, label_width = 1000L)$tree
  col <- function(nm, as = as.numeric) {
    v <- if (nm %in% names(r)) r[[nm]] else rep(NA, nrow(r))
    suppressWarnings(as(v))
  }
  i <- match(d$run, as.character(r$run))
  ofv <- d$ofv
  tab <- data.frame(
    run = d$run,
    description = d$description,
    parent = d$parent,
    ofv = ofv,
    dofv_parent = d$dofv_parent,
    dofv_best = ofv - suppressWarnings(min(ofv, na.rm = TRUE)),
    n_par = d$n_par,
    aic = col("aic")[i],
    bic = col("bic")[i],
    minimization = col("min_status", as.character)[i],
    covariance = col("cov_status", as.character)[i],
    max_eta_shrink = col("max_eta_shrink")[i],
    pct_rse_gt50 = col("pct_rse_gt50")[i],
    status = as.character(d$status),
    stringsAsFactors = FALSE)
  if (all(is.na(ofv))) tab$dofv_best <- NA_real_

  if (sort == "ofv") {
    tab <- tab[order(is.na(tab$ofv), tab$ofv), , drop = FALSE]
  } else if (sort == "run") {
    num <- suppressWarnings(as.numeric(sub("^\\D*(\\d+).*$", "\\1", tab$run)))
    tab <- tab[order(sub("\\d.*$", "", tab$run), ifelse(is.na(num), Inf, num), tab$run), , drop = FALSE]
  }
  rownames(tab) <- NULL
  class(tab) <- c("nm_runs_table", "data.frame")

  if (!is.null(output_file)) {
    .nm_runs_table_docx(tab, output_file, title = title)
    attr(tab, "output_file") <- output_file
    if (verbose) message("Word file written: ", output_file)
  }
  tab
}

# Values shown in the table (print and Word): rounded numbers, empty cells
# for missing values.
.nm_runs_table_display <- function(tab) {
  f <- function(v, d) ifelse(is.na(v), "", formatC(v, format = "f", digits = d, big.mark = ""))
  txt <- function(v) ifelse(is.na(v), "", as.character(v))
  data.frame(
    Run = txt(tab$run),
    Description = txt(tab$description),
    Parent = txt(tab$parent),
    OFV = f(tab$ofv, 2),
    `dOFV (parent)` = f(tab$dofv_parent, 2),
    `dOFV (best)` = f(tab$dofv_best, 2),
    `N par` = txt(tab$n_par),
    AIC = f(tab$aic, 1),
    BIC = f(tab$bic, 1),
    Minimization = txt(tab$minimization),
    Covariance = txt(tab$covariance),
    `Max shrinkage (%)` = f(tab$max_eta_shrink, 1),
    `RSE > 50% (%)` = f(tab$pct_rse_gt50, 1),
    check.names = FALSE, stringsAsFactors = FALSE)
}

.nm_runs_table_docx <- function(tab, output_file, title) {
  for (pkg in c("flextable", "officer")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for the Word export. ",
           "Install it with install.packages(\"", pkg, "\").", call. = FALSE)
  }
  disp <- .nm_runs_table_display(tab)
  ft <- flextable::flextable(disp)
  ft <- flextable::fontsize(ft, size = 8, part = "all")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::bg(ft, bg = "#D9E2F3", part = "header")
  ft <- flextable::align(ft, j = 4:ncol(disp), align = "right", part = "body")
  ft <- flextable::align(ft, align = "center", part = "header")
  best <- which(!is.na(tab$dofv_best) & tab$dofv_best == 0)
  if (length(best)) ft <- flextable::bold(ft, i = best, part = "body")
  bad <- which(tab$status %in% c("not converged", "not completed"))
  if (length(bad)) ft <- flextable::color(ft, i = bad, color = "#C00000", part = "body")
  ft <- flextable::width(ft, j = 1:ncol(disp),
                         width = c(0.6, 1.9, 0.6, 0.75, 0.7, 0.7, 0.4, 0.7, 0.7, 1.2, 0.9, 0.8, 0.6))
  ft <- flextable::border_outer(ft, border = officer::fp_border(color = "black", width = 1))

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, title, style = "heading 1")
  doc <- officer::body_add_par(doc, paste0("Generated on ", format(Sys.Date(), "%d/%m/%Y"),
                                           " - ", nrow(tab), " run(s)."), style = "Normal")
  doc <- flextable::body_add_flextable(doc, ft)
  doc <- officer::body_add_par(doc, paste(
    "dOFV (parent): OFV of the run minus OFV of the run it was derived from;",
    "dOFV (best): minus the lowest OFV of the project (best run in bold).",
    "N par: number of estimated parameters. Max shrinkage: largest ETA shrinkage.",
    "RSE > 50%: share of the parameters estimated with a relative standard error above 50%.",
    "Runs in red did not converge or did not complete."), style = "Normal")
  doc <- officer::body_end_section_landscape(doc)
  d <- dirname(output_file)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  print(doc, target = output_file)
  invisible(output_file)
}

#' @export
print.nm_runs_table <- function(x, ...) {
  cat("Summary of", nrow(x), "run(s)\n")
  disp <- .nm_runs_table_display(x)
  disp$Description <- ifelse(nchar(disp$Description) > 30,
                             paste0(substr(disp$Description, 1, 27), "..."), disp$Description)
  print(disp[, c("Run", "Parent", "OFV", "dOFV (parent)", "dOFV (best)", "N par",
                 "AIC", "BIC", "Covariance", "Max shrinkage (%)", "RSE > 50% (%)")],
        row.names = FALSE, right = FALSE)
  if (!is.null(attr(x, "output_file"))) cat("Word file:", attr(x, "output_file"), "\n")
  invisible(x)
}
