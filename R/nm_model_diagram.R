# Diagram of the models tested in a project (action 7): tree of the runs, each
# run linked to the run it was derived from (its parent), with its OFV, its
# dOFV to the parent and its convergence state.

# Parent written in the control stream as a comment, following the PsN and
# Pirana convention: ";; 1. Based on: 3", ";; Based on: run3", "; based on run3".
# Returns the token written (NA when there is none).
.nm_md_parent_comment <- function(lines) {
  rx <- "^\\s*;+\\s*(?:[0-9]+\\.\\s*)?based\\s+on\\s*:?\\s*([^\\s;,]+)"
  m <- regmatches(lines, regexec(rx, lines, ignore.case = TRUE, perl = TRUE))
  m <- m[lengths(m) == 2L]
  if (!length(m)) return(NA_character_)
  tok <- trimws(m[[1]][2])
  tok <- sub("\\.(ctl|mod)$", "", tok, ignore.case = TRUE)
  if (!nzchar(tok)) NA_character_ else tok
}

# Text of $PROBLEM (first line), used as a short description.
.nm_md_problem <- function(lines) {
  i <- grep("^\\s*\\$PROB", lines, ignore.case = TRUE)
  if (!length(i)) return(NA_character_)
  d <- trimws(sub("^\\s*\\$PROB[A-Z]*", "", lines[i[1]], ignore.case = TRUE))
  if (nzchar(d)) d else NA_character_
}

# Control stream of a run: in its run folder (same stem), else next to the
# run folder (Wings layout: run2.ctl beside run2.reg/).
.nm_md_find_ctl <- function(run, dir) {
  if (is.na(dir) || !nzchar(dir)) return(NA_character_)
  for (d in unique(c(dir, dirname(dir)))) {
    f <- list.files(d, pattern = "(?i)\\.(ctl|mod)$", full.names = TRUE)
    hit <- f[tolower(tools::file_path_sans_ext(basename(f))) == tolower(run)]
    if (length(hit)) return(hit[1])
  }
  NA_character_
}

# Match the token written in "Based on" to a run name: exact, then without
# or with the "run" prefix (PsN writes the number only: "Based on: 3" for
# run3.mod), case ignored.
.nm_md_match_run <- function(tok, runs) {
  if (is.na(tok)) return(NA_character_)
  cand <- unique(c(tok, paste0("run", tok), sub("^run", "", tok, ignore.case = TRUE)))
  hit <- runs[match(tolower(cand), tolower(runs))]
  hit <- hit[!is.na(hit)]
  if (length(hit)) hit[1] else NA_character_
}

.nm_md_status <- function(r) {
  state <- if ("run_state" %in% names(r)) as.character(r$run_state) else rep(NA_character_, nrow(r))
  min_ok <- if ("min_ok" %in% names(r)) r$min_ok else rep(NA, nrow(r))
  cov_ok <- if ("cov_ok" %in% names(r)) r$cov_ok else rep(NA, nrow(r))
  out <- ifelse(!is.na(state) & state == "simulation", "simulation",
         ifelse(!is.na(state) & state != "completed", "not completed",
         ifelse(is.na(r$ofv), "not completed",
         ifelse(.nm_md_true(min_ok) & .nm_md_true(cov_ok), "converged, covariance OK",
         ifelse(.nm_md_true(min_ok), "converged, no covariance", "not converged")))))
  factor(out, levels = c("converged, covariance OK", "converged, no covariance",
                         "not converged", "not completed", "simulation"))
}
.nm_md_true <- function(x) !is.na(x) & as.logical(x)

# Tidy layout of a forest: leaves get consecutive rows, a parent is centred on
# its children, the depth gives the column. Children are sorted by name
# (numbers in the names sorted as numbers: run2 before run10).
.nm_md_layout <- function(run, parent) {
  n <- length(run)
  key <- function(v) {
    num <- suppressWarnings(as.numeric(sub("^\\D*(\\d+).*$", "\\1", v)))
    order(sub("\\d.*$", "", v), ifelse(is.na(num), Inf, num), v)
  }
  kids <- lapply(run, function(r) { k <- run[!is.na(parent) & parent == r]; k[key(k)] })
  names(kids) <- run
  roots <- run[is.na(parent)]
  roots <- roots[key(roots)]
  x <- stats::setNames(rep(NA_real_, n), run)
  y <- x
  nxt <- 0
  visit <- character()
  place <- function(r, depth) {
    visit <<- c(visit, r)
    x[[r]] <<- depth
    k <- kids[[r]]
    if (!length(k)) { nxt <<- nxt + 1; y[[r]] <<- nxt; return(invisible()) }
    for (c in k) place(c, depth + 1)
    y[[r]] <<- mean(y[k])
  }
  for (r in roots) place(r, 0)
  # rows in depth-first order: each run comes right after its parent
  data.frame(run = visit, x = unname(x[visit]), y = unname(y[visit]), stringsAsFactors = FALSE)
}

#' Diagram of the models tested in a project
#'
#' Draws the tree of the runs of a project: each run is linked to the run it
#' was derived from (its parent), and shows its OFV, the change of OFV from
#' its parent (dOFV) and the change in the number of estimated parameters,
#' coloured by convergence state.
#'
#' The parent of a run is read in its control stream, in a comment written
#' with the PsN and Pirana convention, e.g. `;; 1. Based on: run3` or
#' `;; Based on: 3` (the number alone stands for `run3`). The `parents`
#' argument completes or corrects these links. A run whose parent is not
#' found is drawn as a root.
#'
#' @param x Result of [nm_extract_runs()], or the path of an analysis
#'   directory, which is then read with `nm_extract_runs(x)`.
#' @param parents Optional links that complete or replace those read in the
#'   control streams: a named character vector (`c(run2 = "run1")`) or a
#'   data frame with columns `run` and `parent`. `NA` removes a link.
#' @param label_width Maximum width of the description (text of `$PROBLEM`)
#'   in the table; longer texts are cut.
#'
#' @return An object of class `nm_model_diagram`: a list with `tree`, a
#'   data frame with one row per run (`run`, `parent`, `parent_source`:
#'   "control stream", "argument" or "none", `ofv`, `dofv_parent`, `n_par`,
#'   `dpar_parent`, `status`, `description`, and the layout coordinates `x`
#'   and `y`; rows in depth-first order, each run right after its parent),
#'   and `plot`, the ggplot2 tree. `plot()` draws it; it can be
#'   saved with [nm_export_plot()].
#'
#' @examples
#' ex <- system.file("extdata", "example", package = "nmpilot")
#' d <- nm_model_diagram(ex, parents = c(run2 = "run1"))
#' d
#' plot(d)
#' @seealso [nm_extract_runs()], [nm_compare_runs()]
#' @export
nm_model_diagram <- function(x, parents = NULL, label_width = 40) {
  if (is.character(x) && length(x) == 1L) {
    if (!dir.exists(x)) stop("Directory not found: ", x, call. = FALSE)
    x <- nm_extract_runs(x, verbose = FALSE)
  }
  r <- if (is.list(x) && !is.null(x$runs)) x$runs else x
  if (!is.data.frame(r) || !all(c("run", "ofv") %in% names(r)))
    stop("`x` must be the result of nm_extract_runs() or an analysis directory.", call. = FALSE)
  r <- as.data.frame(r, stringsAsFactors = FALSE)
  if (!nrow(r)) stop("No run found.", call. = FALSE)
  runs <- as.character(r$run)
  if (anyDuplicated(runs)) stop("Duplicated run names: ", paste(unique(runs[duplicated(runs)]), collapse = ", "),
                                call. = FALSE)

  # Parents written in the control streams, and description.
  dirs <- if ("dir" %in% names(r)) as.character(r$dir) else rep(NA_character_, nrow(r))
  tok <- desc <- rep(NA_character_, nrow(r))
  for (i in seq_len(nrow(r))) {
    ctl <- .nm_md_find_ctl(runs[i], dirs[i])
    if (is.na(ctl)) next
    L <- readLines(ctl, warn = FALSE)
    tok[i] <- .nm_md_parent_comment(L)
    desc[i] <- .nm_md_problem(L)
  }
  parent <- vapply(tok, .nm_md_match_run, character(1), runs = runs, USE.NAMES = FALSE)
  src <- ifelse(is.na(parent), "none", "control stream")
  unknown <- !is.na(tok) & is.na(parent)
  if (any(unknown))
    warning("Parent written in the control stream but not found among the runs: ",
            paste0(runs[unknown], " (", tok[unknown], ")", collapse = ", "), call. = FALSE)

  # Links given by the user.
  if (!is.null(parents)) {
    if (is.data.frame(parents)) {
      if (!all(c("run", "parent") %in% names(parents)))
        stop("`parents` must have columns `run` and `parent`.", call. = FALSE)
      pv <- stats::setNames(as.character(parents$parent), as.character(parents$run))
    } else {
      if (is.null(names(parents)) || any(!nzchar(names(parents))))
        stop("`parents` must be a named vector, e.g. c(run2 = \"run1\").", call. = FALSE)
      pv <- stats::setNames(as.character(parents), names(parents))
    }
    bad <- setdiff(names(pv), runs)
    if (length(bad)) stop("Unknown run(s) in `parents`: ", paste(bad, collapse = ", "), call. = FALSE)
    badp <- setdiff(stats::na.omit(unname(pv)), runs)
    if (length(badp)) stop("Unknown parent(s) in `parents`: ", paste(badp, collapse = ", "), call. = FALSE)
    for (k in names(pv)) {
      i <- match(k, runs)
      parent[i] <- pv[[k]]
      src[i] <- if (is.na(pv[[k]])) "none" else "argument"
    }
  }

  # A run cannot be its own ancestor: such links are removed.
  self <- !is.na(parent) & parent == runs
  for (i in seq_along(runs)) {
    seen <- runs[i]; p <- parent[i]
    while (!is.na(p) && !self[i]) {
      if (p %in% seen) { self[i] <- TRUE; break }
      seen <- c(seen, p); p <- parent[match(p, runs)]
    }
  }
  if (any(self)) {
    warning("Circular parent links removed for: ", paste(runs[self], collapse = ", "), call. = FALSE)
    parent[self] <- NA_character_; src[self] <- "none"
  }

  ofv <- suppressWarnings(as.numeric(r$ofv))
  npar <- if ("n_par" %in% names(r)) suppressWarnings(as.integer(r$n_par)) else rep(NA_integer_, nrow(r))
  ip <- match(parent, runs)
  desc_cut <- ifelse(!is.na(desc) & nchar(desc) > label_width,
                     paste0(substr(desc, 1L, label_width - 3L), "..."), desc)
  tree <- data.frame(
    run = runs, parent = parent, parent_source = src,
    ofv = ofv, dofv_parent = ofv - ofv[ip],
    n_par = npar, dpar_parent = npar - npar[ip],
    status = .nm_md_status(r), description = desc_cut,
    stringsAsFactors = FALSE)
  lay <- .nm_md_layout(tree$run, tree$parent)
  tree$x <- lay$x[match(tree$run, lay$run)]
  tree$y <- lay$y[match(tree$run, lay$run)]
  tree <- tree[match(lay$run, tree$run), , drop = FALSE]
  rownames(tree) <- NULL

  structure(list(tree = tree, plot = .nm_md_plot(tree)), class = "nm_model_diagram")
}

.nm_md_plot <- function(tree) {
  fmt <- function(v, d = 1) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = d, big.mark = " "))
  lab <- paste0(tree$run, "\nOFV ", fmt(tree$ofv))
  has_p <- !is.na(tree$parent)
  dv <- tree$dofv_parent[has_p]
  lab[has_p] <- paste0(lab[has_p], "\ndOFV ", ifelse(!is.na(dv) & dv > 0, "+", ""),
                       fmt(dv),
                       ifelse(is.na(tree$dpar_parent[has_p]), "",
                              paste0(" (", ifelse(tree$dpar_parent[has_p] >= 0, "+", ""),
                                     tree$dpar_parent[has_p], " par.)")))
  tree$label <- lab
  # Elbow links: parent -> mid-column -> child.
  e <- tree[has_p, , drop = FALSE]
  edges <- do.call(rbind, lapply(seq_len(nrow(e)), function(k) {
    p <- tree[tree$run == e$parent[k], ]
    xm <- p$x + 0.5
    data.frame(g = k, x = c(p$x, xm, xm, e$x[k]), y = c(p$y, p$y, e$y[k], e$y[k]))
  }))
  cols <- c("converged, covariance OK" = "#1b9e77", "converged, no covariance" = "#e6ab02",
            "not converged" = "#d95f02", "not completed" = "#999999", "simulation" = "#7570b3")
  p <- ggplot2::ggplot(tree, ggplot2::aes(x = .data$x, y = .data$y))
  if (!is.null(edges))
    p <- p + ggplot2::geom_path(data = edges, ggplot2::aes(x = .data$x, y = .data$y, group = .data$g),
                                colour = "grey55", linewidth = 0.4)
  p +
    ggplot2::geom_label(ggplot2::aes(label = .data$label, fill = .data$status),
                        size = 3, lineheight = 0.95, key_glyph = "rect",
                        label.padding = ggplot2::unit(0.25, "lines")) +
    ggplot2::scale_fill_manual(values = cols, drop = TRUE, name = NULL) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(add = 0.7)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::labs(title = "Models tested", x = NULL, y = NULL,
                  caption = "dOFV and number of estimated parameters relative to the parent run") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_text(face = "bold"),
                   plot.caption = ggplot2::element_text(colour = "grey40", size = 8))
}

#' @export
print.nm_model_diagram <- function(x, ...) {
  t <- x$tree
  n_root <- sum(is.na(t$parent))
  cat("Diagram of the models tested:", nrow(t), "run(s),", n_root, "root(s)\n\n")
  show <- t[, c("run", "parent", "ofv", "dofv_parent", "dpar_parent", "status", "description")]
  show$ofv <- round(show$ofv, 3); show$dofv_parent <- round(show$dofv_parent, 3)
  print(show, row.names = FALSE)
  if (any(t$parent_source == "none" & !is.na(t$ofv)) && n_root > 1L)
    cat("\nRuns without parent are drawn as roots: write \";; Based on: <run>\" in their",
        "control stream or use `parents =`.\n")
  invisible(x)
}

#' @export
plot.nm_model_diagram <- function(x, ...) {
  print(x$plot)
  invisible(x$plot)
}
