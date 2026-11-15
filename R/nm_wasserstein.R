# =============================================================================
# nm_wasserstein.R
# Distribution-based covariate assessment with the Wasserstein distance.
#
# Implements the method of:
#   Simon N, Hulot JS, von Fabeck K. Distribution-based covariate assessment
#   using Wasserstein distance in population pharmacokinetic models.
#   Front Pharmacol. 2026;17:1804989. doi:10.3389/fphar.2026.1804989
#
# For each run and each covariate, the empirical Bayes estimates (EBE) of one
# ETA are split into two groups of subjects (levels of a categorical
# covariate, or extreme quartiles Q1 / Q4 of a continuous one) and compared
# with the 1-Wasserstein distance W and a permutation test. A covariate that
# is well accounted for by the model leaves no residual distributional signal
# in the ETA. A Wasserstein-based R2 compares the overall size of the ETA
# distribution (its distance to 0) with a reference run.
# =============================================================================


# ---- Wasserstein distance (order 1, one dimension) --------------------------

# W1 between two samples: area between their empirical distribution functions,
#   W1 = integral |F_x(t) - F_y(t)| dt
# equal to the quantile-function form integral_0^1 |Fx^-1(u) - Fy^-1(u)| du.
.nm_w1 <- function(x, y) {
  x <- sort(x[is.finite(x)])
  y <- sort(y[is.finite(y)])
  if (!length(x) || !length(y)) return(NA_real_)
  z <- sort(unique(c(x, y)))
  if (length(z) < 2L) return(0)
  fx <- findInterval(z, x) / length(x)
  fy <- findInterval(z, y) / length(y)
  sum(abs(fx - fy)[-length(z)] * diff(z))
}

# W1 between a sample and a point mass at 0 (mean absolute value).
.nm_w1_to_zero <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(abs(x))
}

# Permutation test of W1 between two groups.
#   p = (number of permuted W >= observed W + 1) / (n_perm + 1)
.nm_w_perm <- function(x, y, n_perm = 10000L, seed = NULL) {
  x <- as.numeric(x); y <- as.numeric(y)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L)
    return(list(W = NA_real_, p = NA_real_))
  w_obs <- .nm_w1(x, y)
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv())
            else assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  pooled <- c(x, y)
  nx <- length(x)
  w_perm <- vapply(seq_len(n_perm), function(i) {
    s <- sample(pooled)
    .nm_w1(s[seq_len(nx)], s[-seq_len(nx)])
  }, numeric(1))
  list(W = w_obs, p = (sum(w_perm >= w_obs) + 1) / (n_perm + 1))
}


# ---- Inputs -------------------------------------------------------------------

# Locate the files of one run: a run directory (.reg/.nm7/any) or a path prefix.
.nm_w_run_files <- function(run) {
  if (dir.exists(run)) {
    d <- run
    phi <- list.files(d, pattern = "(?i)\\.phi$", full.names = TRUE)
    if (!length(phi)) stop("No .phi file in run directory: ", run, call. = FALSE)
    stem <- tools::file_path_sans_ext(basename(phi[1]))
    if (length(phi) > 1L) {
      dstem <- tools::file_path_sans_ext(basename(d))
      hit <- phi[tolower(tools::file_path_sans_ext(basename(phi))) == tolower(dstem)]
      if (length(hit)) { phi <- hit; stem <- dstem }
    }
    phi <- phi[1]
  } else {
    phi <- paste0(run, ".phi")
    if (!file.exists(phi)) stop("Run not found (no directory and no ", phi, ").", call. = FALSE)
    d <- dirname(run); stem <- basename(run)
  }
  pick <- function(ext) {
    f <- file.path(d, paste0(stem, ".", ext))
    if (file.exists(f)) return(f)
    g <- list.files(d, pattern = paste0("(?i)\\.", ext, "$"), full.names = TRUE)
    if (length(g) == 1L) g else NA_character_
  }
  list(dir = d, stem = stem, phi = phi, ext = pick("ext"), ctl = pick("ctl"))
}

# Read a NONMEM .phi file; with several TABLE blocks, the last one is kept.
.nm_w_read_phi <- function(file) {
  lines <- .nm_read_lines(file)
  starts <- grep("^TABLE NO\\.", lines)
  if (length(starts)) lines <- lines[(max(starts) + 1L):length(lines)]
  lines <- lines[nzchar(trimws(lines))]
  df <- utils::read.table(text = lines, header = TRUE, check.names = FALSE,
                          stringsAsFactors = FALSE)
  if (!"ID" %in% names(df)) stop("No ID column in .phi file: ", file, call. = FALSE)
  df
}

# Diagonal OMEGA(k,k) from the final row of a .ext file (NA if unavailable).
.nm_w_omega <- function(ext, k) {
  if (is.na(ext) || !file.exists(ext)) return(NA_real_)
  lines <- .nm_read_lines(ext)
  starts <- grep("^TABLE NO\\.", lines)
  if (length(starts)) lines <- lines[(max(starts) + 1L):length(lines)]
  df <- tryCatch(utils::read.table(text = lines, header = TRUE, check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(df)) return(NA_real_)
  col <- sprintf("OMEGA(%d,%d)", k, k)
  fin <- df[df$ITERATION == -1000000000, , drop = FALSE]
  if (!nrow(fin) || !col %in% names(fin)) return(NA_real_)
  as.numeric(fin[[col]][1])
}

# ETA index from "ETA(2)", "ETA2", 2 or a PK parameter name ("CL").
.nm_w_eta_index <- function(eta, ctl) {
  if (is.numeric(eta)) return(as.integer(eta))
  e <- toupper(gsub("\\s", "", eta))
  if (grepl("^ETA\\(?[0-9]+\\)?$", e)) return(as.integer(gsub("\\D", "", e)))
  if (is.na(ctl)) stop("eta = '", eta, "' needs the run's .ctl to identify the ETA; ",
                       "none was found. Use eta = \"ETA(n)\" instead.", call. = FALSE)
  m <- nm_identify_parameters(ctl, likely_pk_only = FALSE)$mapping_detail_all
  m <- m[m$element_type == "ETA" & toupper(m$pk_parameter) == e, , drop = FALSE]
  if (!nrow(m)) stop("No ETA is linked to parameter '", eta, "' in ", ctl, call. = FALSE)
  d <- if ("depth" %in% names(m)) m$depth else rep(0, nrow(m))
  as.integer(m$index[order(d)][1])
}

# Covariates: one row per subject (first record of each ID).
.nm_w_covariates <- function(data, ctl) {
  if (is.null(data)) {
    if (is.na(ctl)) stop("`data` is needed: no .ctl found to locate the data set.", call. = FALSE)
    data <- .nm_extract_data_path(ctl)
    if (!file.exists(data))
      stop("Data set referenced by $DATA not found: ", data,
           ". Supply it with `data =`.", call. = FALSE)
  }
  if (is.character(data)) {
    data <- utils::read.csv(data, na.strings = c(".", "", "NA"), check.names = FALSE,
                            stringsAsFactors = FALSE, comment.char = "")
    first <- names(data)[1]
    data <- data[!grepl("^[#@;A-Za-z]", as.character(data[[first]])) |
                   suppressWarnings(!is.na(as.numeric(data[[first]]))), , drop = FALSE]
    names(data) <- sub("^[#@;]+", "", names(data))
  }
  data <- as.data.frame(data)
  if (!"ID" %in% names(data)) stop("The covariate data has no ID column.", call. = FALSE)
  data$ID <- suppressWarnings(as.numeric(data$ID))
  data <- data[!is.na(data$ID), , drop = FALSE]
  data[!duplicated(data$ID), , drop = FALSE]
}

# Two groups of subject IDs for one covariate.
# Categorical (non-numeric, or at most `max_levels` distinct values): every
# level is compared with the reference (first) level. Continuous: the lowest
# quartile Q1 is compared with the highest quartile Q4 (quartiles defined once
# on the subjects of the covariate data).
.nm_w_groups <- function(cov, covariate, max_levels) {
  v <- cov[[covariate]]
  ok <- !is.na(v)
  ids <- cov$ID[ok]; v <- v[ok]
  num <- suppressWarnings(as.numeric(v))
  categorical <- any(is.na(num)) || length(unique(v)) <= max_levels
  if (categorical) {
    lv <- sort(unique(v))
    if (length(lv) < 2L) return(list())
    lapply(lv[-1], function(l) list(
      covariate = covariate, type = "categorical",
      comparison = paste0(l, " vs ", lv[1]),
      group_1 = ids[v == l], group_0 = ids[v == lv[1]],
      label_1 = paste0(covariate, " = ", l), label_0 = paste0(covariate, " = ", lv[1])))
  } else {
    q <- dplyr::ntile(num, 4)
    list(list(covariate = covariate, type = "continuous", comparison = "Q4 vs Q1",
              group_1 = ids[q == 4], group_0 = ids[q == 1],
              label_1 = "Q4 (high)", label_0 = "Q1 (low)"))
  }
}


# ---- Main function --------------------------------------------------------------

#' Distribution-based covariate assessment with the Wasserstein distance
#'
#' Implements the method of Simon, Hulot and von Fabeck (2026). For each run
#' and each covariate, the empirical Bayes estimates of one ETA (read from the
#' `.phi` file) are split into two groups of subjects and compared with the
#' 1-Wasserstein distance `W` and a permutation test, together with a Welch
#' t-test for reference. A covariate that the model accounts for leaves no
#' residual distributional signal in the ETA: comparing a base and a final
#' model shows whether the signal has been removed. A Wasserstein-based R2
#' (`R2W`) compares the distance of each run's ETA distribution to 0 with the
#' one of a reference run.
#'
#' @param runs Character vector of run directories (`.reg`, `.nm7` or any
#'   directory holding the `.phi` file) or run path prefixes. Names, if any,
#'   are used as labels; otherwise the run names are used.
#' @param covariates Character vector of covariate columns of the data set.
#' @param eta The ETA to analyse: a PK parameter name (e.g. `"CL"`, identified
#'   from the control stream with [nm_identify_parameters()]), `"ETA(n)"` or
#'   an integer index.
#' @param data Covariate data: a data.frame or the path of a CSV file with an
#'   `ID` column (one row per subject is kept, the first one). If `NULL`
#'   (default), the data set named in `$DATA` of the first run's control
#'   stream is read.
#' @param ref Reference run for `R2W`: index or label (default 1).
#' @param n_perm Number of permutations (default 10000).
#' @param seed Random seed for the permutation tests (default 2025). Each test
#'   uses `seed + k`, so results are reproducible and independent of the
#'   order of the other tests. The user's random state is restored.
#' @param max_levels A numeric covariate with at most this many distinct
#'   values is treated as categorical (default 5).
#'
#' @return An object of class `nm_wasserstein`, a list with:
#'   \describe{
#'     \item{`tests`}{one row per run x covariate comparison: `run`,
#'       `covariate`, `type`, `comparison`, `n_1`, `n_0`, `W`, `p_W`,
#'       `t_stat`, `p_t`.}
#'     \item{`r2w`}{one row per run: `run`, `eta`, `n`, `W_to_zero`, `R2W`
#'       (`1 - W_to_zero / W_to_zero[ref]`), `omega` and `shrinkage_pct`
#'       (eta-shrinkage `1 - sd(EBE)/sqrt(omega)`, when the `.ext` is found).
#'       A large difference in shrinkage between runs can make the ETA
#'       distribution narrower without any real gain.}
#'     \item{`data`}{the EBE merged with the comparison groups (long format),
#'       used by [plot.nm_wasserstein()].}
#'     \item{`settings`}{the arguments used.}
#'   }
#'
#' @references Simon N, Hulot JS, von Fabeck K. Distribution-based covariate
#'   assessment using Wasserstein distance in population pharmacokinetic
#'   models. *Front Pharmacol.* 2026;17:1804989.
#'   \doi{10.3389/fphar.2026.1804989}
#'
#' @examples
#' \dontrun{
#' w <- nm_wasserstein(
#'   runs       = c(BASE = "C:/nm76g64/run/mtx/base.reg",
#'                  FINAL = "C:/nm76g64/run/mtx/final.reg"),
#'   covariates = c("CLCR", "CT01"),
#'   eta        = "CL"
#' )
#' w                      # tests and R2W
#' plot(w)                # ECDF of the ETA by group, per run and covariate
#' plot(w, type = "r2w")  # Wasserstein-based R2
#' }
#' @seealso [nm_identify_parameters()], [nm_compare()]
#' @export
nm_wasserstein <- function(runs, covariates, eta = "CL", data = NULL, ref = 1,
                           n_perm = 10000L, seed = 2025L, max_levels = 5L) {
  if (!length(runs)) stop("`runs` is empty.", call. = FALSE)
  if (!length(covariates)) stop("`covariates` is empty.", call. = FALSE)
  n_perm <- as.integer(n_perm)
  if (is.na(n_perm) || n_perm < 1L) stop("`n_perm` must be a positive integer.", call. = FALSE)

  files <- lapply(runs, .nm_w_run_files)
  labels <- names(runs)
  if (is.null(labels)) labels <- rep("", length(runs))
  labels[!nzchar(labels)] <- vapply(files[!nzchar(labels)], `[[`, "", "stem")
  if (anyDuplicated(labels)) stop("Run labels must be unique: ", paste(labels, collapse = ", "),
                                  ". Name the `runs` vector.", call. = FALSE)

  cov <- .nm_w_covariates(data, files[[1]]$ctl)
  miss <- setdiff(covariates, names(cov))
  if (length(miss)) stop("Covariate(s) not found in the data: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  groups <- unlist(lapply(covariates, function(cv) .nm_w_groups(cov, cv, max_levels)),
                   recursive = FALSE)
  if (!length(groups)) stop("No covariate defines two groups of subjects.", call. = FALSE)

  ebe <- vector("list", length(runs)); r2 <- vector("list", length(runs))
  for (i in seq_along(runs)) {
    f <- files[[i]]
    k <- .nm_w_eta_index(eta, f$ctl)
    phi <- .nm_w_read_phi(f$phi)
    col <- grep(sprintf("^ETA\\(%d\\)$", k), names(phi), value = TRUE)
    if (!length(col)) stop("ETA(", k, ") not found in ", f$phi, call. = FALSE)
    e <- data.frame(run = labels[i], ID = as.numeric(phi$ID),
                    ETA = as.numeric(phi[[col]]), stringsAsFactors = FALSE)
    ebe[[i]] <- e
    om <- .nm_w_omega(f$ext, k)
    r2[[i]] <- data.frame(run = labels[i], eta = sprintf("ETA(%d)", k),
                          n = sum(is.finite(e$ETA)),
                          W_to_zero = .nm_w1_to_zero(e$ETA), omega = om,
                          shrinkage_pct = if (is.finite(om) && om > 0)
                            100 * (1 - stats::sd(e$ETA, na.rm = TRUE) / sqrt(om)) else NA_real_,
                          stringsAsFactors = FALSE)
  }
  r2w <- do.call(rbind, r2)
  ref_i <- if (is.numeric(ref)) as.integer(ref) else match(ref, labels)
  if (is.na(ref_i) || ref_i < 1L || ref_i > length(labels))
    stop("`ref` must be a run index or label.", call. = FALSE)
  r2w$R2W <- 1 - r2w$W_to_zero / r2w$W_to_zero[ref_i]
  r2w <- r2w[, c("run", "eta", "n", "W_to_zero", "R2W", "omega", "shrinkage_pct")]

  tests <- list(); long <- list(); k <- 0L
  for (i in seq_along(runs)) for (g in seq_along(groups)) {
    gr <- groups[[g]]
    e <- ebe[[i]]
    x1 <- e$ETA[e$ID %in% gr$group_1]
    x0 <- e$ETA[e$ID %in% gr$group_0]
    ws <- .nm_w_perm(x1, x0, n_perm = n_perm, seed = seed + (i - 1L) * length(groups) + g)
    tt <- tryCatch(stats::t.test(x1, x0), error = function(err) NULL)
    k <- k + 1L
    tests[[k]] <- data.frame(
      run = labels[i], covariate = gr$covariate, type = gr$type,
      comparison = gr$comparison, n_1 = sum(is.finite(x1)), n_0 = sum(is.finite(x0)),
      W = ws$W, p_W = ws$p,
      t_stat = if (is.null(tt)) NA_real_ else unname(tt$statistic),
      p_t = if (is.null(tt)) NA_real_ else tt$p.value, stringsAsFactors = FALSE)
    long[[k]] <- rbind(
      data.frame(run = labels[i], covariate = gr$covariate, comparison = gr$comparison,
                 group = gr$label_1, ETA = x1, stringsAsFactors = FALSE),
      data.frame(run = labels[i], covariate = gr$covariate, comparison = gr$comparison,
                 group = gr$label_0, ETA = x0, stringsAsFactors = FALSE))
  }
  tests <- do.call(rbind, tests); rownames(tests) <- NULL
  long <- do.call(rbind, long); rownames(long) <- NULL
  long$run <- factor(long$run, levels = labels)

  structure(list(
    tests = tests, r2w = r2w, data = long,
    settings = list(runs = stats::setNames(unname(unlist(runs)), labels),
                    covariates = covariates, eta = eta, ref = labels[ref_i],
                    n_perm = n_perm, seed = seed, max_levels = max_levels)),
    class = "nm_wasserstein")
}


# ---- S3 methods -------------------------------------------------------------------

.nm_w_fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001",
                                                         formatC(p, format = "f", digits = 3)))

#' @param x An `nm_wasserstein` object.
#' @param ... Ignored.
#' @rdname nm_wasserstein
#' @export
print.nm_wasserstein <- function(x, ...) {
  cat("Wasserstein covariate assessment (", x$settings$n_perm, " permutations)\n", sep = "")
  cat("ETA: ", paste(unique(x$r2w$eta), collapse = ", "), " | reference run: ",
      x$settings$ref, "\n\n", sep = "")
  t <- x$tests
  t$W <- signif(t$W, 3); t$t_stat <- round(t$t_stat, 2)
  t$p_W <- .nm_w_fmt_p(t$p_W); t$p_t <- .nm_w_fmt_p(t$p_t)
  print(t, row.names = FALSE)
  cat("\nWasserstein-based R2 (R2W = 1 - W/W_ref, W = distance of the ETA to 0):\n")
  r <- x$r2w
  r$W_to_zero <- signif(r$W_to_zero, 3); r$R2W <- round(r$R2W, 3)
  r$omega <- signif(r$omega, 3); r$shrinkage_pct <- round(r$shrinkage_pct, 1)
  print(r, row.names = FALSE)
  invisible(x)
}

#' Plots for `nm_wasserstein` objects
#'
#' @param x An `nm_wasserstein` object.
#' @param type `"ecdf"` (default): empirical cumulative distribution of the
#'   ETA in the two groups, one panel per run and comparison, annotated with
#'   W, its p-value, t and its p-value. `"r2w"`: bar chart of the
#'   Wasserstein-based R2 per run.
#' @param covariate Optional covariate(s) to show (ECDF plot only).
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.nm_wasserstein <- function(x, type = c("ecdf", "r2w"), covariate = NULL, ...) {
  type <- match.arg(type)
  if (type == "r2w") {
    d <- x$r2w
    d$run <- factor(d$run, levels = d$run)
    return(ggplot2::ggplot(d, ggplot2::aes(x = .data$run, y = .data$R2W)) +
      ggplot2::geom_col(fill = "gray60", colour = "black", linewidth = 0.3) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
      ggplot2::scale_y_continuous(labels = function(v) paste0(round(100 * v), "%")) +
      ggplot2::labs(title = "Wasserstein-based R2",
                    subtitle = paste0("R2W = 1 - W / W_ref, reference run: ", x$settings$ref),
                    x = "Run", y = "R2W") +
      ggplot2::theme_bw(base_size = 11))
  }
  d <- x$data
  tt <- x$tests
  if (!is.null(covariate)) {
    d <- d[d$covariate %in% covariate, , drop = FALSE]
    tt <- tt[tt$covariate %in% covariate, , drop = FALSE]
    if (!nrow(d)) stop("No data for covariate(s): ", paste(covariate, collapse = ", "), call. = FALSE)
  }
  d$panel <- paste0(d$covariate, " (", d$comparison, ")")
  tt$panel <- paste0(tt$covariate, " (", tt$comparison, ")")
  plab <- function(nm, p) ifelse(!is.na(p) & p < 0.001, paste0(nm, " < 0.001"),
                                 paste0(nm, " = ", .nm_w_fmt_p(p)))
  tt$label <- paste0("W = ", signif(tt$W, 3), ", ", plab("p W", tt$p_W),
                     "\nt = ", round(tt$t_stat, 2), ", ", plab("p t", tt$p_t))
  tt$run <- factor(tt$run, levels = levels(d$run))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$ETA, colour = .data$group)) +
    ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
    ggplot2::geom_label(data = tt, ggplot2::aes(x = Inf, y = 0, label = .data$label),
                        inherit.aes = FALSE, hjust = 1.02, vjust = -0.1, size = 2.8,
                        fill = "white") +
    ggplot2::facet_grid(panel ~ run) +
    ggplot2::labs(x = paste0("Individual deviation (", unique(x$r2w$eta)[1], ")"),
                  y = "Empirical cumulative distribution", colour = "Group") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}
