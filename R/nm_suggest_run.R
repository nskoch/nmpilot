# =============================================================================
# nm_suggest_run.R
# Closed loop after nm_suggest_model(): write the control streams of the best
# candidates, run them with nm_run(), and compare the NONMEM ranking with the
# data-driven score. Covariate screening on the individual estimates.
# =============================================================================

#' Write, run and rank the candidates of a model suggestion
#'
#' Writes the control streams of the chosen candidates of
#' [nm_suggest_model()] with [nm_write_ctl()] (`from = x`), launches them with
#' [nm_run()], reads their objective function value and ranks them by BIC
#' (OFV + number of estimated parameters x log(number of observations)). The
#' rank by BIC is compared with the rank of the data-driven score (Spearman
#' correlation), which shows whether the suggestion anticipated the NONMEM
#' result. A lag time of a candidate is added after `$PK` as `ALAG1 =
#' THETA(n)` with its initial value, since [nm_write_ctl()] writes none.
#'
#' Each candidate gets a dataset that matches its compartment numbering:
#' doses in the depot (CMT 1) and observations in the central compartment
#' (CMT 2) for first-order and transit absorption (ADVAN2/4/12 or ODE);
#' doses and observations in the central compartment (CMT 1) for IV and
#' zero-order input (ADVAN1/3/11), with `RATE = -2` on the doses for a
#' modelled zero-order duration (`D1`). When the original dataset does not
#' match, a copy `<data>_<prefix><n>.csv` is written in `analysis_dir`
#' (records copied verbatim, only CMT and RATE changed); the original file
#' is never modified. Datasets with several observed compartments (parent
#' and metabolite, PK and PD) or doses by several routes are not handled:
#' the function stops.
#'
#' @param x An `nm_model_suggestion` object.
#' @param analysis_dir Analysis directory where the control streams are
#'   written and run. It must contain the dataset, or `data` must point to
#'   it.
#' @param nmfe_path Path to the nmfe launcher (see [nm_run()]).
#' @param candidates Candidate numbers (default: the first 3).
#' @param prefix Prefix of the control stream names (`<prefix><n>.ctl`).
#' @param data Dataset file name written in `$DATA` (default: the file
#'   analysed by [nm_suggest_model()]), relative to `analysis_dir` or
#'   absolute.
#' @param run If `FALSE`, only the control streams are written.
#' @param overwrite Replace existing control streams and adapted datasets,
#'   and run the candidates again even when a previous run of the same name
#'   has terminated (otherwise [nm_run()] skips it and its old results are
#'   read).
#' @param verbose Print progress messages.
#' @param n_cores Number of candidates run at the same time, passed to
#'   [nm_run()] (default 1).
#' @param ... Further arguments passed to [nm_write_ctl()] (e.g. `method`).
#' @return An object of class `nm_suggest_run`: a data.frame with one row
#'   per candidate (`candidate`, `model`, `score`, `rank_score`, `ctl`,
#'   `data` (dataset used), `ofv`, `n_par`, `n_obs`, `bic`, `rank_bic`,
#'   `min_ok`), with the
#'   Spearman correlation between score and BIC ranks in attribute
#'   `spearman` and the [nm_run()] log in attribute `log`.
#' @seealso [nm_suggest_model()], [nm_write_ctl()], [nm_run()]
#' @examples
#' \dontrun{
#' s <- nm_suggest_model("C:/nm76g64/run/etude/data.csv")
#' r <- nm_suggest_run(s, analysis_dir = "C:/nm76g64/run/etude",
#'                     nmfe_path = "C:/nm76g64/run/nmfe76.bat", candidates = 1:4)
#' r
#' }
#' @export
nm_suggest_run <- function(x, analysis_dir, nmfe_path = NULL, candidates = 1:3,
                           prefix = "sugg", data = NULL, run = TRUE,
                           overwrite = FALSE, verbose = TRUE, n_cores = 1L, ...) {
  if (!inherits(x, "nm_model_suggestion")) stop("`x` must be an nm_model_suggestion object.")
  if (!dir.exists(analysis_dir)) stop("`analysis_dir` does not exist: ", analysis_dir)
  candidates <- as.integer(candidates)
  if (any(!candidates %in% seq_along(x$specs)))
    stop("`candidates` must be between 1 and ", length(x$specs), ".")
  if (run && is.null(nmfe_path)) stop("`nmfe_path` is needed when `run = TRUE`.")
  if (is.null(data)) {
    if (is.null(x$data_summary$data_file)) stop("`data` is needed: the suggestion was made from a data.frame.")
    data <- basename(x$data_summary$data_file)
    if (!file.exists(file.path(analysis_dir, data)))
      stop("Dataset ", data, " not found in `analysis_dir`; give its path in `data`.")
  }
  data_path <- if (grepl("^([A-Za-z]:|/|\\\\)", data)) data else file.path(analysis_dir, data)
  if (!file.exists(data_path)) stop("Dataset not found: ", data_path)
  ctl_files <- character(0); data_used <- character(0)
  for (i in candidates) {
    f <- file.path(analysis_dir, sprintf("%s%d.ctl", prefix, i))
    if (file.exists(f) && !overwrite) stop("File exists: ", f, " (use overwrite = TRUE).")
    dp <- .nm_sr_prepare_data(data_path, x$specs[[i]], analysis_dir,
                              tag = sprintf("%s%d", prefix, i), overwrite = overwrite, verbose = verbose)
    d_rel <- if (identical(dp, data_path)) data else basename(dp)
    ctl <- nm_write_ctl(file = NULL, data = d_rel, input = dp, from = x, candidate = i,
                        verbose = FALSE, ...)
    ctl <- .nm_sr_add_alag(unclass(ctl), x$specs[[i]]$alag)
    writeLines(ctl, f)
    ctl_files <- c(ctl_files, f); data_used <- c(data_used, dp)
    if (verbose) message("Written: ", f)
  }
  n_obs <- x$data_summary$n_obs
  cand <- x$candidates[candidates, ]
  out <- data.frame(candidate = candidates,
                    model = sprintf("%d cpt | %s | %s | %s", cand$cpts, cand$input, cand$elimination, cand$error),
                    score = cand$score, rank_score = rank(-cand$score, ties.method = "min"),
                    ctl = ctl_files, data = data_used, stringsAsFactors = FALSE)
  out$n_par <- vapply(candidates, function(i) {
    sp <- x$specs[[i]]
    length(sp$theta) + length(sp$alag) + length(sp$iiv) + length(sp$sigma_used)
  }, numeric(1))
  out$n_obs <- n_obs
  out$ofv <- NA_real_; out$min_ok <- NA
  log <- NULL
  if (run) {
    # with overwrite, the control streams have just been rewritten: previous
    # results of the same runs must not be reused
    log <- nm_run(ctl_files, nmfe_path, analysis_dir, skip_terminated = !overwrite,
                  verbose = verbose, n_cores = n_cores)
    m <- match(tools::file_path_sans_ext(basename(ctl_files)), log$run)
    out$ofv <- log$ofv[m]
    out$min_ok <- log$min_ok[m]
  }
  out$bic <- out$ofv + out$n_par * log(n_obs)
  out$rank_bic <- ifelse(is.finite(out$bic), rank(out$bic, ties.method = "min", na.last = "keep"), NA)
  ok <- is.finite(out$bic)
  rho <- if (sum(ok) >= 3) suppressWarnings(stats::cor(out$rank_score[ok], out$rank_bic[ok], method = "spearman")) else NA
  attr(out, "spearman") <- rho
  attr(out, "log") <- log
  class(out) <- c("nm_suggest_run", "data.frame")
  out
}

#' @rdname nm_suggest_run
#' @export
print.nm_suggest_run <- function(x, ...) {
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  d <- d[, c("candidate", "model", "score", "rank_score", "ofv", "n_par", "bic", "rank_bic", "min_ok")]
  cat("Candidates run with NONMEM\n")
  print(d, row.names = FALSE)
  rho <- attr(x, "spearman")
  cat(sprintf("\nSpearman correlation between score rank and BIC rank: %s\n",
              if (is.finite(rho)) sprintf("%.2f", rho) else "NA (fewer than 3 successful runs)"))
  invisible(x)
}

# Dataset matching the compartment numbering of a candidate (see the
# documentation of nm_suggest_run). Returns the path of the dataset to use:
# the original one when nothing has to change, otherwise a copy where only
# CMT (and RATE for a modelled zero-order duration) are changed.
.nm_sr_prepare_data <- function(data_path, sp, analysis_dir, tag, sep = ",",
                                overwrite = FALSE, verbose = TRUE) {
  st <- .nm_ctl_struct(sp$cpts, sp$route, sp$absorption, sp$elimination, sp$transit_n)
  target_obs  <- st$central_cmt
  target_dose <- if (st$first_ord || st$transit) 1L else st$central_cmt
  need_rate   <- isTRUE(st$zero_d1)
  df <- utils::read.csv(data_path, sep = sep, colClasses = "character", check.names = FALSE,
                        na.strings = character(0), quote = "", comment.char = "")
  up <- toupper(sub("^[#@]", "", trimws(names(df))))
  col <- function(n) { j <- which(up == n); if (length(j)) j[1] else NA_integer_ }
  i_evid <- col("EVID"); i_amt <- col("AMT"); i_cmt <- col("CMT"); i_rate <- col("RATE")
  if (is.na(i_evid) && is.na(i_amt)) stop("No EVID or AMT column in ", data_path)
  first <- trimws(df[[1]])
  data_row <- !substr(first, 1, 1) %in% c("#", "@") & !grepl("^[A-Za-z]", first)
  num <- function(j) suppressWarnings(as.numeric(df[[j]]))
  is_dose <- data_row & if (!is.na(i_evid)) num(i_evid) %in% c(1, 4) else (!is.na(num(i_amt)) & num(i_amt) > 0)
  is_obs  <- data_row & if (!is.na(i_evid)) num(i_evid) %in% 0 else !is_dose
  changed <- FALSE
  if (!is.na(i_cmt)) {
    cm <- num(i_cmt)
    obs_c <- unique(cm[is_obs]); obs_c <- obs_c[is.finite(obs_c)]
    if (length(obs_c) > 1)
      stop("Several observed compartments in the dataset (CMT ", paste(sort(obs_c), collapse = ", "),
           "): nm_suggest_run() handles a single observed compartment.")
    dose_c <- unique(cm[is_dose]); dose_c <- dose_c[is.finite(dose_c)]
    if (length(dose_c) > 1)
      stop("Doses in several compartments (CMT ", paste(sort(dose_c), collapse = ", "),
           "): nm_suggest_run() handles a single route of administration.")
    newc <- df[[i_cmt]]
    if (length(dose_c) && dose_c != target_dose) newc[is_dose & cm %in% dose_c] <- as.character(target_dose)
    # every other record in the observed compartment (EVID 0, 2...) follows it
    if (length(obs_c) && obs_c != target_obs) newc[data_row & !is_dose & cm %in% obs_c] <- as.character(target_obs)
    if (!identical(newc, df[[i_cmt]])) { df[[i_cmt]] <- newc; changed <- TRUE }
  }
  if (need_rate) {
    if (is.na(i_rate)) {
      rate <- ifelse(is_dose, "-2", ifelse(data_row, "0", "."))
      pos <- if (!is.na(i_amt)) i_amt else ncol(df)
      df <- cbind(df[, seq_len(pos), drop = FALSE], RATE = rate,
                  df[, setdiff(seq_len(ncol(df)), seq_len(pos)), drop = FALSE],
                  stringsAsFactors = FALSE)
      changed <- TRUE
    } else {
      r <- df[[i_rate]]
      if (any(num(i_rate)[is_dose] != -2, na.rm = TRUE) || any(is.na(num(i_rate)[is_dose]))) {
        r[is_dose] <- "-2"; df[[i_rate]] <- r; changed <- TRUE
      }
    }
  }
  if (!changed) return(data_path)
  out <- file.path(analysis_dir, sprintf("%s_%s.csv", tools::file_path_sans_ext(basename(data_path)), tag))
  if (file.exists(out) && !overwrite) stop("File exists: ", out, " (use overwrite = TRUE).")
  utils::write.table(df, out, sep = sep, quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
  if (verbose) message(sprintf("Dataset adapted to candidate %s (CMT %d -> doses, %d -> observations%s): %s",
                               tag, target_dose, target_obs, if (need_rate) ", RATE = -2 on doses" else "", out))
  out
}

# Add ALAG1 = THETA(n) after the $PK line and its $THETA (lag-time candidates).
.nm_sr_add_alag <- function(ctl, alag) {
  if (is.null(alag)) return(ctl)
  n_theta <- sum(lengths(regmatches(ctl, gregexpr("THETA\\(\\d+\\)", ctl))) > 0)
  used <- as.integer(unlist(regmatches(ctl, gregexpr("(?<=THETA\\()\\d+(?=\\))", ctl, perl = TRUE))))
  k <- if (length(used)) max(used) + 1L else n_theta + 1L
  ipk <- grep("^\\s*\\$PK", ctl)[1]
  if (is.na(ipk)) stop("No $PK record: the lag time cannot be added.")
  ctl <- append(ctl, sprintf("  ALAG1 = THETA(%d)   ; lag time (nm_suggest_model)", k), after = ipk)
  last_theta <- max(grep("^\\s*\\$THETA", ctl))
  j <- last_theta
  while (j < length(ctl) && !grepl("^\\s*\\$", ctl[j + 1])) j <- j + 1L
  while (j > last_theta && !nzchar(trimws(ctl[j]))) j <- j - 1L
  append(ctl, sprintf(" (0, %s) ; ALAG1", format(alag)), after = j)
}


#' Covariate screening on individual PK estimates before modelling
#'
#' Applies the Wasserstein-distance comparison of [nm_wasserstein()] to the
#' individual estimates obtained by [nm_suggest_model()] (NCA clearance,
#' central volume and absorption parameter of the rich intervals, log
#' scale), before any NONMEM run. For each covariate, subjects are split into
#' groups (each level against the first for a categorical covariate; Q4
#' against Q1, or all pairs of quartiles, for a continuous one) and the
#' distributions of the log estimates are compared with the 1-Wasserstein
#' distance, a permutation test and a Welch t-test. It points to covariates
#' worth testing in the model; it is no substitute for the covariate model
#' (individual NCA estimates carry their own estimation noise, and only
#' subjects with a rich interval are used).
#'
#' @param x An `nm_model_suggestion` object (individual estimates are needed,
#'   i.e. at least 3 rich intervals).
#' @param covariates Character vector of covariate columns.
#' @param data Covariate data (data.frame or CSV path with an `ID` column).
#'   Default: the dataset analysed by [nm_suggest_model()].
#' @param parameters Individual estimates to analyse among `"CL"`, `"Vc"`,
#'   `"KA"`, `"DUR"`, `"MTT"` (default `c("CL", "Vc")`).
#' @param n_perm Number of permutations (default 2000).
#' @param seed Random seed of the permutation tests; the user's random state
#'   is restored.
#' @param max_levels A numeric covariate with at most this many distinct
#'   values is treated as categorical (default 5).
#' @param continuous `"extremes"` (Q4 vs Q1, default) or `"pairwise"`.
#' @return A data.frame of class `nm_suggest_covariates`: `parameter`,
#'   `covariate`, `type`, `comparison`, `n_1`, `n_0`, `ratio_median`
#'   (median of group 1 / median of group 0), `W`, `p_W`, `p_t`, `p_W_adj`
#'   (Bonferroni over the comparisons of a covariate and parameter).
#' @references Simon N, Hulot JS, von Fabeck K. Distribution-based covariate
#'   assessment using Wasserstein distance in population pharmacokinetic
#'   models. *Front Pharmacol.* 2026;17:1804989.
#'   \doi{10.3389/fphar.2026.1804989}
#' @seealso [nm_suggest_model()], [nm_wasserstein()]
#' @examples
#' ex <- system.file("extdata", "example", "example_data.csv", package = "nmpilot")
#' s <- nm_suggest_model(ex, verbose = FALSE, design_eval = FALSE, max_profiles = 15)
#' nm_suggest_covariates(s, covariates = c("WT", "GENO"), n_perm = 500)
#' @export
nm_suggest_covariates <- function(x, covariates, data = NULL, parameters = c("CL", "Vc"),
                                  n_perm = 2000L, seed = 2025L, max_levels = 5L,
                                  continuous = c("extremes", "pairwise")) {
  if (!inherits(x, "nm_model_suggestion")) stop("`x` must be an nm_model_suggestion object.")
  continuous <- match.arg(continuous)
  ind <- x$individual
  if (is.null(ind) || !nrow(ind))
    stop("No individual estimates: at least 3 rich intervals are needed (sparse data).")
  parameters <- intersect(parameters, c("CL", "Vc", "KA", "DUR", "MTT"))
  if (!length(parameters)) stop("`parameters` must be among CL, Vc, KA, DUR, MTT.")
  if (is.null(data)) {
    data <- x$data_summary$data_file
    if (is.null(data)) stop("`data` is needed: the suggestion was made from a data.frame.")
  }
  cov <- .nm_w_covariates(data, NA)
  miss <- setdiff(covariates, names(cov))
  if (length(miss)) stop("Covariate(s) not found in the data: ", paste(miss, collapse = ", "))
  per <- stats::aggregate(ind[, parameters, drop = FALSE], by = list(ID = ind$ID),
                          FUN = function(v) { v <- v[is.finite(v) & v > 0]; if (length(v)) stats::median(v) else NA })
  per$ID <- suppressWarnings(as.numeric(per$ID))
  cov <- cov[cov$ID %in% per$ID, , drop = FALSE]
  rows <- list(); k <- 0L
  for (cv in covariates) {
    groups <- .nm_w_groups(cov, cv, max_levels, continuous)
    for (p in parameters) {
      v <- stats::setNames(log(per[[p]]), per$ID)
      res_p <- list()
      for (g in groups) {
        k <- k + 1L
        a <- v[as.character(g$group_1)]; b <- v[as.character(g$group_0)]
        a <- a[is.finite(a)]; b <- b[is.finite(b)]
        pw <- .nm_w_perm(a, b, n_perm = n_perm, seed = seed + k)
        pt <- if (length(a) >= 2 && length(b) >= 2)
          tryCatch(stats::t.test(a, b)$p.value, error = function(e) NA_real_) else NA_real_
        res_p[[length(res_p) + 1]] <- data.frame(
          parameter = p, covariate = cv, type = g$type, comparison = g$comparison,
          n_1 = length(a), n_0 = length(b),
          ratio_median = if (length(a) && length(b)) round(exp(stats::median(a) - stats::median(b)), 3) else NA,
          W = pw$W, p_W = pw$p, p_t = pt, stringsAsFactors = FALSE)
      }
      if (length(res_p)) {
        r <- do.call(rbind, res_p)
        r$p_W_adj <- pmin(1, r$p_W * nrow(r))
        rows[[length(rows) + 1]] <- r
      }
    }
  }
  if (!length(rows)) stop("No comparison could be made (covariates with a single level?).")
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("nm_suggest_covariates", "data.frame")
  out
}

#' @rdname nm_suggest_covariates
#' @param ... Ignored.
#' @export
print.nm_suggest_covariates <- function(x, ...) {
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  d$W <- signif(d$W, 3)
  for (v in c("p_W", "p_t", "p_W_adj")) d[[v]] <- .nm_w_fmt_p(d[[v]])
  cat("Covariate screening on individual estimates (log scale, Wasserstein distance)\n")
  print(d, row.names = FALSE)
  invisible(x)
}
