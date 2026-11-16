# =============================================================================
# nm_compare(): statistics and plots for the direct comparison of two runs
# =============================================================================
#
#   - Vuong (1989) and Clarke (2007) tests on the individual contributions to
#     the objective function (iOFV of the .phi files): non-nested models
#   - likelihood ratio test when the user declares the models nested
#   - subject bootstrap of the delta iOFV (no re-estimation)
#   - influence: number of subjects whose removal reverses the sign of dOFV
#   - delta iOFV versus covariates
#   - parameters and ETAs of the two runs side by side
#   - run-level quality criteria (AIC, BIC, RSE, condition number, shrinkage)
#   - smoothers with subject-bootstrap bands for the gain map and the
#     superimposed diagnostics
# =============================================================================


# ---- Run information ----------------------------------------------------------

# Run metadata and parameters through nm_extract_runs() (NULL when unavailable).
.nm_cmp_runinfo <- function(path, paths) {
  stem <- NULL
  if (!is.null(paths$phi) && file.exists(paths$phi))
    stem <- tools::file_path_sans_ext(basename(paths$phi))
  if (dir.exists(path)) {
    d <- path
  } else {
    d <- dirname(path); stem <- tools::file_path_sans_ext(basename(path))
  }
  r <- tryCatch(suppressWarnings(suppressMessages(
    nm_extract_runs(d, runs = stem, verbose = FALSE))), error = function(e) NULL)
  if (is.null(r) || !nrow(r$runs)) return(NULL)
  run <- as.data.frame(r$runs[1, , drop = FALSE])
  pr <- as.data.frame(r$params[r$params$run == run$run, , drop = FALSE])
  cor_f <- list.files(run$dir, pattern = "(?i)\\.cor$", full.names = TRUE)
  cor_f <- cor_f[tools::file_path_sans_ext(basename(cor_f)) == run$run]
  list(run = run, params = pr,
       cond = if (length(cor_f)) .nm_cmp_condition(cor_f[1]) else NA_real_,
       ctl = { f <- list.files(run$dir, pattern = "(?i)\\.ctl$", full.names = TRUE)
               f <- f[tools::file_path_sans_ext(basename(f)) == run$run]
               if (length(f)) f[1] else NA_character_ })
}

# Condition number of the correlation matrix of the estimates (.cor file;
# NONMEM writes the standard errors on the diagonal).
.nm_cmp_condition <- function(file) {
  lines <- .nm_read_lines(file)
  tb <- grep("^TABLE NO", lines)
  if (length(tb)) lines <- lines[(max(tb) + 1L):length(lines)]
  m <- tryCatch(utils::read.table(text = lines, header = TRUE, check.names = FALSE,
                                  stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(m) || ncol(m) < 2L) return(NA_real_)
  x <- as.matrix(m[, -1, drop = FALSE])
  keep <- which(diag(x) != 0)
  if (length(keep) < 2L) return(NA_real_)
  x <- x[keep, keep, drop = FALSE]
  diag(x) <- 1
  ev <- eigen((x + t(x)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) return(Inf)
  max(ev) / min(ev)
}


# ---- Tests on the individual contributions --------------------------------------

# d_i = (iOFV_1 - iOFV_2) / 2 = log L_2,i - log L_1,i : positive favours run 2.
.nm_cmp_tests <- function(subj, k, ofv, n_obs, nested, n_boot, seed, labels) {
  ok <- is.finite(subj$iOFV_1) & is.finite(subj$iOFV_2)
  d <- (subj$iOFV_1[ok] - subj$iOFV_2[ok]) / 2
  N <- length(d)
  if (N < 5L) return(NULL)
  dk <- k[2] - k[1]
  fav <- function(x) ifelse(is.na(x), "-", ifelse(x > 0, labels[2], labels[1]))
  rows <- list()
  omega <- sqrt(mean((d - mean(d))^2))
  corr <- c(none = 0, AIC = dk, BIC = dk * log(N) / 2)
  if (!is.finite(dk)) corr <- corr["none"]
  for (cn in names(corr)) {
    z <- (sum(d) - corr[[cn]]) / (sqrt(N) * omega)
    rows[[length(rows) + 1L]] <- data.frame(
      test = "Vuong", correction = cn, statistic = round(z, 3),
      p = 2 * stats::pnorm(-abs(z)), favours = fav(z), stringsAsFactors = FALSE)
    dd <- d - corr[[cn]] / N
    m <- sum(dd != 0); b <- sum(dd > 0)
    p <- if (m > 0) stats::binom.test(b, m, 0.5)$p.value else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(
      test = "Clarke", correction = cn, statistic = b,
      p = p, favours = fav(b - m / 2), stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows)
  tab$statistic_is <- ifelse(tab$test == "Vuong", "z", sprintf("subjects favouring %s (of %d)", labels[2], N))
  tab <- tab[order(tab$test == "Clarke", match(tab$correction, c("none", "AIC", "BIC"))), ]
  if (isTRUE(nested) && is.finite(dk) && dk != 0 && all(is.finite(ofv))) {
    big <- if (dk > 0) 2L else 1L
    stat <- ofv[3L - big] - ofv[big]
    tab <- rbind(data.frame(test = "LRT (nested)", correction = sprintf("df = %d", abs(dk)),
                            statistic = round(stat, 3),
                            p = if (stat > 0) stats::pchisq(stat, abs(dk), lower.tail = FALSE) else 1,
                            favours = if (stat > 0 && stats::pchisq(stat, abs(dk), lower.tail = FALSE) < 0.05)
                              labels[big] else labels[3L - big],
                            statistic_is = "chi-square", stringsAsFactors = FALSE), tab)
  }
  rownames(tab) <- NULL
  tab$favours[tab$p >= 0.05 & !grepl("^LRT", tab$test)] <- "none (p >= 0.05)"

  # Subject bootstrap of the delta iOFV (iOFV_2 - iOFV_1), estimates not refitted
  di <- -2 * d
  set.seed(seed)
  bs <- vapply(seq_len(n_boot), function(i) sum(sample(di, N, replace = TRUE)), numeric(1))
  boot <- list(dofv = bs, ci = stats::quantile(bs, c(0.025, 0.5, 0.975), names = FALSE),
               p_run2_better = mean(bs < 0),
               p_run2_better_aic = if (is.finite(dk)) mean(bs + 2 * dk < 0) else NA_real_,
               n_boot = n_boot)

  # Influence: subjects to remove (largest contributions first) to reverse dOFV
  tot <- sum(di)
  s <- if (tot < 0) sort(di) else sort(di, decreasing = TRUE)
  cum <- tot - cumsum(s)
  n_rev <- if (tot == 0) 0L else which(sign(cum) != sign(tot) | cum == 0)[1]
  infl <- list(n_reverse = if (is.na(n_rev)) NA_integer_ else as.integer(n_rev),
               share_top5 = if (tot != 0) sum(s[seq_len(min(5L, N))]) / tot else NA_real_,
               n_subjects = N, omega = omega)
  list(table = tab, bootstrap = boot, influence = infl, k = k)
}


# ---- Delta iOFV versus covariates -------------------------------------------------

.nm_cmp_cov_tests <- function(subj, cov, covariates, max_levels = 5L) {
  if (!"diOFV" %in% names(subj)) return(NULL)
  m <- match(subj$ID, cov$ID)
  rows <- list()
  for (cv in covariates) {
    v <- cov[[cv]][m]
    num <- suppressWarnings(as.numeric(as.character(v)))
    is_num <- all(is.na(v) | !is.na(num))
    ok <- !is.na(v) & is.finite(subj$diOFV)
    if (sum(ok) < 6L || length(unique(v[ok])) < 2L) next
    cont <- is_num && length(unique(v[ok])) > max_levels
    if (cont) {
      ct <- suppressWarnings(stats::cor.test(subj$diOFV[ok], num[ok], method = "spearman", exact = FALSE))
      st <- unname(ct$estimate); p <- ct$p.value; what <- "Spearman rho"
    } else {
      kt <- stats::kruskal.test(subj$diOFV[ok], factor(v[ok]))
      g <- tapply(subj$diOFV[ok], factor(v[ok]), sum)
      st <- max(g) - min(g); p <- kt$p.value; what <- "range of summed diOFV by level"
    }
    rows[[length(rows) + 1L]] <- data.frame(covariate = cv, type = if (cont) "continuous" else "categorical",
                                            statistic = round(st, 3), statistic_is = what, p = p,
                                            stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$p_adj <- pmin(1, out$p * nrow(out))
  out[order(out$p), , drop = FALSE]
}


# ---- Parameters and ETAs -------------------------------------------------------------

.nm_cmp_params <- function(p1, p2) {
  if (is.null(p1) || is.null(p2)) return(NULL)
  keep <- function(p) {
    p <- p[p$type %in% c("THETA", "OMEGA", "SIGMA", "ETA_SD", "EPS_SD"), , drop = FALSE]
    p$key <- paste(p$type, p$name)
    p[!duplicated(p$key), c("key", "type", "name", "value", "se", "rse_pct", "fixed", "shrinkage")]
  }
  a <- keep(p1); b <- keep(p2)
  out <- merge(a, b, by = c("key", "type", "name"), all = TRUE, suffixes = c("_1", "_2"))
  out$ratio <- ifelse(is.finite(out$value_1) & out$value_1 != 0, out$value_2 / out$value_1, NA_real_)
  ord <- match(out$type, c("THETA", "OMEGA", "ETA_SD", "SIGMA", "EPS_SD"))
  out <- out[order(ord, out$name), setdiff(names(out), "key"), drop = FALSE]
  rownames(out) <- NULL
  out
}

# Paired EBEs for the ETAs with the same parameter name in both runs.
.nm_cmp_etas <- function(phi1, phi2, p1, p2) {
  if (is.null(phi1) || is.null(phi2) || is.null(p1) || is.null(p2)) return(NULL)
  nm_of <- function(p) {
    o <- p[p$type == "OMEGA", , drop = FALSE]
    stats::setNames(o$name, o$index)
  }
  n1 <- nm_of(p1); n2 <- nm_of(p2)
  common <- intersect(n1, n2)
  common <- common[!is.na(common) & nzchar(common)]
  if (!length(common)) return(NULL)
  out <- list()
  for (nm in common) {
    k1 <- names(n1)[n1 == nm][1]; k2 <- names(n2)[n2 == nm][1]
    c1 <- sprintf("ETA(%s)", k1); c2 <- sprintf("ETA(%s)", k2)
    if (!c1 %in% names(phi1) || !c2 %in% names(phi2)) next
    d <- merge(phi1[, c("ID", c1)], phi2[, c("ID", c2)], by = "ID")
    names(d) <- c("ID", "eta_1", "eta_2")
    d$parameter <- nm
    out[[nm]] <- d
  }
  if (!length(out)) return(NULL)
  out <- do.call(rbind, out); rownames(out) <- NULL
  out
}


# ---- Smoother with a subject-bootstrap band -------------------------------------------

.nm_cmp_band <- function(x, y, id, n_boot = 200L, n_grid = 60L, span = 0.75, log_x = FALSE, seed = 1L) {
  ok <- is.finite(x) & is.finite(y) & (!log_x | x > 0)
  x <- x[ok]; y <- y[ok]; id <- id[ok]
  if (length(unique(x)) < 4L) return(NULL)
  tx <- if (log_x) log(x) else x
  rng <- stats::quantile(tx, c(0.02, 0.98), names = FALSE)
  grid <- seq(rng[1], rng[2], length.out = n_grid)
  fit1 <- function(xx, yy) {
    f <- tryCatch(suppressWarnings(stats::loess(yy ~ xx, span = span, degree = 1,
                                                family = "gaussian",
                                                control = stats::loess.control(surface = "direct"))), error = function(e) NULL)
    if (is.null(f)) return(rep(NA_real_, n_grid))
    suppressWarnings(stats::predict(f, newdata = data.frame(xx = grid)))
  }
  est <- fit1(tx, y)
  ids <- unique(id); rows <- split(seq_along(id), id)
  set.seed(seed)
  B <- vapply(seq_len(n_boot), function(b) {
    s <- unlist(rows[as.character(sample(ids, length(ids), replace = TRUE))], use.names = FALSE)
    fit1(tx[s], y[s])
  }, numeric(n_grid))
  data.frame(x = if (log_x) exp(grid) else grid, fit = est,
             lo = apply(B, 1, stats::quantile, 0.025, na.rm = TRUE),
             hi = apply(B, 1, stats::quantile, 0.975, na.rm = TRUE))
}


# ---- Summary ------------------------------------------------------------------------------

#' @rdname nm_compare
#' @param object An `nm_compare` object.
#' @return `summary()`: a data.frame of decision criteria (class
#'   `nm_compare_summary`), one row per criterion with the value for each run
#'   and the run it favours.
#' @export
summary.nm_compare <- function(object, ...) {
  x <- object; lab <- x$meta$labels
  rows <- list()
  add <- function(criterion, v1, v2, favours, note = "") {
    rows[[length(rows) + 1L]] <<- data.frame(criterion = criterion, run1 = v1, run2 = v2,
                                             favours = favours, note = note, stringsAsFactors = FALSE)
  }
  low <- function(a, b, digits = 3) {
    if (!is.finite(a) || !is.finite(b)) "-" else if (signif(a, digits) == signif(b, digits)) "tie"
    else if (b < a) lab[2] else lab[1]
  }
  f <- function(v, d = 2) ifelse(is.finite(v), formatC(v, format = "f", digits = d), "NA")
  ofv <- x$meta$ofv; k <- x$meta$n_par; nobs <- x$meta$n_obs
  add("OFV", f(ofv[1]), f(ofv[2]), low(ofv[1], ofv[2]), sprintf("dOFV = %s", f(ofv[2] - ofv[1])))
  add("Estimated parameters", f(k[1], 0), f(k[2], 0), "-")
  aic <- ofv + 2 * k; bic <- ofv + k * log(nobs)
  add("AIC", f(aic[1]), f(aic[2]), low(aic[1], aic[2]))
  add("BIC (log n_obs)", f(bic[1]), f(bic[2]), low(bic[1], bic[2]), sprintf("n_obs = %d", nobs))
  t <- x$tests
  if (!is.null(t)) {
    for (i in seq_len(nrow(t$table))) {
      r <- t$table[i, ]
      add(paste0(r$test, if (!grepl("^LRT", r$test)) paste0(" (", r$correction, ")") else ""),
          "", "", r$favours, sprintf("%s = %s, p = %s", r$statistic_is, format(r$statistic),
                                     format.pval(r$p, digits = 3)))
    }
    b <- t$bootstrap
    add("Bootstrap P(dOFV < 0)", "", f(b$p_run2_better, 3),
        if (b$p_run2_better >= 0.95) lab[2] else if (b$p_run2_better <= 0.05) lab[1] else "undecided",
        sprintf("95%% CI of dOFV [%s; %s], %d subject resamples, no refit",
                f(b$ci[1], 1), f(b$ci[3], 1), b$n_boot))
    if (is.finite(b$p_run2_better_aic))
      add("Bootstrap P(dAIC < 0)", "", f(b$p_run2_better_aic, 3),
          if (b$p_run2_better_aic >= 0.95) lab[2] else if (b$p_run2_better_aic <= 0.05) lab[1] else "undecided")
    inf <- t$influence
    add("Subjects to remove to reverse dOFV", "", f(inf$n_reverse, 0), "-",
        sprintf("of %d; top 5 subjects carry %s %% of dOFV", inf$n_subjects, f(100 * inf$share_top5, 0)))
  }
  m <- x$metrics
  for (mt in c("RMSE(CWRES)", "Prop |CWRES|>2")) if (mt %in% m$metric) {
    v <- m[m$metric == mt, ]; add(mt, f(v$run1, 3), f(v$run2, 3), low(v$run1, v$run2))
  }
  p <- attr(m, "wilcox_absCWRES_p")
  if (!is.null(p)) add("Wilcoxon paired |CWRES|", "", "", if (p < 0.05) low(m$run1[m$metric == "Mean |CWRES|"], m$run2[m$metric == "Mean |CWRES|"]) else "none (p >= 0.05)",
                       sprintf("p = %s", format.pval(p, digits = 3)))
  ri <- x$run_info
  if (!is.null(ri)) {
    g <- function(field) vapply(ri, function(r) if (is.null(r) || is.null(r[[field]])) NA_real_ else as.numeric(r[[field]]), 0)
    rse <- vapply(ri, function(r) if (is.null(r)) NA_real_ else suppressWarnings(max(r$rse_max, na.rm = TRUE)), 0)
    add("Max RSE (%)", f(rse[1], 1), f(rse[2], 1), low(rse[1], rse[2]))
    cn <- g("cond"); add("Condition number", f(cn[1], 1), f(cn[2], 1), low(cn[1], cn[2]), "from the .cor file")
    sh <- g("max_eta_shrink"); add("Max ETA shrinkage (%)", f(sh[1], 1), f(sh[2], 1), low(sh[1], sh[2]))
    st <- vapply(ri, function(r) if (is.null(r)) NA_character_ else as.character(r$min_status), "")
    cs <- vapply(ri, function(r) if (is.null(r)) NA_character_ else as.character(r$cov_status), "")
    short <- function(v) sub("^MINIMIZATION ", "", v)
    add("Minimization", short(st[1]), short(st[2]), "-"); add("Covariance step", cs[1], cs[2], "-")
  }
  out <- do.call(rbind, rows)
  names(out)[2:3] <- lab
  class(out) <- c("nm_compare_summary", "data.frame")
  out
}

#' @export
print.nm_compare_summary <- function(x, ...) {
  cat("Decision criteria (favours: run preferred by the criterion)\n\n")
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  print(d[, 1:4], row.names = FALSE, right = FALSE)
  nt <- d[nzchar(d$note), , drop = FALSE]
  if (nrow(nt)) {
    cat("\nDetails:\n")
    cat(sprintf("  %-36s %s", nt$criterion, nt$note), sep = "\n")
  }
  cat("\nVuong/Clarke: individual contributions from the .phi files; corrections for the",
      "\nnumber of parameters (AIC: k, BIC: k*log(N subjects)/2). Bootstrap without refit:",
      "\nconditional on the estimates, optimistic for the larger model.\n")
  invisible(x)
}


# ---- Plots --------------------------------------------------------------------------------

.nm_cmp_plot_gain <- function(x, n_band) {
  d <- x$data_obs
  if (!"dAbsCWRES" %in% names(d)) return(NULL)
  lab <- x$meta$labels
  panels <- list()
  if ("TIME" %in% names(d)) panels$TIME <- list(x = d$TIME, log = FALSE)
  pc <- intersect(c("PRED_1", "PRED"), names(d))[1]
  if (!is.na(pc)) panels$PRED <- list(x = d[[pc]], log = all(d[[pc]] > 0, na.rm = TRUE))
  if (!length(panels)) return(NULL)
  pts <- do.call(rbind, lapply(names(panels), function(nm)
    data.frame(panel = paste("versus", nm, if (panels[[nm]]$log) "(log scale)" else ""),
               x = if (panels[[nm]]$log) log10(panels[[nm]]$x) else panels[[nm]]$x, y = d$dAbsCWRES)))
  bands <- do.call(rbind, lapply(names(panels), function(nm) {
    b <- .nm_cmp_band(panels[[nm]]$x, d$dAbsCWRES, d$ID, n_boot = n_band, log_x = panels[[nm]]$log)
    if (is.null(b)) return(NULL)
    b$panel <- paste("versus", nm, if (panels[[nm]]$log) "(log scale)" else "")
    if (panels[[nm]]$log) b$x <- log10(b$x)
    b
  }))
  p <- ggplot2::ggplot(pts, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey30") +
    ggplot2::geom_point(alpha = 0.25, size = 0.8)
  if (!is.null(bands)) p <- p +
    ggplot2::geom_ribbon(data = bands, ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi),
                         inherit.aes = FALSE, fill = "#2c7fb8", alpha = 0.3) +
    ggplot2::geom_line(data = bands, ggplot2::aes(x = .data$x, y = .data$fit), colour = "#08519c", linewidth = 0.9)
  p + ggplot2::facet_wrap(~panel, scales = "free_x") +
    ggplot2::labs(title = "Where does the fit change?",
                  subtitle = sprintf("|CWRES %s| - |CWRES %s|; below 0: %s better. Band: 95%% subject bootstrap",
                                     lab[2], lab[1], lab[2]),
                  x = NULL, y = "Delta |CWRES|") +
    ggplot2::theme_bw()
}

.nm_cmp_plot_overlay <- function(x, n_band) {
  d <- x$data_obs; lab <- x$meta$labels
  spec <- list(c("CWRES", "TIME"), c("CWRES", "PRED"), c("DV", "PRED"), c("DV", "IPRED"))
  pts <- list(); bands <- list()
  for (s in spec) {
    for (r in 1:2) {
      yv <- if (s[1] == "DV") (if ("DV_1" %in% names(d)) d[[paste0("DV_", r)]] else d$DV) else d[[paste0(s[1], "_", r)]]
      xv <- if (s[2] == "TIME") d$TIME else d[[paste0(s[2], "_", r)]]
      if (is.null(xv) || is.null(yv)) next
      pn <- paste(s[1], "vs", s[2])
      pts[[length(pts) + 1L]] <- data.frame(panel = pn, run = lab[r], x = xv, y = yv)
      b <- .nm_cmp_band(xv, yv, d$ID, n_boot = n_band, seed = r)
      if (!is.null(b)) { b$panel <- pn; b$run <- lab[r]; bands[[length(bands) + 1L]] <- b }
    }
  }
  if (!length(pts)) return(NULL)
  pts <- do.call(rbind, pts)
  pts$panel <- factor(pts$panel, levels = unique(pts$panel))
  ref <- data.frame(panel = factor(levels(pts$panel), levels = levels(pts$panel)))
  ref$slope <- ifelse(grepl("^DV", ref$panel), 1, 0)
  p <- ggplot2::ggplot(pts, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$run)) +
    ggplot2::geom_point(alpha = 0.15, size = 0.7) +
    ggplot2::geom_abline(data = ref, ggplot2::aes(slope = .data$slope, intercept = 0),
                         colour = "grey30", linetype = "dashed")
  if (length(bands)) {
    bands <- do.call(rbind, bands)
    bands$panel <- factor(bands$panel, levels = levels(pts$panel))
    p <- p + ggplot2::geom_ribbon(data = bands, ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$hi,
                                                            fill = .data$run),
                                  inherit.aes = FALSE, alpha = 0.25) +
      ggplot2::geom_line(data = bands, ggplot2::aes(x = .data$x, y = .data$fit, colour = .data$run), linewidth = 0.9)
  }
  p + ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::scale_colour_manual(values = stats::setNames(c("#d95f0e", "#2c7fb8"), lab)) +
    ggplot2::scale_fill_manual(values = stats::setNames(c("#d95f0e", "#2c7fb8"), lab)) +
    ggplot2::labs(title = "Superimposed diagnostics", subtitle = "Loess with 95% subject-bootstrap band",
                  x = NULL, y = NULL, colour = NULL, fill = NULL) +
    ggplot2::theme_bw()
}

.nm_cmp_plot_bootstrap <- function(x) {
  b <- x$tests$bootstrap
  if (is.null(b)) return(NULL)
  lab <- x$meta$labels
  ggplot2::ggplot(data.frame(d = b$dofv), ggplot2::aes(x = .data$d)) +
    ggplot2::geom_histogram(bins = 40, fill = "grey70", colour = "white") +
    ggplot2::geom_vline(xintercept = 0, colour = "red") +
    ggplot2::geom_vline(xintercept = b$ci[c(1, 3)], linetype = "dashed") +
    ggplot2::labs(title = "Subject bootstrap of dOFV (no refit)",
                  subtitle = sprintf("P(%s better) = %.3f; 95%% CI [%.1f; %.1f]",
                                     lab[2], b$p_run2_better, b$ci[1], b$ci[3]),
                  x = sprintf("dOFV (%s - %s)", lab[2], lab[1]), y = "Resamples") +
    ggplot2::theme_bw()
}

.nm_cmp_plot_cov <- function(x) {
  s <- x$data_subj; cov <- x$covariates
  if (is.null(cov) || !"diOFV" %in% names(s)) return(NULL)
  lab <- x$meta$labels
  m <- match(s$ID, cov$data$ID)
  dd <- do.call(rbind, lapply(cov$names, function(cv) {
    v <- cov$data[[cv]][m]
    num <- suppressWarnings(as.numeric(as.character(v)))
    data.frame(covariate = cv, value = as.character(v), num = num, diOFV = s$diOFV,
               cont = all(is.na(v) | !is.na(num)) && length(unique(v[!is.na(v)])) > 5L)
  }))
  plots <- lapply(split(dd, dd$covariate), function(z) {
    base <- ggplot2::ggplot(z) + ggplot2::geom_hline(yintercept = 0, colour = "grey30") + ggplot2::theme_bw() +
      ggplot2::labs(title = z$covariate[1], y = sprintf("diOFV (%s - %s)", lab[2], lab[1]), x = NULL)
    if (z$cont[1]) base + ggplot2::geom_point(ggplot2::aes(x = .data$num, y = .data$diOFV), alpha = 0.6) +
      ggplot2::geom_smooth(ggplot2::aes(x = .data$num, y = .data$diOFV), method = "loess", formula = y ~ x, se = FALSE)
    else base + ggplot2::geom_boxplot(ggplot2::aes(x = .data$value, y = .data$diOFV), outlier.shape = NA) +
      ggplot2::geom_jitter(ggplot2::aes(x = .data$value, y = .data$diOFV), width = 0.15, alpha = 0.6)
  })
  if (length(plots) == 1L) plots[[1]] else plots
}

.nm_cmp_plot_params <- function(x) {
  p <- x$parameters
  if (is.null(p)) return(NULL)
  lab <- x$meta$labels
  p <- p[p$type %in% c("THETA", "OMEGA", "SIGMA") & !(p$fixed_1 %in% TRUE & p$fixed_2 %in% TRUE), , drop = FALSE]
  if (!nrow(p)) return(NULL)
  long <- rbind(data.frame(par = paste(p$type, p$name), run = lab[1], est = p$value_1, se = p$se_1),
                data.frame(par = paste(p$type, p$name), run = lab[2], est = p$value_2, se = p$se_2))
  long <- long[is.finite(long$est), , drop = FALSE]
  long$run <- factor(long$run, levels = lab)
  ggplot2::ggplot(long, ggplot2::aes(x = .data$run, y = .data$est, colour = .data$run)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$est - 1.96 * .data$se, ymax = .data$est + 1.96 * .data$se),
                           width = 0.2, na.rm = TRUE) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~par, scales = "free_y") +
    ggplot2::scale_colour_manual(values = stats::setNames(c("#d95f0e", "#2c7fb8"), lab)) +
    ggplot2::labs(title = "Parameter estimates (95% Wald CI)", x = NULL, y = NULL, colour = NULL) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "none")
}

.nm_cmp_plot_etas <- function(x) {
  e <- x$etas
  if (is.null(e)) return(NULL)
  lab <- x$meta$labels
  ggplot2::ggplot(e, ggplot2::aes(x = .data$eta_1, y = .data$eta_2)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") + ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red") +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::facet_wrap(~parameter, scales = "free") +
    ggplot2::labs(title = "EBEs of the two runs", subtitle = "Red: identity. A flatter cloud: variability explained in run 2 (e.g. by a covariate)",
                  x = paste("ETA", lab[1]), y = paste("ETA", lab[2])) +
    ggplot2::theme_bw()
}
