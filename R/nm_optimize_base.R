# =============================================================================
# nm_optimize_base(): proposals to improve a base model without iterative runs
# =============================================================================
#
# Level 1 (outputs of the base run only, no NONMEM run):
#   - residual error : resmod-like models fitted in R on CWRES
#   - random effects : EBE correlations, distribution, shrinkage
#   - covariates     : ETA-covariate screening (Spearman / Kruskal-Wallis)
#
# Level 2 (one extra NONMEM evaluation run, MAXEVAL=0, written by
# nm_derivatives_ctl()): the table of the derivatives of Y with respect to
# the ETAs (G) and EPSs (H, and the cross derivatives EPS x ETA) gives a
# first-order linearization of the model around the EBEs (Khandelwal et al.,
# AAPS J 2011; Svensson & Karlsson, 2017). Every candidate extension
# (residual model, OMEGA structure, IOV, covariates on the ETAs) is then
# re-estimated in R on the linearized model, under the FOCE-I objective.
#
#   y_ij = IPRED0_ij + G_ij (eta_i - eta0_i) + eps_ij (H_ij + D_ij (eta_i - eta0_i))
#
# A covariate effect is a shift of the mean of an ETA: exact for a parameter
# written TV * EXP(ETA) with a power (log-covariate) or exponential effect,
# first-order otherwise.
# =============================================================================


# ---- Derivatives control stream ---------------------------------------------

# Block of an abbreviated-code record ($ERROR or $PRED): line range.
.nm_ob_record_range <- function(lines, rec) {
  s <- grep(paste0("^\\s*\\$", rec), lines, ignore.case = TRUE)[1]
  if (is.na(s)) return(NULL)
  all_rec <- grep("^\\s*\\$[A-Za-z]", lines)
  e <- all_rec[all_rec > s][1]
  e <- if (is.na(e)) length(lines) else e - 1L
  while (e > s && !nzchar(trimws(lines[e]))) e <- e - 1L
  c(s, e)
}

# Matrix (lower triangle) from the OMEGA(i,j) / SIGMA(i,j) entries of a .ext row.
.nm_ob_ext_matrix <- function(final_row, rec) {
  nm <- grep(paste0("^", rec, "\\("), names(final_row), value = TRUE)
  if (!length(nm)) return(matrix(numeric(0), 0, 0))
  ij <- do.call(rbind, lapply(regmatches(nm, gregexpr("[0-9]+", nm)), as.integer))
  n <- max(ij)
  m <- matrix(0, n, n)
  for (r in seq_len(nrow(ij))) {
    v <- as.numeric(final_row[[nm[r]]])
    m[ij[r, 1], ij[r, 2]] <- v; m[ij[r, 2], ij[r, 1]] <- v
  }
  m
}

# Groups of indices linked by non-zero off-diagonal elements (connected components).
.nm_ob_blocks <- function(m, tol = 0) {
  n <- nrow(m); grp <- seq_len(n)
  if (n > 1L) for (i in 2:n) for (j in seq_len(i - 1L)) {
    if (abs(m[i, j]) > tol) grp[grp == grp[i]] <- grp[j]
  }
  unname(split(seq_len(n), grp))
}

# $OMEGA / $SIGMA records, every value FIX, blocks written only where needed.
.nm_ob_write_random <- function(rec, m) {
  if (!length(m)) return(character(0))
  out <- character(0)
  for (b in .nm_ob_blocks(m)) {
    if (length(b) == 1L) {
      out <- c(out, sprintf("$%s %s FIX", rec, format(m[b, b], digits = 10)))
    } else {
      if (any(diff(b) != 1L))
        stop("Non-contiguous ", rec, " block (", paste(b, collapse = ", "),
             "): cannot be written back.", call. = FALSE)
      out <- c(out, sprintf("$%s BLOCK(%d) FIX", rec, length(b)))
      for (i in seq_along(b))
        out <- c(out, paste0("  ", paste(format(m[b[i], b[seq_len(i)]], digits = 10), collapse = " ")))
    }
  }
  out
}

# Estimation options of the base run, reused for the MAXEVAL=0 evaluation.
.nm_ob_est_options <- function(lines) {
  est <- toupper(paste(.nm_ctl_records_text(lines, "EST[A-Z]*"), collapse = " "))
  laplace <- grepl("LAPLAC", est)
  inter <- grepl("INTER", est) && !grepl("NOINTER", est)
  fo <- grepl("METHOD\\s*=\\s*(0|ZERO)\\b", est) || !nzchar(est)
  other <- grepl("METHOD\\s*=\\s*(SAEM|IMP|IMPMAP|BAYES|ITS|DIRECT|NUTS|CHAIN)", est)
  list(text = paste(c("$ESTIMATION", if (fo && !other) "METHOD=0" else "METHOD=1",
                      if (inter || other) "INTERACTION", if (laplace) "LAPLACIAN",
                      "MAXEVAL=0 POSTHOC NOABORT"), collapse = " "),
       inter = inter || other, fo = fo && !other, other = other)
}

#' Control stream of derivatives for the linearization of a base model
#'
#' Writes a copy of a base control stream that evaluates the model once at the
#' final estimates (`MAXEVAL=0`, all parameters `FIX`) and tabulates, for every
#' record, the individual prediction and the derivatives of `Y` with respect to
#' the ETAs (`G_ETAk`), the EPSs (`H_EPSj`) and the cross derivatives
#' (`D_EPSj_ETAk`, interaction). The derivatives are read through verbatim code
#' (`G(k,1)`, `HH(j,1)`, `HH(j,k+1)` in `$ERROR`; `G` and `H` in `$PRED`)
#' placed at the end of the block (`"LAST`), as in the PsN `linearize` tool.
#' The table is the input of level 2 of [nm_optimize_base()].
#'
#' @param run Run directory of the base model (e.g. `"run1.reg"`) or path
#'   prefix (`"dir/run1"`) whose `.ctl` and `.ext` files are used.
#' @param out_dir Directory of the new control stream. Default: the analysis
#'   directory (parent of a `.reg`/`.nm7` run directory), so that the `$DATA`
#'   path keeps the convention of [nm_run()].
#' @param name Name of the new run. Default `"<run>_deriv"`.
#' @param extra Additional data items to tabulate (e.g. an occasion column
#'   `"OCC"` for inter-occasion variability). Items absent from `$INPUT` are
#'   dropped with a warning.
#' @return Invisibly, the path of the control stream written. The table file
#'   is `<name>.dtab`.
#' @details Not testable without NONMEM: the verbatim code follows the
#'   NONMEM help (ERROR: `G(i,1)` is the partial derivative of F with respect
#'   to ETA(i); `HH(j,1)` with respect to EPS(j); `HH(j,i+1)` the second
#'   derivative with respect to EPS(j) and ETA(i)). Models with `F_FLAG`
#'   (M3 likelihood) or simulation-only records are not supported.
#' @seealso [nm_optimize_base()]
#' @export
nm_derivatives_ctl <- function(run, out_dir = NULL, name = NULL, extra = NULL) {
  f <- .nm_w_run_files(run)
  if (is.na(f$ctl)) stop("No .ctl found for run: ", run, call. = FALSE)
  if (is.na(f$ext)) stop("No .ext found for run: ", run, call. = FALSE)
  lines <- .nm_read_lines(f$ctl)
  if (any(grepl("\\bF_FLAG\\b", .nm_ctl_strip_comments(lines), ignore.case = TRUE)))
    stop("F_FLAG (M3 / categorical likelihood) is not supported by the linearization.",
         call. = FALSE)
  fin <- .nm_vpc_read_ext_final(f$ext)
  th <- as.numeric(fin[1, grep("^THETA", names(fin))])
  om <- .nm_ob_ext_matrix(fin, "OMEGA"); sg <- .nm_ob_ext_matrix(fin, "SIGMA")
  n_eta <- nrow(om); n_eps <- nrow(sg)
  if (!n_eta) stop("No OMEGA in the .ext file: a model without ETA cannot be linearized.",
                   call. = FALSE)
  est <- .nm_ob_est_options(lines)
  if (est$other)
    warning("The base run does not use FO/FOCE: the evaluation uses METHOD=1 INTERACTION.",
            call. = FALSE)

  rec <- if (!is.null(.nm_ob_record_range(lines, "ERROR"))) "ERROR" else "PRED"
  rg <- .nm_ob_record_range(lines, rec)
  if (is.null(rg)) stop("Neither $ERROR nor $PRED found in ", f$ctl, call. = FALSE)
  hh <- if (rec == "ERROR") "HH" else "H"
  g_nm <- sprintf("G_ETA%d", seq_len(n_eta))
  h_nm <- sprintf("H_EPS%d", seq_len(n_eps))
  d_nm <- as.vector(outer(seq_len(n_eps), seq_len(n_eta),
                          function(j, k) sprintf("D_EPS%d_ETA%d", j, k)))
  d_j <- rep(seq_len(n_eps), n_eta); d_k <- rep(seq_len(n_eta), each = n_eps)
  code <- c("; ---- nmpilot: derivatives of Y for the linearization ----",
            paste0("  ", c(g_nm, h_nm, d_nm), " = 0"),
            "\"LAST",
            sprintf("\"  %s = G(%d,1)", g_nm, seq_len(n_eta)),
            sprintf("\"  %s = %s(%d,1)", h_nm, hh, seq_len(n_eps)),
            sprintf("\"  %s = %s(%d,%d)", d_nm, hh, d_j, d_k + 1L))
  lines <- append(lines, code, after = rg[2])

  # Data items for the table
  inp <- .nm_ctl_input_items(lines)
  if (!"ID" %in% inp) stop("No ID item in $INPUT.", call. = FALSE)
  items <- intersect(c("ID", "TIME", "EVID", "MDV"), inp)
  if (!is.null(extra)) {
    miss <- setdiff(extra, inp)
    if (length(miss)) warning("Item(s) not in $INPUT, dropped: ", paste(miss, collapse = ", "),
                              call. = FALSE)
    items <- unique(c(items, intersect(extra, inp)))
  }
  if (!"MDV" %in% items) items <- c(items, "MDV")
  pred_item <- if (est$fo) "PRED" else if (est$inter) "CIPREDI" else "CIPRED"

  stem <- f$stem
  if (is.null(name)) name <- paste0(stem, "_deriv")
  if (is.null(out_dir)) {
    out_dir <- if (grepl("\\.(reg|nm7)$", f$dir, ignore.case = TRUE)) dirname(f$dir) else f$dir
  }
  body <- lines[!grepl("^\\s*\\$PROB", lines, ignore.case = TRUE)]
  body <- .nm_vpc_remove_records(body)
  for (r in c("SIML", "SIMU", "MSFI", "NONP", "SCAT", "PRIOR"))
    body <- .nm_strip_record(body, r)
  out <- c(sprintf("$PROBLEM %s - derivatives of %s for nm_optimize_base (MAXEVAL=0)", name, stem),
           body, "",
           "; ---- final estimates of the base run, fixed ----",
           paste0("$THETA ", vapply(th, format, "", digits = 10), " FIX"),
           .nm_ob_write_random("OMEGA", om),
           .nm_ob_write_random("SIGMA", sg),
           "", est$text, "",
           paste("$TABLE", paste(c(items, "DV", pred_item, sprintf("ETA%d", seq_len(n_eta))),
                                 collapse = " ")),
           paste0("  ", paste(c(g_nm, h_nm, d_nm), collapse = " ")),
           paste0("  ONEHEADER NOPRINT NOAPPEND FORMAT=s1PE16.8 FILE=", name, ".dtab"), "")
  path <- file.path(out_dir, paste0(name, ".ctl"))
  writeLines(out, path)
  invisible(path)
}


# ---- Inputs -------------------------------------------------------------------

# First NONMEM table of a run directory holding all the requested columns.
.nm_ob_find_table <- function(dir, need) {
  cand <- list.files(dir, full.names = TRUE)
  cand <- cand[!dir.exists(cand)]
  cand <- cand[!grepl("\\.(ext|phi|lst|ctl|cov|cor|coi|shk|grd|xml|csv|mod|dtab|phm|smr|cpu)$",
                      cand, ignore.case = TRUE)]
  for (f in cand) {
    l1 <- tryCatch(readLines(f, n = 2L, warn = FALSE), error = function(e) character(0))
    if (length(l1) < 2L || !grepl("^\\s*TABLE NO", l1[1])) next
    if (all(need %in% strsplit(trimws(l1[2]), "\\s+")[[1]])) return(f)
  }
  NA_character_
}

.nm_ob_final_ofv <- function(ext) {
  fin <- .nm_vpc_read_ext_final(ext)
  as.numeric(fin$OBJ[1])
}


# ---- Level 1: residual error on CWRES (resmod-like) --------------------------

# -2 log-likelihood of one subject for a Gaussian with covariance C.
.nm_ob_mvn <- function(r, C) {
  L <- tryCatch(chol(C), error = function(e) NULL)
  if (is.null(L)) return(1e10)
  z <- backsolve(L, r, transpose = TRUE)
  2 * sum(log(diag(L))) + sum(z^2)
}

.nm_ob_l1_residual <- function(tab) {
  d <- tab[is.finite(tab$CWRES) & is.finite(tab$IPRED), , drop = FALSE]
  if (nrow(d) < 20L) return(NULL)
  d <- d[order(d$ID, d$TIME), , drop = FALSE]
  ids <- split(seq_len(nrow(d)), d$ID)
  ip <- d$IPRED
  gm <- exp(mean(log(ip[ip > 0])))
  lip <- log(pmax(ip, gm * 1e-3) / gm)
  tcut <- stats::median(d$TIME)
  late <- d$TIME > tcut
  gh <- .nm_ob_gauss_hermite(15L)

  nll <- function(p, type) {
    mu <- p[1]; w2 <- exp(2 * p[2]); s2 <- exp(2 * p[3])
    tot <- 0
    for (ix in ids) {
      r <- d$CWRES[ix] - mu; n <- length(ix)
      J <- matrix(w2, n, n)
      if (type == "base") {
        tot <- tot + .nm_ob_mvn(r, J + diag(s2, n))
      } else if (type == "power") {
        tot <- tot + .nm_ob_mvn(r, J + diag(s2 * exp(2 * p[4] * lip[ix]), n))
      } else if (type == "time_varying") {
        tot <- tot + .nm_ob_mvn(r, J + diag(s2 * ifelse(late[ix], exp(2 * p[4]), 1), n))
      } else if (type == "autocorrelation") {
        tt <- d$TIME[ix]
        hl <- exp(p[4])
        tot <- tot + .nm_ob_mvn(r, J + s2 * exp(-abs(outer(tt, tt, "-")) * log(2) / hl))
      } else if (type == "iiv_on_eps") {
        we <- exp(p[4])
        li <- vapply(seq_along(gh$x), function(q) {
          -0.5 * .nm_ob_mvn(r, J + diag(s2 * exp(sqrt(2) * we * gh$x[q]), n))
        }, numeric(1))
        mx <- max(li)
        tot <- tot - 2 * (mx + log(sum(gh$w * exp(li - mx)) / sqrt(pi)))
      }
    }
    tot
  }
  v0 <- stats::var(d$CWRES)
  start <- c(mean(d$CWRES), log(sqrt(0.1 * v0)), log(sqrt(0.9 * v0)))
  fit <- function(type, extra) {
    o <- tryCatch(stats::optim(c(start, extra), nll, type = type, method = "BFGS",
                               control = list(maxit = 300)), error = function(e) NULL)
    if (is.null(o)) return(list(ofv = NA_real_, par = rep(NA_real_, 4)))
    list(ofv = o$value, par = o$par)
  }
  b <- fit("base", numeric(0))
  types <- c("power", "time_varying", "autocorrelation", "iiv_on_eps")
  ex0 <- list(power = 0, time_varying = 0,
              autocorrelation = log(max(diff(range(d$TIME)) / 10, 1e-3)),
              iiv_on_eps = log(0.3))
  res <- lapply(types, function(ty) fit(ty, ex0[[ty]]))
  est <- vapply(seq_along(types), function(i) {
    p <- res[[i]]$par[4]
    switch(types[i], power = p, time_varying = exp(p), autocorrelation = exp(p),
           iiv_on_eps = exp(p))
  }, numeric(1))
  lab <- c(power = "delta (variance ~ IPRED^(2*delta))",
           time_varying = sprintf("SD ratio TIME > %g vs <= %g", signif(tcut, 3), signif(tcut, 3)),
           autocorrelation = "half-life of the correlation (time units)",
           iiv_on_eps = "SD of the ETA on EPS (log scale)")
  out <- data.frame(model = types, df = 1L,
                    dOFV = round(vapply(res, `[[`, numeric(1), "ofv") - b$ofv, 2),
                    estimate = signif(est, 3), meaning = unname(lab[types]),
                    stringsAsFactors = FALSE)
  out$suggestion <- ifelse(
    out$dOFV > -10, "-",
    c(power = "change the residual model: proportional/combined/power",
      time_varying = "time-varying residual error (or structural misfit at late/early times)",
      autocorrelation = "autocorrelated residuals: structural misfit or AR(1) error",
      iiv_on_eps = "individual-specific residual variance: ETA on EPS")[out$model])
  attr(out, "base") <- list(ofv = b$ofv, mean = b$par[1], sd_iiv = exp(b$par[2]),
                            sd = exp(b$par[3]), n_obs = nrow(d), tcut = tcut)
  out
}

.nm_ob_gauss_hermite <- function(n) {
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- sqrt(i / 2); J[cbind(i + 1L, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}


# ---- Level 1: random effects ---------------------------------------------------

.nm_ob_l1_random <- function(eta, omega) {
  K <- ncol(eta)
  nm <- colnames(eta)
  sd_w <- sqrt(pmax(diag(omega)[seq_len(K)], 0))
  shr <- ifelse(sd_w > 0, 1 - apply(eta, 2, stats::sd) / sd_w, NA_real_)
  skew <- apply(eta, 2, function(v) mean((v - mean(v))^3) / stats::sd(v)^3)
  sw <- apply(eta, 2, function(v) if (length(unique(v)) > 2 && length(v) <= 5000)
    stats::shapiro.test(v)$p.value else NA_real_)
  eta_tab <- data.frame(eta = nm, omega = signif(diag(omega)[seq_len(K)], 3),
                        shrinkage_pct = round(100 * shr, 1), skewness = round(skew, 2),
                        p_shapiro = signif(sw, 2), stringsAsFactors = FALSE)
  eta_tab$suggestion <- ifelse(
    !is.na(shr) & shr > 0.4, "high shrinkage: ETA poorly informed (consider removing/fixing it)",
    ifelse(eta_tab$omega < 0.005, "very small variance: consider removing the ETA",
           ifelse(abs(skew) > 0.5 & !is.na(sw) & sw < 0.01,
                  "skewed ETA: consider a Box-Cox / logit transformation", "-")))
  pairs <- NULL
  if (K > 1L) {
    pr <- utils::combn(K, 2)
    pairs <- do.call(rbind, lapply(seq_len(ncol(pr)), function(i) {
      a <- pr[1, i]; b <- pr[2, i]
      ct <- stats::cor.test(eta[, a], eta[, b])
      data.frame(eta_1 = nm[a], eta_2 = nm[b], r = round(unname(ct$estimate), 2),
                 p = signif(ct$p.value, 2), max_shrinkage_pct = round(100 * max(shr[a], shr[b]), 1),
                 omega_cov = signif(omega[a, b], 3), stringsAsFactors = FALSE)
    }))
    pairs$suggestion <- ifelse(
      pairs$omega_cov != 0, "already in an OMEGA block",
      ifelse(abs(pairs$r) >= 0.3 & pairs$p < 0.05,
             ifelse(pairs$max_shrinkage_pct < 30, "consider an OMEGA BLOCK",
                    "correlation possibly induced by shrinkage"), "-"))
  }
  list(eta = eta_tab, pairs = pairs)
}


# ---- Level 1: covariates on EBEs ----------------------------------------------

# Covariate types and design columns (one row per subject).
.nm_ob_cov_info <- function(cov, covariates, max_levels) {
  lapply(stats::setNames(covariates, covariates), function(cv) {
    v <- cov[[cv]]
    num <- suppressWarnings(as.numeric(as.character(v)))
    is_num <- all(is.na(v) | !is.na(num))
    nlev <- length(unique(v[!is.na(v)]))
    list(type = if (is_num && nlev > max_levels) "continuous" else "categorical",
         values = if (is_num) num else as.character(v), n_levels = nlev)
  })
}

.nm_ob_l1_covariates <- function(eta, ids, cov, info, shrink) {
  m <- match(ids, cov$ID)
  rows <- list()
  for (cv in names(info)) {
    v <- info[[cv]]$values[m]
    for (k in seq_len(ncol(eta))) {
      e <- eta[, k]; ok <- !is.na(v)
      if (sum(ok) < 6L || length(unique(v[ok])) < 2L) next
      if (info[[cv]]$type == "continuous") {
        ct <- suppressWarnings(stats::cor.test(e[ok], v[ok], method = "spearman", exact = FALSE))
        stat <- unname(ct$estimate); p <- ct$p.value; what <- "Spearman rho"
      } else {
        kt <- stats::kruskal.test(e[ok], factor(v[ok]))
        g <- tapply(e[ok], factor(v[ok]), stats::median)
        stat <- max(g) - min(g); p <- kt$p.value; what <- "range of medians (Kruskal-Wallis)"
      }
      rows[[length(rows) + 1L]] <- data.frame(
        eta = colnames(eta)[k], covariate = cv, type = info[[cv]]$type, statistic = round(stat, 3),
        statistic_is = what, p = p, shrinkage_pct = round(100 * shrink[k], 1),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$p_adj <- pmin(1, out$p * nrow(out))
  out$suggestion <- ifelse(out$p_adj < 0.05,
                           ifelse(out$shrinkage_pct < 30, "candidate covariate",
                                  "candidate (EBE shrinkage > 30 %: check at level 2)"), "-")
  out$p <- signif(out$p, 2); out$p_adj <- signif(out$p_adj, 2)
  out[order(out$p), , drop = FALSE]
}


# ---- Level 2: linearized model ---------------------------------------------------

# Data of the linearization, one element per subject.
.nm_ob_prepare <- function(dt, omega, sigma, X, occ_col) {
  K <- nrow(omega); J <- nrow(sigma)
  g_nm <- sprintf("G_ETA%d", seq_len(K)); h_nm <- sprintf("H_EPS%d", seq_len(J))
  pred_nm <- intersect(c("CIPREDI", "CIPRED", "PRED", "IPRED"), names(dt))[1]
  need <- c("ID", "DV", g_nm, h_nm, sprintf("ETA%d", seq_len(K)))
  miss <- setdiff(need, names(dt))
  if (length(miss) || is.na(pred_nm))
    stop("Columns missing from the derivatives table: ",
         paste(c(miss, if (is.na(pred_nm)) "CIPREDI"), collapse = ", "), call. = FALSE)
  if ("MDV" %in% names(dt)) dt <- dt[dt$MDV == 0, , drop = FALSE]
  if ("EVID" %in% names(dt)) dt <- dt[dt$EVID == 0, , drop = FALSE]
  if (!"TIME" %in% names(dt)) dt$TIME <- stats::ave(dt$ID, dt$ID, FUN = seq_along)
  ids <- unique(dt$ID)
  subj <- lapply(ids, function(i) {
    s <- dt[dt$ID == i, , drop = FALSE]
    G <- as.matrix(s[, g_nm, drop = FALSE]); H <- as.matrix(s[, h_nm, drop = FALSE])
    eta0 <- as.numeric(s[1, sprintf("ETA%d", seq_len(K))])
    D <- lapply(seq_len(J), function(j) {
      nmj <- sprintf("D_EPS%d_ETA%d", j, seq_len(K))
      if (all(nmj %in% names(s))) as.matrix(s[, nmj, drop = FALSE]) else matrix(0, nrow(s), K)
    })
    ip <- s[[pred_nm]]
    occ <- if (!is.null(occ_col)) as.integer(factor(s[[occ_col]])) else rep(1L, nrow(s))
    list(id = i, y = s$DV, ipred0 = ip, eta0 = eta0, G = G, H = H, D = D,
         f0 = ip - drop(G %*% eta0), r0 = s$DV - ip + drop(G %*% eta0),
         time = s$TIME, occ = occ, n = nrow(s),
         x = if (is.null(X)) numeric(0) else X[match(i, X$ID), -1, drop = TRUE])
  })
  list(subj = subj, K = K, J = J, omega0 = omega, sigma0 = sigma, n_obs = nrow(dt),
       cache = new.env(),
       tcut = stats::median(dt$TIME), ids = ids)
}

# Covariate design: one numeric column per effect, one row per subject.
.nm_ob_cov_design <- function(cov, ids, info, form) {
  m <- match(ids, cov$ID)
  cols <- list(ID = ids); terms <- list()
  for (cv in names(info)) {
    v <- info[[cv]]$values[m]
    if (info[[cv]]$type == "continuous") {
      med <- stats::median(v, na.rm = TRUE)
      if (anyNA(v)) { warning("Missing ", cv, " imputed by the median.", call. = FALSE); v[is.na(v)] <- med }
      if (form == "power" && all(v > 0)) {
        cols[[cv]] <- log(v / med); lab <- sprintf("(%s/%g)^theta", cv, signif(med, 4))
      } else {
        cols[[cv]] <- v - med; lab <- sprintf("exp(theta*(%s-%g))", cv, signif(med, 4))
      }
      terms[[cv]] <- list(cols = cv, label = lab, df = 1L)
    } else {
      tb <- sort(table(v), decreasing = TRUE)
      if (anyNA(v)) { warning("Missing ", cv, " set to its most frequent level.", call. = FALSE); v[is.na(v)] <- names(tb)[1] }
      lev <- sort(unique(as.character(v)))
      ref <- names(tb)[1]
      cn <- character(0)
      for (l in setdiff(lev, ref)) {
        nmc <- paste0(cv, "_", l); cols[[nmc]] <- as.numeric(as.character(v) == l); cn <- c(cn, nmc)
      }
      terms[[cv]] <- list(cols = cn, label = sprintf("theta^(%s != %s)", cv, ref), df = length(cn))
    }
  }
  list(X = as.data.frame(cols, check.names = FALSE), terms = terms)
}

# Number of parameters of each part of a candidate.
.nm_ob_npar <- function(spec, prep) {
  n_om <- sum(vapply(spec$blocks, function(b) length(b) * (length(b) + 1L) / 2L, numeric(1)))
  n_res <- switch(spec$res, base = 1L, add = 1L, prop = 1L, comb = 2L, power = 2L, timevar = 2L)
  c(omega = n_om, mu = prep$K, beta = length(spec$beta_cols), res = n_res,
    iov = length(spec$iov))
}

.nm_ob_unpack <- function(p, spec, prep) {
  np <- .nm_ob_npar(spec, prep); K <- prep$K
  i <- 0L; take <- function(n) { v <- p[i + seq_len(n)]; i <<- i + n; v }
  Om <- matrix(0, K, K)
  for (b in spec$blocks) {
    nb <- length(b); L <- matrix(0, nb, nb)
    L[lower.tri(L, diag = TRUE)] <- take(nb * (nb + 1L) / 2L)
    diag(L) <- exp(diag(L))
    Om[b, b] <- L %*% t(L)
  }
  mu <- take(K)
  beta <- matrix(0, K, length(spec$beta_cols))
  bt <- take(np[["beta"]])
  if (length(bt)) beta[cbind(spec$beta_eta, seq_along(bt))] <- bt
  res <- take(np[["res"]])
  iov <- exp(2 * take(np[["iov"]]))
  list(Om = Om, mu = mu, beta = beta, res = res, iov = iov)
}

# Residual variance of each observation and its derivatives with respect to
# the ETAs (n x K), given the linearized prediction and the ETA shifts.
.nm_ob_resvar <- function(type, res, s, fhat, deta, sigma0, tcut) {
  n <- s$n; K <- ncol(s$G)
  if (type %in% c("base", "timevar")) {
    # dY/dEPS around the EBEs, h_j(eta) = H_j * exp(D_j (eta - eta0) / H_j):
    # same value and first derivative as the linear expansion at eta0, but it
    # cannot cross zero (a linear h reaching 0 makes the variance vanish and
    # the individual objective unbounded).
    hv <- s$H; dh <- vector("list", length(s$D))
    for (j in seq_along(s$D)) {
      lin <- rowSums(s$D[[j]] * deta)
      Hj <- s$H[, j]
      pos <- Hj != 0
      if (all(pos)) {
        ex <- exp(pmin(lin / Hj, 30))
        hv[, j] <- Hj * ex
      } else {
        ex <- rep(1, n)
        ex[pos] <- exp(pmin(lin[pos] / Hj[pos], 30))
        hv[, j] <- Hj * ex + (!pos) * lin
      }
      dh[[j]] <- s$D[[j]] * ex
    }
    hs <- hv %*% sigma0
    b <- rowSums(hs * hv)
    db <- 2 * hs[, 1] * dh[[1]]
    if (length(dh) > 1L) for (j in 2:length(dh)) db <- db + 2 * hs[, j] * dh[[j]]
    sc <- if (type == "base") rep(exp(2 * res[1]), n)
          else ifelse(s$time > tcut, exp(2 * res[2]), exp(2 * res[1]))
    R <- sc * b; dR <- sc * db
  } else if (type == "add") {
    R <- rep(exp(2 * res[1]), n); dR <- matrix(0, n, K)
  } else if (type == "prop") {
    R <- exp(2 * res[1]) * fhat^2; dR <- 2 * exp(2 * res[1]) * fhat * s$G
  } else if (type == "comb") {
    R <- exp(2 * res[1]) + exp(2 * res[2]) * fhat^2; dR <- 2 * exp(2 * res[2]) * fhat * s$G
  } else if (type == "power") {
    af <- pmax(abs(fhat), 1e-12)
    R <- exp(2 * res[1]) * af^(2 * res[2])
    dR <- exp(2 * res[1]) * 2 * res[2] * af^(2 * res[2] - 1) * sign(fhat) * s$G
  }
  low <- R < 1e-12
  R[low] <- 1e-12; dR[low, ] <- 0
  list(R = R, dR = dR)
}

# FOCE-I objective of the linearized model (without the n*log(2*pi) constant),
# as in NONMEM with INTERACTION:
#   OFV_i = sum(log R + e^2/R) + u' Psi^-1 u + log|Psi|
#           + log|Psi^-1 + Z' R^-1 Z + 1/2 dR' R^-2 dR|
# at the mode u of the individual objective (Fisher scoring), where u holds
# the ETAs with a variance (IIV) and the occasion ETAs (IOV).
.nm_ob_ofv <- function(p, spec, prep, max_iter = 10L, details = FALSE) {
  u <- .nm_ob_unpack(p, spec, prep)
  if (any(!is.finite(u$Om))) return(1e10)
  K <- prep$K
  kept <- sort(unlist(spec$blocks))
  Ok <- u$Om[kept, kept, drop = FALSE]
  Oinv <- if (length(kept)) tryCatch(solve(Ok), error = function(e) NULL) else matrix(0, 0, 0)
  if (is.null(Oinv)) return(1e10)
  ld_O <- if (length(kept)) as.numeric(determinant(Ok, logarithm = TRUE)$modulus) else 0
  tot <- 0
  eb <- if (details) matrix(NA_real_, length(prep$subj), K) else NULL
  for (si in seq_along(prep$subj)) {
    s <- prep$subj[[si]]
    mu_i <- u$mu + if (length(s$x)) drop(u$beta %*% unlist(s$x)[spec$beta_cols_idx]) else 0
    # random-effect components: (eta index, rows)
    comp_k <- kept; comp_rows <- rep(list(rep(TRUE, s$n)), length(kept))
    psi_inv <- Oinv; ld_psi <- ld_O
    if (length(spec$iov)) {
      occs <- sort(unique(s$occ))
      for (o in occs) for (ii in seq_along(spec$iov)) {
        comp_k <- c(comp_k, spec$iov[ii]); comp_rows <- c(comp_rows, list(s$occ == o))
      }
      d_iov <- rep(u$iov, times = length(occs))
      q0 <- nrow(psi_inv); q <- q0 + length(d_iov)
      P <- matrix(0, q, q)
      if (q0) P[seq_len(q0), seq_len(q0)] <- psi_inv
      P[cbind(q0 + seq_along(d_iov), q0 + seq_along(d_iov))] <- 1 / d_iov
      psi_inv <- P; ld_psi <- ld_psi + sum(log(d_iov))
    }
    q <- length(comp_k)
    M <- if (q) vapply(seq_len(q), function(l) as.numeric(comp_rows[[l]]), numeric(s$n)) else matrix(0, s$n, 0)
    if (!is.matrix(M)) M <- matrix(M, s$n)
    Z <- if (q) s$G[, comp_k, drop = FALSE] * M else matrix(0, s$n, 0)
    mu_mat <- matrix(mu_i, s$n, K, byrow = TRUE)
    eta0_mat <- matrix(s$eta0, s$n, K, byrow = TRUE)
    state <- function(uu) {
      eo <- mu_mat
      if (q) for (l in seq_len(q)) eo[, comp_k[l]] <- eo[, comp_k[l]] + uu[l] * M[, l]
      fh <- s$f0 + rowSums(s$G * eo)
      rv <- .nm_ob_resvar(spec$res, u$res, s, fh, eo - eta0_mat, prep$sigma0, prep$tcut)
      e <- s$y - fh
      dRu <- if (q) rv$dR[, comp_k, drop = FALSE] * M else matrix(0, s$n, 0)
      list(e = e, R = rv$R, dRu = dRu,
           obj = sum(log(rv$R) + e^2 / rv$R) + if (q) sum(uu * (psi_inv %*% uu)) else 0)
    }
    key <- paste(si, q)
    u_def <- if (q) c(s$eta0[kept] - mu_i[kept], rep(0, q - length(kept))) else numeric(0)
    mode_find <- function(uu, n_it) {
      st <- state(uu); conv <- !q
      if (q) for (it in seq_len(n_it)) {
        g <- drop(crossprod(st$dRu, 1 / st$R - st$e^2 / st$R^2) - 2 * crossprod(Z, st$e / st$R)) +
          2 * drop(psi_inv %*% uu)
        # observed Hessian without the second derivatives of R; Fisher
        # scoring when it is not positive definite
        w1 <- 2 * st$e^2 / st$R^3 - 1 / st$R^2
        ZdR <- crossprod(Z, st$dRu * (2 * st$e / st$R^2))
        H <- 2 * crossprod(Z, Z / st$R) + crossprod(st$dRu, st$dRu * w1) + ZdR + t(ZdR) + 2 * psi_inv
        cH <- tryCatch(chol(H), error = function(e) NULL)
        if (is.null(cH)) {
          H <- 2 * crossprod(Z, Z / st$R) + crossprod(st$dRu, st$dRu / st$R^2) + 2 * psi_inv
          cH <- tryCatch(chol(H), error = function(e) NULL)
          if (is.null(cH)) break
        }
        step <- backsolve(cH, backsolve(cH, g, transpose = TRUE))
        if (any(!is.finite(step))) break
        if (max(abs(step)) < 1e-6) { conv <- TRUE; break }
        lam <- 1; sc <- NULL
        repeat {
          cand <- uu - lam * step
          tr <- state(cand)
          if (is.finite(tr$obj) && tr$obj <= st$obj + 1e-10) { sc <- tr; break }
          lam <- lam / 2
          if (lam < 1e-6) break
        }
        if (is.null(sc)) { conv <- max(abs(step)) < 1e-4; break }
        dobj <- st$obj - sc$obj
        uu <- cand; st <- sc
        if (lam == 1 && dobj < 1e-10) { conv <- TRUE; break }
      }
      list(uu = uu, st = st, conv = conv)
    }
    warm <- if (!is.null(prep$cache)) prep$cache[[key]] else NULL
    mf <- if (!is.null(warm) && length(warm) == q) mode_find(warm, 25L) else NULL
    if (is.null(mf) || !mf$conv) {
      mf2 <- mode_find(u_def, 80L)
      if (is.null(mf) || mf2$st$obj < mf$st$obj) mf <- mf2
    }
    uu <- mf$uu; st <- mf$st
    if (!is.null(prep$cache) && q && mf$conv) assign(key, uu, envir = prep$cache)
    ofv_i <- st$obj + ld_psi
    if (q) {
      A <- psi_inv + crossprod(Z, Z / st$R) + 0.5 * crossprod(st$dRu, st$dRu / st$R^2)
      ld <- determinant(A, logarithm = TRUE)
      if (ld$sign <= 0) return(1e10)
      ofv_i <- ofv_i + as.numeric(ld$modulus)
    }
    tot <- tot + ofv_i
    if (details) { e_i <- mu_i; if (length(kept)) e_i[kept] <- e_i[kept] + uu[seq_along(kept)]; eb[si, ] <- e_i }
  }
  if (details) return(list(ofv = tot, eta = eb))
  if (!is.finite(tot)) 1e10 else tot
}

# Residual variance of each observation.
.nm_ob_resvar_fast <- function(type, res, s, fhat, deta, sigma0) {
  v <- switch(type,
    base = , timevar = {
      hv <- s$H
      if (any(vapply(s$D, function(d) any(d != 0), logical(1))))
        for (j in seq_along(s$D)) hv[, j] <- hv[, j] + rowSums(s$D[[j]] * deta)
      b <- rowSums((hv %*% sigma0) * hv)
      if (type == "base") exp(2 * res[1]) * b
      else b * ifelse(s$time > attr(s, "tcut"), exp(2 * res[2]), exp(2 * res[1]))
    },
    add = rep(exp(2 * res[1]), s$n),
    prop = exp(2 * res[1]) * fhat^2,
    comb = exp(2 * res[1]) + exp(2 * res[2]) * fhat^2,
    power = exp(2 * res[1]) * abs(fhat)^(2 * res[2]))
  pmax(v, 1e-12)
}

# Fast approximation of the objective (objective = "fast"): modes given by the
# linear predictor with the residual variance held fixed, no residual-variance
# term in the determinant. Close dOFV, OFV a few units away from NONMEM's
# when the residual error depends on the ETAs.
.nm_ob_ofv_fast <- function(p, spec, prep, n_iter = 2L) {
  u <- .nm_ob_unpack(p, spec, prep)
  if (any(!is.finite(u$Om))) return(1e10)
  tot <- 0
  for (si in seq_along(prep$subj)) {
    s <- prep$subj[[si]]
    attr(s, "tcut") <- prep$tcut
    mu_i <- u$mu + if (length(s$x)) drop(u$beta %*% unlist(s$x)[spec$beta_cols_idx]) else 0
    r <- s$r0 - drop(s$G %*% mu_i)
    GOG <- s$G %*% u$Om %*% t(s$G)
    Giov <- NULL
    if (length(spec$iov)) {
      Giov <- lapply(sort(unique(s$occ)), function(o) {
        Go <- s$G[, spec$iov, drop = FALSE] * (s$occ == o); Go
      })
      for (Go in Giov) GOG <- GOG + Go %*% (u$iov * t(Go))
    }
    deta <- matrix(0, s$n, prep$K)
    fhat <- s$ipred0
    for (it in seq_len(n_iter)) {
      R <- .nm_ob_resvar_fast(spec$res, u$res, s, fhat, deta, prep$sigma0)
      V <- GOG + diag(R, s$n)
      L <- tryCatch(chol(V), error = function(e) NULL)
      if (is.null(L)) return(1e10)
      w <- backsolve(L, backsolve(L, r, transpose = TRUE))
      b <- drop(u$Om %*% crossprod(s$G, w))
      eta_obs <- matrix(mu_i + b, s$n, prep$K, byrow = TRUE)
      if (length(spec$iov)) for (Go in Giov) {
        kap <- u$iov * drop(crossprod(Go, w))
        rows <- rowSums(Go != 0) > 0
        eta_obs[rows, spec$iov] <- eta_obs[rows, spec$iov] + matrix(kap, sum(rows), length(kap), byrow = TRUE)
      }
      deta <- eta_obs - matrix(s$eta0, s$n, prep$K, byrow = TRUE)
      fhat <- s$f0 + rowSums(s$G * eta_obs)
    }
    z <- backsolve(L, r, transpose = TRUE)
    tot <- tot + 2 * sum(log(diag(L))) + sum(z^2)
  }
  if (!is.finite(tot)) 1e10 else tot
}


# Starting values of a candidate from the base estimates.
.nm_ob_start <- function(spec, prep, from = NULL) {
  K <- prep$K
  Om <- prep$omega0
  if (!is.null(from)) Om <- from$Om
  p <- numeric(0)
  for (b in spec$blocks) {
    S <- Om[b, b, drop = FALSE]
    d <- pmax(diag(S), 1e-4)
    S[S == 0 & row(S) != col(S)] <- 0
    diag(S) <- d
    L <- tryCatch(t(chol(S)), error = function(e) diag(sqrt(d), length(b)))
    diag(L) <- log(diag(L))
    p <- c(p, L[lower.tri(L, diag = TRUE)])
  }
  p <- c(p, if (is.null(from)) rep(0, K) else from$mu)
  p <- c(p, rep(0, length(spec$beta_cols)))
  # residual parameters started from the residual variance of the base run
  # at the EBEs (R0 = h' SIGMA h), projected on each variance model
  f <- unlist(lapply(prep$subj, `[[`, "ipred0"))
  R0 <- unlist(lapply(prep$subj, function(s) rowSums((s$H %*% prep$sigma0) * s$H)))
  ok <- is.finite(R0) & R0 > 0 & is.finite(f)
  f <- f[ok]; R0 <- R0[ok]; f2 <- f^2
  comb0 <- tryCatch({
    X <- cbind(1, f2); cf <- stats::lm.fit(X, R0)$coefficients
    if (any(!is.finite(cf)) || any(cf <= 0)) stop("neg")
    log(sqrt(cf))
  }, error = function(e) c(log(sqrt(min(R0))), log(sqrt(stats::median(R0 / pmax(f2, 1e-12))))))
  pw <- tryCatch({
    k <- abs(f) > 0
    cf <- stats::lm.fit(cbind(1, log(abs(f[k]))), log(R0[k]))$coefficients
    c(cf[1] / 2, cf[2] / 2)
  }, error = function(e) c(log(sqrt(stats::median(R0))), 0.5))
  p <- c(p, switch(spec$res,
                   base = 0, add = log(sqrt(mean(R0))),
                   prop = log(sqrt(stats::median(R0 / pmax(f2, 1e-12)))),
                   comb = comb0, power = unname(pw),
                   timevar = c(0, 0)))
  c(p, rep(log(0.1), length(spec$iov)))
}

.nm_ob_fit <- function(spec, prep, from = NULL) {
  rm(list = ls(prep$cache), envir = prep$cache)   # warm starts of the modes: per candidate
  st <- .nm_ob_start(spec, prep, from)
  fn <- if (identical(prep$objective, "fast")) .nm_ob_ofv_fast else .nm_ob_ofv
  o <- tryCatch(stats::optim(st, fn, spec = spec, prep = prep, method = "BFGS",
                             control = list(maxit = 1000, reltol = 1e-8)),
                error = function(e) NULL)
  if (is.null(o)) return(list(ofv = NA_real_, par = st, conv = 99L, u = NULL))
  list(ofv = o$value, par = o$par, conv = o$convergence, u = .nm_ob_unpack(o$par, spec, prep))
}

# Candidate specification.
.nm_ob_spec <- function(blocks, res = "base", beta = NULL, iov = integer(0), design = NULL) {
  bc <- character(0); be <- integer(0)
  if (length(beta)) for (b in beta) {
    cols <- design$terms[[b$covariate]]$cols
    bc <- c(bc, cols); be <- c(be, rep(b$eta, length(cols)))
  }
  xcols <- if (is.null(design)) character(0) else names(design$X)[-1]
  list(blocks = blocks, res = res, beta_cols = bc, beta_eta = be,
       beta_cols_idx = match(bc, xcols), iov = iov, beta = beta)
}

.nm_ob_lbl_beta <- function(beta, etas) {
  if (!length(beta)) return("-")
  paste(vapply(beta, function(b) paste0(b$covariate, " on ", etas[b$eta]), ""), collapse = " + ")
}


# ---- Main function --------------------------------------------------------------

#' Proposals to improve a base model without iterative NONMEM runs
#'
#' Examines a base run and proposes changes to the residual error model, the
#' random effects and the covariate model.
#'
#' **Level 1** uses only the outputs of the base run (`.phi`, `.ext`, a table
#' with `CWRES`, `IPRED` and `TIME`, and the data set): resmod-like models on
#' the CWRES (power, time-varying, autocorrelated residuals, individual
#' residual variance), EBE correlations, distribution and shrinkage, and an
#' ETA-covariate screening.
#'
#' **Level 2** needs one extra NONMEM evaluation run written by
#' [nm_derivatives_ctl()]. Its table of derivatives gives a first-order
#' linearization of the model around the EBEs; every candidate is re-estimated
#' in R under the FOCE-I objective of the linearized model: residual models
#' (additive, proportional, combined, power, time-varying), OMEGA structure
#' (full block, pairwise blocks, removal of an ETA, inter-occasion variability),
#' covariates on the ETAs (univariate tests, then forward inclusion and
#' backward elimination), and finally the combined proposal.
#'
#' The changes in objective function (`dOFV`) are those of the linearized
#' model: they approximate the changes of the real model, and the proposal
#' must be confirmed by a NONMEM run of the modified model. The quality of the
#' approximation is reported in `check` (linearized vs NONMEM OFV of the base
#' model).
#'
#' @param run Run directory (e.g. `"run1.reg"`) or path prefix of the base run.
#' @param covariates Covariates to test (columns of the data set). `NULL`: no
#'   covariate analysis.
#' @param deriv Table written by the run of [nm_derivatives_ctl()]
#'   (`<run>_deriv.dtab`). `NULL`: searched for next to the run, otherwise the
#'   function runs level 1 only (unless `nmfe_path` is given).
#' @param nmfe_path Path of nmfe: if given and no derivatives table exists, the
#'   derivatives control stream is written and run with [nm_run()].
#' @param data Data set (path or data.frame) holding the covariates. Default:
#'   the `$DATA` file of the run.
#' @param occ Occasion column for inter-occasion variability (level 2; it must
#'   be in the derivatives table: `nm_derivatives_ctl(extra = "OCC")`).
#' @param level `1` or `2` (2 runs both levels when the table is available).
#' @param cov_form Form of continuous covariate effects: `"power"`
#'   (`(COV/median)^theta`) or `"exponential"` (`exp(theta*(COV-median))`).
#' @param p_forward,p_backward Significance levels of the forward inclusion
#'   and backward elimination of covariates (linearized model).
#' @param max_levels A numeric covariate with at most this number of distinct
#'   values is treated as categorical.
#' @param objective Objective of the linearized model. `"focei"` (default):
#'   the FOCE-I objective of NONMEM with interaction (individual modes of the
#'   full individual objective, residual-variance term in the determinant);
#'   at the NONMEM estimates it reproduces the OFV of the base run. `"fast"`:
#'   an approximation (modes given by the linear predictor with the residual
#'   variance held fixed), about three times faster, whose OFV differs from
#'   NONMEM's by a few units when the residual error depends on the ETAs;
#'   the dOFV are close.
#' @param verbose Print progress.
#' @return An object of class `nm_optimize_base`: lists `level1` (`residual`,
#'   `random`, `pairs`, `covariates`) and `level2` (`check`, `residual`,
#'   `random`, `covariates`, `scm`, `proposal`), and `proposal_text`.
#' @references Khandelwal A, Harling K, Jonsson EN, Hooker AC, Karlsson MO.
#'   A fast method for testing covariates in population PK/PD models.
#'   AAPS J 2011;13:464-72.
#'   Ibrahim MMA, Nordgren R, Kjellsson MC, Karlsson MO. Model-based residual
#'   post-processing for residual model identification. AAPS J 2018;20:81.
#' @seealso [nm_derivatives_ctl()], [nm_wasserstein()]
#' @examples
#' ex <- system.file("extdata", "example", "run1.reg", package = "nmpilot")
#' o <- nm_optimize_base(ex, covariates = c("WT", "GENO"), level = 1, verbose = FALSE)
#' o
#' @export
nm_optimize_base <- function(run, covariates = NULL, deriv = NULL, nmfe_path = NULL,
                             data = NULL, occ = NULL, level = 2L,
                             cov_form = c("power", "exponential"),
                             p_forward = 0.05, p_backward = 0.01, max_levels = 5L,
                             objective = c("focei", "fast"), verbose = TRUE) {
  cov_form <- match.arg(cov_form)
  objective <- match.arg(objective)
  f <- .nm_w_run_files(run)
  if (is.na(f$ext)) stop("No .ext file for run: ", run, call. = FALSE)
  say <- function(...) if (verbose) message(...)
  fin <- .nm_vpc_read_ext_final(f$ext)
  omega <- .nm_ob_ext_matrix(fin, "OMEGA"); sigma <- .nm_ob_ext_matrix(fin, "SIGMA")
  ofv_nm <- as.numeric(fin$OBJ[1])
  phi <- .nm_w_read_phi(f$phi)
  eta_cols <- grep("^ETA\\(", names(phi), value = TRUE)
  eta <- as.matrix(phi[, eta_cols, drop = FALSE])
  K <- ncol(eta)
  etas <- paste0("ETA", seq_len(K))
  if (!is.na(f$ctl)) {
    pm <- tryCatch(nm_identify_parameters(f$ctl, likely_pk_only = FALSE)$mapping_detail_all,
                   error = function(e) NULL)
    if (!is.null(pm)) {
      pm <- pm[pm$element_type == "ETA", , drop = FALSE]
      for (k in seq_len(K)) {
        h <- pm[pm$index == k, , drop = FALSE]
        h <- h$pk_parameter[order(if ("depth" %in% names(h)) h$depth else seq_len(nrow(h)))]
        if (length(h) && !is.na(h[1]) && nzchar(h[1])) etas[k] <- paste0("ETA", k, " (", h[1], ")")
      }
    }
  }
  colnames(eta) <- etas
  shrink <- 1 - apply(eta, 2, stats::sd) / sqrt(pmax(diag(omega)[seq_len(K)], 1e-12))

  cov <- info <- NULL
  if (length(covariates)) {
    cov <- .nm_w_covariates(data, f$ctl)
    miss <- setdiff(covariates, names(cov))
    if (length(miss)) stop("Covariate(s) not in the data: ", paste(miss, collapse = ", "), call. = FALSE)
    info <- .nm_ob_cov_info(cov, covariates, max_levels)
  }

  # ---- Level 1
  say("Level 1: diagnostics from the outputs of ", f$stem)
  tabf <- .nm_ob_find_table(f$dir, c("ID", "CWRES"))
  l1_res <- NULL
  if (!is.na(tabf)) {
    tab <- .nm_read_table(tabf)
    if (!"IPRED" %in% names(tab)) tab$IPRED <- tab[[intersect(c("CIPREDI", "CIPRED", "PRED"), names(tab))[1]]]
    if ("MDV" %in% names(tab)) tab <- tab[tab$MDV == 0, , drop = FALSE]
    if ("EVID" %in% names(tab)) tab <- tab[tab$EVID == 0, , drop = FALSE]
    if (!"TIME" %in% names(tab)) tab$TIME <- stats::ave(tab$ID, tab$ID, FUN = seq_along)
    l1_res <- .nm_ob_l1_residual(tab)
  } else say("  no table with CWRES found: residual diagnostics skipped")
  l1_ran <- .nm_ob_l1_random(eta, omega)
  l1_cov <- if (length(covariates)) .nm_ob_l1_covariates(eta, phi$ID, cov, info, shrink) else NULL
  out <- list(run = f$stem, ofv_nonmem = ofv_nm, etas = etas,
              level1 = list(residual = l1_res, random = l1_ran$eta, pairs = l1_ran$pairs,
                            covariates = l1_cov),
              level2 = NULL)

  # ---- Level 2
  if (level >= 2L) {
    if (is.null(deriv)) {
      ad <- if (grepl("\\.(reg|nm7)$", f$dir, ignore.case = TRUE)) dirname(f$dir) else f$dir
      nm_d <- paste0(f$stem, "_deriv")
      cand <- c(file.path(ad, paste0(nm_d, c(".reg", ".nm7", "")), paste0(nm_d, ".dtab")),
                file.path(f$dir, paste0(nm_d, ".dtab")))
      deriv <- cand[file.exists(cand)][1]
      if (is.na(deriv) && !is.null(nmfe_path)) {
        say("Level 2: writing and running ", nm_d, ".ctl")
        ctl_d <- nm_derivatives_ctl(run, extra = occ)
        nm_run(nmfe_path, ctl_files = basename(ctl_d), analysis_dir = ad,
               skip_terminated = FALSE, verbose = verbose)
        deriv <- cand[file.exists(cand)][1]
      }
    }
    if (is.null(deriv) || is.na(deriv) || !file.exists(deriv)) {
      say("Level 2 skipped: no derivatives table. Write it with nm_derivatives_ctl(\"",
          run, "\") and run it, or give `nmfe_path`.")
    } else {
      say("Level 2: linearized model from ", basename(deriv))
      out$level2 <- .nm_ob_level2(.nm_read_table(deriv), omega, sigma, ofv_nm, etas, cov, info,
                                  occ, cov_form, p_forward, p_backward, say,
                                  objective)
    }
  }
  out$proposal_text <- .nm_ob_proposal_text(out)
  class(out) <- "nm_optimize_base"
  out
}

.nm_ob_level2 <- function(dt, omega, sigma, ofv_nm, etas, cov, info, occ, cov_form,
                          p_forward, p_backward, say, objective = "focei") {
  K <- nrow(omega)
  design <- NULL
  if (!is.null(info)) design <- .nm_ob_cov_design(cov, unique(dt$ID), info, cov_form)
  if (!is.null(occ) && !occ %in% names(dt)) {
    warning("Occasion column '", occ, "' not in the derivatives table: IOV not tested.", call. = FALSE)
    occ <- NULL
  }
  prep <- .nm_ob_prepare(dt, omega, sigma, if (is.null(design)) NULL else design$X, occ)
  prep$objective <- objective
  base_blocks <- .nm_ob_blocks(omega)
  base_blocks <- base_blocks[vapply(base_blocks, function(b) any(diag(omega)[b] > 0), logical(1))]

  # Check: base model at the NONMEM estimates, then re-estimated
  sp0 <- .nm_ob_spec(base_blocks, "base", design = design)
  st0 <- .nm_ob_start(sp0, prep)
  ofv_at_nm <- if (objective == "fast") .nm_ob_ofv_fast(st0, sp0, prep) else .nm_ob_ofv(st0, sp0, prep)
  fb <- .nm_ob_fit(sp0, prep)
  check <- data.frame(item = c("OFV NONMEM (base run)",
                               sprintf("OFV linearized at the NONMEM estimates (%s)",
                                       if (objective == "fast") "fast approximation" else "FOCE-I"),
                               "OFV linearized re-estimated (reference of the dOFV)",
                               "subjects", "observations"),
                      value = c(round(ofv_nm, 2), round(ofv_at_nm, 2), round(fb$ofv, 2),
                                length(prep$subj), prep$n_obs), stringsAsFactors = FALSE)
  if (abs(ofv_at_nm - ofv_nm) > max(10, 0.01 * abs(ofv_nm)))
    warning("Linearized OFV differs from the NONMEM OFV by ", round(ofv_at_nm - ofv_nm, 1),
            ": the linearization may be poor (check the derivatives table).", call. = FALSE)
  row <- function(area, cand, fit, df) data.frame(
    area = area, candidate = cand, df = df, ofv = round(fit$ofv, 2),
    dOFV = round(fit$ofv - fb$ofv, 2),
    p = signif(ifelse(df > 0, stats::pchisq(pmax(fb$ofv - fit$ofv, 0), df, lower.tail = FALSE),
                      ifelse(df < 0, stats::pchisq(pmax(fit$ofv - fb$ofv, 0), -df, lower.tail = FALSE), NA)), 3),
    stringsAsFactors = FALSE)

  # Residual error
  say("  residual error models")
  res_types <- c(add = "additive", prop = "proportional", comb = "combined",
                 power = "power", timevar = sprintf("time-varying (cut at %g)", signif(prep$tcut, 3)))
  res_fit <- lapply(names(res_types), function(t) .nm_ob_fit(.nm_ob_spec(base_blocks, t, design = design), prep, fb$u))
  names(res_fit) <- names(res_types)
  rtab <- do.call(rbind, lapply(names(res_types), function(t)
    row("residual", res_types[[t]], res_fit[[t]], .nm_ob_npar(.nm_ob_spec(base_blocks, t), prep)[["res"]] - 1L)))
  rtab$estimates <- vapply(names(res_types), function(t) {
    r <- res_fit[[t]]$u$res
    if (is.null(r)) return(NA_character_)
    switch(t, add = sprintf("SD=%.3g", exp(r[1])), prop = sprintf("CV=%.3g", exp(r[1])),
           comb = sprintf("SD=%.3g, CV=%.3g", exp(r[1]), exp(r[2])),
           power = sprintf("SD=%.3g*F^%.3g", exp(r[1]), r[2]),
           timevar = sprintf("SD ratio late/early=%.3g", exp(r[2] - r[1])))
  }, character(1))
  rtab$note <- "re-estimated in linear space; base = residual model of the run (scale re-estimated)"

  # Random effects
  say("  random effects")
  kept <- sort(unlist(base_blocks))
  rand <- list()
  if (length(kept) > 1L) {
    full <- list(kept)
    if (!identical(base_blocks, full)) rand[["full OMEGA block"]] <- full
    for (pr in utils::combn(kept, 2, simplify = FALSE)) {
      bb <- base_blocks
      if (any(vapply(bb, function(b) all(pr %in% b), logical(1)))) next
      bb <- bb[!vapply(bb, function(b) any(b %in% pr), logical(1))]
      merged <- sort(unique(unlist(base_blocks[vapply(base_blocks, function(b) any(b %in% pr), logical(1))])))
      bb <- c(bb, list(merged))
      rand[[sprintf("block %s", paste(etas[merged], collapse = "-"))]] <- bb
    }
  }
  for (k in kept) {
    bb <- lapply(base_blocks, setdiff, k); bb <- bb[lengths(bb) > 0]
    rand[[sprintf("remove %s", etas[k])]] <- bb
  }
  ran_rows <- list()
  for (nm in names(rand)) {
    sp <- .nm_ob_spec(rand[[nm]], "base", design = design)
    fit <- .nm_ob_fit(sp, prep, fb$u)
    df <- .nm_ob_npar(sp, prep)[["omega"]] - .nm_ob_npar(sp0, prep)[["omega"]]
    r <- row("random", nm, fit, df)
    r$estimates <- if (grepl("^remove", nm)) "-" else {
      Om <- fit$u$Om; cc <- stats::cov2cor(Om[kept, kept, drop = FALSE])
      idx <- which(lower.tri(cc) & abs(cc) > 0, arr.ind = TRUE)
      if (!nrow(idx)) "-" else paste(sprintf("r(%s,%s)=%.2f", etas[kept][idx[, 2]], etas[kept][idx[, 1]],
                                             cc[idx]), collapse = "; ")
    }
    ran_rows[[nm]] <- r
  }
  if (!is.null(occ)) for (k in kept) {
    sp <- .nm_ob_spec(base_blocks, "base", iov = k, design = design)
    fit <- .nm_ob_fit(sp, prep, fb$u)
    r <- row("random", sprintf("IOV on %s", etas[k]), fit, 1L)
    r$estimates <- sprintf("omega_IOV=%.3g", fit$u$iov); ran_rows[[r$candidate]] <- r
  }
  rantab <- if (length(ran_rows)) do.call(rbind, ran_rows) else NULL
  if (!is.null(rantab)) rownames(rantab) <- NULL

  # Covariates
  covtab <- scm <- NULL; sel <- list()
  if (!is.null(design)) {
    say("  covariates: univariate tests")
    pairs <- expand.grid(k = kept, cv = names(design$terms), stringsAsFactors = FALSE)
    uni_fit <- lapply(seq_len(nrow(pairs)), function(i) {
      b <- list(list(eta = pairs$k[i], covariate = pairs$cv[i]))
      .nm_ob_fit(.nm_ob_spec(base_blocks, "base", beta = b, design = design), prep, fb$u)
    })
    cov_rows <- lapply(seq_len(nrow(pairs)), function(i) {
      b <- list(list(eta = pairs$k[i], covariate = pairs$cv[i]))
      sp <- .nm_ob_spec(base_blocks, "base", beta = b, design = design)
      fit <- uni_fit[[i]]
      r <- row("covariate", sprintf("%s on %s", pairs$cv[i], etas[pairs$k[i]]), fit,
               design$terms[[pairs$cv[i]]]$df)
      r$form <- design$terms[[pairs$cv[i]]]$label
      r$theta <- paste(signif(fit$u$beta[pairs$k[i], seq_along(sp$beta_cols)], 3), collapse = ", ")
      r
    })
    covtab <- do.call(rbind, cov_rows)
    covtab <- covtab[order(covtab$dOFV), , drop = FALSE]; rownames(covtab) <- NULL

    say("  covariates: forward inclusion / backward elimination")
    cur_ofv <- fb$ofv; cur_u <- fb$u; steps <- list()
    pool <- lapply(seq_len(nrow(pairs)), function(i) list(eta = pairs$k[i], covariate = pairs$cv[i]))
    repeat {
      rest <- pool[!vapply(pool, function(b) any(vapply(sel, identical, logical(1), b)), logical(1))]
      if (!length(rest)) break
      tr <- if (!length(sel)) uni_fit[match(vapply(rest, function(b) paste(b$eta, b$covariate), ""),
                                            vapply(pool, function(b) paste(b$eta, b$covariate), ""))]
            else lapply(rest, function(b) {
              sp <- .nm_ob_spec(base_blocks, "base", beta = c(sel, list(b)), design = design)
              .nm_ob_fit(sp, prep, cur_u)
            })
      d <- vapply(tr, `[[`, numeric(1), "ofv") - cur_ofv
      dfs <- vapply(rest, function(b) design$terms[[b$covariate]]$df, numeric(1))
      crit <- stats::qchisq(1 - p_forward, dfs)
      ok <- which(-d > crit)
      if (!length(ok)) break
      best <- ok[which.min(d[ok] + crit[ok])]
      sel <- c(sel, list(rest[[best]]))
      steps[[length(steps) + 1L]] <- data.frame(step = "forward", change = paste("+", .nm_ob_lbl_beta(list(rest[[best]]), etas)),
                                                dOFV = round(d[best], 2), model = .nm_ob_lbl_beta(sel, etas),
                                                stringsAsFactors = FALSE)
      cur_ofv <- tr[[best]]$ofv; cur_u <- tr[[best]]$u
    }
    repeat {
      if (!length(sel)) break
      tr <- lapply(seq_along(sel), function(i) .nm_ob_fit(.nm_ob_spec(base_blocks, "base", beta = sel[-i], design = design), prep, cur_u))
      d <- vapply(tr, `[[`, numeric(1), "ofv") - cur_ofv
      dfs <- vapply(sel, function(b) design$terms[[b$covariate]]$df, numeric(1))
      crit <- stats::qchisq(1 - p_backward, dfs)
      ok <- which(d < crit)
      if (!length(ok)) break
      worst <- ok[which.min(d[ok] - crit[ok])]
      steps[[length(steps) + 1L]] <- data.frame(step = "backward", change = paste("-", .nm_ob_lbl_beta(list(sel[[worst]]), etas)),
                                                dOFV = round(d[worst], 2), model = .nm_ob_lbl_beta(sel[-worst], etas),
                                                stringsAsFactors = FALSE)
      sel <- sel[-worst]; cur_ofv <- tr[[worst]]$ofv; cur_u <- tr[[worst]]$u
    }
    scm <- if (length(steps)) do.call(rbind, steps) else
      data.frame(step = "-", change = "no covariate retained", dOFV = NA, model = "-")
  }

  # Combined proposal
  say("  combined proposal")
  crit1 <- stats::qchisq(0.99, 1)
  best_res <- "base"
  if (any(rtab$dOFV < -crit1 * pmax(rtab$df, 1), na.rm = TRUE)) {
    adj <- rtab$dOFV + crit1 * pmax(rtab$df, 1)
    best_res <- names(res_types)[which.min(adj)]
  }
  blocks <- base_blocks; iov <- integer(0); ran_sel <- character(0)
  if (!is.null(rantab)) {
    blk <- rantab[grepl("block", rantab$candidate) & rantab$dOFV < -stats::qchisq(0.99, pmax(rantab$df, 1)), , drop = FALSE]
    if (nrow(blk)) {
      pick <- blk$candidate[which.min(blk$dOFV + stats::qchisq(0.99, pmax(blk$df, 1)))]
      blocks <- rand[[pick]]; ran_sel <- c(ran_sel, pick)
    }
    rm <- rantab[grepl("^remove", rantab$candidate) & rantab$dOFV < stats::qchisq(0.95, 1), , drop = FALSE]
    for (cnd in rm$candidate) {
      k <- kept[etas[kept] == sub("^remove ", "", cnd)]
      if (!length(unlist(lapply(blocks, setdiff, k)))) next
      blocks <- lapply(blocks, setdiff, k); blocks <- blocks[lengths(blocks) > 0]
      ran_sel <- c(ran_sel, cnd)
    }
    io <- rantab[grepl("^IOV", rantab$candidate) & rantab$dOFV < -crit1, , drop = FALSE]
    for (cnd in io$candidate) {
      k <- kept[etas[kept] == sub("^IOV on ", "", cnd)]
      if (k %in% unlist(blocks)) { iov <- c(iov, k); ran_sel <- c(ran_sel, cnd) }
    }
  }
  sel <- Filter(function(b) b$eta %in% unlist(blocks), sel)
  spf <- .nm_ob_spec(blocks, best_res, beta = sel, iov = iov, design = design)
  ff <- .nm_ob_fit(spf, prep, NULL)
  dff <- sum(.nm_ob_npar(spf, prep)) - sum(.nm_ob_npar(sp0, prep))
  proposal <- data.frame(
    part = c("residual error", "random effects", "covariates", "joint (linearized)"),
    proposal = c(if (best_res == "base") "keep the residual model" else res_types[[best_res]],
                 if (length(ran_sel)) paste(ran_sel, collapse = "; ") else "keep the OMEGA structure",
                 if (length(sel)) .nm_ob_lbl_beta(sel, etas) else "no covariate",
                 sprintf("dOFV = %.1f for %+d parameters", ff$ofv - fb$ofv, dff)),
    stringsAsFactors = FALSE)

  list(check = check, residual = rtab, random = rantab, covariates = covtab, scm = scm,
       proposal = proposal, fit = list(base = fb, proposal = ff, spec = spf), etas = etas)
}

.nm_ob_proposal_text <- function(x) {
  txt <- character(0)
  r <- x$level1$residual
  if (!is.null(r)) for (i in which(r$suggestion != "-"))
    txt <- c(txt, sprintf("[level 1, residual] %s (dOFV %.1f)", r$suggestion[i], r$dOFV[i]))
  e <- x$level1$random
  for (i in which(e$suggestion != "-")) txt <- c(txt, sprintf("[level 1, random] %s: %s", e$eta[i], e$suggestion[i]))
  p <- x$level1$pairs
  if (!is.null(p)) for (i in which(p$suggestion != "-" & p$suggestion != "already in an OMEGA block"))
    txt <- c(txt, sprintf("[level 1, random] %s-%s (r = %.2f): %s", p$eta_1[i], p$eta_2[i], p$r[i], p$suggestion[i]))
  cv <- x$level1$covariates
  if (!is.null(cv)) for (i in which(cv$suggestion != "-"))
    txt <- c(txt, sprintf("[level 1, covariate] %s on %s: %s (adjusted p %s)", cv$covariate[i], cv$eta[i], cv$suggestion[i], format(cv$p_adj[i])))
  if (!is.null(x$level2)) {
    pr <- x$level2$proposal
    txt <- c(txt, sprintf("[level 2] %s: %s", pr$part, pr$proposal))
  }
  if (!length(txt)) txt <- "No change suggested."
  txt
}

#' @rdname nm_optimize_base
#' @param x An `nm_optimize_base` object.
#' @param ... Ignored.
#' @export
print.nm_optimize_base <- function(x, ...) {
  cat(sprintf("Optimization of base run %s (NONMEM OFV %.2f)\n", x$run, x$ofv_nonmem))
  cat("\n== Level 1: outputs of the base run ==\n")
  if (!is.null(x$level1$residual)) {
    cat("\nResidual error (models on CWRES, dOFV vs CWRES ~ N(mu, sigma) with random intercept):\n")
    print(x$level1$residual[, c("model", "dOFV", "estimate", "meaning", "suggestion")], row.names = FALSE)
  }
  cat("\nRandom effects:\n"); print(x$level1$random, row.names = FALSE)
  if (!is.null(x$level1$pairs)) { cat("\nEBE correlations:\n"); print(x$level1$pairs, row.names = FALSE) }
  if (!is.null(x$level1$covariates)) {
    cat("\nETA-covariate screening (Bonferroni-adjusted):\n")
    print(x$level1$covariates[, c("eta", "covariate", "statistic", "p", "p_adj", "shrinkage_pct", "suggestion")], row.names = FALSE)
  }
  if (!is.null(x$level2)) {
    l2 <- x$level2
    cat("\n== Level 2: linearized model ==\n")
    print(l2$check, row.names = FALSE)
    cat("\nResidual error:\n"); print(l2$residual[, c("candidate", "df", "dOFV", "p", "estimates")], row.names = FALSE)
    if (!is.null(l2$random)) { cat("\nRandom effects:\n"); print(l2$random[, c("candidate", "df", "dOFV", "p", "estimates")], row.names = FALSE) }
    if (!is.null(l2$covariates)) {
      cat("\nCovariates (univariate):\n"); print(l2$covariates[, c("candidate", "df", "dOFV", "p", "form", "theta")], row.names = FALSE)
      cat("\nStepwise selection:\n"); print(l2$scm, row.names = FALSE)
    }
    cat("\nProposal:\n"); print(l2$proposal, row.names = FALSE)
    cat("\ndOFV are those of the linearized model: confirm with a NONMEM run of the modified model.\n")
  }
  cat("\nSummary:\n"); cat(paste0("  - ", x$proposal_text), sep = "\n")
  invisible(x)
}

#' @rdname nm_optimize_base
#' @param area For `plot()`: `"all"`, `"residual"`, `"random"` or `"covariate"`.
#' @export
plot.nm_optimize_base <- function(x, area = c("all", "residual", "random", "covariate"), ...) {
  area <- match.arg(area)
  if (is.null(x$level2)) {
    r <- x$level1$residual
    if (is.null(r)) stop("Nothing to plot: no level 2 and no level-1 residual diagnostics.")
    d <- data.frame(area = "residual (CWRES, level 1)", candidate = r$model, dOFV = r$dOFV, df = 1)
  } else {
    l2 <- x$level2
    d <- rbind(l2$residual[, c("area", "candidate", "dOFV", "df")],
               if (!is.null(l2$random)) l2$random[, c("area", "candidate", "dOFV", "df")],
               if (!is.null(l2$covariates)) l2$covariates[, c("area", "candidate", "dOFV", "df")])
    if (area != "all") d <- d[d$area == area, , drop = FALSE]
  }
  # Effect class: better / worse at p < 0.01 (chi-square with |df| degrees of
  # freedom, 1 for non-nested alternatives with df = 0), otherwise unchanged.
  # A removal (df < 0) whose dOFV exceeds the threshold is a worse model:
  # the random effect must be kept.
  crit <- stats::qchisq(0.99, pmax(abs(d$df), 1))
  d$effect <- ifelse(is.na(d$dOFV), "not estimated",
                     ifelse(d$dOFV < -crit, "better (p < 0.01)",
                            ifelse(d$dOFV > crit, "worse (p < 0.01)", "no significant change")))
  d$effect <- factor(d$effect, levels = c("better (p < 0.01)", "no significant change",
                                          "worse (p < 0.01)", "not estimated"))
  d <- d[order(d$area, d$dOFV), , drop = FALSE]
  d$candidate <- factor(paste(d$area, d$candidate, sep = "\r"),
                        levels = rev(unique(paste(d$area, d$candidate, sep = "\r"))))
  d$label <- ifelse(is.na(d$dOFV), "NA", formatC(round(d$dOFV, 1) + 0, format = "f", digits = 1))
  d$label[d$label == "-0.0"] <- "0.0"
  d$hjust <- ifelse(!is.na(d$dOFV) & d$dOFV < 0, 1.1, -0.1)
  d$dOFV0 <- ifelse(is.na(d$dOFV), 0, d$dOFV)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$dOFV0, y = .data$candidate, fill = .data$effect)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, hjust = .data$hjust), size = 3) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    ggplot2::scale_y_discrete(labels = function(z) sub("^.*\r", "", z)) +
    # one panel per area, each with its own dOFV scale
    ggplot2::facet_wrap(ggplot2::vars(.data$area), ncol = 1, scales = "free") +
    ggplot2::scale_fill_manual(values = c(`better (p < 0.01)` = "#2a6f97",
                                          `no significant change` = "grey70",
                                          `worse (p < 0.01)` = "#c0504d",
                                          `not estimated` = "white"),
                               drop = TRUE, name = NULL) +
    ggplot2::labs(x = "dOFV (linearized model)", y = NULL,
                  title = paste("Candidate changes to", x$run),
                  subtitle = "Removal of an ETA in red: the ETA is needed") +
    ggplot2::theme_bw()
}
