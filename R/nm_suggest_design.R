# =============================================================================
# nm_suggest_design.R
# Expected precision of a candidate model for a sampling design: population
# Fisher information matrix (first-order approximation), relative standard
# errors and expected eta-shrinkage, computed before any NONMEM run.
# =============================================================================

# ---------------------------------------------------------------------------
# Analytic predictions (linear models, superposition of doses)
# ---------------------------------------------------------------------------

# Macro constants of a unit bolus in the central compartment:
# C(t) = sum_i A_i exp(-lam_i t), from the micro parameters (CL, volumes, Q).
.nm_sd_macro <- function(theta, st) {
  v1 <- theta[[st$vc]]
  k10 <- theta[["CL"]] / v1
  n <- 1L + length(st$peri)
  if (n == 1L) return(list(A = 1 / v1, lam = k10))
  K <- matrix(0, n, n)
  K[1, 1] <- -k10
  for (j in seq_along(st$peri)) {
    q <- theta[[st$peri[[j]]$q]]; vp <- theta[[st$peri[[j]]$vol]]
    K[1, 1] <- K[1, 1] - q / v1
    K[1, j + 1] <- q / vp
    K[j + 1, 1] <- q / v1
    K[j + 1, j + 1] <- -q / vp
  }
  e <- eigen(K)
  vec <- Re(e$vectors)
  A <- vec[1, ] * solve(vec)[, 1] / v1
  list(A = A, lam = -Re(e$values))
}

# Predicted concentrations (DV scale) of one subject.
.nm_sd_predict <- function(theta, sp, st, doses, times, factor) {
  m <- .nm_sd_macro(theta, st)
  out <- numeric(length(times))
  for (d in seq_len(nrow(doses))) {
    tt <- times - doses$TIME[d]
    pos <- tt > 0
    if (!any(pos)) next
    t1 <- tt[pos]
    contrib <- numeric(length(t1))
    for (i in seq_along(m$lam)) {
      lam <- m$lam[i]
      lg <- if (sp$route == "IV") {
        if (doses$TINF[d] > 0) .nm_sm_g_zo(t1, lam, doses$TINF[d]) else -lam * t1
      } else switch(sp$input,
        `first-order` = .nm_sm_g_fo(t1, lam, theta[["KA"]]),
        `first-order+lag` = .nm_sm_g_fo(t1 - theta[["ALAG1"]], lam, theta[["KA"]]),
        `zero-order-d1` = .nm_sm_g_zo(t1, lam, theta[["DUR"]]),
        transit = {
          a <- sp$transit_n + 1; ktr <- a / theta[["MTT"]]
          if (ktr <= lam) stop("Transit rate below a disposition rate: analytic solution not available.")
          .nm_sm_g_tr(t1, lam, ktr, a)
        })
      contrib <- contrib + m$A[i] * exp(lg)
    }
    out[pos] <- out[pos] + doses$AMT[d] * contrib
  }
  out * factor
}

# ---------------------------------------------------------------------------
# Population FIM (first-order approximation, block-diagonal)
# ---------------------------------------------------------------------------
.nm_sd_evaluate <- function(sp, design, times = NULL, n_subjects = NULL) {
  if (sp$elimination != "linear")
    stop("Design evaluation is available for linear elimination only.")
  st <- .nm_ctl_struct(sp$cpts, sp$route, sp$absorption, sp$elimination, sp$transit_n)
  theta <- sp$theta
  if (!is.null(sp$alag)) theta <- c(theta, ALAG1 = sp$alag)
  eta_p <- names(sp$iiv)
  omega <- .nm_ctl_cv_to_omega(sp$iiv)
  names(omega) <- eta_p
  sig <- sp$sigma_used
  s_add <- if ("add" %in% names(sig)) sig[["add"]] else 0
  s_prop <- if ("prop" %in% names(sig)) sig[["prop"]] else 0
  factor <- design$factor
  ids <- unique(design$obs$ID)
  np <- length(theta); ne <- length(eta_p); ns <- length(sig)
  Fb <- matrix(0, np, np); Fv <- matrix(0, ne + ns, ne + ns)
  shr <- numeric(ne); n_used <- 0L
  h <- 1e-4
  Om <- diag(omega, nrow = ne)
  for (id in ids) {
    dd <- design$doses[design$doses$ID == id, , drop = FALSE]
    tt <- if (is.null(times)) design$obs$TIME[design$obs$ID == id] else min(dd$TIME) + times
    if (!nrow(dd) || !length(tt)) next
    f <- .nm_sd_predict(theta, sp, st, dd, tt, factor)
    keep <- is.finite(f) & f > 1e-10 * max(f, na.rm = TRUE)
    if (sum(keep) < 1) next
    f <- f[keep]; tt <- tt[keep]
    J <- sapply(seq_len(np), function(k) {
      up <- theta; dn <- theta
      up[k] <- theta[k] * (1 + h); dn[k] <- theta[k] * (1 - h)
      (.nm_sd_predict(up, sp, st, dd, tt, factor) - .nm_sd_predict(dn, sp, st, dd, tt, factor)) /
        (2 * h * theta[k])
    })
    if (is.null(dim(J))) J <- matrix(J, nrow = length(f))
    H <- J[, match(eta_p, names(theta)), drop = FALSE] %*% diag(theta[eta_p], nrow = ne)
    rv <- s_add + s_prop * f^2
    V <- H %*% Om %*% t(H) + diag(rv, nrow = length(f))
    Vi <- tryCatch(solve(V), error = function(e) NULL)
    if (is.null(Vi)) next
    Fb <- Fb + t(J) %*% Vi %*% J
    dV <- c(lapply(seq_len(ne), function(k) H[, k] %*% t(H[, k])),
            lapply(names(sig), function(s) if (s == "add") diag(1, length(f)) else diag(f^2, nrow = length(f))))
    VdV <- lapply(dV, function(d) Vi %*% d)
    for (a in seq_along(dV)) for (b in a:length(dV)) {
      v <- 0.5 * sum(VdV[[a]] * t(VdV[[b]]))
      Fv[a, b] <- Fv[a, b] + v
      if (a != b) Fv[b, a] <- Fv[b, a] + v
    }
    # posterior variance of the ETAs (linearised model) -> shrinkage
    if (ne) {
      C <- tryCatch(solve(t(H) %*% diag(1 / rv, nrow = length(f)) %*% H + solve(Om)), error = function(e) NULL)
      if (!is.null(C)) shr <- shr + diag(C) / omega
    }
    n_used <- n_used + 1L
  }
  if (!n_used) stop("No subject with a usable design.")
  scale <- if (!is.null(n_subjects)) n_subjects / n_used else 1
  Fb <- Fb * scale; Fv <- Fv * scale
  inv <- function(M) tryCatch(solve(M), error = function(e) matrix(NA_real_, nrow(M), ncol(M)))
  se_b <- sqrt(pmax(diag(inv(Fb)), 0))
  se_v <- if (nrow(Fv)) sqrt(pmax(diag(inv(Fv)), 0)) else numeric(0)
  val_v <- c(omega, sig)
  pars <- data.frame(
    parameter = c(names(theta), if (ne) paste0("OMEGA(", eta_p, ")"), if (ns) paste0("SIGMA(", names(sig), ")")),
    type = c(rep("theta", np), rep("omega", ne), rep("sigma", ns)),
    value = c(unname(theta), unname(val_v)),
    se = c(se_b, se_v), stringsAsFactors = FALSE)
  pars$rse_pct <- round(100 * pars$se / pars$value, 1)
  shrink <- data.frame(parameter = eta_p,
                       shrinkage_pct = if (ne) round(100 * (1 - sqrt(pmax(1 - shr / n_used, 0))), 1) else numeric(0),
                       stringsAsFactors = FALSE)
  structure(list(parameters = pars, shrinkage = shrink,
                 n_subjects = if (is.null(n_subjects)) n_used else n_subjects,
                 n_obs = nrow(design$obs), spec = sp,
                 times = times),
            class = "nm_suggest_design")
}

#' Expected precision of a suggested model for a sampling design
#'
#' Computes, before any NONMEM run, the expected relative standard errors
#' (RSE) of the parameters of a candidate model from [nm_suggest_model()] and
#' the expected eta-shrinkage, for the sampling design of the dataset or for
#' another design. The population Fisher information matrix is computed with
#' the first-order (FO) approximation, block-diagonal between fixed effects
#' and variance parameters (as in PopED's reduced FIM), with numerical
#' derivatives; predictions use the analytic solution of the linear model,
#' with superposition of all doses (`ADDL`/`II` expanded). The expected
#' shrinkage of each ETA is `1 - sqrt(1 - mean(C_kk / omega_k))`, where
#' `C` is the posterior covariance of the ETAs of each subject in the
#' linearised model.
#'
#' The values depend on the initial values of the candidate (THETA, IIV,
#' SIGMA), which are approximate: they rank designs and parameters rather
#' than predict the RSE reported by NONMEM. Michaelis-Menten candidates are
#' not supported (no analytic solution).
#'
#' @param x An `nm_model_suggestion` object.
#' @param candidate Candidate number (default 1).
#' @param times Optional sampling times after the first dose, applied to every
#'   subject in place of the observed times (to compare designs).
#' @param n_subjects Optional number of subjects: the information is scaled
#'   from the subjects of the dataset to this number.
#' @return An object of class `nm_suggest_design`: `parameters` (parameter,
#'   type, value, se, rse_pct), `shrinkage` (parameter, shrinkage_pct),
#'   `n_subjects`, `n_obs`, `spec`, `times`.
#' @seealso [nm_suggest_model()]
#' @examples
#' ex <- system.file("extdata", "example", "example_data.csv", package = "nmtools")
#' s <- nm_suggest_model(ex, verbose = FALSE, design_eval = FALSE, max_profiles = 15)
#' nm_suggest_design(s)
#' # a sparse design with 3 samples per subject
#' nm_suggest_design(s, times = c(1, 4, 12))
#' @export
nm_suggest_design <- function(x, candidate = 1L, times = NULL, n_subjects = NULL) {
  if (!inherits(x, "nm_model_suggestion")) stop("`x` must be an nm_model_suggestion object.")
  if (candidate < 1 || candidate > length(x$specs)) stop("`candidate` must be between 1 and ", length(x$specs), ".")
  .nm_sd_evaluate(x$specs[[candidate]], x$design_data, times = times, n_subjects = n_subjects)
}

#' @rdname nm_suggest_design
#' @param ... Ignored.
#' @export
print.nm_suggest_design <- function(x, ...) {
  sp <- x$spec
  cat(sprintf("Expected precision: %d cpt | %s | %s error | %d subjects%s\n", sp$cpts, sp$input,
              sp$error, x$n_subjects, if (!is.null(x$times)) " | design: times given" else " | design of the dataset"))
  p <- x$parameters
  p$value <- signif(p$value, 3); p$se <- signif(p$se, 3)
  print(p, row.names = FALSE)
  if (nrow(x$shrinkage)) {
    cat("\nExpected eta-shrinkage (%)\n")
    print(x$shrinkage, row.names = FALSE)
  }
  invisible(x)
}

#' @rdname nm_suggest_design
#' @export
plot.nm_suggest_design <- function(x, ...) {
  p <- x$parameters
  p$what <- "RSE (%)"; p$y <- p$rse_pct
  s <- x$shrinkage
  d <- p[, c("parameter", "what", "y")]
  if (nrow(s)) d <- rbind(d, data.frame(parameter = paste0("ETA(", s$parameter, ")"),
                                         what = "eta-shrinkage (%)", y = s$shrinkage_pct))
  d$parameter <- factor(d$parameter, levels = rev(unique(d$parameter)))
  lim <- data.frame(what = c("RSE (%)", "eta-shrinkage (%)"), y = c(50, 30))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$y, y = .data$parameter)) +
    ggplot2::geom_col(fill = "#4393C3") +
    ggplot2::geom_vline(data = lim, ggplot2::aes(xintercept = .data$y), linetype = 2, colour = "firebrick") +
    ggplot2::facet_wrap(~what, scales = "free") +
    ggplot2::labs(x = NULL, y = NULL, title = "Expected precision for the sampling design",
                  subtitle = "Dashed lines: RSE 50%, shrinkage 30%") +
    ggplot2::theme_bw()
}
