# =============================================================================
# nm_suggest_model.R
# Data-driven suggestion of structural PK models before any NONMEM run:
# number of compartments, input model, elimination, residual error model,
# initial values (THETA, IIV, SIGMA), parameters likely to be hard to
# estimate, and alerts on features outside the tested structures.
# =============================================================================

#' Suggest structural PK models from a NONMEM dataset
#'
#' Analyses a NONMEM dataset (rich or sparse data) without running NONMEM and
#' ranks candidate structural models: number of compartments (1 to 3), input
#' model (IV bolus or infusion, first-order, first-order with lag time,
#' zero-order with estimated duration, transit compartments), elimination
#' (linear or Michaelis-Menten) and residual error model (additive,
#' proportional or combined). Each candidate comes with initial values for
#' THETA, IIV and SIGMA in the parameter names of [nm_write_ctl()], the
#' parameters likely to be hard to estimate and, for linear models, the
#' precision expected from the sampling design ([nm_suggest_design()]).
#'
#' Method (heuristic; it chooses a starting model and does not replace model
#' comparison in NONMEM):
#' 1. Records are split into dosing intervals; `ADDL`/`II` doses are
#'    expanded. An interval is "rich" when it holds at least `min_obs_rich`
#'    quantified concentrations, with at least 3 points after the peak for
#'    extravascular data.
#' 2. NCA on each rich interval (linear-up log-down AUC, `lambda_z` by best
#'    adjusted R2, CL = Dose/AUCinf after a first dose, Dose/AUCtau at steady
#'    state).
#' 3. Number of compartments: sums of 1 to 3 exponentials fitted to the
#'    post-peak phase (log scale); BIC weights averaged over intervals. A
#'    phase is kept only when it is identifiable (rates at least 1.5-fold
#'    apart, each amplitude at least 1% of the total, fast phase decaying
#'    within the sampling range).
#' 4. Input (extravascular data): first-order, first-order + lag,
#'    zero-order (duration) and gamma (transit) inputs, each convolved with
#'    the disposition retained for the interval and fitted to the whole
#'    profile; BLQ samples before the first quantified one enter as a
#'    censoring penalty (M3-like). BIC weights averaged.
#' 5. Elimination: dose proportionality (power model on AUC, or on
#'    dose-normalised concentrations by time bin for sparse data) and
#'    curvature of the log-concentration decline (concave = saturable).
#' 6. Residual error: the best structural fit of each interval is refitted
#'    by iteratively reweighted least squares under a combined variance
#'    `a^2 + b^2 f^2`; on these common, leverage-corrected (PRESS) residuals,
#'    the variance models `a^2` (additive), `b^2 f^2` (proportional) and
#'    `a^2 + b^2 f^2` (combined) are fitted by maximum likelihood and compared
#'    by BIC, the combined model being kept only when its BIC is better by
#'    more than 15 (threshold chosen by simulation; 10 until version 0.9.0).
#'    Points of leverage above 0.7 and predictions below twice
#'    the LLOQ are left out. The estimates give the initial SIGMA values
#'    (variances).
#' 7. IIV: variance of the log of the individual estimates of CL, central
#'    volume and absorption parameter (median per subject), given as CV%.
#'    It includes the estimation noise of each individual value and tends to
#'    overestimate the true IIV, which is the safe side for initial values.
#' 8. With fewer than 3 rich intervals, all observations are pooled by time
#'    after dose (naive pooled, dose-normalised, binned) and steps 2 to 5
#'    are run on the binned geometric means (indicative only); the residual
#'    error model and the IIV cannot be separated in that case and default
#'    values are used.
#'
#' The candidate score is the product of the relative supports of the four
#' components; it sums to 1 over all combinations and is a ranking aid, not
#' a probability.
#'
#' Alerts (`$alerts`) flag features that none of the tested structures
#' describes: double peaks (enterohepatic recycling, several absorption
#' sites), a clearance that differs between first dose and steady state
#' (auto-induction, inhibition, time-varying CL), an accumulation different
#' from the one predicted by the first-dose terminal slope, and a
#' less-than-proportional exposure (saturable absorption or bioavailability).
#'
#' @param data Path to a NONMEM CSV file, or a data.frame.
#' @param cols Optional named list mapping roles to column names: `id`,
#'   `time`, `dv`, `amt`, `evid`, `mdv`, `cmt`, `rate`, `addl`, `ii`, `ss`,
#'   `blq`. Missing roles are guessed from standard names.
#' @param route `"auto"` (default: from the CMT of doses and observations,
#'   then from the shape of the profiles), `"EV"` or `"IV"`.
#' @param lloq Optional LLOQ (DV scale). Concentrations below it are BLQ.
#'   Without it, the lowest quantified concentration is used as LLOQ for the
#'   BLQ samples.
#' @param comment_char Characters marking ignored lines (first field).
#' @param sep Field separator of the CSV file.
#' @param dose_conc_factor Unit factor between amount/volume and the DV
#'   units, such that CL = factor * Dose / AUC: 1 (default) when DV is in
#'   dose units per volume unit (mg and mg/L), 1000 for doses in mg and DV in
#'   ng/mL (volumes then in L). It is passed to [nm_write_ctl()] by `from`,
#'   which writes the scale `S = V/factor`, so that the initial values match
#'   the control stream.
#' @param min_obs_rich Minimum number of quantified concentrations in a
#'   dosing interval for it to be analysed individually (default 6).
#' @param n_bins Number of time-after-dose bins for the pooled analysis.
#' @param max_cpts Maximum number of compartments tested (1 to 3).
#' @param n_candidates Number of candidate models returned.
#' @param max_profiles Maximum number of rich intervals analysed
#'   individually (a reproducible random subset is drawn above it).
#' @param design_eval If `TRUE` (default), the expected precision (RSE) and
#'   eta-shrinkage of each linear candidate are computed for the sampling
#'   design of the dataset with [nm_suggest_design()].
#' @param seed Seed of the subset of intervals; the user's random state is
#'   restored.
#' @param verbose Print progress messages.
#'
#' @return An object of class `nm_model_suggestion`, a list with:
#'   \describe{
#'     \item{`candidates`}{one row per candidate: `rank`, `cpts`, `route`,
#'       `input`, `elimination`, `error`, `score`, `parameters`,
#'       `theta_init`, `iiv_cv`, `sigma_init`, `max_rse_pct`,
#'       `max_shrinkage_pct`, `hard_to_estimate`, `notes`.}
#'     \item{`specs`}{one list per candidate with the arguments used by
#'       `nm_write_ctl(from = , candidate = )`.}
#'     \item{`calls`}{the matching [nm_write_ctl()] calls, as text.}
#'     \item{`support`}{relative support of each level of the four
#'       components.}
#'     \item{`design`}{expected RSE and shrinkage of each candidate
#'       (`nm_suggest_design` objects, `NULL` for Michaelis-Menten).}
#'     \item{`alerts`}{data.frame `type`, `message`.}
#'     \item{`nca`}{one row per analysed interval.}
#'     \item{`individual`}{individual estimates (CL, central volume,
#'       absorption parameter) used for the IIV and by
#'       [nm_suggest_covariates()].}
#'     \item{`evidence`, `profiles`, `data_summary`, `flags`, `notes`,
#'       `settings`}{details of the analysis.}
#'   }
#'
#' @seealso [nm_write_ctl()], [nm_suggest_design()], [nm_suggest_run()],
#'   [nm_suggest_covariates()], [plot.nm_model_suggestion()]
#'
#' @examples
#' ex <- system.file("extdata", "example", "example_data.csv", package = "nmpilot")
#' s <- nm_suggest_model(ex, verbose = FALSE, max_profiles = 15)
#' s
#' s$candidates[, c("rank", "cpts", "input", "elimination", "error", "score")]
#' nm_write_ctl(from = s, candidate = 1, verbose = FALSE)
#' @export
nm_suggest_model <- function(data,
                             cols = NULL,
                             route = c("auto", "EV", "IV"),
                             lloq = NULL,
                             comment_char = c("#", "@"),
                             sep = ",",
                             dose_conc_factor = 1,
                             min_obs_rich = 6L,
                             n_bins = 10L,
                             max_cpts = 3L,
                             n_candidates = 6L,
                             max_profiles = 60L,
                             design_eval = TRUE,
                             seed = 1L,
                             verbose = TRUE) {
  route <- match.arg(route)
  max_cpts <- as.integer(max_cpts)
  if (length(max_cpts) != 1L || !max_cpts %in% 1:3) stop("`max_cpts` must be 1, 2 or 3.")
  if (!is.numeric(dose_conc_factor) || length(dose_conc_factor) != 1L ||
      dose_conc_factor <= 0) stop("`dose_conc_factor` must be a positive number.")

  msg <- function(...) if (verbose) message(...)
  notes <- character(0)
  add_note <- function(x) notes <<- c(notes, x)
  data_file <- if (is.character(data) && length(data) == 1L) normalizePath(data, winslash = "/", mustWork = FALSE) else NULL

  # 1. Read and map ----------------------------------------------------------
  df  <- .nm_sm_read(data, sep = sep, comment_char = comment_char)
  map <- .nm_sm_map(names(df), cols)
  for (r in c("id", "time", "dv", "amt"))
    if (is.na(map[[r]])) stop(sprintf("Column for role '%s' not found; give it in `cols`.", r))
  msg(sprintf("Columns: %s",
              paste(sprintf("%s=%s", names(map)[!is.na(map)], unlist(map[!is.na(map)])),
                    collapse = ", ")))
  rec <- .nm_sm_records(df, map, lloq)
  if (!any(rec$is_dose)) stop("No dosing record found (EVID 1/4 or AMT > 0).")
  if (!any(rec$is_obs))  stop("No observation record found.")
  doses <- .nm_sm_expand_doses(rec)
  obs   <- rec[rec$is_obs, c("ID", "TIME", "DV", "BLQ", "CMT")]
  lloq_dv <- if (!is.null(lloq)) lloq else if (any(obs$BLQ)) min(obs$DV[!obs$BLQ & obs$DV > 0]) else NA
  if (is.null(lloq) && any(obs$BLQ))
    add_note(sprintf("No `lloq` given: the lowest quantified concentration (%s) is used as LLOQ for BLQ samples.",
                     format(signif(lloq_dv, 3))))

  # 2. Route -----------------------------------------------------------------
  rt <- .nm_sm_route(doses, obs, route)
  route_used <- rt$route
  infusion   <- rt$infusion
  if (!is.null(rt$note)) add_note(rt$note)
  msg(sprintf("Route: %s%s", route_used, if (infusion) " (infusion)" else ""))

  # 3. Dosing intervals --------------------------------------------------------
  prof  <- .nm_sm_profiles(doses, obs)
  pinfo <- prof$info
  pobs  <- prof$obs
  rich_all <- pinfo$pid[pinfo$n_quant >= min_obs_rich]
  rich_ids <- rich_all
  if (length(rich_ids) > max_profiles) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    set.seed(seed)
    rich_ids <- sort(sample(rich_ids, max_profiles))
    if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv())
    add_note(sprintf("%d rich intervals analysed (random subset of the %d available, `max_profiles`).",
                     max_profiles, length(rich_all)))
  }

  # 4. Analysis of each rich interval -------------------------------------------
  nca_rows <- list(); cpt_rows <- list(); abs_rows <- list(); curv_rows <- list()
  res_rows <- list(); grid_rows <- list(); obs_rows <- list(); ind_rows <- list()
  fit_rows <- list()
  for (pid in rich_ids) {
    pi <- pinfo[pinfo$pid == pid, ]
    po <- pobs[pobs$pid == pid, ]
    res <- .nm_sm_analyse_profile(po$TAD, po$DV, po$BLQ, dose = pi$AMT,
                                  type = pi$type, tau = pi$tau,
                                  tinf = if (infusion) pi$TINF else 0,
                                  route = route_used, max_cpts = max_cpts,
                                  factor = dose_conc_factor, lloq_dv = lloq_dv)
    if (is.null(res)) next
    nca_rows[[pid]]  <- cbind(data.frame(pid = pid, ID = pi$ID, dose = pi$AMT,
                                         type = pi$type, tau = pi$tau,
                                         stringsAsFactors = FALSE), res$nca)
    cpt_rows[[pid]]  <- cbind(data.frame(pid = pid, stringsAsFactors = FALSE), res$cpts)
    if (!is.null(res$abs)) abs_rows[[pid]] <- cbind(data.frame(pid = pid, stringsAsFactors = FALSE), res$abs)
    curv_rows[[pid]] <- data.frame(pid = pid, curv = res$curv$coef, p = res$curv$p)
    res_rows[[pid]]  <- cbind(pid = pid, res$residuals, stringsAsFactors = FALSE)
    grid_rows[[pid]] <- cbind(pid = pid, res$grid, stringsAsFactors = FALSE)
    obs_rows[[pid]]  <- cbind(pid = pid, res$obs, stringsAsFactors = FALSE)
    ind_rows[[pid]]  <- .nm_sm_individual(pid, pi$ID, res, route_used, dose_conc_factor)
    fit_rows[[pid]]  <- res$best_fit
  }
  nca  <- .nm_sm_bind(nca_rows)
  cpt  <- .nm_sm_bind(cpt_rows)
  absf <- .nm_sm_bind(abs_rows)
  curv <- .nm_sm_bind(curv_rows)
  resid <- .nm_sm_bind(res_rows)
  indiv <- .nm_sm_bind(ind_rows)
  n_rich <- if (is.null(nca)) 0L else nrow(nca)
  n_subj <- length(unique(obs$ID))
  subj_rich <- if (n_rich) length(unique(nca$ID)) else 0L
  msg(sprintf("%d subjects, %d rich interval(s) from %d subject(s).", n_subj, n_rich, subj_rich))

  # 5. Pooled analysis -------------------------------------------------------------
  pooled <- .nm_sm_pooled(pinfo, pobs, n_bins = n_bins, route = route_used,
                          infusion = infusion, max_cpts = max_cpts,
                          factor = dose_conc_factor)
  source_ev <- if (n_rich >= 3L) "individual" else if (!is.null(pooled$res)) "pooled" else "none"
  if (source_ev == "none")
    stop("Not enough data: fewer than 3 rich intervals and too few observations to pool.")
  if (source_ev == "pooled") {
    add_note(sprintf("Fewer than 3 rich intervals: evidence comes from naive pooled data (%s); indicative only.",
                     pooled$basis))
    msg("Evidence from naive pooled data.")
  }
  design <- if (subj_rich >= 0.5 * n_subj) "rich" else if (subj_rich == 0) "sparse" else "mixed"

  # 6. Compartments ----------------------------------------------------------------
  wc <- if (source_ev == "individual")
    .nm_sm_weights_table(cpt, prefix = "bic_cpt", levels = 1:max_cpts) else
    .nm_sm_weights_table(cbind(pid = "pooled", pooled$res$cpts), prefix = "bic_cpt", levels = 1:max_cpts)
  p_cpt <- stats::setNames(wc$mean_w, as.character(1:max_cpts))

  # 7. Input ---------------------------------------------------------------------------
  abs_models <- c("first-order", "first-order+lag", "zero-order-d1", "transit")
  wa <- NULL
  if (route_used == "IV") {
    p_abs <- c(IV = 1)
  } else {
    a_src <- if (source_ev == "individual") absf else
      if (!is.null(pooled$res$abs)) cbind(pid = "pooled", pooled$res$abs) else NULL
    if (!is.null(a_src) && nrow(a_src) > 0 &&
        any(is.finite(as.matrix(a_src[, paste0("bic_", .nm_sm_absid(abs_models)), drop = FALSE])))) {
      wa <- .nm_sm_weights_table(a_src, prefix = "bic_", levels = .nm_sm_absid(abs_models))
      p_abs <- wa$mean_w
    } else {
      p_abs <- c(0.55, 0.15, 0.15, 0.15)
      add_note("Absorption phase not characterisable (fewer than 2 points before the peak): first-order absorption chosen by parsimony.")
    }
    names(p_abs) <- abs_models
  }

  # 8. Elimination -----------------------------------------------------------------
  el <- .nm_sm_elimination(nca, curv, pinfo, pobs, n_bins = n_bins)
  p_el <- c(linear = 1 - el$p_mm, `michaelis-menten` = el$p_mm)
  if (!is.null(el$note)) add_note(el$note)

  # 9. Residual error ----------------------------------------------------------------
  err <- .nm_sm_error_models(Filter(Negate(is.null), fit_rows), source_ev, obs, lloq_dv)
  p_err <- err$support
  if (!is.null(err$note)) add_note(err$note)

  # 10. Initial values -----------------------------------------------------------------
  iv <- .nm_sm_initial_values(nca, cpt, absf, pooled, source_ev, route_used,
                              dose_conc_factor, obs, max_cpts)
  if (!is.null(iv$note)) add_note(iv$note)
  iiv <- .nm_sm_iiv(indiv, source_ev)
  if (!is.null(iiv$note)) add_note(iiv$note)

  # 11. Flags and alerts -------------------------------------------------------------
  flags <- .nm_sm_flags(nca, cpt, absf, pooled, source_ev, route_used, el, obs, n_subj)
  alerts <- .nm_sm_alerts(nca, obs_rows, el)

  # 12. Candidates -----------------------------------------------------------------------
  grid <- expand.grid(cpts = 1:max_cpts, input = names(p_abs), elimination = names(p_el),
                      error = names(p_err), stringsAsFactors = FALSE)
  grid$score <- p_cpt[as.character(grid$cpts)] * p_abs[grid$input] *
    p_el[grid$elimination] * p_err[grid$error]
  grid$score <- grid$score / sum(grid$score)
  grid <- utils::head(grid[order(-grid$score), ], n_candidates)
  rownames(grid) <- NULL

  design_data <- list(doses = doses[, c("ID", "TIME", "AMT", "TINF")],
                      obs = obs[!obs$BLQ, c("ID", "TIME")], factor = dose_conc_factor)
  specs <- list(); des <- list()
  cand <- data.frame(rank = seq_len(nrow(grid)), cpts = grid$cpts, route = route_used,
                     input = if (route_used == "IV") ifelse(infusion, "IV infusion", "IV bolus") else grid$input,
                     elimination = grid$elimination, error = grid$error,
                     score = round(grid$score, 3), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(grid))) {
    sp <- .nm_sm_spec(grid$cpts[i], route_used, grid$input[i], grid$elimination[i],
                      grid$error[i], iv, iiv, err)
    specs[[i]] <- sp
    des[i] <- list(NULL)
    if (design_eval && sp$elimination == "linear") {
      d <- tryCatch(.nm_sd_evaluate(sp, design_data), error = function(e) NULL)
      if (!is.null(d)) des[[i]] <- d
    }
    hd <- .nm_sm_hard_for(sp, flags, des[[i]])
    cand$parameters[i] <- paste(names(sp$theta), collapse = ", ")
    cand$theta_init[i] <- paste(sprintf("%s=%s", names(sp$theta), .nm_sm_fmt(sp$theta)), collapse = ", ")
    cand$iiv_cv[i] <- if (length(sp$iiv)) paste(sprintf("%s=%s", names(sp$iiv), .nm_sm_fmt(sp$iiv)), collapse = ", ") else ""
    cand$sigma_init[i] <- paste(sprintf("%s=%s", names(sp$sigma_used), .nm_sm_fmt(sp$sigma_used)), collapse = ", ")
    cand$max_rse_pct[i] <- if (!is.null(des[[i]])) round(max(des[[i]]$parameters$rse_pct, na.rm = TRUE), 1) else NA
    cand$max_shrinkage_pct[i] <- if (!is.null(des[[i]]) && nrow(des[[i]]$shrinkage))
      round(max(des[[i]]$shrinkage$shrinkage_pct, na.rm = TRUE), 1) else NA
    cand$hard_to_estimate[i] <- if (length(hd)) paste(sprintf("%s (%s)", names(hd), hd), collapse = "; ") else ""
    cand$notes[i] <- sp$note
  }
  calls <- vapply(specs, .nm_sm_call, "")
  names(specs) <- names(calls) <- names(des) <- paste0("candidate_", seq_len(nrow(grid)))

  out <- list(
    candidates = cand,
    specs = specs,
    calls = calls,
    support = list(compartments = p_cpt, input = p_abs, elimination = p_el, error = p_err),
    design = des,
    alerts = alerts,
    nca = nca,
    individual = indiv,
    evidence = list(source = source_ev, compartments = wc$table,
                    input = if (!is.null(wa)) wa$table else NULL,
                    elimination = el, error = err$table, curvature = curv, pooled = pooled),
    profiles = list(obs = .nm_sm_bind(obs_rows), grid = .nm_sm_bind(grid_rows), residuals = resid,
                    best = if (n_rich) stats::setNames(nca$best_fit, nca$pid) else NULL),
    data_summary = list(n_subjects = n_subj, n_obs = sum(!obs$BLQ), n_blq = sum(obs$BLQ),
                        n_profiles = nrow(pinfo), n_rich_profiles = n_rich,
                        n_subjects_rich = subj_rich, design = design,
                        route = route_used, infusion = infusion,
                        dose_levels = sort(unique(signif(doses$AMT, 4))),
                        obs_per_subject = stats::median(table(obs$ID[!obs$BLQ])),
                        data_file = data_file, columns = attr(df, "orig_names")),
    design_data = design_data,
    flags = flags,
    notes = unique(notes),
    settings = list(dose_conc_factor = dose_conc_factor, min_obs_rich = min_obs_rich,
                    n_bins = n_bins, max_cpts = max_cpts, lloq = lloq_dv)
  )
  class(out) <- "nm_model_suggestion"
  out
}


#' @rdname nm_suggest_model
#' @param x An `nm_model_suggestion` object.
#' @param ... Ignored.
#' @export
print.nm_model_suggestion <- function(x, ...) {
  ds <- x$data_summary
  sup <- function(v) paste(sprintf("%s = %.2f", names(v), v), collapse = ", ")
  cat("Structural model suggestion\n")
  cat("===========================\n")
  cat(sprintf("Subjects: %d | quantified obs: %d | BLQ: %d | median obs/subject: %s\n",
              ds$n_subjects, ds$n_obs, ds$n_blq, format(ds$obs_per_subject)))
  cat(sprintf("Design: %s (%d rich intervals from %d subjects) | evidence: %s\n",
              ds$design, ds$n_rich_profiles, ds$n_subjects_rich, x$evidence$source))
  cat(sprintf("Route: %s%s | dose levels: %s\n", ds$route, if (ds$infusion) " (infusion)" else "",
              paste(utils::head(ds$dose_levels, 8), collapse = ", ")))
  cat("\nRelative support\n")
  cat("  Compartments: ", sup(x$support$compartments), "\n")
  cat("  Input:        ", sup(x$support$input), "\n")
  cat("  Elimination:  ", sup(x$support$elimination), "\n")
  cat("  Error model:  ", sup(x$support$error), "\n")
  cat("\nCandidate models\n")
  cc <- x$candidates
  for (i in seq_len(nrow(cc))) {
    cat(sprintf("%d. %d cpt | %s | %s | %s error | score %.3f\n", cc$rank[i], cc$cpts[i],
                cc$input[i], cc$elimination[i], cc$error[i], cc$score[i]))
    cat(sprintf("   THETA: %s\n", cc$theta_init[i]))
    if (nzchar(cc$iiv_cv[i])) cat(sprintf("   IIV (CV%%): %s | SIGMA: %s\n", cc$iiv_cv[i], cc$sigma_init[i]))
    if (is.finite(cc$max_rse_pct[i]))
      cat(sprintf("   expected precision: max RSE %.0f%%, max eta-shrinkage %s%%\n", cc$max_rse_pct[i],
                  format(cc$max_shrinkage_pct[i])))
    if (nzchar(cc$hard_to_estimate[i])) cat(sprintf("   hard to estimate: %s\n", cc$hard_to_estimate[i]))
    if (nzchar(cc$notes[i])) cat(sprintf("   note: %s\n", cc$notes[i]))
  }
  if (!is.null(x$alerts) && nrow(x$alerts)) {
    cat("\nAlerts\n")
    for (i in seq_len(nrow(x$alerts))) cat(sprintf(" ! [%s] %s\n", x$alerts$type[i], x$alerts$message[i]))
  }
  if (length(x$notes)) {
    cat("\nNotes\n")
    for (n in x$notes) cat(" -", n, "\n")
  }
  cat("\nControl stream of candidate 1: nm_write_ctl(from = x, candidate = 1)\n")
  invisible(x)
}


#' Plot a model suggestion
#'
#' @param x An `nm_model_suggestion` object.
#' @param type `"profiles"` (observed concentrations of the analysed
#'   intervals with the best fit of each, log scale), `"support"` (relative
#'   support of each component), `"pooled"` (binned pooled data with their
#'   fit) or `"design"` (expected RSE and eta-shrinkage of a candidate).
#' @param candidate Candidate shown by `type = "design"`.
#' @param n_max Maximum number of intervals shown by `type = "profiles"`.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.nm_model_suggestion <- function(x, type = c("profiles", "support", "pooled", "design"),
                                     candidate = 1L, n_max = 12L, ...) {
  type <- match.arg(type)
  if (type == "profiles") {
    o <- x$profiles$obs; g <- x$profiles$grid
    if (is.null(o)) stop("No rich interval analysed individually: use type = \"pooled\".")
    ids <- utils::head(unique(o$pid), n_max)
    o <- o[o$pid %in% ids & o$time > 0, ]; g <- g[g$pid %in% ids & g$pred > 0, ]
    lab <- x$profiles$best[ids]
    lab <- stats::setNames(sprintf("%s (%s)", ids, lab), ids)
    o$panel <- lab[o$pid]; g$panel <- lab[g$pid]
    lq <- x$settings$lloq
    o$dv_plot <- ifelse(o$blq, if (is.finite(lq)) lq / 2 else NA, o$dv)
    return(ggplot2::ggplot() +
             ggplot2::geom_line(data = g, ggplot2::aes(x = .data$time, y = .data$pred), colour = "#2166AC") +
             ggplot2::geom_point(data = o[!o$blq, ], ggplot2::aes(x = .data$time, y = .data$dv_plot)) +
             ggplot2::geom_point(data = o[o$blq & is.finite(o$dv_plot), ],
                                 ggplot2::aes(x = .data$time, y = .data$dv_plot), shape = 4, colour = "grey50") +
             ggplot2::scale_y_log10() +
             ggplot2::facet_wrap(~panel, scales = "free") +
             ggplot2::labs(x = "Time after dose", y = "Concentration (log scale)",
                           title = "Analysed intervals and best fit",
                           subtitle = "Crosses: BLQ samples, drawn at LLOQ/2") +
             ggplot2::theme_bw())
  }
  if (type == "support") {
    s <- x$support
    d <- do.call(rbind, lapply(names(s), function(n)
      data.frame(component = n, level = names(s[[n]]), support = as.numeric(s[[n]]))))
    d$component <- factor(d$component, levels = names(s))
    return(ggplot2::ggplot(d, ggplot2::aes(x = .data$level, y = .data$support)) +
             ggplot2::geom_col(fill = "#4393C3") +
             ggplot2::facet_wrap(~component, scales = "free_x") +
             ggplot2::coord_cartesian(ylim = c(0, 1)) +
             ggplot2::labs(x = NULL, y = "Relative support", title = "Support of each model component") +
             ggplot2::theme_bw() +
             ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)))
  }
  if (type == "pooled") {
    p <- x$evidence$pooled
    if (is.null(p$res)) stop("No pooled analysis available (too few observations).")
    b <- p$bins; g <- p$res$grid
    g <- g[g$pred > 0, ]
    return(ggplot2::ggplot() +
             ggplot2::geom_line(data = g, ggplot2::aes(x = .data$time, y = .data$pred), colour = "#2166AC") +
             ggplot2::geom_point(data = b, ggplot2::aes(x = .data$tad, y = .data$gm, size = .data$n)) +
             ggplot2::scale_y_log10() +
             ggplot2::labs(x = "Time after dose", y = "Dose-normalised concentration (geometric mean, log)",
                           size = "n obs", title = sprintf("Naive pooled data (%s)", p$basis)) +
             ggplot2::theme_bw())
  }
  d <- x$design[[candidate]]
  if (is.null(d)) stop("No design evaluation for this candidate (Michaelis-Menten or failure).")
  plot(d)
}


# ---------------------------------------------------------------------------
# Reading and mapping
# ---------------------------------------------------------------------------
.nm_sm_read <- function(data, sep, comment_char) {
  if (is.data.frame(data)) {
    df <- as.data.frame(data, stringsAsFactors = FALSE)
    for (n in names(df)) if (is.character(df[[n]]) || is.factor(df[[n]])) {
      v <- trimws(as.character(df[[n]])); v[v %in% c(".", "", "NA")] <- NA; df[[n]] <- v
    }
  } else {
    if (!file.exists(data)) stop("File not found: ", data)
    lines <- readLines(data, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    header <- lines[1]
    for (cc in comment_char) header <- sub(paste0("^\\s*\\", cc), "", header)
    body <- lines[-1]
    first <- substr(trimws(body), 1, 1)
    body <- body[!first %in% comment_char]
    df <- utils::read.table(text = c(header, body), sep = sep, header = TRUE,
                            colClasses = "character", na.strings = c(".", "", "NA"),
                            check.names = FALSE, strip.white = TRUE, quote = "\"")
  }
  orig <- trimws(names(df))
  names(df) <- toupper(orig)
  attr(df, "orig_names") <- orig
  for (n in names(df)) {
    v <- suppressWarnings(as.numeric(as.character(df[[n]])))
    if (all(is.na(v) == is.na(df[[n]]))) df[[n]] <- v
  }
  df
}

.nm_sm_map <- function(nms, cols) {
  guess <- list(id = c("ID"), time = c("TIME", "TIM"), dv = c("DV", "CONC", "LNDV"),
                amt = c("AMT", "DOSE"), evid = "EVID", mdv = "MDV", cmt = "CMT",
                rate = "RATE", addl = "ADDL", ii = "II", ss = "SS",
                blq = c("BLQ", "BQL", "CENS", "LLQ"))
  map <- lapply(names(guess), function(r) {
    if (!is.null(cols[[r]])) {
      v <- toupper(cols[[r]])
      if (!v %in% nms) stop(sprintf("Column '%s' (role '%s') not found.", cols[[r]], r))
      return(v)
    }
    hit <- intersect(guess[[r]], nms)
    if (length(hit)) hit[1] else NA_character_
  })
  names(map) <- names(guess)
  if (!is.na(map$dv) && map$dv == "LNDV")
    warning("DV column 'LNDV' found: log-transformed DV is not supported; give an untransformed DV in `cols`.")
  map
}

.nm_sm_records <- function(df, map, lloq) {
  g <- function(r, default = NA) if (is.na(map[[r]])) rep(default, nrow(df)) else df[[map[[r]]]]
  rec <- data.frame(ID = g("id"), TIME = g("time"), DV = g("dv"), AMT = g("amt", 0),
                    EVID = g("evid"), MDV = g("mdv"), CMT = g("cmt"), RATE = g("rate", 0),
                    ADDL = g("addl", 0), II = g("ii", 0), SS = g("ss", 0), BLQF = g("blq", 0))
  rec$AMT[is.na(rec$AMT)] <- 0
  for (v in c("RATE", "ADDL", "II", "SS", "BLQF")) rec[[v]][is.na(rec[[v]])] <- 0
  if (any(is.na(rec$ID)) || any(is.na(rec$TIME)))
    stop("Missing ID or TIME values: TIME must be numeric (relative time).")
  if (!is.na(map$evid)) {
    rec$is_dose <- rec$EVID %in% c(1, 4)
    rec$is_obs  <- rec$EVID %in% 0 & !is.na(rec$DV)
    if (!is.na(map$mdv)) rec$is_obs <- rec$is_obs & rec$MDV %in% 0
    rec$reset <- rec$EVID %in% c(3, 4)
  } else {
    rec$is_dose <- rec$AMT > 0
    rec$is_obs  <- !rec$is_dose & !is.na(rec$DV)
    if (!is.na(map$mdv)) rec$is_obs <- rec$is_obs & rec$MDV %in% 0
    rec$reset <- FALSE
  }
  rec$BLQ <- rec$is_obs & (rec$BLQF == 1 | is.na(rec$DV) | rec$DV <= 0 |
                             (if (!is.null(lloq)) rec$DV < lloq else FALSE))
  rec <- rec[order(rec$ID, rec$TIME, !rec$is_dose), ]
  rec
}

.nm_sm_expand_doses <- function(rec) {
  d <- rec[rec$is_dose & rec$AMT > 0, c("ID", "TIME", "AMT", "CMT", "RATE", "ADDL", "II", "SS", "reset")]
  if (any(d$ADDL > 0 & d$II > 0)) {
    extra <- do.call(rbind, lapply(which(d$ADDL > 0 & d$II > 0), function(i) {
      k <- seq_len(d$ADDL[i])
      e <- d[rep(i, length(k)), ]
      e$TIME <- d$TIME[i] + k * d$II[i]; e$ADDL <- 0; e$SS <- 0; e$reset <- FALSE
      e
    }))
    d <- rbind(d, extra)
  }
  d$TINF <- ifelse(d$RATE > 0, d$AMT / d$RATE, 0)
  d <- d[order(d$ID, d$TIME), ]
  rownames(d) <- NULL
  d
}

# ---------------------------------------------------------------------------
# Route
# ---------------------------------------------------------------------------
.nm_sm_route <- function(doses, obs, route) {
  infusion <- any(doses$RATE > 0)
  note <- NULL
  if (route != "auto") return(list(route = route, infusion = infusion && route == "IV", note = NULL))
  if (any(doses$RATE < 0)) {
    note <- "RATE = -1/-2 found on doses: modelled zero-order input; analysed as extravascular input."
    return(list(route = "EV", infusion = FALSE, note = note))
  }
  if (!all(is.na(doses$CMT)) && !all(is.na(obs$CMT))) {
    obs_cmt <- as.numeric(names(which.max(table(obs$CMT))))
    ev <- doses$CMT != obs_cmt
    if (any(ev) && any(!ev))
      note <- sprintf("Doses in several compartments (%d extravascular, %d intravascular records): the main route is analysed.",
                      sum(ev, na.rm = TRUE), sum(!ev, na.rm = TRUE))
    r <- if (mean(ev, na.rm = TRUE) >= 0.5) "EV" else "IV"
    return(list(route = r, infusion = infusion && r == "IV", note = note))
  }
  if (infusion) return(list(route = "IV", infusion = TRUE,
                            note = "No CMT column: RATE > 0 on doses taken as IV infusion."))
  # shape: peak at first sample after dose -> IV bolus
  first_peak <- vapply(split(seq_len(nrow(doses)), doses$ID), function(ii) {
    id <- doses$ID[ii[1]]; t0 <- doses$TIME[ii[1]]
    t1 <- if (length(ii) > 1) doses$TIME[ii[2]] else Inf
    o <- obs[obs$ID == id & obs$TIME > t0 & obs$TIME < t1 & !obs$BLQ, ]
    if (nrow(o) < 3) return(NA)
    which.max(o$DV) == 1
  }, logical(1))
  r <- if (mean(first_peak, na.rm = TRUE) >= 0.8) "IV" else "EV"
  list(route = r, infusion = FALSE,
       note = sprintf("No CMT column: route guessed as %s from the profile shapes; give `route` to override.", r))
}

# ---------------------------------------------------------------------------
# Profiles (dosing intervals)
# ---------------------------------------------------------------------------
.nm_sm_profiles <- function(doses, obs) {
  info <- list(); po <- list()
  for (id in unique(doses$ID)) {
    d <- doses[doses$ID == id, ]
    o <- obs[obs$ID == id, ]
    nd <- nrow(d)
    reset_idx <- cumsum(d$reset | seq_len(nd) == 1)
    for (j in seq_len(nd)) {
      t0 <- d$TIME[j]
      t1 <- if (j < nd) d$TIME[j + 1] else Inf
      # a pre-dose sample recorded at the time of the next dose closes the interval
      sel <- o$TIME > t0 & (o$TIME < t1 | (j < nd & o$TIME == t1))
      oj <- o[sel, ]
      if (!nrow(oj)) next
      prev <- which(reset_idx == reset_idx[j] & seq_len(nd) < j)
      dose_rank <- length(prev) + 1L
      regular <- dose_rank >= 5 && {
        tt <- d$TIME[c(prev, j)]; ii <- diff(tt); aa <- d$AMT[c(prev, j)]
        length(ii) >= 4 && stats::sd(utils::tail(ii, 4)) / mean(utils::tail(ii, 4)) < 0.1 &&
          length(unique(utils::tail(aa, 5))) == 1
      }
      type <- if (d$SS[j] > 0 || regular) "ss" else if (dose_rank == 1) "single" else "other"
      tau <- if (d$SS[j] > 0 && d$II[j] > 0) d$II[j] else if (j < nd) t1 - t0 else
        if (dose_rank > 1) t0 - d$TIME[j - 1] else Inf
      if (type == "single" && j < nd && is.finite(t1)) type <- "single_interval"
      pid <- sprintf("%s_%d", id, j)
      info[[pid]] <- data.frame(pid = pid, ID = id, dose_index = j, AMT = d$AMT[j],
                                TINF = d$TINF[j], type = type, tau = tau,
                                n_quant = sum(!oj$BLQ), n_obs = nrow(oj),
                                stringsAsFactors = FALSE)
      po[[pid]] <- data.frame(pid = pid, ID = id, TAD = oj$TIME - t0, DV = oj$DV,
                              BLQ = oj$BLQ, AMT = d$AMT[j], stringsAsFactors = FALSE)
    }
  }
  list(info = .nm_sm_bind(info), obs = .nm_sm_bind(po))
}

# ---------------------------------------------------------------------------
# Analysis of one profile: NCA + exponential fits + absorption fits
# ---------------------------------------------------------------------------
.nm_sm_analyse_profile <- function(tad, dv, blq, dose, type, tau, tinf, route,
                                   max_cpts, factor, lloq_dv = NA) {
  ok <- !blq & is.finite(dv) & dv > 0 & tad > 0
  t <- tad[ok]; c <- dv[ok]
  if (length(t) < 4) return(NULL)
  o <- order(t); t <- t[o]; c <- c[o]
  y <- c / dose
  imax <- which.max(c)
  tmax <- t[imax]
  n_pre <- sum(t < tmax)
  if (route == "EV" && (length(t) - imax) < 3) return(NULL)
  # BLQ samples before the first quantified one (lag evidence)
  first_q <- min(t)
  n_blq_early <- sum(blq & tad > 0 & tad < first_q)

  # --- NCA ---
  nca <- .nm_sm_nca(t, c, dose, type, tau, tinf, route, imax, factor)

  # --- exponential fits on the post-peak phase ---
  if (route == "IV" && tinf == 0) {
    sel <- seq_along(t)
    tt <- t; yy <- y
  } else if (route == "IV") {
    sel <- which(t > tinf + 1e-9)
    tt <- t[sel] - tinf; yy <- y[sel]
  } else {
    sel <- which(t > tmax)
    tt <- t[sel]; yy <- y[sel]
  }
  cpts <- .nm_sm_fit_exps(tt, log(yy), max_cpts, lz = nca$lambda_z)
  kb_exp <- which.min(ifelse(is.finite(cpts$bic), cpts$bic, Inf))
  if (!length(kb_exp) || !is.finite(cpts$bic[kb_exp])) kb_exp <- 1L
  exp_pred <- function(times) {
    p <- cpts$par[[kb_exp]]
    if (is.null(p)) return(rep(NA_real_, length(times)))
    log(colSums(p$A * exp(-outer(p$lam, times))))
  }
  # predictions of the multi-exponential fit (before any amplitude conversion)
  shift <- if (route == "IV" && tinf > 0) tinf else 0
  fit_exp <- list(time = t[sel], pred = exp(exp_pred(tt)) * dose)
  tg_e <- seq(min(tt), max(tt), length.out = 150)
  grid_exp <- data.frame(time = tg_e + shift, pred = exp(exp_pred(tg_e)) * dose)
  # convert post-infusion amplitudes to bolus amplitudes
  if (route == "IV" && tinf > 0) cpts$par <- lapply(cpts$par, function(p) {
    if (is.null(p)) return(p)
    p$A_raw <- p$A
    p$A <- p$A * p$lam * tinf / (1 - exp(-p$lam * tinf)); p
  })
  # --- curvature of the post-peak log decline (MM vs multi-exponential) ---
  # extravascular: points after 2 x tmax only, the end of absorption also
  # gives a concave shape just after the peak
  cs <- if (route == "EV") tt > 2 * tmax else rep(TRUE, length(tt))
  curv <- .nm_sm_curvature(tt[cs], log(yy[cs]))

  # --- absorption fits ---
  # full-profile fits: each input model convolved with the disposition
  # (sum of exponentials) retained for this profile; early BLQ samples enter
  # as a censoring penalty (M3-like) when an LLOQ is known
  ab <- NULL
  if (route == "EV" && n_pre >= 2) {
    kb <- which.min(ifelse(is.finite(cpts$bic), cpts$bic, Inf))
    if (!length(kb) || !is.finite(cpts$bic[kb])) kb <- 1L
    blq_t <- tad[blq & tad > 0 & tad < first_q]
    ab <- .nm_sm_fit_absorption(t, log(y), k = kb, disp = cpts$par[[kb]],
                                tmax = tmax, cmax = y[imax], lz = nca$lambda_z,
                                t_first = first_q, blq_t = blq_t,
                                lloq_norm = if (is.finite(lloq_dv)) lloq_dv / dose else NA)
  }

  cpt_df <- data.frame(base::t(stats::setNames(cpts$bic, paste0("bic_cpt", seq_len(max_cpts)))),
                       n_post = length(tt))
  for (k in seq_len(max_cpts)) {
    p <- cpts$par[[k]]
    if (!is.null(p) && !isTRUE(p$valid)) p <- NULL   # non-identifiable fit: not used for initial values
    for (i in seq_len(k)) {
      cpt_df[[sprintf("A%d_%d", k, i)]]   <- if (is.null(p)) NA else p$A[i]
      cpt_df[[sprintf("lam%d_%d", k, i)]] <- if (is.null(p)) NA else p$lam[i]
    }
  }
  nca$n_pre_tmax <- n_pre
  nca$n_blq_early <- n_blq_early
  nca$t_first <- first_q
  # best structural fit of the profile: residuals (error models) and plots
  best <- paste0(kb_exp, "-exp")
  res_df <- data.frame(time = fit_exp$time, dv = c[sel], pred = fit_exp$pred)
  grid <- grid_exp
  pe <- cpts$par[[kb_exp]]
  best_fit <- NULL
  if (!is.null(pe)) {
    k <- length(pe$lam)
    p0 <- c(log(pe$A_raw %||% pe$A), log(pe$lam[k]))
    if (k > 1) for (i in 1:(k - 1)) p0 <- c(p0, log(pe$lam[i] / pe$lam[i + 1] - 1))
    best_fit <- list(fun = .nm_sm_exp_fun(k), par = p0, t = tt, y = yy, dose = dose)
  }
  if (!is.null(ab)) {
    b <- c(ab$bic_fo, ab$bic_lag, ab$bic_zo, ab$bic_tr)
    if (any(is.finite(b))) {
      j <- which.min(ifelse(is.finite(b), b, Inf))
      best <- c("first-order", "first-order+lag", "zero-order-d1", "transit")[j]
      res_df <- data.frame(time = t, dv = c, pred = exp(attr(ab, "logpred")[, j]) * dose)
      g <- attr(ab, "grid")
      grid <- data.frame(time = g$time, pred = exp(g[[j + 1]]) * dose)
      fb <- attr(ab, "fits")[[j]]
      if (!is.null(fb)) best_fit <- list(fun = fb$fun, par = fb$par, t = t, y = y, dose = dose)
    }
  }
  nca$best_fit <- best
  list(nca = nca, cpts = cpt_df, abs = ab, curv = curv, residuals = res_df, grid = grid,
       obs = data.frame(time = tad, dv = dv, blq = blq), best_fit = best_fit)
}

.nm_sm_nca <- function(t, c, dose, type, tau, tinf, route, imax, factor) {
  n <- length(t)
  # terminal slope: best adjusted R2 on the last k points (k >= 3), after the peak
  start_min <- if (route == "IV" && tinf == 0) 1 else imax + 1
  lz <- NA; r2 <- NA; nlz <- NA; tlz <- NA
  if (n - start_min + 1 >= 3) {
    best <- -Inf
    for (s in seq(n - 2, start_min)) {
      tt <- t[s:n]; lc <- log(c[s:n])
      f <- stats::lm(lc ~ tt)
      b <- stats::coef(f)[2]
      if (!is.finite(b) || b >= 0) next
      k <- length(tt)
      r2a <- 1 - (1 - summary(f)$r.squared) * (k - 1) / (k - 2)
      if (!is.finite(r2a)) r2a <- 1
      if (r2a > best + 1e-4 || (abs(r2a - best) <= 1e-4 && k > nlz)) {
        best <- r2a; lz <- -b; r2 <- r2a; nlz <- k; tlz <- t[s]
      }
    }
  }
  # AUC (linear-up log-down)
  if (route == "IV" && tinf == 0) {
    c0 <- if (n >= 2 && c[2] < c[1]) c[1] * exp(log(c[1] / c[2]) * t[1] / (t[2] - t[1])) else c[1]
  } else if (type == "ss") {
    c0 <- c[n]
  } else c0 <- 0
  tt <- c(0, t); cc <- c(c0, c)
  auc <- 0; aumc <- 0
  for (i in 2:length(tt)) {
    dt <- tt[i] - tt[i - 1]
    if (cc[i] < cc[i - 1] && cc[i] > 0 && cc[i - 1] > 0) {
      k <- log(cc[i - 1] / cc[i]) / dt
      auc  <- auc + (cc[i - 1] - cc[i]) / k
      aumc <- aumc + (cc[i - 1] * tt[i - 1] - cc[i] * tt[i]) / k + (cc[i - 1] - cc[i]) / k^2
    } else {
      auc  <- auc + dt * (cc[i] + cc[i - 1]) / 2
      aumc <- aumc + dt * (cc[i] * tt[i] + cc[i - 1] * tt[i - 1]) / 2
    }
  }
  tlast <- t[n]; clast <- c[n]
  auc_inf <- if (is.finite(lz)) auc + clast / lz else NA
  aumc_inf <- if (is.finite(lz)) aumc + clast * tlast / lz + clast / lz^2 else NA
  pct_extrap <- if (is.finite(auc_inf)) 100 * (auc_inf - auc) / auc_inf else NA
  cl <- NA; auc_tau <- NA
  if (type %in% c("single", "single_interval") && is.finite(auc_inf)) {
    cl <- factor * dose / auc_inf
  } else if (type == "ss" && is.finite(tau)) {
    auc_tau <- auc
    if (tlast < tau && is.finite(lz)) auc_tau <- auc + clast / lz * (1 - exp(-lz * (tau - tlast)))
    cl <- factor * dose / auc_tau
  }
  auc_0tau <- NA
  if (type != "ss" && is.finite(tau)) {
    auc_0tau <- if (tlast < tau && is.finite(lz)) auc + clast / lz * (1 - exp(-lz * (tau - tlast))) else auc
  }
  mrt <- if (is.finite(aumc_inf) && type != "ss") aumc_inf / auc_inf else NA
  thalf <- if (is.finite(lz)) log(2) / lz else NA
  span <- if (is.finite(lz)) (tlast - tlz) / thalf else NA
  data.frame(cmax = c[imax], tmax = t[imax], tlast = tlast, n_quant = n,
             lambda_z = lz, r2_adj = r2, n_lz = nlz, t_half = thalf, span_hl = span,
             auc_last = auc, auc_inf = auc_inf, pct_extrap = pct_extrap, auc_tau = auc_tau, auc_0tau = auc_0tau,
             cl = cl, vz = if (is.finite(cl) && is.finite(lz)) cl / lz else NA,
             mrt = mrt, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Optimisation helper (least squares on log concentrations)
# ---------------------------------------------------------------------------
.nm_sm_optim <- function(fun, starts, t, y) {
  obj <- function(p) {
    pr <- fun(p, t)
    if (any(!is.finite(pr))) return(1e10)
    sum((y - pr)^2)
  }
  .nm_sm_multistart(obj, starts)
}

# coarse Nelder-Mead from every start, then refinement of the best one
.nm_sm_multistart <- function(obj, starts) {
  best <- NULL
  for (s in starts) {
    v0 <- obj(s)
    if (!is.finite(v0) || v0 >= 1e10) next
    f <- tryCatch(stats::optim(s, obj, method = "Nelder-Mead",
                               control = list(maxit = 800, reltol = 1e-8)),
                  error = function(e) NULL)
    if (!is.null(f) && f$value < 1e9 && (is.null(best) || f$value < best$value)) best <- f
  }
  if (is.null(best)) return(NULL)
  f <- tryCatch(stats::optim(best$par, obj, method = "Nelder-Mead",
                             control = list(maxit = 6000, reltol = 1e-11)),
                error = function(e) NULL)
  if (!is.null(f) && f$value <= best$value) best <- f
  best
}

.nm_sm_bic <- function(sse, n, k) n * log(max(sse, 1e-12) / n) + k * log(n)

.nm_sm_lse <- function(m) {  # row-wise log-sum-exp
  mx <- apply(m, 1, max)
  mx + log(rowSums(exp(m - mx)))
}

# Sums of exponentials: log C = log(sum A_i exp(-lam_i t)), lam_1 > lam_2 > ...
.nm_sm_exp_fun <- function(k) {
  function(p, t) {
    a <- p[1:k]
    lam <- numeric(k); lam[k] <- exp(p[k + 1])
    if (k > 1) for (i in (k - 1):1) lam[i] <- lam[i + 1] * (1 + exp(p[k + 1 + i]))
    m <- sapply(seq_len(k), function(i) a[i] - lam[i] * t)
    if (is.null(dim(m))) m <- matrix(m, nrow = length(t))
    .nm_sm_lse(m)
  }
}

.nm_sm_fit_exps <- function(t, y, max_cpts, lz = NA) {
  n <- length(t)
  min_n <- c(3, 6, 9)
  bic <- rep(NA_real_, max_cpts); par <- vector("list", max_cpts)
  if (!is.finite(lz) || lz <= 0) {
    f <- stats::lm(y ~ t); lz <- max(-stats::coef(f)[2], 1e-3, na.rm = TRUE)
  }
  for (k in seq_len(max_cpts)) {
    if (n < min_n[k]) next
    fun <- .nm_sm_exp_fun(k)
    starts <- list()
    for (mult in list(c(1), c(5, 1), c(15, 3, 1), c(3, 1), c(10, 1), c(30, 6, 1), c(8, 2, 1))) {
      if (length(mult) != k) next
      lam <- lz * mult
      a0 <- rep(log(mean(exp(y))) , k) + log(1 / k) + log(mult)
      p <- c(a0, log(lam[k]))
      if (k > 1) for (i in 1:(k - 1)) p <- c(p, log(lam[i] / lam[i + 1] - 1))
      starts[[length(starts) + 1]] <- p
    }
    f <- .nm_sm_optim(fun, starts, t, y)
    if (is.null(f)) next
    p <- f$par
    lam <- numeric(k); lam[k] <- exp(p[k + 1])
    if (k > 1) for (i in (k - 1):1) lam[i] <- lam[i + 1] * (1 + exp(p[k + 1 + i]))
    A <- exp(p[1:k])
    frac <- A / sum(A)
    # identifiability: phases well separated and each phase visible
    valid <- TRUE
    if (k > 1) {
      if (any(lam[-k] / lam[-1] < 1.5)) valid <- FALSE
      if (any(frac < 0.01)) valid <- FALSE
      # each extra phase must be supported by data: fast phase decays within sampling range
      if (any(1 / lam[-k] > diff(range(t)))) valid <- FALSE
    }
    bic[k] <- if (valid) .nm_sm_bic(f$value, n, 2 * k) else NA
    par[[k]] <- list(A = A, lam = lam, sse = f$value, valid = valid)
  }
  list(bic = bic, par = par)
}

.nm_sm_curvature <- function(t, y) {
  if (length(t) < 5) return(list(coef = NA, p = NA))
  ts <- (t - mean(t)) / stats::sd(t)
  f <- stats::lm(y ~ ts + I(ts^2))
  s <- summary(f)$coefficients
  if (nrow(s) < 3) return(list(coef = NA, p = NA))
  list(coef = s[3, 1], p = s[3, 4])
}

# ---------------------------------------------------------------------------
# Absorption models (local one-compartment fits, dose-normalised data)
# ---------------------------------------------------------------------------
.nm_sm_absid <- function(x) c(`first-order` = "fo", `first-order+lag` = "lag",
                              `zero-order-d1` = "zo", transit = "tr")[x]

# log of the response of exp(-lam*t) to each unit input, vectorised over t
.nm_sm_g_fo <- function(t, lam, ka) {
  out <- rep(-Inf, length(t)); pos <- t > 0
  tt <- t[pos]
  out[pos] <- if (abs(ka - lam) / lam < 1e-6) log(ka * tt) - lam * tt else
    log(ka / abs(ka - lam)) + log(abs(exp(-lam * tt) - exp(-ka * tt)))
  out
}
.nm_sm_g_zo <- function(t, lam, D) {
  ifelse(t <= D, log((1 - exp(-lam * t)) / (lam * D)),
         log((1 - exp(-lam * D)) / (lam * D)) - lam * (t - D))
}
.nm_sm_g_tr <- function(t, lam, ktr, a) {
  a * log(ktr) - a * log(ktr - lam) - lam * t +
    stats::pgamma(t, shape = a, rate = ktr - lam, log.p = TRUE)
}

.nm_sm_fit_absorption <- function(t, y, k, disp, tmax, cmax, lz, t_first,
                                  blq_t = numeric(0), lloq_norm = NA) {
  n_obs <- length(t)
  use_blq <- length(blq_t) > 0 && is.finite(lloq_norm)
  lgl <- if (use_blq) log(lloq_norm) else NA
  k0 <- if (is.finite(lz) && lz > 0) lz else log(2) / max(t)
  if (is.null(disp)) disp <- list(A = exp(max(y)), lam = k0)
  # disposition: p[1:k] = log A_i, p[k+1] = log lam_k, p[k+1+i] = log(lam_i/lam_(i+1) - 1)
  dl <- function(p) {
    lam <- numeric(k); lam[k] <- exp(p[k + 1])
    if (k > 1) for (i in (k - 1):1) lam[i] <- lam[i + 1] * (1 + exp(p[k + 1 + i]))
    list(a = p[1:k], lam = lam)
  }
  np_d <- 2 * k
  disp_start <- c(log(disp$A), log(disp$lam[k]))
  if (k > 1) for (i in 1:(k - 1)) disp_start <- c(disp_start, log(disp$lam[i] / disp$lam[i + 1] - 1))
  # predicted log C from a list of per-exponential log responses
  comb <- function(dd, glist) {
    m <- sapply(seq_len(k), function(i) dd$a[i] + glist[[i]])
    if (is.null(dim(m))) m <- matrix(m, ncol = k)
    .nm_sm_lse(m)
  }
  # objective with BLQ censoring penalty: predictions above LLOQ at BLQ times
  fit <- function(pred_fun, starts) {
    tt_all <- if (use_blq) c(t, blq_t) else t
    obj <- function(p) {
      pr <- pred_fun(p, tt_all)
      po <- pr[seq_len(n_obs)]
      if (any(!is.finite(po))) return(1e10)
      sse <- sum((y - po)^2)
      if (use_blq) {
        pb <- pr[-seq_len(n_obs)]
        pb[!is.finite(pb)] <- -Inf
        sse <- sse + sum(pmax(pb - lgl, 0)^2)
      }
      sse
    }
    .nm_sm_multistart(obj, starts)
  }
  n_eff <- n_obs + if (use_blq) length(blq_t) else 0
  bic <- function(f, np) {
    if (is.null(f) || n_eff < np + 3) return(NA)
    .nm_sm_bic(f$value, n_eff, np)
  }
  # rescale amplitudes so that the start matches Cmax
  rescale <- function(pred_fun, p) {
    pm <- pred_fun(p, tmax)
    if (is.finite(pm)) p[1:k] <- p[1:k] + log(cmax) - pm
    p
  }
  # first-order
  pf_fo <- function(p, tt) { dd <- dl(p); ka <- exp(p[np_d + 1])
    comb(dd, lapply(dd$lam, function(l) .nm_sm_g_fo(tt, l, ka))) }
  st <- lapply(c(2, 5, 15, 50), function(m) rescale(pf_fo, c(disp_start, log(k0 * m))))
  f_fo <- fit(pf_fo, st)
  # first-order + lag (0 < tlag < first quantified time)
  pf_lg <- function(p, tt) { dd <- dl(p); ka <- exp(p[np_d + 1]); tl <- t_first * stats::plogis(p[np_d + 2])
    comb(dd, lapply(dd$lam, function(l) .nm_sm_g_fo(tt - tl, l, ka))) }
  st <- list()
  for (m in c(4, 20)) for (q in stats::qlogis(c(0.4, 0.9)))
    st[[length(st) + 1]] <- rescale(pf_lg, c(disp_start, log(k0 * m), q))
  if (!is.null(f_fo)) st[[length(st) + 1]] <- c(f_fo$par, stats::qlogis(0.3))
  f_lg <- fit(pf_lg, st)
  # zero-order, estimated duration D
  pf_zo <- function(p, tt) { dd <- dl(p); D <- exp(p[np_d + 1])
    comb(dd, lapply(dd$lam, function(l) .nm_sm_g_zo(tt, l, D))) }
  st <- lapply(c(0.5, 0.8, 1, 1.2), function(m) rescale(pf_zo, c(disp_start, log(tmax * m))))
  f_zo <- fit(pf_zo, st)
  # gamma (transit) input: shape a > 1, ktr > fastest disposition rate
  pf_tr <- function(p, tt) { dd <- dl(p); ktr <- dd$lam[1] + exp(p[np_d + 1]); a <- 1 + exp(p[np_d + 2])
    comb(dd, lapply(dd$lam, function(l) .nm_sm_g_tr(tt, l, ktr, a))) }
  st <- list()
  for (a in c(2, 4, 8)) {
    ktr <- max(a / tmax, disp$lam[1] * 1.5)
    st[[length(st) + 1]] <- rescale(pf_tr, c(disp_start, log(ktr - disp$lam[1]), log(a - 1)))
  }
  f_tr <- fit(pf_tr, st)

  out <- data.frame(bic_fo = bic(f_fo, np_d + 1), bic_lag = bic(f_lg, np_d + 2),
                    bic_zo = bic(f_zo, np_d + 1), bic_tr = bic(f_tr, np_d + 2),
                    n_abs = n_eff, k_disp = k)
  dd <- if (!is.null(f_fo)) dl(f_fo$par) else NULL
  out$fo_V  <- if (!is.null(dd) && k == 1) exp(-dd$a[1]) else NA
  out$fo_k  <- if (!is.null(dd)) dd$lam[k] else NA
  out$fo_ka <- if (!is.null(f_fo)) exp(f_fo$par[np_d + 1]) else NA
  for (i in seq_len(k)) {
    out[[sprintf("fo_A%d_%d", k, i)]]   <- if (!is.null(dd)) exp(dd$a[i]) else NA
    out[[sprintf("fo_lam%d_%d", k, i)]] <- if (!is.null(dd)) dd$lam[i] else NA
  }
  out$lag_ka <- if (!is.null(f_lg)) exp(f_lg$par[np_d + 1]) else NA
  out$lag_tlag <- if (!is.null(f_lg)) t_first * stats::plogis(f_lg$par[np_d + 2]) else NA
  out$zo_D  <- if (!is.null(f_zo)) exp(f_zo$par[np_d + 1]) else NA
  if (!is.null(f_tr)) {
    d2 <- dl(f_tr$par); ktr <- d2$lam[1] + exp(f_tr$par[np_d + 1]); a <- 1 + exp(f_tr$par[np_d + 2])
    out$tr_mtt <- a / ktr; out$tr_shape <- a
  } else { out$tr_mtt <- NA; out$tr_shape <- NA }
  # predictions (log, dose-normalised) at the observed times and on a grid,
  # used for the residual error models and the plots
  tg <- seq(max(t) / 200, max(t), length.out = 150)
  fits <- list(fo = list(f_fo, pf_fo), lag = list(f_lg, pf_lg),
               zo = list(f_zo, pf_zo), tr = list(f_tr, pf_tr))
  lp <- sapply(fits, function(z) if (is.null(z[[1]])) rep(NA_real_, n_obs) else z[[2]](z[[1]]$par, t))
  lg <- sapply(fits, function(z) if (is.null(z[[1]])) rep(NA_real_, length(tg)) else z[[2]](z[[1]]$par, tg))
  if (is.null(dim(lp))) lp <- matrix(lp, nrow = n_obs, dimnames = list(NULL, names(fits)))
  attr(out, "logpred") <- lp
  attr(out, "grid") <- data.frame(time = tg, lg, check.names = FALSE)
  attr(out, "fits") <- lapply(fits, function(z) if (is.null(z[[1]])) NULL else list(fun = z[[2]], par = z[[1]]$par))
  out
}

# ---------------------------------------------------------------------------
# BIC weights
# ---------------------------------------------------------------------------
.nm_sm_weights_table <- function(tab, prefix, levels) {
  # columns: <prefix><level>, e.g. bic_cpt1 or bic_fo
  cols <- paste0(prefix, levels)
  m <- as.matrix(tab[, cols, drop = FALSE])
  w <- t(apply(m, 1, function(b) {
    if (all(!is.finite(b))) return(rep(NA_real_, length(b)))
    b[!is.finite(b)] <- Inf
    e <- exp(-0.5 * (b - min(b)))
    e / sum(e)
  }))
  if (length(levels) == 1) w <- matrix(w, ncol = 1)
  colnames(w) <- paste0("w_", levels)
  keep <- stats::complete.cases(w)
  mean_w <- if (any(keep)) colMeans(w[keep, , drop = FALSE]) else rep(1 / length(levels), length(levels))
  best <- apply(w, 1, function(r) if (all(is.na(r))) NA else levels[which.max(r)])
  table_out <- cbind(tab[, "pid", drop = FALSE], round(w, 3), best = best)
  list(table = table_out, mean_w = mean_w)
}

# ---------------------------------------------------------------------------
# Pooled (naive) analysis for sparse data
# ---------------------------------------------------------------------------
.nm_sm_pooled <- function(pinfo, pobs, n_bins, route, infusion, max_cpts, factor) {
  first_types <- c("single", "single_interval")
  fp <- pinfo$pid[pinfo$type %in% first_types]
  q <- pobs[!pobs$BLQ & pobs$DV > 0 & pobs$TAD > 0, ]
  basis <- "first-dose intervals"
  qf <- q[q$pid %in% fp, ]
  if (nrow(qf) >= 3 * n_bins) {
    q <- qf; ptype <- "single"
  } else {
    qs <- q[q$pid %in% pinfo$pid[pinfo$type == "ss"], ]
    if (nrow(qs) >= 3 * n_bins) { q <- qs; ptype <- "ss"; basis <- "steady-state intervals" }
    else { ptype <- "other"; basis <- "all intervals (superposition ignored)" }
  }
  if (nrow(q) < max(12, 2 * n_bins)) return(list(res = NULL, basis = basis, bins = NULL))
  q$y <- q$DV / q$AMT
  br <- unique(stats::quantile(q$TAD, probs = seq(0, 1, length.out = n_bins + 1)))
  if (length(br) < 4) return(list(res = NULL, basis = basis, bins = NULL))
  q$bin <- cut(q$TAD, br, include.lowest = TRUE)
  bins <- do.call(rbind, lapply(split(q, q$bin), function(b) {
    if (nrow(b) < 3) return(NULL)
    data.frame(tad = stats::median(b$TAD), gm = exp(mean(log(b$y))), n = nrow(b))
  }))
  if (is.null(bins) || nrow(bins) < 4) return(list(res = NULL, basis = basis, bins = bins))
  tau <- if (ptype == "ss") stats::median(pinfo$tau[pinfo$type == "ss"]) else Inf
  tinf <- if (infusion) stats::median(pinfo$TINF) else 0
  res <- .nm_sm_analyse_profile(bins$tad, bins$gm, rep(FALSE, nrow(bins)), dose = 1,
                                type = ptype, tau = tau, tinf = tinf, route = route,
                                max_cpts = max_cpts, factor = factor)
  list(res = res, basis = basis, bins = bins, type = ptype)
}

# ---------------------------------------------------------------------------
# Elimination: dose proportionality and curvature
# ---------------------------------------------------------------------------
.nm_sm_elimination <- function(nca, curv, pinfo, pobs, n_bins) {
  out <- list(method = NA, slope = NA, lower = NA, upper = NA, n_levels = NA,
              n_concave = NA, n_convex = NA, p_mm = 0.05, note = NULL, assessable = FALSE)
  # curvature
  if (!is.null(curv) && nrow(curv)) {
    out$n_concave <- sum(curv$curv < 0 & curv$p < 0.05, na.rm = TRUE)
    out$n_convex  <- sum(curv$curv > 0 & curv$p < 0.05, na.rm = TRUE)
  }
  # rich: power model on AUC
  fit <- NULL
  if (!is.null(nca) && nrow(nca)) {
    a <- nca
    a$auc <- ifelse(a$type == "ss", a$auc_tau,
                    ifelse(is.finite(a$pct_extrap) & a$pct_extrap < 30, a$auc_inf, NA))
    a <- a[is.finite(a$auc) & a$auc > 0, ]
    lev <- unique(a$dose)
    if (nrow(a) >= 4 && length(lev) >= 2 && max(lev) / min(lev) >= 2) {
      fit <- if (length(unique(a$type)) > 1)
        stats::lm(log(auc) ~ log(dose) + type, data = a) else stats::lm(log(auc) ~ log(dose), data = a)
      ci <- stats::confint(fit)["log(dose)", ]
      out$method <- "power model on AUC"; out$slope <- stats::coef(fit)["log(dose)"]
      out$lower <- ci[1]; out$upper <- ci[2]; out$n_levels <- length(lev)
    }
  }
  # sparse: dose-normalised concentrations by time bin
  if (is.null(fit)) {
    q <- pobs[!pobs$BLQ & pobs$DV > 0 & pobs$TAD > 0, ]
    lev <- unique(q$AMT)
    if (nrow(q) >= 20 && length(lev) >= 2 && max(lev) / min(lev) >= 2) {
      br <- unique(stats::quantile(q$TAD, probs = seq(0, 1, length.out = n_bins + 1)))
      q$bin <- cut(q$TAD, br, include.lowest = TRUE)
      q$type <- pinfo$type[match(q$pid, pinfo$pid)]
      f <- tryCatch(stats::lm(log(DV / AMT) ~ bin + log(AMT), data = q), error = function(e) NULL)
      if (!is.null(f) && "log(AMT)" %in% names(stats::coef(f))) {
        ci <- stats::confint(f)["log(AMT)", ]
        out$method <- "dose-normalised concentrations by time bin"
        out$slope <- 1 + stats::coef(f)["log(AMT)"]
        out$lower <- 1 + ci[1]; out$upper <- 1 + ci[2]; out$n_levels <- length(lev)
        fit <- f
      }
    }
  }
  if (!is.null(fit)) {
    out$assessable <- TRUE
    # one-sided p-value for slope > 1
    se <- (out$upper - out$lower) / (2 * stats::qt(0.975, stats::df.residual(fit)))
    p1 <- stats::pt((out$slope - 1) / se, stats::df.residual(fit), lower.tail = FALSE)
    out$p_greater <- p1
    out$p_mm <- if (p1 < 0.025) 0.85 else if (p1 < 0.15) 0.3 else 0.08
    if (out$upper < 1)
      out$note <- sprintf("Less than dose-proportional exposure (slope %.2f [%.2f; %.2f]): consider saturable absorption or bioavailability rather than saturable elimination.",
                          out$slope, out$lower, out$upper)
  } else {
    out$note <- "Single dose level (or range < 2-fold): saturable elimination cannot be assessed from dose proportionality."
  }
  if (is.finite(out$n_concave) && out$n_concave >= 2 && out$n_concave > 2 * out$n_convex) {
    out$p_mm <- max(out$p_mm, 0.5)
    out$note <- c(out$note, sprintf("%d profile(s) with a concave log-concentration decline (typical of saturable elimination).",
                                     out$n_concave))
    out$note <- paste(out$note, collapse = " ")
  }
  out
}

# ---------------------------------------------------------------------------
# Initial values
# ---------------------------------------------------------------------------
.nm_sm_med <- function(x) { x <- x[is.finite(x) & x > 0]; if (length(x)) stats::median(x) else NA }

.nm_sm_micro <- function(A, lam) {
  # unit-dose macro constants -> micro constants (central volume = 1/sum(A))
  k <- length(lam); a <- A / sum(A)
  if (k == 1) return(list(k10 = lam))
  if (k == 2) {
    k21 <- a[1] * lam[2] + a[2] * lam[1]
    k10 <- lam[1] * lam[2] / k21
    k12 <- lam[1] + lam[2] - k21 - k10
    if (any(c(k21, k10, k12) <= 0)) return(NULL)
    return(list(k10 = k10, k12 = k12, k21 = k21))
  }
  b <- a[1] * (lam[2] + lam[3]) + a[2] * (lam[1] + lam[3]) + a[3] * (lam[1] + lam[2])
  cc <- a[1] * lam[2] * lam[3] + a[2] * lam[1] * lam[3] + a[3] * lam[1] * lam[2]
  disc <- b^2 - 4 * cc
  if (disc <= 0) return(NULL)
  k21 <- (b + sqrt(disc)) / 2; k31 <- (b - sqrt(disc)) / 2
  k10 <- prod(lam) / (k21 * k31)
  s1 <- sum(lam) - k10 - k21 - k31
  s2 <- lam[1] * lam[2] + lam[1] * lam[3] + lam[2] * lam[3] - k10 * (k21 + k31) - k21 * k31
  # k12 + k13 = s1 ; k12*k31 + k13*k21 = s2
  k12 <- (s2 - s1 * k21) / (k31 - k21); k13 <- s1 - k12
  if (any(c(k21, k31, k10, k12, k13) <= 0)) return(NULL)
  list(k10 = k10, k12 = k12, k21 = k21, k13 = k13, k31 = k31)
}

.nm_sm_initial_values <- function(nca, cpt, absf, pooled, source_ev, route,
                                  factor, obs, max_cpts) {
  iv <- list(); note <- NULL
  if (source_ev == "individual") {
    cl <- .nm_sm_med(nca$cl); vz <- .nm_sm_med(nca$vz); lz <- .nm_sm_med(nca$lambda_z)
    cp <- cpt; ab <- absf
  } else {
    r <- pooled$res
    cl <- .nm_sm_med(r$nca$cl); vz <- .nm_sm_med(r$nca$vz); lz <- .nm_sm_med(r$nca$lambda_z)
    cp <- r$cpts; ab <- r$abs
  }
  if (!is.finite(cl)) {
    note <- "CL could not be derived (no first-dose or steady-state AUC): default value used."
    cl <- 5
  }
  if (!is.finite(vz)) vz <- if (is.finite(lz)) cl / lz else 50
  iv$CL <- cl; iv$Vz <- vz
  # one compartment central volume
  v1c <- if (route == "EV" && !is.null(ab) && nrow(ab)) .nm_sm_med(ab$fo_V) * factor else NA
  if (route == "IV" && !is.null(cp)) {
    A1 <- .nm_sm_med(cp$A1_1); if (is.finite(A1)) v1c <- factor / A1
  }
  iv$V_1cpt <- if (is.finite(v1c)) v1c else vz
  # multi-compartment: median macro constants -> micro constants scaled to CL
  for (k in 2:3) {
    if (k > max_cpts || is.null(cp)) next
    A <- sapply(seq_len(k), function(i) .nm_sm_med(cp[[sprintf("A%d_%d", k, i)]]))
    lam <- sapply(seq_len(k), function(i) .nm_sm_med(cp[[sprintf("lam%d_%d", k, i)]]))
    # extravascular: prefer the amplitudes of the full first-order fits, which
    # are not distorted by absorption
    if (route == "EV" && !is.null(ab) && sprintf("fo_A%d_1", k) %in% names(ab)) {
      A2 <- sapply(seq_len(k), function(i) .nm_sm_med(ab[[sprintf("fo_A%d_%d", k, i)]]))
      l2 <- sapply(seq_len(k), function(i) .nm_sm_med(ab[[sprintf("fo_lam%d_%d", k, i)]]))
      if (all(is.finite(c(A2, l2)))) { A <- A2; lam <- l2 }
    }
    mc <- if (all(is.finite(c(A, lam)))) .nm_sm_micro(A, lam) else NULL
    if (is.null(mc)) {
      # fallback heuristic
      v1 <- 0.3 * vz
      iv[[paste0("cpt", k)]] <- if (k == 2) list(V1 = v1, Q = 2 * cl, V2 = vz - v1, heuristic = TRUE) else
        list(V1 = v1, Q2 = 2 * cl, V2 = 0.4 * vz, Q3 = 0.3 * cl, V3 = 0.3 * vz, heuristic = TRUE)
      next
    }
    V1 <- cl / mc$k10
    vals <- c(V1, unlist(mc) * V1)
    if (any(!is.finite(vals)) || any(vals > 1e3 * max(vz, cl)) || any(vals < 1e-3 * min(vz, cl))) {
      v1 <- 0.3 * vz
      iv[[paste0("cpt", k)]] <- if (k == 2) list(V1 = v1, Q = 2 * cl, V2 = vz - v1, heuristic = TRUE) else
        list(V1 = v1, Q2 = 2 * cl, V2 = 0.4 * vz, Q3 = 0.3 * cl, V3 = 0.3 * vz, heuristic = TRUE)
      next
    }
    if (k == 2) iv$cpt2 <- list(V1 = V1, Q = mc$k12 * V1, V2 = mc$k12 * V1 / mc$k21, heuristic = FALSE)
    else iv$cpt3 <- list(V1 = V1, Q2 = mc$k12 * V1, V2 = mc$k12 * V1 / mc$k21,
                         Q3 = mc$k13 * V1, V3 = mc$k13 * V1 / mc$k31, heuristic = FALSE)
  }
  # absorption
  if (route == "EV") {
    iv$KA   <- if (!is.null(ab)) .nm_sm_med(ab$fo_ka) else NA
    iv$KA_lag <- if (!is.null(ab)) .nm_sm_med(ab$lag_ka) else NA
    iv$ALAG <- if (!is.null(ab)) .nm_sm_med(ab$lag_tlag) else NA
    iv$DUR  <- if (!is.null(ab)) .nm_sm_med(ab$zo_D) else NA
    iv$MTT  <- if (!is.null(ab)) .nm_sm_med(ab$tr_mtt) else NA
    iv$NTR  <- if (!is.null(ab)) .nm_sm_med(ab$tr_shape) else NA
    tmax <- if (source_ev == "individual") .nm_sm_med(nca$tmax) else .nm_sm_med(pooled$res$nca$tmax)
    if (!is.finite(iv$KA)) {   # from tmax with k = lambda_z
      k <- if (is.finite(lz)) lz else cl / vz
      f <- function(ka) log(ka / k) / (ka - k) - tmax
      iv$KA <- tryCatch(stats::uniroot(f, c(k * 1.001, k * 1000))$root, error = function(e) 1)
    }
    if (!is.finite(iv$KA_lag)) iv$KA_lag <- iv$KA
    if (!is.finite(iv$ALAG)) iv$ALAG <- if (is.finite(tmax)) 0.2 * tmax else 0.5
    if (!is.finite(iv$DUR))  iv$DUR  <- if (is.finite(tmax)) tmax else 1
    if (!is.finite(iv$MTT))  iv$MTT  <- if (is.finite(tmax)) 0.8 * tmax else 1
    if (!is.finite(iv$NTR))  iv$NTR  <- 3
  }
  # Michaelis-Menten: KM ~ median concentration (dose units / volume), VMAX = CL * KM
  q <- obs$DV[!obs$BLQ]
  iv$KM <- stats::median(q) / factor
  iv$VMAX <- cl * iv$KM
  iv$note <- note
  iv
}

# ---------------------------------------------------------------------------
# Flags: parameters likely to be hard to estimate
# ---------------------------------------------------------------------------
.nm_sm_flags <- function(nca, cpt, absf, pooled, source_ev, route, el, obs, n_subj) {
  f <- list()
  q <- obs[!obs$BLQ, ]
  per_subj <- as.numeric(table(q$ID))
  f$median_obs_subject <- stats::median(per_subj)
  # in the pooled case, the same quantities are read on the binned profile
  src <- if (source_ev == "individual") nca else pooled$res$nca
  if (source_ev != "individual") { cpt <- pooled$res$cpts; absf <- pooled$res$abs }
  n_pre  <- stats::median(src$n_pre_tmax)
  span   <- stats::median(src$span_hl, na.rm = TRUE)
  extrap <- stats::median(src$pct_extrap, na.rm = TRUE)
  tfirst <- stats::median(src$t_first); tmax <- stats::median(src$tmax)
  lagev  <- if (source_ev == "individual") mean(nca$n_blq_early > 0) else NA
  f$n_pre_tmax <- n_pre; f$span_hl <- span; f$pct_extrap <- extrap
  f$lag_evidence <- lagev
  # absorption
  f$abs_poor <- route == "EV" && (is.na(n_pre) || n_pre < 2)
  f$abs_shape_poor <- route == "EV" && (is.na(n_pre) || n_pre < 3)
  f$lag_late <- route == "EV" && is.finite(tfirst) && is.finite(tmax) && tfirst > 0.5 * tmax
  f$flipflop <- FALSE
  if (route == "EV" && !is.null(absf) && nrow(absf)) {
    ratio <- absf$fo_ka / absf$fo_k
    f$flipflop <- mean(ratio < 2, na.rm = TRUE) >= 0.5
  }
  # terminal phase
  f$terminal_short <- (is.finite(span) && span < 2) || (is.finite(extrap) && extrap > 20)
  # distribution
  f$dist_poor <- FALSE; f$lam_ratio <- NA
  if (!is.null(cpt) && nrow(cpt) && all(c("lam2_1", "lam2_2") %in% names(cpt))) {
    lr <- stats::median(cpt$lam2_1 / cpt$lam2_2, na.rm = TRUE)
    f$lam_ratio <- lr
    f$dist_poor <- !is.finite(lr) || lr < 3
  }
  f$sparse <- source_ev == "pooled"
  f$few_obs <- f$median_obs_subject < 4
  f$mm_not_assessable <- !el$assessable
  f$p_mm <- el$p_mm
  f$n_subj <- n_subj
  f
}

.nm_sm_fmt <- function(x) vapply(x, function(v) format(signif(v, 3), scientific = FALSE,
                                                         drop0trailing = TRUE), "")

.nm_sm_bind <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (!length(lst)) return(NULL)
  allc <- unique(unlist(lapply(lst, names)))
  lst <- lapply(lst, function(d) { for (m in setdiff(allc, names(d))) d[[m]] <- NA; d[, allc] })
  out <- do.call(rbind, lst)
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Individual estimates (IIV initial values, covariate screening)
# ---------------------------------------------------------------------------
.nm_sm_individual <- function(pid, id, res, route, factor) {
  vc <- NA; ka <- NA; dur <- NA; mtt <- NA
  ab <- res$abs
  if (!is.null(ab)) {
    k <- ab$k_disp
    A <- unlist(ab[1, sprintf("fo_A%d_%d", k, seq_len(k))])
    if (all(is.finite(A))) vc <- factor / sum(A)
    ka <- ab$fo_ka; dur <- ab$zo_D; mtt <- ab$tr_mtt
  } else {
    b <- unlist(res$cpts[1, grep("^bic_cpt", names(res$cpts))])
    if (any(is.finite(b))) {
      k <- which.min(ifelse(is.finite(b), b, Inf))
      A <- unlist(res$cpts[1, sprintf("A%d_%d", k, seq_len(k))])
      if (all(is.finite(A))) vc <- factor / sum(A)
    }
  }
  data.frame(pid = pid, ID = id, CL = res$nca$cl, Vc = vc, KA = ka, DUR = dur, MTT = mtt,
             stringsAsFactors = FALSE)
}

.nm_sm_iiv <- function(indiv, source_ev) {
  default <- c(CL = 30, Vc = 30, KA = 50, DUR = 50, MTT = 50)
  cv <- default; src <- stats::setNames(rep("default", 5), names(default))
  note <- NULL
  if (source_ev == "individual" && !is.null(indiv) && nrow(indiv)) {
    per <- stats::aggregate(indiv[, names(default)], by = list(ID = indiv$ID),
                            FUN = function(v) { v <- v[is.finite(v) & v > 0]; if (length(v)) stats::median(v) else NA })
    for (p in names(default)) {
      v <- per[[p]]; v <- v[is.finite(v) & v > 0]
      if (length(v) >= 5) {
        om <- stats::var(log(v))
        cv[p] <- min(100, max(10, 100 * sqrt(exp(om) - 1)))
        src[p] <- "data"
      }
    }
  } else {
    note <- "IIV initial values are defaults (CV 30% on CL and volume, 50% on absorption): individual estimates are not available from pooled data."
  }
  list(cv = round(cv), source = src, note = note)
}

# ---------------------------------------------------------------------------
# Residual error models
# ---------------------------------------------------------------------------
.nm_sm_error_models <- function(fits, source_ev, obs, lloq_dv) {
  models <- c("additive", "proportional", "combined")
  q <- obs$DV[!obs$BLQ & obs$DV > 0]
  add_default <- if (is.finite(lloq_dv)) (lloq_dv / 2)^2 else (0.05 * stats::median(q))^2
  default <- list(additive = c(prop = 0.04, add = add_default),
                  proportional = c(prop = 0.04, add = add_default),
                  combined = c(prop = 0.04, add = add_default))
  n_tot <- sum(vapply(fits, function(f) length(f$t), numeric(1)))
  if (source_ev != "individual" || length(fits) < 3 || n_tot < 20) {
    return(list(support = stats::setNames(c(0.2, 0.4, 0.4), models), sigma = default, table = NULL,
                note = "Residual error model not assessable without rich individual profiles: proportional or combined error assumed by default; SIGMA initial values are defaults."))
  }
  # Extended least squares: the best structural model of each interval is
  # refitted under the weighting of each error model, so that the residuals
  # reflect that model (a fit on the log scale assumes a proportional error).
  refit <- function(f, wfun) {
    y <- f$y * f$dose
    obj <- function(p) {
      pr <- exp(f$fun(p, f$t)) * f$dose
      if (any(!is.finite(pr))) return(1e10)
      sum((y - pr)^2 * wfun(pr))
    }
    r <- tryCatch(stats::optim(f$par, obj, method = "Nelder-Mead",
                               control = list(maxit = 3000, reltol = 1e-10)), error = function(e) NULL)
    par <- if (is.null(r) || !is.finite(r$value) || r$value >= 1e10) f$par else r$par
    press(f, par, y, wfun)
  }
  # PRESS residuals e_i / (1 - h_ii): leverage-corrected (leave-one-out)
  # residuals, since each interval has few points for several parameters
  press <- function(f, par, y, wfun) {
    fp <- exp(f$fun(par, f$t)) * f$dose
    J <- sapply(seq_along(par), function(j) {
      d <- rep(0, length(par)); d[j] <- 1e-5
      (exp(f$fun(par + d, f$t)) - exp(f$fun(par - d, f$t))) * f$dose / 2e-5
    })
    if (is.null(dim(J))) J <- matrix(J, nrow = length(fp))
    w <- wfun(fp)
    h <- tryCatch({
      M <- solve(t(J) %*% (J * w))
      rowSums((J %*% M) * J) * w
    }, error = function(e) rep(NA_real_, length(fp)))
    # high-leverage points (last samples of a profile) carry little
    # information on the residual error and 1/(1-h) would inflate them
    ok <- is.finite(h) & h < 0.7 & is.finite(fp) & fp > 0
    list(y = fp[ok] + (y[ok] - fp[ok]) / (1 - h[ok]), f = fp[ok], np = length(par))
  }
  # observations predicted below 2 x LLOQ are left out of the variance
  # models: near the LLOQ, censoring keeps only the values pushed upwards by
  # the error, which mimics an additive floor
  f_min <- if (is.finite(lloq_dv)) 2 * lloq_dv else 0
  ml <- function(res, s2fun, start) {
    y <- unlist(lapply(res, `[[`, "y")); fp <- unlist(lapply(res, `[[`, "f"))
    ok <- is.finite(fp) & fp > 0 & fp >= f_min
    y <- y[ok]; fp <- fp[ok]
    nll <- function(p) { s2 <- s2fun(p, fp); 0.5 * sum(log(2 * pi * s2) + (y - fp)^2 / s2) }
    r <- tryCatch(stats::optim(start, nll, method = if (length(start) == 1) "BFGS" else "Nelder-Mead"),
                  error = function(e) list(par = start, value = nll(start)))
    list(par = r$par, value = r$value, n = length(y))
  }
  # Starting weights: additive variance from an unweighted refit,
  # proportional variance from the log-scale fits.
  r_add <- lapply(fits, refit, wfun = function(pr) 1)
  r_prop <- lapply(fits, function(f) press(f, f$par, f$y * f$dose, function(pr) 1 / pr^2))
  m_add0 <- ml(r_add, function(p, f) rep(exp(2 * p[1]), length(f)), log(stats::sd(unlist(lapply(r_add, function(r) r$y - r$f)))))
  m_prop0 <- ml(r_prop, function(p, f) exp(2 * p[1]) * f^2,
                log(stats::sd(unlist(lapply(r_prop, function(r) (r$y - r$f) / r$f)))))
  a2 <- exp(2 * m_add0$par[1]) / 2; b2 <- exp(2 * m_prop0$par[1]) / 2
  fmed <- stats::median(unlist(lapply(r_prop, `[[`, "f")))
  # Common structural fit: iteratively reweighted least squares under the
  # combined variance a2 + b2 f^2, so that the three variance models are
  # compared on the same (leverage-corrected) residuals and are nested.
  for (it in 1:3) {
    res <- lapply(fits, refit, wfun = function(pr) 1 / (a2 + b2 * pr^2))
    mc <- ml(res, function(p, f) exp(2 * p[1]) + exp(2 * p[2]) * f^2,
             c(log(sqrt(a2)), log(sqrt(b2))))
    a2 <- max(exp(2 * mc$par[1]), 1e-6 * fmed^2); b2 <- max(exp(2 * mc$par[2]), 1e-6)
  }
  res <- lapply(fits, refit, wfun = function(pr) 1 / (a2 + b2 * pr^2))
  f_all <- unlist(lapply(res, `[[`, "f")); e_all <- unlist(lapply(res, function(r) r$y - r$f))
  keep <- is.finite(f_all) & f_all >= f_min
  e_all <- e_all[keep]; f_all <- f_all[keep]
  m_add <- ml(res, function(p, f) rep(exp(2 * p[1]), length(f)), log(stats::sd(e_all)))
  m_prop <- ml(res, function(p, f) exp(2 * p[1]) * f^2, log(stats::sd(e_all / f_all)))
  m_comb <- ml(res, function(p, f) exp(2 * p[1]) + exp(2 * p[2]) * f^2, c(log(sqrt(a2)), log(sqrt(b2))))
  # the combined model nests the two others: keep its -2LL at most theirs
  m_comb$value <- min(m_comb$value, m_add$value, m_prop$value)
  n <- m_add$n
  # the combined model must be strongly supported (BIC better by more than
  # 15): misfit of individual profiles mimics an extra error component. The
  # threshold (10 until 0.9.0, "very strong" for Kass and Raftery 1995) was
  # set to 15 after the simulation validation (point E28): correct error
  # model in rich designs 87% on average over the three true error models,
  # as with 10, with the proportional model found more often (86% vs 82%)
  # and the combined model a little less often (90% vs 95%).
  bic <- c(2 * m_add$value + log(n), 2 * m_prop$value + log(n), 2 * m_comb$value + 2 * log(n) + 15)
  w <- exp(-0.5 * (bic - min(bic))); w <- w / sum(w)
  # degrees of freedom: individual parameters were estimated on the same data
  dfc <- 1
  est <- list(additive = c(add = exp(2 * m_add$par[1]) * dfc),
              proportional = c(prop = exp(2 * m_prop$par[1]) * dfc),
              combined = c(prop = exp(2 * m_comb$par[2]) * dfc, add = exp(2 * m_comb$par[1]) * dfc))
  # initial values: a component estimated near 0 would start NONMEM on a
  # boundary; floors of 10% of the default additive variance and a 5% CV
  sigma <- list(additive = c(prop = 0.04, add = max(unname(est$additive["add"]), 0.1 * add_default)),
                proportional = c(prop = max(unname(est$proportional["prop"]), 0.0025), add = add_default),
                combined = c(prop = max(unname(est$combined["prop"]), 0.0025),
                             add = max(unname(est$combined["add"]), 0.1 * add_default)))
  tab <- data.frame(model = models,
                    sd_add = c(sqrt(est$additive["add"]), NA, sqrt(est$combined["add"])),
                    cv_prop_pct = 100 * c(NA, sqrt(est$proportional["prop"]), sqrt(est$combined["prop"])),
                    bic = bic, weight = round(w, 3), row.names = NULL)
  list(support = stats::setNames(w, models), sigma = sigma, table = tab,
       note = "SIGMA initial values come from the residuals of individual fits (common reweighted fit, leverage-corrected residuals): they include part of the structural misfit.")
}

# ---------------------------------------------------------------------------
# Alerts: features outside the tested structures
# ---------------------------------------------------------------------------
.nm_sm_alerts <- function(nca, obs_rows, el) {
  al <- list()
  add <- function(type, m) al[[length(al) + 1]] <<- data.frame(type = type, message = m, stringsAsFactors = FALSE)
  # double peaks
  if (length(obs_rows)) {
    dp <- vapply(obs_rows, function(o) {
      o <- o[!o$blq & o$dv > 0 & o$time > 0, ]
      o <- o[order(o$time), ]
      cc <- o$dv; n <- length(cc); cmax <- max(cc)
      if (n < 4) return(FALSE)
      for (i in 1:(n - 2)) for (k in (i + 2):n) {
        lo <- min(cc[i], cc[k]); m <- min(cc[(i + 1):(k - 1)])
        if (lo > 1.3 * m && lo - m > 0.15 * cmax && lo > 0.3 * cmax) return(TRUE)
      }
      FALSE
    }, logical(1))
    if (sum(dp) >= 2 && mean(dp) >= 0.2)
      add("double peak", sprintf("%d of %d intervals show a second peak after the maximum: enterohepatic recycling or several absorption sites are not described by the candidate models.",
                                 sum(dp), length(dp)))
  }
  if (!is.null(nca) && nrow(nca)) {
    first <- nca[nca$type %in% c("single", "single_interval"), ]
    ss <- nca[nca$type == "ss", ]
    # clearance at first dose vs steady state
    a <- first[is.finite(first$cl), c("ID", "cl")]; b <- ss[is.finite(ss$cl), c("ID", "cl")]
    m <- data.frame()
    if (nrow(a) && nrow(b)) {
      a <- stats::aggregate(cl ~ ID, a, stats::median); b <- stats::aggregate(cl ~ ID, b, stats::median)
      m <- merge(a, b, by = "ID", suffixes = c("_first", "_ss"))
    }
    if (nrow(m) >= 3) {
      r <- m$cl_ss / m$cl_first
      p <- tryCatch(suppressWarnings(stats::wilcox.test(log(r))$p.value), error = function(e) NA)
      if ((stats::median(r) > 1.25 || stats::median(r) < 0.8) && is.finite(p) && p < 0.1)
        add("time-varying clearance", sprintf("Steady-state CL is %.2f times the first-dose CL (median of %d subjects, p = %.3f): auto-induction, inhibition or time-varying CL.",
                                              stats::median(r), nrow(m), p))
    }
    # accumulation: observed vs predicted from the first-dose terminal slope
    a <- first[is.finite(first$auc_0tau) & is.finite(first$lambda_z), c("ID", "auc_0tau", "lambda_z", "dose")]
    b <- ss[is.finite(ss$auc_tau) & is.finite(ss$tau), c("ID", "auc_tau", "tau", "dose")]
    m <- merge(a, b, by = "ID", suffixes = c("_first", "_ss"))
    if (nrow(m) >= 3) {
      rac_obs <- (m$auc_tau / m$dose_ss) / (m$auc_0tau / m$dose_first)
      rac_pred <- 1 / (1 - exp(-m$lambda_z * m$tau))
      r <- stats::median(rac_obs / rac_pred)
      if (r > 1.33 || r < 0.75)
        add("accumulation", sprintf("Observed accumulation is %.2f times the one predicted from the first-dose terminal slope (median of %d subjects): non-linear or time-varying kinetics, or a long terminal phase not seen after the first dose.",
                                    r, nrow(m)))
    }
  }
  if (isTRUE(el$assessable) && is.finite(el$upper) && el$upper < 1)
    add("less than dose-proportional", sprintf("Exposure increases less than the dose (slope %.2f [%.2f; %.2f]): saturable absorption or bioavailability (dose effect on F or KA) is not among the candidate models.",
                                                el$slope, el$lower, el$upper))
  if (!length(al)) return(data.frame(type = character(0), message = character(0), stringsAsFactors = FALSE))
  do.call(rbind, al)
}

# ---------------------------------------------------------------------------
# Candidate specification (arguments of nm_write_ctl) and its call
# ---------------------------------------------------------------------------
.nm_sm_spec <- function(cpts, route, input, elimination, error, iv, iiv, err) {
  absorption <- switch(input, `first-order` = "first-order", `first-order+lag` = "first-order",
                       `zero-order-d1` = "zero-order-d1", transit = "transit", "first-order")
  transit_n <- if (input == "transit") max(1L, as.integer(round(iv$NTR - 1))) else 2L
  st <- .nm_ctl_struct(cpts, route, absorption, elimination, transit_n)
  dist <- if (cpts > 1) iv[[paste0("cpt", cpts)]] else NULL
  peri_q <- if (cpts == 2) "Q" else c("Q2", "Q3")
  peri_v <- if (cpts == 2) "V2" else c("V2", "V3")
  th <- stats::setNames(numeric(length(st$params)), st$params)
  for (p in st$params) {
    th[p] <- if (p == "CL") iv$CL else if (p == "VMAX") iv$VMAX else if (p == "KM") iv$KM
    else if (p == st$vc) { if (cpts == 1) iv$V_1cpt else dist$V1 }
    else if (p == "KA") { if (input == "first-order+lag") iv$KA_lag else iv$KA }
    else if (p == "DUR") iv$DUR else if (p == "MTT") iv$MTT
    else {
      j <- which(vapply(st$peri, function(z) z$q == p || z$vol == p, logical(1)))
      if (st$peri[[j]]$q == p) dist[[peri_q[j]]] else dist[[peri_v[j]]]
    }
  }
  th <- signif(th, 3)
  alag <- if (input == "first-order+lag") signif(iv$ALAG, 3) else NULL
  el_p <- if (elimination == "linear") "CL" else "VMAX"
  iiv_v <- stats::setNames(iiv$cv[c("CL", "Vc")], c(el_p, st$vc))
  ap <- intersect(st$params, c("KA", "DUR", "MTT"))
  if (length(ap)) iiv_v <- c(iiv_v, stats::setNames(iiv$cv[ap], ap))
  sig <- err$sigma[[error]]
  sig_full <- c(prop = unname(sig["prop"]), add = unname(sig["add"]))
  used <- switch(error, additive = "add", proportional = "prop", combined = c("prop", "add"))
  note <- character(0)
  if (!is.null(alag))
    note <- c(note, sprintf("ALAG1 = %s to add by hand (nm_write_ctl writes no lag time).", format(alag)))
  if (cpts > 1 && isTRUE(dist$heuristic))
    note <- c(note, "Distribution values are heuristic (macro constants not identifiable).")
  if (input == "zero-order-d1")
    note <- c(note, "Zero-order duration (ADVAN1/3/11): doses and observations in CMT 1 and RATE = -2 on doses; nm_suggest_run() adapts the dataset.")
  list(cpts = cpts, route = route, input = input, absorption = absorption, transit_n = transit_n,
       elimination = elimination, error = error, theta = th, alag = alag,
       iiv = iiv_v, sigma = signif(sig_full, 3), sigma_used = signif(sig_full[used], 3),
       note = paste(note, collapse = " "))
}

.nm_sm_call <- function(sp) {
  vec <- function(x) paste(sprintf("%s = %s", names(x), .nm_sm_fmt(x)), collapse = ", ")
  args <- sprintf('cpts = %d, route = "%s"', sp$cpts, sp$route)
  if (sp$route == "EV") args <- paste0(args, sprintf(', absorption = "%s"', sp$absorption))
  if (sp$input == "transit") args <- paste0(args, sprintf(", transit_n = %dL", sp$transit_n))
  args <- paste0(args, sprintf(', elimination = "%s", error = "%s"', sp$elimination, sp$error))
  args <- paste0(args, sprintf(", theta_init = c(%s)", vec(sp$theta)))
  args <- paste0(args, sprintf(", iiv = c(%s)", vec(sp$iiv)))
  args <- paste0(args, sprintf(", sigma_init = c(%s)", vec(sp$sigma)))
  sprintf("nm_write_ctl(%s)", args)
}

# ---------------------------------------------------------------------------
# Parameters likely to be hard to estimate
# ---------------------------------------------------------------------------
.nm_sm_hard_for <- function(sp, fl, des = NULL) {
  params <- names(sp$theta); cpts <- sp$cpts; input <- sp$input
  h <- list()
  add <- function(p, why) h[[p]] <<- c(h[[p]], why)
  if (sp$route == "EV") {
    for (p in intersect(params, c("KA", "DUR", "MTT"))) {
      if (fl$abs_poor) add(p, "fewer than 2 samples before the peak")
      else if (p %in% c("MTT", "DUR") && fl$abs_shape_poor) add(p, "fewer than 3 samples before the peak")
      if (fl$sparse) add(p, "sparse data: absorption informed by pooled data only")
    }
    if ("KA" %in% params && fl$flipflop) add("KA", "KA close to lambda_z: possible flip-flop, KA and CL/V poorly separated")
    if (input == "first-order+lag") {
      if (fl$lag_late) add("ALAG1", "first sample later than half of tmax")
      if (fl$abs_poor) add("ALAG1", "absorption phase poorly sampled")
    }
    if (input == "transit") add("transit_n", "number of transit compartments rarely identifiable: test 2 to 5, or estimate it (Savic 2007)")
  }
  if ("CL" %in% params && fl$terminal_short) add("CL", "terminal phase covers < 2 half-lives or AUC extrapolation > 20%")
  if (cpts >= 2) {
    st <- .nm_ctl_struct(cpts, sp$route, sp$absorption, sp$elimination, sp$transit_n)
    dist_p <- unlist(lapply(st$peri, function(z) c(z$q, z$vol)))
    for (p in dist_p) {
      if (fl$dist_poor) add(p, if (is.finite(fl$lam_ratio)) sprintf("distribution and terminal phases poorly separated (lambda ratio %.1f)", fl$lam_ratio) else "distribution phase not identified")
      if (fl$sparse) add(p, "sparse data: distribution phase informed by pooled data only")
    }
    if (cpts == 3) for (p in utils::tail(dist_p, 2)) add(p, "third compartment: requires rich and long sampling")
  }
  if (sp$elimination == "michaelis-menten") {
    for (p in c("VMAX", "KM")) {
      if (fl$mm_not_assessable) add(p, "single dose level: VMAX and KM not separately identifiable")
      else if (fl$p_mm < 0.5) add(p, "little evidence of saturation in the data: VMAX and KM strongly correlated")
    }
  }
  # expected precision for the design of the dataset
  if (!is.null(des)) {
    pr <- des$parameters
    for (i in which(is.finite(pr$rse_pct) & pr$rse_pct > 50))
      add(pr$parameter[i], sprintf("expected RSE %.0f%%", pr$rse_pct[i]))
    for (i in which(!is.finite(pr$rse_pct)))
      add(pr$parameter[i], "not identifiable with this design (singular information matrix)")
    sh <- des$shrinkage
    for (i in which(is.finite(sh$shrinkage_pct) & sh$shrinkage_pct > 30))
      add(sprintf("ETA(%s)", sh$parameter[i]),
          sprintf("expected eta-shrinkage %.0f%%: consider removing this ETA", sh$shrinkage_pct[i]))
  }
  if (fl$few_obs)
    add("IIV", sprintf("median %.0f quantified obs/subject: IIV beyond CL and V likely to show high shrinkage",
                       fl$median_obs_subject))
  vapply(h, function(x) paste(unique(x), collapse = ", "), "")
}
