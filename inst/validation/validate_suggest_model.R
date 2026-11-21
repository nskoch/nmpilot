# =============================================================================
# Validation of nm_suggest_model() by simulation
#
# For each scenario (known structure) and design (rich or sparse), datasets
# are simulated and nm_suggest_model() is run. The script reports how often
# the top-ranked candidate has the true number of compartments, input,
# elimination and residual error model, and how often the true model is
# among the first 3 candidates.
#
# Usage (from R, nmpilot installed):
#   source(system.file("validation", "validate_suggest_model.R", package = "nmpilot"))
#   res <- validate_suggest_model(n_rep = 20, out_dir = "validation_results")
#   res$summary
#
# Linear models are simulated with the analytic solution used by
# nm_suggest_design(); Michaelis-Menten elimination with a fixed-step
# Runge-Kutta (RK4) solver. IIV: log-normal, CV about 30%. Rich design:
# 24 subjects, 11 to 14 samples over 24-48 h. Sparse design: 120 subjects,
# 3 samples at random times.
# =============================================================================

.val_sim_mm <- function(th, dose, times, dt = 0.01) {
  # one compartment, first-order absorption, Michaelis-Menten elimination
  a <- c(dose, 0); out <- numeric(length(times)); tg <- seq(0, max(times), by = dt)
  der <- function(a) c(-th[["KA"]] * a[1],
                       th[["KA"]] * a[1] - th[["VMAX"]] * (a[2] / th[["V"]]) / (th[["KM"]] + a[2] / th[["V"]]))
  for (t in tg) {
    m <- which(abs(times - t) < dt / 2)
    if (length(m)) out[m] <- a[2] / th[["V"]]
    k1 <- der(a); k2 <- der(a + dt / 2 * k1); k3 <- der(a + dt / 2 * k2); k4 <- der(a + dt * k3)
    a <- a + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
  }
  out
}

.val_scenarios <- function() {
  list(
    S1 = list(label = "1 cpt, first-order", cpts = 1, route = "EV", input = "first-order",
              elimination = "linear", error = "proportional",
              theta = c(CL = 5, V = 50, KA = 1.2), doses = 100),
    S2 = list(label = "2 cpt, first-order + lag", cpts = 2, route = "EV", input = "first-order+lag",
              elimination = "linear", error = "proportional",
              theta = c(CL = 5, V2 = 20, Q = 10, V3 = 80, KA = 1.5), alag = 0.6, doses = c(50, 100, 200)),
    S3 = list(label = "2 cpt, IV bolus", cpts = 2, route = "IV", input = "IV bolus",
              elimination = "linear", error = "proportional",
              theta = c(CL = 5, V1 = 20, Q = 15, V2 = 60), doses = 100),
    S4 = list(label = "1 cpt, transit", cpts = 1, route = "EV", input = "transit",
              elimination = "linear", error = "proportional",
              theta = c(CL = 5, V = 50, MTT = 2), transit_n = 4L, doses = 100),
    S5 = list(label = "1 cpt, zero-order", cpts = 1, route = "EV", input = "zero-order-d1",
              elimination = "linear", error = "proportional",
              theta = c(CL = 5, V = 50, DUR = 2), doses = 100),
    S6 = list(label = "1 cpt, Michaelis-Menten", cpts = 1, route = "EV", input = "first-order",
              elimination = "michaelis-menten", error = "proportional",
              theta = c(VMAX = 20, KM = 2, V = 20, KA = 1.5), doses = c(10, 100, 1000)),
    S7 = list(label = "1 cpt, additive error", cpts = 1, route = "EV", input = "first-order",
              elimination = "linear", error = "additive",
              theta = c(CL = 5, V = 50, KA = 1.2), doses = 100),
    S8 = list(label = "1 cpt, combined error", cpts = 1, route = "EV", input = "first-order",
              elimination = "linear", error = "combined",
              theta = c(CL = 5, V = 50, KA = 1.2), doses = 100)
  )
}

.val_simulate <- function(sc, design, seed) {
  set.seed(seed)
  n <- if (design == "rich") 24 else 120
  rich_times <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 24, 36, 48)
  if (sc$route == "IV") rich_times <- c(0.1, rich_times[-1])
  absorption <- switch(sc$input, `first-order+lag` = "first-order", `IV bolus` = "first-order", sc$input)
  sp <- list(cpts = sc$cpts, route = sc$route, input = if (sc$route == "IV") "IV" else sc$input,
             absorption = absorption, transit_n = if (is.null(sc$transit_n)) 2L else sc$transit_n, elimination = "linear")
  st <- if (sc$elimination == "linear")
    nmpilot:::.nm_ctl_struct(sc$cpts, sc$route, absorption, "linear", sp$transit_n) else NULL
  err <- switch(sc$error, proportional = c(0.15, 0), additive = c(0, 0.08), combined = c(0.15, 0.03))
  rows <- list()
  for (i in seq_len(n)) {
    th <- sc$theta * exp(stats::rnorm(length(sc$theta), 0, 0.3))
    if (!is.null(sc$alag)) th <- c(th, ALAG1 = sc$alag)
    amt <- sc$doses[(i - 1) %% length(sc$doses) + 1]
    times <- if (design == "rich") rich_times else sort(round(stats::runif(3, 0.25, 24), 2))
    f <- if (sc$elimination == "linear")
      nmpilot:::.nm_sd_predict(th, sp, st, data.frame(ID = i, TIME = 0, AMT = amt, TINF = 0), times, 1) else
      .val_sim_mm(th, amt, times)
    dv <- f * (1 + stats::rnorm(length(f), 0, err[1])) + stats::rnorm(length(f), 0, err[2])
    blq <- dv < 0.02
    rows[[i]] <- rbind(
      data.frame(ID = i, TIME = 0, AMT = amt, DV = ".", EVID = 1, MDV = 1,
                 CMT = if (sc$route == "IV") 2 else 1, BLQ = 0),
      data.frame(ID = i, TIME = times, AMT = ".", DV = ifelse(blq, 0, round(dv, 5)), EVID = 0,
                 MDV = 0, CMT = 2, BLQ = as.integer(blq)))
  }
  f <- tempfile(fileext = ".csv")
  utils::write.csv(do.call(rbind, rows), f, row.names = FALSE, quote = FALSE)
  f
}

validate_suggest_model <- function(n_rep = 20, designs = c("rich", "sparse"),
                                   scenarios = names(.val_scenarios()), out_dir = NULL,
                                   seed = 1, verbose = TRUE) {
  scs <- .val_scenarios()[scenarios]
  res <- list()
  for (sn in names(scs)) for (design in designs) for (r in seq_len(n_rep)) {
    sc <- scs[[sn]]
    f <- .val_simulate(sc, design, seed = seed + 1000 * r + match(sn, names(.val_scenarios())))
    t0 <- Sys.time()
    s <- tryCatch(nm_suggest_model(f, verbose = FALSE, design_eval = FALSE, lloq = 0.02,
                                   route = if (sc$input == "zero-order-d1") "EV" else "auto"),
                  error = function(e) NULL)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    unlink(f)
    if (is.null(s)) {
      res[[length(res) + 1]] <- data.frame(scenario = sn, label = sc$label, design = design, rep = r,
                                           ok_run = FALSE, cpts = NA, input = NA, elimination = NA,
                                           error = NA, all = NA, top3 = NA, seconds = el)
      next
    }
    c1 <- s$candidates[1, ]
    top3 <- utils::head(s$candidates, 3)
    is_true <- function(d) d$cpts == sc$cpts & d$input == sc$input & d$elimination == sc$elimination &
      d$error == sc$error
    res[[length(res) + 1]] <- data.frame(
      scenario = sn, label = sc$label, design = design, rep = r, ok_run = TRUE,
      cpts = c1$cpts == sc$cpts, input = c1$input == sc$input,
      elimination = c1$elimination == sc$elimination, error = c1$error == sc$error,
      all = is_true(c1), top3 = any(is_true(top3)), seconds = el)
    if (verbose) message(sprintf("%s %s rep %d: %s (%.0f s)", sn, design, r,
                                 if (is_true(c1)) "correct" else "different", el))
  }
  res <- do.call(rbind, res)
  summ <- stats::aggregate(cbind(cpts, input, elimination, error, all, top3) ~ scenario + label + design,
                           data = res[res$ok_run, ], FUN = function(v) round(100 * mean(v)))
  summ <- summ[order(summ$design, summ$scenario), ]
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(res, file.path(out_dir, "validation_runs.csv"), row.names = FALSE)
    utils::write.csv(summ, file.path(out_dir, "validation_summary.csv"), row.names = FALSE)
  }
  list(runs = res, summary = summ)
}
