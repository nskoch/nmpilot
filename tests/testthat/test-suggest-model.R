# Tests of nm_suggest_model() and related functions on simulated datasets
# with a known structure. Data are simulated with the analytic solution used
# by nm_suggest_design() (linear models, superposition of doses).

sim_spec <- function(cpts, route, input, theta, alag = NULL, transit_n = 2L) {
  absorption <- switch(input, `first-order` = "first-order", `first-order+lag` = "first-order",
                       `zero-order-d1` = "zero-order-d1", transit = "transit", "first-order")
  list(cpts = cpts, route = route, input = input, absorption = absorption,
       transit_n = transit_n, elimination = "linear", theta = theta, alag = alag)
}

# n subjects; doses: data.frame(TIME, AMT) applied to every subject (AMT may be
# scaled by subject through dose_levels); IIV 25% log-normal on all THETAs.
sim_data <- function(sp, n, times, dose_levels = 100, doses = data.frame(TIME = 0),
                     prop = 0.12, add = 0, iiv = 0.25, seed = 1, theta_fun = NULL,
                     extra_input = NULL, lloq = NULL) {
  set.seed(seed)
  st <- nmtools:::.nm_ctl_struct(sp$cpts, sp$route, sp$absorption, sp$elimination, sp$transit_n)
  rows <- list()
  for (i in seq_len(n)) {
    th <- sp$theta * exp(stats::rnorm(length(sp$theta), 0, iiv))
    if (!is.null(sp$alag)) th <- c(th, ALAG1 = sp$alag)
    amt <- dose_levels[(i - 1) %% length(dose_levels) + 1]
    dd <- data.frame(ID = i, TIME = doses$TIME, AMT = amt, TINF = 0)
    f <- if (is.null(theta_fun)) nmtools:::.nm_sd_predict(th, sp, st, dd, times, 1) else
      theta_fun(th, sp, st, dd, times, amt)
    dv <- f * (1 + stats::rnorm(length(f), 0, prop)) + stats::rnorm(length(f), 0, add)
    blq <- if (!is.null(lloq)) dv < lloq else dv <= 0
    dv <- ifelse(blq, 0, round(dv, 5))
    dcmt <- if (sp$route == "IV") 2 else 1
    rows[[i]] <- rbind(
      data.frame(ID = i, TIME = doses$TIME, AMT = amt, DV = ".", EVID = 1, MDV = 1, CMT = dcmt,
                 BLQ = 0, COV = i %% 2, RND = stats::rnorm(1)),
      data.frame(ID = i, TIME = times, AMT = ".", DV = dv, EVID = 0, MDV = 0, CMT = 2,
                 BLQ = as.integer(blq), COV = i %% 2, RND = 0))
  }
  d <- do.call(rbind, rows)
  d <- d[order(d$ID, d$TIME, -d$EVID), ]
  f <- tempfile(fileext = ".csv")
  utils::write.csv(d, f, row.names = FALSE, quote = FALSE)
  f
}

rich_t <- c(0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 24)

test_that("rich one-compartment oral data: structure, values and outputs", {
  f <- sim_data(sim_spec(1, "EV", "first-order", c(CL = 5, V = 50, KA = 1.2)), 16, rich_t)
  s <- nm_suggest_model(f, verbose = FALSE)
  expect_s3_class(s, "nm_model_suggestion")
  c1 <- s$candidates[1, ]
  expect_equal(c1$cpts, 1L)
  expect_equal(c1$input, "first-order")
  expect_equal(c1$elimination, "linear")
  expect_equal(c1$error, "proportional")
  expect_equal(sum(unlist(s$support$compartments)), 1, tolerance = 1e-8)
  th <- s$specs[[1]]$theta
  expect_equal(unname(th["CL"]), 5, tolerance = 0.3)
  expect_equal(unname(th["V"]), 50, tolerance = 0.3)
  expect_equal(unname(th["KA"]), 1.2, tolerance = 0.4)
  expect_equal(sqrt(unname(s$specs[[1]]$sigma_used["prop"])), 0.12, tolerance = 0.3)
  expect_true(all(c("CL", "V", "KA") %in% names(s$specs[[1]]$iiv)))
  expect_s3_class(s$design[[1]], "nm_suggest_design")
  expect_output(print(s), "Candidate models")
  for (ty in c("profiles", "support", "pooled", "design"))
    expect_s3_class(plot(s, type = ty), "ggplot")
  # every candidate gives a valid control stream
  for (i in seq_along(s$specs)) {
    ctl <- nm_write_ctl(from = s, candidate = i, verbose = FALSE)
    chk <- nm_ctl_check(ctl)
    expect_false(any(chk$level == "ERROR"), info = paste("candidate", i))
  }
  # the text call matches the from = form
  ctl_call <- eval(parse(text = sub("nm_write_ctl\\(", "nm_write_ctl(verbose = FALSE, input = s$data_summary$data_file, data = basename(s$data_summary$data_file), ",
                                    s$calls[[1]])))
  expect_identical(unclass(ctl_call), unclass(nm_write_ctl(from = s, candidate = 1, verbose = FALSE)))
})

test_that("two-compartment IV bolus data are recognised", {
  f <- sim_data(sim_spec(2, "IV", "IV", c(CL = 5, V1 = 20, Q = 15, V2 = 60)), 14,
                c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 24))
  s <- nm_suggest_model(f, verbose = FALSE, design_eval = FALSE)
  expect_equal(s$data_summary$route, "IV")
  expect_equal(s$candidates$cpts[1], 2L)
  expect_equal(s$candidates$input[1], "IV bolus")
  expect_equal(unname(s$specs[[1]]$theta["V1"]), 20, tolerance = 0.35)
})

test_that("lag time and zero-order input are recognised", {
  f <- sim_data(sim_spec(1, "EV", "first-order+lag", c(CL = 5, V = 50, KA = 1.5), alag = 0.6), 16,
                c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 24), lloq = 0.05)
  s <- nm_suggest_model(f, verbose = FALSE, design_eval = FALSE, lloq = 0.05)
  expect_equal(s$candidates$input[1], "first-order+lag")
  expect_equal(s$specs[[1]]$alag, 0.6, tolerance = 0.3)
  expect_message(nm_write_ctl(from = s, candidate = 1), "ALAG1")
  f2 <- sim_data(sim_spec(1, "EV", "zero-order-d1", c(CL = 5, V = 50, DUR = 2)), 16,
                 c(0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 24))
  s2 <- nm_suggest_model(f2, route = "EV", verbose = FALSE, design_eval = FALSE)
  expect_equal(s2$candidates$input[1], "zero-order-d1")
  expect_equal(unname(s2$specs[[1]]$theta["DUR"]), 2, tolerance = 0.25)
})

test_that("the residual error model follows the simulated error", {
  sp <- sim_spec(1, "EV", "first-order", c(CL = 5, V = 50, KA = 1.2))
  fa <- sim_data(sp, 16, c(rich_t, 36, 48), prop = 0, add = 0.08, seed = 2)
  sa <- nm_suggest_model(fa, verbose = FALSE, design_eval = FALSE)
  expect_gt(sa$support$error[["additive"]], 0.5)
  expect_equal(sqrt(sa$specs[[which(sa$candidates$error == "additive")[1]]]$sigma_used[["add"]]),
               0.08, tolerance = 0.35)
  fc <- sim_data(sp, 16, c(rich_t, 36, 48), prop = 0.15, add = 0.03, seed = 3)
  sc <- nm_suggest_model(fc, verbose = FALSE, design_eval = FALSE)
  expect_gt(sc$support$error[["combined"]], 0.5)
  expect_false(is.null(sc$evidence$error))
})

test_that("sparse data use the pooled analysis and flag poorly informed ETAs", {
  set.seed(5)
  sp <- sim_spec(1, "EV", "first-order", c(CL = 5, V = 50, KA = 1.2))
  times <- function() sort(round(stats::runif(3, 0.25, 24), 2))
  # different random times per subject: build the dataset subject by subject
  files <- lapply(1:120, function(i) utils::read.csv(sim_data(sp, 1, times(), seed = i)))
  d <- do.call(rbind, Map(function(x, i) { x$ID <- i; x }, files, seq_along(files)))
  f <- tempfile(fileext = ".csv"); utils::write.csv(d, f, row.names = FALSE, quote = FALSE)
  s <- nm_suggest_model(f, verbose = FALSE)
  expect_equal(s$evidence$source, "pooled")
  expect_equal(s$data_summary$design, "sparse")
  expect_equal(s$candidates$cpts[1], 1L)
  expect_true(any(grepl("pooled", s$notes)))
  sh <- s$design[[1]]$shrinkage
  expect_gt(sh$shrinkage_pct[sh$parameter == "KA"], sh$shrinkage_pct[sh$parameter == "CL"])
  expect_match(s$candidates$hard_to_estimate[1], "KA")
})

test_that("nm_suggest_design scales with the number of subjects and the design", {
  f <- sim_data(sim_spec(1, "EV", "first-order", c(CL = 5, V = 50, KA = 1.2)), 12, rich_t)
  s <- nm_suggest_model(f, verbose = FALSE, design_eval = FALSE)
  d1 <- nm_suggest_design(s, n_subjects = 12)
  d2 <- nm_suggest_design(s, n_subjects = 48)
  expect_equal(d2$parameters$rse_pct, d1$parameters$rse_pct / 2, tolerance = 0.02)
  rich <- nm_suggest_design(s, times = rich_t)
  sparse <- nm_suggest_design(s, times = c(6, 12, 24))
  ka <- function(d) d$parameters$rse_pct[d$parameters$parameter == "KA"]
  expect_gt(ka(sparse), ka(rich))
  expect_output(print(rich), "Expected precision")
  expect_s3_class(plot(rich), "ggplot")
  expect_error(nm_suggest_design(s, candidate = 99), "between")
})

test_that("macro constants of the analytic solution match the compartmental model", {
  st <- nmtools:::.nm_ctl_struct(2L, "IV", "first-order", "linear")
  m <- nmtools:::.nm_sd_macro(c(CL = 5, V1 = 20, Q = 15, V2 = 60), st)
  # sum of amplitudes = 1/V1; AUC of a unit dose = 1/CL
  expect_equal(sum(m$A), 1 / 20)
  expect_equal(sum(m$A / m$lam), 1 / 5)
})

test_that("alerts flag double peaks, time-varying CL and less-than-proportional exposure", {
  # double peak: half of the dose absorbed 5 h later
  sp_dp <- sim_spec(1, "EV", "first-order", c(CL = 10, V = 30, KA = 3))
  dp <- function(th, sp, st, dd, times, amt) {
    d2 <- rbind(transform(dd, AMT = amt / 2), transform(dd, AMT = amt / 2, TIME = dd$TIME + 5))
    nmtools:::.nm_sd_predict(th, sp, st, d2, times, 1)
  }
  f <- sim_data(sp_dp, 12, c(0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 5.5, 6, 7, 8, 10, 12, 24),
                prop = 0.05, iiv = 0.15, theta_fun = dp)
  s <- nm_suggest_model(f, verbose = FALSE, design_eval = FALSE)
  expect_true("double peak" %in% s$alerts$type)
  sp <- sim_spec(1, "EV", "first-order", c(CL = 5, V = 50, KA = 1.5))
  # steady-state CL 1.8-fold the first-dose CL (q12h, 8 doses)
  doses <- data.frame(TIME = seq(0, 84, by = 12))
  tv <- function(th, sp, st, dd, times, amt) {
    first <- times < 12
    out <- nmtools:::.nm_sd_predict(th, sp, st, dd, times, 1)
    th2 <- th; th2["CL"] <- th["CL"] * 1.8
    out[!first] <- nmtools:::.nm_sd_predict(th2, sp, st, dd, times[!first], 1)
    out
  }
  ft <- sim_data(sp, 10, c(0.5, 1, 2, 3, 4, 6, 8, 11.9, 84.5, 85, 86, 87, 88, 90, 92, 95.9),
                 doses = doses, prop = 0.05, theta_fun = tv)
  st <- nm_suggest_model(ft, verbose = FALSE, design_eval = FALSE)
  expect_true("time-varying clearance" %in% st$alerts$type)
  # exposure proportional to dose^0.6
  lp <- function(th, sp, st, dd, times, amt)
    nmtools:::.nm_sd_predict(th, sp, st, dd, times, 1) * (amt / 10)^-0.4
  fl <- sim_data(sp, 18, rich_t, dose_levels = c(10, 100, 1000), prop = 0.08, theta_fun = lp)
  sl <- nm_suggest_model(fl, verbose = FALSE, design_eval = FALSE)
  expect_true("less than dose-proportional" %in% sl$alerts$type)
  expect_lt(sl$support$elimination[["michaelis-menten"]], 0.2)
})

test_that("nm_write_ctl(from =) respects explicit arguments", {
  f <- sim_data(sim_spec(1, "EV", "first-order", c(CL = 5, V = 50, KA = 1.2)), 10, rich_t)
  s <- nm_suggest_model(f, verbose = FALSE, design_eval = FALSE)
  ctl <- nm_write_ctl(from = s, candidate = 1, cpts = 2, verbose = FALSE)
  expect_true(any(grepl("ADVAN4", ctl)))
  expect_true(any(grepl(sprintf("\\(0, %s\\)\\s*; CL", format(s$specs[[1]]$theta[["CL"]])), ctl)))
  ctl2 <- nm_write_ctl(from = s, candidate = 1, error = "additive", verbose = FALSE)
  expect_true(any(grepl("\\$SIGMA 0 FIX ; ERR\\(1\\) proportional", ctl2)))
  expect_error(nm_write_ctl(from = list()), "nm_model_suggestion")
  expect_error(nm_write_ctl(from = s, candidate = 50), "between")
})

test_that("nm_suggest_run writes, runs and ranks the candidates", {
  root <- tempfile("sr"); dir.create(root)
  f <- sim_data(sim_spec(1, "EV", "first-order+lag", c(CL = 5, V = 50, KA = 1.5), alag = 0.6), 12,
                c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 24), lloq = 0.05)
  file.copy(f, file.path(root, "data.csv"))
  s <- nm_suggest_model(file.path(root, "data.csv"), verbose = FALSE, design_eval = FALSE, lloq = 0.05)
  # fake launcher: the OFV is the number of lines of the control stream
  if (.Platform$OS.type == "windows") {
    nmfe <- file.path(root, "fake.bat")
    # line count in pure cmd: no external "find" (Rtools/Git may shadow it)
    writeLines(c("@echo off", "set N=0",
                 "for /f \"usebackq delims=\" %%l in (\"%~1\") do set /a N+=1",
                 "echo MINIMIZATION SUCCESSFUL> %2",
                 "echo #OBJV:  OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT: %N%>> %2"), nmfe)
  } else {
    nmfe <- file.path(root, "fake.sh")
    writeLines(c("#!/bin/sh", "N=$(wc -l < \"$1\")", "echo \"MINIMIZATION SUCCESSFUL\" > \"$2\"",
                 "echo \"#OBJV:  OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT: $N\" >> \"$2\""), nmfe)
    Sys.chmod(nmfe, "755")
  }
  r <- nm_suggest_run(s, root, nmfe, candidates = 1:3, verbose = FALSE)
  expect_s3_class(r, "nm_suggest_run")
  expect_equal(nrow(r), 3L)
  expect_true(all(file.exists(r$ctl)))
  expect_true(all(is.finite(r$ofv)))
  expect_equal(r$bic, r$ofv + r$n_par * log(s$data_summary$n_obs))
  expect_true(is.finite(attr(r, "spearman")))
  expect_output(print(r), "Spearman")
  # the lag-time candidate gets ALAG1 and stays consistent
  lag_i <- which(s$candidates$input[1:3] == "first-order+lag")
  if (length(lag_i)) {
    l <- readLines(r$ctl[lag_i[1]])
    expect_true(any(grepl("ALAG1 = THETA\\(", l)))
    expect_false(any(nm_ctl_check(l)$level == "ERROR"))
  }
  # a second call without overwrite stops; with overwrite, the runs are launched again
  expect_error(nm_suggest_run(s, root, nmfe, candidates = 1:3, verbose = FALSE), "exists")
  r_again <- nm_suggest_run(s, root, nmfe, candidates = 1:3, overwrite = TRUE, verbose = FALSE)
  expect_true(all(attr(r_again, "log")$action == "launched"))
  # written but not run
  r2 <- nm_suggest_run(s, root, candidates = 1, prefix = "dry", run = FALSE, verbose = FALSE)
  expect_true(is.na(r2$ofv))
  expect_error(nm_suggest_run(s, root, nmfe, prefix = "dry", candidates = 1, verbose = FALSE), "exists")
})

test_that("nm_suggest_covariates detects a covariate effect on CL", {
  ex <- system.file("extdata", "example", "example_data.csv", package = "nmtools")
  s <- nm_suggest_model(ex, verbose = FALSE, design_eval = FALSE, max_profiles = 60)
  cv <- nm_suggest_covariates(s, c("GENO", "WT"), n_perm = 500)
  expect_s3_class(cv, "nm_suggest_covariates")
  g <- cv[cv$parameter == "CL" & cv$covariate == "GENO", ]
  expect_lt(g$p_W, 0.05)
  expect_gt(g$ratio_median, 1.3)
  expect_output(print(cv), "Wasserstein")
  expect_error(nm_suggest_covariates(s, "NOPE"), "not found")
})

test_that("input checks", {
  expect_error(nm_suggest_model(tempfile()), "not found")
  bad <- data.frame(ID = 1, TIME = 0, DV = 1)
  expect_error(nm_suggest_model(bad, verbose = FALSE), "amt")
  expect_error(nm_suggest_model(bad, max_cpts = 4, verbose = FALSE), "max_cpts")
})

test_that("conc_factor writes the unit conversion in the scale", {
  ctl <- nm_write_ctl(conc_factor = 1000, verbose = FALSE)
  expect_true(any(grepl("^S2 = V/1000", ctl)))
  expect_true(any(grepl("^S2 = V ", nm_write_ctl(verbose = FALSE))))
  expect_true(any(grepl("^S1 = V1/1000", nm_write_ctl(route = "IV", cpts = 2, conc_factor = 1000, verbose = FALSE))))
  expect_error(nm_write_ctl(conc_factor = 0, verbose = FALSE), "positive")
})

test_that("dose_conc_factor gives the same volumes in L whatever the DV units", {
  f <- sim_data(sim_spec(1, "EV", "first-order", c(CL = 5, V = 50, KA = 1.2)), 12, rich_t)
  d <- utils::read.csv(f, colClasses = "character")
  d$DV[d$EVID == "0"] <- as.character(as.numeric(d$DV[d$EVID == "0"]) * 1000)   # mg/L -> ng/mL
  f_ng <- tempfile(fileext = ".csv"); utils::write.csv(d, f_ng, row.names = FALSE, quote = FALSE)
  s1 <- nm_suggest_model(f, verbose = FALSE, design_eval = FALSE)
  s2 <- nm_suggest_model(f_ng, verbose = FALSE, design_eval = FALSE, dose_conc_factor = 1000)
  expect_equal(s2$specs[[1]]$theta, s1$specs[[1]]$theta, tolerance = 0.01)
  ctl <- nm_write_ctl(from = s2, candidate = 1, verbose = FALSE)
  expect_true(any(grepl("^S2 = V/1000", ctl)))
  expect_true(any(grepl("^S2 = V ", nm_write_ctl(from = s1, candidate = 1, verbose = FALSE))))
  # SIGMA: the additive variance follows the DV units, the proportional one does not
  i <- which(s1$candidates$error == "combined")[1]
  j <- which(s2$candidates$error == "combined")[1]
  if (!is.na(i) && !is.na(j)) {
    expect_equal(s2$specs[[j]]$sigma[["prop"]], s1$specs[[i]]$sigma[["prop"]], tolerance = 0.05)
    expect_equal(s2$specs[[j]]$sigma[["add"]] / 1e6, s1$specs[[i]]$sigma[["add"]], tolerance = 0.05)
  }
})

test_that("datasets are adapted to the compartment numbering of each candidate", {
  root <- tempfile("prep"); dir.create(root)
  src <- file.path(root, "d.csv")
  writeLines(c("#ID,TIME,AMT,DV,EVID,MDV,CMT", "# comment line,,,,,,",
               "1,0,100,.,1,1,1", "1,1,.,5.2,0,0,2", "1,4,.,3.1,0,0,2"), src)
  spec <- function(cpts, route, input, absorption) list(cpts = cpts, route = route, input = input,
    absorption = absorption, transit_n = 2L, elimination = "linear")
  # first-order: depot 1, central 2 -> original file
  expect_identical(nmtools:::.nm_sr_prepare_data(src, spec(1, "EV", "first-order", "first-order"),
                                                 root, "a", verbose = FALSE), src)
  expect_identical(nmtools:::.nm_sr_prepare_data(src, spec(2, "EV", "transit", "transit"),
                                                 root, "b", verbose = FALSE), src)
  # zero-order duration: everything in CMT 1, RATE = -2 on doses, rest untouched
  out <- nmtools:::.nm_sr_prepare_data(src, spec(1, "EV", "zero-order-d1", "zero-order-d1"),
                                       root, "c", verbose = FALSE)
  expect_false(identical(out, src))
  l <- readLines(out)
  expect_equal(l[1], "#ID,TIME,AMT,RATE,DV,EVID,MDV,CMT")
  expect_equal(l[2], "# comment line,,,.,,,,")
  expect_equal(l[3:5], c("1,0,100,-2,.,1,1,1", "1,1,.,0,5.2,0,0,1", "1,4,.,0,3.1,0,0,1"))
  expect_identical(readLines(src)[3], "1,0,100,.,1,1,1")   # original unchanged
  # IV with observations coded CMT 2: moved to the central compartment
  out_iv <- nmtools:::.nm_sr_prepare_data(src, spec(2, "IV", "IV bolus", "first-order"),
                                          root, "d", verbose = FALSE)
  expect_equal(readLines(out_iv)[4], "1,1,.,5.2,0,0,1")
  # existing RATE column: -2 set on doses only
  src2 <- file.path(root, "r.csv")
  writeLines(c("ID,TIME,AMT,RATE,DV,EVID,MDV,CMT", "1,0,100,0,.,1,1,1", "1,1,.,0,5.2,0,0,1"), src2)
  l2 <- readLines(nmtools:::.nm_sr_prepare_data(src2, spec(1, "EV", "zero-order-d1", "zero-order-d1"),
                                                root, "e", verbose = FALSE))
  expect_equal(l2[2:3], c("1,0,100,-2,.,1,1,1", "1,1,.,0,5.2,0,0,1"))
  # several observed compartments or several routes: explicit stop
  src3 <- file.path(root, "m.csv")
  writeLines(c("ID,TIME,AMT,DV,EVID,MDV,CMT", "1,0,100,.,1,1,1", "1,1,.,5,0,0,2", "1,1,.,2,0,0,3"), src3)
  expect_error(nmtools:::.nm_sr_prepare_data(src3, spec(1, "EV", "zero-order-d1", "zero-order-d1"),
                                             root, "f", verbose = FALSE), "observed compartments")
  src4 <- file.path(root, "v.csv")
  writeLines(c("ID,TIME,AMT,DV,EVID,MDV,CMT", "1,0,100,.,1,1,1", "1,0,50,.,1,1,2", "1,1,.,5,0,0,2"), src4)
  expect_error(nmtools:::.nm_sr_prepare_data(src4, spec(1, "EV", "zero-order-d1", "zero-order-d1"),
                                             root, "g", verbose = FALSE), "route")
})

test_that("nm_suggest_run runs a zero-order candidate on an adapted dataset", {
  root <- tempfile("srz"); dir.create(root)
  f <- sim_data(sim_spec(1, "EV", "first-order", c(CL = 5, V = 50, KA = 1.2)), 10, rich_t)
  file.copy(f, file.path(root, "data.csv"))
  s <- nm_suggest_model(file.path(root, "data.csv"), verbose = FALSE, design_eval = FALSE,
                        n_candidates = 24)
  i_zo <- which(s$candidates$input == "zero-order-d1")[1]
  i_fo <- which(s$candidates$input == "first-order")[1]
  skip_if(is.na(i_zo))
  r <- nm_suggest_run(s, root, candidates = c(i_fo, i_zo), run = FALSE, verbose = FALSE)
  expect_equal(normalizePath(r$data[1]), normalizePath(file.path(root, "data.csv")))
  expect_false(identical(normalizePath(r$data[2]), normalizePath(file.path(root, "data.csv"))))
  l <- readLines(r$ctl[2])
  expect_true(any(grepl("^\\$INPUT .*RATE", l)))
  expect_true(any(grepl(basename(r$data[2]), l, fixed = TRUE)))
  expect_false(any(nm_ctl_check(l)$level == "ERROR"))
  d <- utils::read.csv(r$data[2])
  expect_true(all(d$CMT == 1))
  expect_true(all(d$RATE[d$EVID == 1] == -2))
})
