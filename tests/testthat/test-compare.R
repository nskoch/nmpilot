compare_dirs <- function() {
  root <- tempfile("cmp")
  a <- file.path(root, "runA"); b <- file.path(root, "runB")
  dir.create(a, recursive = TRUE); dir.create(b)
  write_sdtab(file.path(a, "sdtab1"))
  write_sdtab(file.path(b, "sdtab1"), shift = 0.1)
  c(a, b)
}

test_that("nm_compare pairs observations and computes deltas", {
  d <- compare_dirs()
  cmp <- nm_compare(d[1], d[2], labels = c("A", "B"))
  expect_s3_class(cmp, "nm_compare")
  expect_named(cmp, c("data_obs", "data_subj", "metrics", "meta", "tests", "parameters",
                      "etas", "covariate_tests", "covariates", "run_info"), ignore.order = TRUE)
  expect_null(cmp$tests)          # no .phi in these fixtures
  expect_equal(nrow(cmp$data_obs), 6L)
  expect_equal(cmp$data_obs$dCWRES, rep(0.1, 6), tolerance = 1e-6)
  expect_equal(nrow(cmp$data_subj), 2L)
  m <- cmp$metrics
  expect_equal(m$run1[m$metric == "Mean |CWRES|"], 0.8 / 3, tolerance = 1e-6)
})

test_that("plot.nm_compare returns plots", {
  d <- compare_dirs()
  cmp <- nm_compare(d[1], d[2])
  p <- plot(cmp)
  expect_true(inherits(p, "ggplot") || is.list(p))
  expect_output(print(cmp))
})

test_that("Vuong and Clarke tests on individual contributions", {
  set.seed(1)
  subj <- data.frame(ID = 1:40, iOFV_1 = rnorm(40, 10, 3))
  subj$iOFV_2 <- subj$iOFV_1 - rnorm(40, 1, 2)          # run 2 better on average
  t <- nmpilot:::.nm_cmp_tests(subj, k = c(5, 6), ofv = c(sum(subj$iOFV_1), sum(subj$iOFV_2)),
                               n_obs = 200, nested = TRUE, n_boot = 500, seed = 1,
                               labels = c("A", "B"))
  d <- (subj$iOFV_1 - subj$iOFV_2) / 2
  z_aic <- (sum(d) - 1) / (sqrt(40) * sqrt(mean((d - mean(d))^2)))
  v <- t$table[t$table$test == "Vuong" & t$table$correction == "AIC", ]
  expect_equal(v$statistic, round(z_aic, 3))
  cl <- t$table[t$table$test == "Clarke" & t$table$correction == "BIC", ]
  expect_equal(cl$statistic, sum(d - log(40) / 2 / 40 > 0))
  expect_true("LRT (nested)" %in% t$table$test)
  expect_equal(t$bootstrap$n_boot, 500)
  expect_true(t$bootstrap$p_run2_better > 0.9)
  # influence: removing the subjects with the largest gains reverses dOFV
  di <- subj$iOFV_2 - subj$iOFV_1
  n <- t$influence$n_reverse
  s <- sort(di)
  expect_gte(sum(di) - sum(s[seq_len(n)]), 0)
  expect_lt(sum(di) - sum(s[seq_len(n - 1)]), 0)
})

test_that("nm_compare on the example runs (base vs final)", {
  ex <- system.file("extdata", "example", package = "nmpilot")
  cmp <- nm_compare(file.path(ex, "run1.reg"), file.path(ex, "run2.reg"),
                    labels = c("base", "final"), nested = TRUE,
                    covariates = c("WT", "GENO"), n_boot = 300)
  expect_equal(cmp$meta$n_par, c(8, 10))
  tt <- cmp$tests$table
  expect_true(all(tt$favours[tt$test %in% c("Vuong", "Clarke")] == "final"))
  expect_equal(nrow(cmp$covariate_tests), 2L)
  p <- cmp$parameters
  expect_true(all(c("CL_GENO", "CL_WT") %in% p$name))
  expect_lt(p$ratio[p$type == "OMEGA" & p$name == "CL"], 0.6)
  expect_setequal(unique(cmp$etas$parameter), c("CL", "V", "KA"))
  expect_true(is.finite(cmp$run_info[[1]]$cond) && cmp$run_info[[1]]$cond > 1)
  s <- summary(cmp)
  expect_s3_class(s, "nm_compare_summary")
  expect_equal(s$favours[s$criterion == "AIC"], "final")
  expect_output(print(s), "Decision criteria")
  expect_output(print(cmp), "Vuong")
  pl <- plot(cmp, which = c("gain_map", "overlay_gof", "bootstrap", "parameters",
                            "eta_scatter", "iofv_covariate"), n_band = 20)
  expect_true(all(vapply(pl, inherits, logical(1), "ggplot")))
  expect_true(all(c("gain_map", "overlay_gof", "iofv_covariate_WT") %in% names(pl)))
})
