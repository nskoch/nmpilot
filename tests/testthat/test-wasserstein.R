# Independent W1: integral of |Qx(u) - Qy(u)| over (0, 1), the quantile form
# used by transport::wasserstein1d for p = 1.
w1_quantile <- function(x, y) {
  x <- sort(x); y <- sort(y); n <- length(x); m <- length(y)
  u <- sort(unique(c(0, seq_len(n) / n, seq_len(m) / m)))
  mid <- (utils::head(u, -1) + utils::tail(u, -1)) / 2
  sum(abs(x[pmin(n, ceiling(mid * n))] - y[pmin(m, ceiling(mid * m))]) * diff(u))
}

wass_study <- function() {
  set.seed(11)
  root <- tempfile("wass")
  dir.create(file.path(root, "base.reg"), recursive = TRUE)
  dir.create(file.path(root, "final.reg"))
  n <- 60
  cov <- data.frame(ID = 1:n, CLCR = round(stats::rnorm(n, 90, 30)),
                    CT01 = stats::rbinom(n, 1, 0.4))
  eta_base <- 0.006 * (cov$CLCR - 90) + 0.5 * cov$CT01 + stats::rnorm(n, 0, 0.2)
  eta_final <- stats::rnorm(n, 0, 0.2)
  write_phi <- function(e, f) writeLines(c(
    "TABLE NO.     1: First Order Conditional Estimation with Interaction",
    " SUBJECT_NO   ID           ETA(1)       ETA(2)       OBJ",
    sprintf("%12d%13.4E%13.5E%13.5E%13.5E", 1:n, cov$ID, e, stats::rnorm(n, 0, 0.1), 10)), f)
  write_phi(eta_base, file.path(root, "base.reg", "base.phi"))
  write_phi(eta_final, file.path(root, "final.reg", "final.phi"))
  # .ext with OMEGA(1,1) for the shrinkage, and a .ctl naming ETA(1) as CL
  for (r in c("base", "final")) {
    writeLines(c("TABLE NO.     1: First Order",
                 " ITERATION    THETA1       OMEGA(1,1)   OMEGA(2,1)   OMEGA(2,2)   OBJ",
                 ext_row(-1000000000, c(5, 0.09, 0, 0.01, 100))),
               file.path(root, paste0(r, ".reg"), paste0(r, ".ext")))
    writeLines(c("$INPUT ID TIME DV CLCR CT01", "$DATA ../cov.csv IGNORE=@",
                 "$THETA (0,5)", "$OMEGA 0.09", " 0.01", "$PK",
                 "  CL = THETA(1)*EXP(ETA(1))", "  V = 50*EXP(ETA(2))"),
               file.path(root, paste0(r, ".reg"), paste0(r, ".ctl")))
  }
  utils::write.csv(cov, file.path(root, "cov.csv"), row.names = FALSE)
  list(root = root, cov = cov, eta_base = eta_base, eta_final = eta_final)
}

test_that("W1 matches exact formulas", {
  set.seed(3)
  for (i in 1:50) {
    x <- stats::rnorm(sample(3:30, 1)); y <- stats::rexp(sample(3:30, 1))
    expect_equal(nmtools:::.nm_w1(x, y), w1_quantile(x, y), tolerance = 1e-12)
  }
  x <- stats::rnorm(20); y <- stats::rnorm(20, 1)
  expect_equal(nmtools:::.nm_w1(x, y), mean(abs(sort(x) - sort(y))), tolerance = 1e-12)
  expect_equal(nmtools:::.nm_w1_to_zero(x), nmtools:::.nm_w1(x, rep(0, 20)), tolerance = 1e-12)
  expect_equal(nmtools:::.nm_w1(x, x), 0)
})

test_that("W1 matches transport::wasserstein1d when available", {
  skip_if_not_installed("transport")
  set.seed(4)
  x <- stats::rnorm(17); y <- stats::rgamma(23, 2)
  expect_equal(nmtools:::.nm_w1(x, y), transport::wasserstein1d(x, y, p = 1), tolerance = 1e-10)
})

test_that("the permutation test is reproducible and keeps the user's random state", {
  set.seed(99); before <- stats::runif(1); set.seed(99)
  a <- nmtools:::.nm_w_perm(1:10 + 0.5, 1:12, n_perm = 200, seed = 1)
  after <- stats::runif(1)
  expect_identical(before, after)
  b <- nmtools:::.nm_w_perm(1:10 + 0.5, 1:12, n_perm = 200, seed = 1)
  expect_identical(a, b)
  expect_gte(a$p, 1 / 201); expect_lte(a$p, 1)
})

test_that("nm_wasserstein reproduces the method of the published script", {
  s <- wass_study()
  w <- nm_wasserstein(c(BASE = file.path(s$root, "base.reg"), FINAL = file.path(s$root, "final.reg")),
                      covariates = c("CLCR", "CT01"), eta = "ETA(1)",
                      data = file.path(s$root, "cov.csv"), n_perm = 500)
  expect_s3_class(w, "nm_wasserstein")
  expect_equal(nrow(w$tests), 4L)

  # Groups as in the script: CT01 1 vs 0; CLCR Q4 vs Q1 (quartiles on subjects)
  q <- dplyr::ntile(s$cov$CLCR, 4)
  b <- w$tests[w$tests$run == "BASE", ]
  x1 <- s$eta_base[s$cov$CT01 == 1]; x0 <- s$eta_base[s$cov$CT01 == 0]
  expect_equal(b$W[b$covariate == "CT01"], w1_quantile(x1, x0), tolerance = 1e-4)
  expect_equal(b$t_stat[b$covariate == "CT01"], unname(stats::t.test(x1, x0)$statistic), tolerance = 1e-3)
  x4 <- s$eta_base[q == 4]; xq1 <- s$eta_base[q == 1]
  expect_equal(b$W[b$covariate == "CLCR"], w1_quantile(x4, xq1), tolerance = 1e-4)
  expect_equal(b$comparison[b$covariate == "CLCR"], "Q4 vs Q1")

  # Strong residual signal in BASE, none in FINAL
  expect_lt(max(b$p_W), 0.01)

  # R2W = 1 - mean|ETA_final| / mean|ETA_base|
  expect_equal(w$r2w$R2W[1], 0)
  expect_equal(w$r2w$R2W[2], 1 - mean(abs(s$eta_final)) / mean(abs(s$eta_base)), tolerance = 1e-4)
  # Shrinkage from OMEGA(1,1) = 0.09
  expect_equal(w$r2w$shrinkage_pct[2], 100 * (1 - stats::sd(s$eta_final) / 0.3), tolerance = 1e-3)
})

test_that("ETA by parameter name, data from $DATA, labels from run names", {
  s <- wass_study()
  w <- nm_wasserstein(file.path(s$root, c("base.reg", "final.reg")),
                      covariates = "CT01", eta = "CL", n_perm = 100)
  expect_equal(unique(w$r2w$eta), "ETA(1)")
  expect_equal(w$r2w$run, c("base", "final"))
  w2 <- nm_wasserstein(file.path(s$root, c("base.reg", "final.reg")),
                       covariates = "CT01", eta = "V", n_perm = 100, ref = "final")
  expect_equal(unique(w2$r2w$eta), "ETA(2)")
  expect_equal(w2$r2w$R2W[2], 0)
})

test_that("print and plot methods", {
  s <- wass_study()
  w <- nm_wasserstein(file.path(s$root, c("base.reg", "final.reg")),
                      covariates = c("CLCR", "CT01"), eta = "CL", n_perm = 100)
  expect_output(print(w), "Wasserstein")
  expect_s3_class(plot(w), "ggplot")
  expect_s3_class(plot(w, covariate = "CT01"), "ggplot")
  expect_s3_class(plot(w, type = "r2w"), "ggplot")
})

test_that("nm_wasserstein checks its inputs", {
  s <- wass_study()
  runs <- file.path(s$root, c("base.reg", "final.reg"))
  expect_error(nm_wasserstein(runs, "AGE", n_perm = 10), "not found in the data")
  expect_error(nm_wasserstein(runs, "CT01", eta = "KA", n_perm = 10), "No ETA is linked")
  expect_error(nm_wasserstein(runs, "CT01", n_perm = 0), "n_perm")
  expect_error(nm_wasserstein(file.path(s$root, "none"), "CT01"), "Run not found")
  expect_error(nm_wasserstein(runs, "CT01", ref = "X", n_perm = 10), "ref")
})

test_that("pairwise quartile comparisons with Bonferroni adjustment", {
  s <- wass_study()
  runs <- file.path(s$root, c("base.reg", "final.reg"))
  w <- nm_wasserstein(runs, covariates = c("CLCR", "CT01"), eta = "CL",
                      continuous = "pairwise", n_perm = 200)
  cl <- w$tests[w$tests$covariate == "CLCR" & w$tests$run == "base", ]
  expect_equal(cl$comparison, c("Q2 vs Q1", "Q3 vs Q1", "Q4 vs Q1",
                                "Q3 vs Q2", "Q4 vs Q2", "Q4 vs Q3"))
  expect_equal(cl$p_W_adj, pmin(1, cl$p_W * 6))
  expect_equal(cl$p_t_adj, pmin(1, cl$p_t * 6))
  ct <- w$tests[w$tests$covariate == "CT01", ]
  expect_equal(ct$p_W_adj, ct$p_W)                # single comparison: no adjustment
  # Q4 vs Q1 identical to the default analysis
  q <- dplyr::ntile(s$cov$CLCR, 4)
  expect_equal(cl$W[cl$comparison == "Q4 vs Q1"],
               w1_quantile(s$eta_base[q == 4], s$eta_base[q == 1]), tolerance = 1e-4)
  expect_s3_class(plot(w, covariate = "CLCR"), "ggplot")
  expect_output(print(w), "Bonferroni")
  expect_error(nm_wasserstein(runs, "CLCR", continuous = "deciles", n_perm = 10))
})
