# nm_optimize_base / nm_derivatives_ctl. NONMEM is not needed: the table of
# derivatives of run1 is emulated (helper-optimize.R).

ex_copy <- function() {
  d <- file.path(tempfile("optb_"))
  dir.create(d)
  file.copy(system.file("extdata", "example", package = "nmtools"), d, recursive = TRUE)
  file.path(d, "example")
}

test_that("nm_derivatives_ctl writes the evaluation run", {
  ex <- ex_copy()
  expect_warning(p <- nm_derivatives_ctl(file.path(ex, "run1.reg"), extra = c("WT", "NOTTHERE")),
                 "NOTTHERE")
  expect_equal(normalizePath(dirname(p)), normalizePath(ex))
  x <- readLines(p)
  expect_true(any(grepl("^\\$ESTIMATION METHOD=1 INTERACTION MAXEVAL=0", x)))
  expect_true(any(x == "\"  G_ETA2 = G(2,1)"))
  expect_true(any(x == "\"  D_EPS1_ETA3 = HH(1,4)"))
  expect_equal(sum(grepl("^\\$THETA .* FIX$", x)), 5L)
  expect_equal(sum(grepl("^\\$OMEGA .* FIX$", x)), 3L)
  expect_false(any(grepl("^\\$COV", x)))
  expect_true(any(grepl("CIPREDI", x)) && any(grepl("FILE=run1_deriv.dtab", x)))
  expect_true(any(grepl("^\\$TABLE.* WT ", x)))
  ctl <- file.path(ex, "run1.reg", "run1.ctl")
  writeLines(sub("Y = IPRED", "F_FLAG = 0\n  Y = IPRED", readLines(ctl)), ctl)
  expect_error(nm_derivatives_ctl(file.path(ex, "run1.reg")), "F_FLAG")
})

test_that("OMEGA blocks are written back", {
  m <- matrix(c(0.1, 0.02, 0, 0.02, 0.2, 0, 0, 0, 0), 3)
  w <- nmtools:::.nm_ob_write_random("OMEGA", m)
  expect_equal(w[1], "$OMEGA BLOCK(2) FIX")
  expect_equal(w[4], "$OMEGA 0 FIX")
})

test_that("level 1 diagnostics find the simulated covariates", {
  ex <- ex_copy()
  o <- nm_optimize_base(file.path(ex, "run1.reg"), covariates = c("WT", "GENO"),
                        level = 1, verbose = FALSE)
  expect_s3_class(o, "nm_optimize_base")
  expect_null(o$level2)
  cv <- o$level1$covariates
  expect_setequal(paste(cv$covariate, cv$eta)[cv$suggestion != "-"],
                  c("WT ETA2 (V)", "GENO ETA1 (CL)"))
  expect_equal(nrow(o$level1$residual), 4L)
  expect_true(all(o$level1$residual$dOFV > -10))   # true residual model = combined, as in run1
  expect_output(print(o), "Level 1")
})

test_that("linearized model reproduces the NONMEM OFV and detects a covariate", {
  ex <- ex_copy()
  rd <- file.path(ex, "run1.reg")
  .emul_deriv(rd, file.path(ex, "example_data.csv"), file.path(rd, "run1_deriv.dtab"))
  fin <- nmtools:::.nm_vpc_read_ext_final(file.path(rd, "run1.ext"))
  om <- nmtools:::.nm_ob_ext_matrix(fin, "OMEGA"); sg <- nmtools:::.nm_ob_ext_matrix(fin, "SIGMA")
  dt <- nmtools:::.nm_read_table(file.path(rd, "run1_deriv.dtab"))
  cov <- nmtools:::.nm_w_covariates(file.path(ex, "example_data.csv"), NA)
  info <- nmtools:::.nm_ob_cov_info(cov, "GENO", 5L)
  des <- nmtools:::.nm_ob_cov_design(cov, unique(dt$ID), info, "power")
  prep <- nmtools:::.nm_ob_prepare(dt, om, sg, des$X, NULL)
  sp0 <- nmtools:::.nm_ob_spec(list(1L, 2L, 3L), "base", design = des)
  st0 <- nmtools:::.nm_ob_start(sp0, prep)
  # FOCE-I objective of the linearized model = NONMEM OFV at the estimates
  ofv0 <- nmtools:::.nm_ob_ofv(st0, sp0, prep)
  expect_lt(abs(ofv0 - (-1147.341)), 0.05)
  # fast approximation: a few units away
  expect_lt(abs(nmtools:::.nm_ob_ofv_fast(st0, sp0, prep) - (-1147.341)), 10)
  f0 <- nmtools:::.nm_ob_fit(sp0, prep)
  sp1 <- nmtools:::.nm_ob_spec(list(1L, 2L, 3L), "base",
                               beta = list(list(eta = 1L, covariate = "GENO")), design = des)
  f1 <- nmtools:::.nm_ob_fit(sp1, prep, f0$u)
  expect_lt(f1$ofv - f0$ofv, -20)
  expect_gt(f1$u$beta[1, 1], log(1.2))
  # wrong (additive) residual model in the base: combined is proposed
  dt$H_EPS1 <- 0.15; dt[grep("^D_EPS", names(dt))] <- 0
  prep2 <- nmtools:::.nm_ob_prepare(dt, om, sg, NULL, NULL)
  b <- nmtools:::.nm_ob_fit(nmtools:::.nm_ob_spec(list(1L, 2L, 3L), "base"), prep2)
  cmb <- nmtools:::.nm_ob_fit(nmtools:::.nm_ob_spec(list(1L, 2L, 3L), "comb"), prep2, b$u)
  expect_lt(cmb$ofv - b$ofv, -100)
})

test_that("full run with both levels (long)", {
  skip_on_cran()
  ex <- ex_copy()
  rd <- file.path(ex, "run1.reg")
  .emul_deriv(rd, file.path(ex, "example_data.csv"), file.path(rd, "run1_deriv.dtab"))
  o <- nm_optimize_base(rd, covariates = c("WT", "GENO"), objective = "fast", verbose = FALSE)
  expect_false(is.null(o$level2))
  pr <- o$level2$proposal
  expect_match(pr$proposal[pr$part == "covariates"], "GENO on ETA1 \\(CL\\)")
  expect_match(pr$proposal[pr$part == "covariates"], "WT on ETA2 \\(V\\)")
  expect_equal(pr$proposal[pr$part == "residual error"], "keep the residual model")
  p <- plot(o)
  expect_s3_class(p, "ggplot")
  eff <- setNames(as.character(p$data$effect), sub("^.*\r", "", as.character(p$data$candidate)))
  expect_equal(unname(eff["remove ETA1 (CL)"]), "worse (p < 0.01)")
  expect_equal(unname(eff["GENO on ETA1 (CL)"]), "better (p < 0.01)")
})

test_that("residual starts follow the base residual model", {
  ex <- ex_copy()
  rd <- file.path(ex, "run1.reg")
  .emul_deriv(rd, file.path(ex, "example_data.csv"), file.path(rd, "run1_deriv.dtab"))
  fin <- nmtools:::.nm_vpc_read_ext_final(file.path(rd, "run1.ext"))
  om <- nmtools:::.nm_ob_ext_matrix(fin, "OMEGA"); sg <- nmtools:::.nm_ob_ext_matrix(fin, "SIGMA")
  dt <- nmtools:::.nm_read_table(file.path(rd, "run1_deriv.dtab"))
  prep <- nmtools:::.nm_ob_prepare(dt, om, sg, NULL, NULL)
  st <- nmtools:::.nm_ob_start(nmtools:::.nm_ob_spec(list(1L, 2L, 3L), "comb"), prep)
  # run1: W = sqrt((THETA(4)*IPRED)^2 + THETA(5)^2), THETA(4) = 0.168, THETA(5) = 0.0138
  expect_equal(unname(exp(st[7:8])), c(0.0138296, 0.168218), tolerance = 1e-3)
})
