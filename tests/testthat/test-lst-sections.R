# Outputs with several $EST, several $PROB, empty or inconsistent .ext,
# failed or interrupted runs, IOV (SAME) and $MIX. The test runs are built
# from the example run2 (simulated data shipped with the package).

ex_lst <- function() readLines(system.file("extdata", "example", "run2.reg", "run2.lst", package = "nmtools"), warn = FALSE)
ex_ext <- function() readLines(system.file("extdata", "example", "run2.reg", "run2.ext", package = "nmtools"), warn = FALSE)

make_run <- function(root, name, lst = NULL, ext = NULL, ext_name = NULL, ctl = NULL) {
  d <- file.path(root, paste0(name, ".nm7")); dir.create(d, recursive = TRUE)
  if (!is.null(lst)) writeLines(lst, file.path(d, paste0(name, ".lst")))
  if (!is.null(ext)) writeLines(ext, file.path(d, ext_name %||% paste0(name, ".ext")))
  if (!is.null(ctl)) writeLines(ctl, file.path(d, paste0(name, ".ctl")))
  d
}
`%||%` <- function(a, b) if (is.null(a)) b else a
row_of <- function(r, run) as.data.frame(r$runs)[r$runs$run == run, , drop = FALSE]

test_that("several $EST: OFV and status come from the last step", {
  L <- ex_lst(); i <- grep("#METH:", L)[1]
  its <- c(" #METH: Iterative Two Stage", "", " OPTIMIZATION WAS NOT COMPLETED", "",
           " #OBJV:********************************************     -1100.000       **********")
  root <- tempfile(); make_run(root, "multi", lst = c(L[seq_len(i - 1)], its, L[i:length(L)]))
  r <- suppressWarnings(nm_extract_runs(root, verbose = FALSE))
  x <- row_of(r, "multi")
  expect_equal(x$lst_ofv, -1224.035, tolerance = 1e-6)
  expect_match(x$min_status, "MINIMIZATION SUCCESSFUL")
  expect_true(x$min_ok)
  expect_equal(x$n_est, 2L)
  expect_equal(x$est_methods, "ITS > FOCEI")
  expect_equal(x$run_state, "completed")
  # no .ext: the estimates are read in the .lst
  expect_equal(x$params_source, "lst")
  expect_gt(x$n_theta, 0)
})

test_that("an empty .ext or an .ext from another run is replaced by the .lst", {
  root <- tempfile()
  make_run(root, "ok", lst = ex_lst(), ext = ex_ext())
  make_run(root, "empty", lst = ex_lst(), ext = character(0))
  E <- ex_ext(); k <- grep("^\\s*-1000000000", E)
  E[k] <- sub("[-0-9.E+]+\\s*$", "99999.123", E[k])
  make_run(root, "stale", lst = ex_lst(), ext = E)
  expect_warning(r <- nm_extract_runs(root, verbose = FALSE), "differs from the .lst")
  ok <- row_of(r, "ok"); em <- row_of(r, "empty"); st <- row_of(r, "stale")
  expect_equal(ok$params_source, "ext")
  expect_equal(em$params_source, "lst")
  expect_equal(st$params_source, "lst")
  expect_equal(st$ofv, -1224.035, tolerance = 1e-6)
  expect_match(st$warnings, "EXT/LST MISMATCH")
  p <- as.data.frame(r$params)
  th <- function(run) p$value[p$run == run & p$type == "THETA"]
  expect_equal(th("empty"), th("ok"), tolerance = 0.01)
  expect_equal(em$n_omega, ok$n_omega)
})

test_that("several $PROB: the last problem with an estimation is read", {
  L <- c(ex_lst(), "1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.6.0",
         " PROBLEM NO.:         2", " SIMULATION", "0SIMULATION STEP PERFORMED", "Stop Time:")
  root <- tempfile(); make_run(root, "twoprob", lst = L)
  x <- row_of(suppressWarnings(nm_extract_runs(root, verbose = FALSE)), "twoprob")
  expect_equal(x$n_problems, 2L)
  expect_equal(x$ofv, -1224.035, tolerance = 1e-6)
  expect_equal(x$run_state, "completed")
})

test_that("interrupted, terminated and NM-TRAN failed runs are recognised", {
  L <- ex_lst(); i <- grep("#METH:", L)[1]
  root <- tempfile()
  make_run(root, "cut", lst = L[seq_len(i + 20)])
  make_run(root, "term", lst = c(L[seq_len(i + 20)], " PROGRAM TERMINATED BY OBJ"))
  make_run(root, "nmtran", lst = c(" AN ERROR WAS FOUND IN THE CONTROL STATEMENTS.",
                                   " AT LINE 12 THE CHARACTERS IN ERROR ARE: XX"))
  r <- suppressWarnings(nm_extract_runs(root, verbose = FALSE))
  expect_equal(row_of(r, "cut")$run_state, "interrupted")
  expect_false(row_of(r, "cut")$min_ok)
  expect_true(is.na(row_of(r, "cut")$ofv))
  expect_equal(row_of(r, "term")$run_state, "terminated")
  expect_equal(row_of(r, "nmtran")$run_state, "nmtran_error")
  expect_equal(row_of(r, "nmtran")$min_status, "NMTRAN_ERROR")
})

test_that("the run name ignores the case of RUN.EXT (NONMEM 7.1)", {
  root <- tempfile(); make_run(root, "caseA", lst = ex_lst(), ext = ex_ext(), ext_name = "CASEA.EXT")
  r <- nm_extract_runs(root, verbose = FALSE)
  expect_equal(r$runs$run, "caseA")
  expect_equal(r$runs$params_source, "ext")
})

test_that("IOV: ETAs declared SAME are found (BLOCK, BLOC, SAME(k))", {
  ctl <- c("$OMEGA 0.1 ; CL", "$OMEGA 0 FIX", "$OMEGA BLOCK(1) 0.1 ; IOV1",
           "$OMEGA BLOCK(1) SAME", "$OMEGA BLOC(2) 0.1 0.01 0.1", "$OMEGA BLOCK(2) SAME(2)",
           ";$OMEGA BLOCK(1) SAME", "$SIGMA 0.1")
  s <- nmtools:::.nm_ctl_same_etas(ctl)
  expect_equal(s$etas, c(4L, 7L, 8L, 9L, 10L))
  expect_equal(s$n_records, 2L)
})

test_that("$MIX: subpopulations, proportion THETAs and .phm assignment", {
  ctl <- c("$THETA 1 2 (0,0.5,1)", "$MIX", "  P(1)=THETA(3)", "  P(2)=1.-THETA(3)", "  NSPOP=2", "$EST")
  m <- nmtools:::.nm_ctl_mixture(ctl)
  expect_true(m$mix); expect_equal(m$n_subpop, 2L); expect_equal(m$p_theta, 3L)
  expect_false(nmtools:::.nm_ctl_mixture("$THETA 1")$mix)

  phm <- tempfile(fileext = ".phm")
  writeLines(c("TABLE NO.     1: First Order Conditional Estimation with Interaction: Problem=1",
               " SUBJECT_NO   ID           SUBPOP       PMIX         OBJ",
               "            1            1            1  8.00000E-01  10.0",
               "            1            1            2  2.00000E-01  12.0",
               "            2            2            1  1.00000E-01  9.0",
               "            2            2            2  9.00000E-01  8.0"), phm)
  d <- nmtools:::.nm_parse_phm(phm)
  expect_equal(d$SUBPOP[d$assigned], c(1L, 2L))
})

test_that("mixture shrinkage: the smallest value over subpopulations is kept", {
  L <- c(" ETAshrink(%):   1.7689E+00  5.0761E+00  1.0000E+02  1.0000E+02",
         " EBVshrink(%):   1.9973E+00  5.0717E+00  1.0000E+02  1.0000E+02",
         " ETAshrink(%):   1.0000E+02  1.0000E+02  1.4253E+00  2.0000E+01")
  expect_equal(nmtools:::.nm_lst_shrink(L, "ETA.?SHRINK"), c(1.7689, 5.0761, 1.4253, 20))
  # without $MIX, SD then VR, as before
  L2 <- c(" ETASHRINKSD(%)  4.2618E+00  1.0118E+01", " ETASHRINKVR(%)  8.3420E+00  1.9213E+01")
  expect_equal(nmtools:::.nm_lst_shrink(L2, "ETA.?SHRINK"), c(4.2618, 10.118, 8.342, 19.213))
})

test_that("nm_run, nm_profile, nm_bootstrap, nm_jackknife and nm_compare read the last step", {
  L <- ex_lst(); i <- grep("#METH:", L)[1]
  its <- c(" #METH: Iterative Two Stage", "", " OPTIMIZATION WAS NOT COMPLETED", "",
           " THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********", "",
           "         TH 1      TH 2", "         1.00E+00  2.00E+00", "",
           " OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********",
           " #OBJV:********************************************     -1100.000       **********")
  f1 <- tempfile(fileext = ".lst"); writeLines(L, f1)
  f2 <- tempfile(fileext = ".lst"); writeLines(c(L[seq_len(i - 1)], its, L[i:length(L)]), f2)
  # nm_run
  expect_equal(nmtools:::.nm_lst_status(f2)$minimization_status, "successful")
  expect_equal(nmtools:::.nm_lst_ofv(f2), nmtools:::.nm_lst_ofv(f1))
  # nm_profile / nm_jackknife
  expect_true(nmtools:::.nm_parse_status(f2)$min_ok)
  expect_equal(nmtools:::.nm_extract_ofv(f2), nmtools:::.nm_extract_ofv(f1))
  expect_equal(nmtools:::.nm_jk_theta_lst(readLines(f2)), nmtools:::.nm_jk_theta_lst(readLines(f1)))
  # nm_bootstrap
  expect_equal(nmtools:::.nm_min_status(f2), "successful")
  # nm_compare
  expect_equal(nmtools:::.nm_cmp_extract_ofv(f2), nmtools:::.nm_cmp_extract_ofv(f1))

  # last step = SAEM: its "COMPLETED" status counts as success
  saem <- c(L[seq_len(i - 1)], " #METH: Stochastic Approximation Expectation-Maximization", "",
            " STOCHASTIC PORTION WAS COMPLETED", " REDUCED STOCHASTIC PORTION WAS COMPLETED",
            " #OBJV:********************************************     -1200.000       **********")
  f3 <- tempfile(fileext = ".lst"); writeLines(saem, f3)
  expect_equal(nmtools:::.nm_lst_step_status(saem)$code, "successful")
  expect_equal(nmtools:::.nm_min_status(f3), "successful")
  # earlier step successful but last one terminated: terminated
  term <- c(L, " #METH: First Order Conditional Estimation with Interaction", "",
            "0MINIMIZATION TERMINATED", " DUE TO ROUNDING ERRORS (ERROR=134)")
  expect_equal(nmtools:::.nm_lst_step_status(term)$code, "rounding_errors")
})

test_that("nm_run: the covariance status ignores 'AFTER CHECKING THAT THE COVARIANCE STEP...'", {
  L <- ex_lst(); k <- grep("MINIMIZATION SUCCESSFUL", L)[1]
  L2 <- append(L, c(" HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.",
                    " REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY",
                    " AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT."), after = k)
  f <- tempfile(fileext = ".lst"); writeLines(L2, f)
  s <- nmtools:::.nm_lst_status(f)
  expect_equal(s$minimization_status, "successful")
  expect_equal(s$cov_status, "successful")
  expect_true(s$run_ok)
})
