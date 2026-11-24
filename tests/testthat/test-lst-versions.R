# Output layouts met in real runs of NONMEM 7.1 to 7.6 (and V/VI) for the
# point A.3: labels ET10/EP10, simulation-only run of NONMEM 7.1, tables
# printed in the .lst, OFV of NONMEM V/VI, and the static control stream check.
# The test runs are built from the example run2 (simulated data shipped with
# the package) or written by hand.

ex_lst2 <- function() readLines(system.file("extdata", "example", "run2.reg", "run2.lst", package = "nmpilot"), warn = FALSE)

make_run2 <- function(root, name, lst) {
  d <- file.path(root, paste0(name, ".nm7")); dir.create(d, recursive = TRUE)
  writeLines(lst, file.path(d, paste0(name, ".lst")))
  d
}
row_of2 <- function(r, run) as.data.frame(r$runs)[r$runs$run == run, , drop = FALSE]

# OMEGA matrix printed as NONMEM does, with labels ET10, ET11... from the 10th ETA
om_block <- function(d) {
  n <- length(d)
  lab <- ifelse(seq_len(n) < 10, paste0("ETA", seq_len(n)), paste0("ET", seq_len(n)))
  out <- c(" OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********", "", "",
           paste0("         ", paste(sprintf("%-8s", lab), collapse = "  ")))
  for (i in seq_len(n)) {
    v <- c(rep(0, i - 1L), d[i])
    out <- c(out, " ", paste0(" ", lab[i]),
             paste0("+        ", paste(sprintf("%8.2E", v), collapse = "  ")))
  }
  c(out, " ", "")
}

test_that("ETA labels ET10, ET11 in the .lst are read (run without .ext)", {
  L <- ex_lst2()
  d <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2)
  i_fin <- grep("FINAL PARAMETER ESTIMATE", L)[1]
  i_om <- grep("OMEGA - COV MATRIX FOR RANDOM EFFECTS", L); i_om <- i_om[i_om > i_fin][1]
  i_si <- grep("SIGMA - COV MATRIX FOR RANDOM EFFECTS", L); i_si <- i_si[i_si > i_om][1]
  L2 <- c(L[seq_len(i_om - 1L)], om_block(d), L[i_si:length(L)])

  v <- .nm_lst_block_values(L2[i_fin:length(L2)])
  expect_equal(unname(v["OMEGA(10,10)"]), 1.1)
  expect_equal(unname(v["OMEGA(11,11)"]), 1.2)
  expect_false(any(grepl("^OMEGA\\(9,1[0-9]\\)$", names(v))))

  root <- tempfile(); make_run2(root, "eta11", L2)
  r <- suppressWarnings(nm_extract_runs(root, verbose = FALSE))
  x <- row_of2(r, "eta11")
  expect_equal(x$params_source, "lst")
  expect_equal(x$n_eta, 11L)
  expect_equal(x$ofv, -1224.035, tolerance = 1e-6)
  p <- as.data.frame(r$params)
  expect_equal(p$value[p$type == "ETA_SD" & p$index == 10], sqrt(1.1), tolerance = 1e-6)
})

test_that("EPS labels EP10 in the .lst are read", {
  B <- c(" SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****", "",
         " EPS9", "+        0.00E+00  4.00E-02",
         " EP10", "+        0.00E+00  0.00E+00  9.00E-02")
  v <- .nm_lst_block_values(B)
  expect_equal(unname(v["SIGMA(10,3)"]), 0.09)
  expect_equal(unname(v["SIGMA(9,2)"]), 0.04)
})

test_that("a correlation whose variance is missing gives NA instead of an error", {
  ext <- list(raw = tibble::tibble(
    raw_name = c("THETA1", "OMEGA(1,1)", "OMEGA(2,1)", "SIGMA(1,1)"),
    value = c(1, 0.1, 0.02, 0.05), fixed = FALSE, se = NA_real_, rse_pct = NA_real_),
    ofv = 10, last_iteration = NA_real_, final_row_available = TRUE,
    table = NULL, source = "lst")
  expect_no_error(p <- .nm_ext_params(ext))
  p <- as.data.frame(p)
  expect_true(all(is.na(p$value[p$type == "ETA_COR"])))
})

test_that("a simulation-only run of NONMEM 7.1 (empty #METH: line) is a simulation", {
  L <- c("$PROBLEM sim", "$SIMULATION (20102) (330333 UNIFORM) ONLYSIM SUBPROBLEMS=1", "",
         "1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.1.0",
         " PROBLEM NO.:         1", " sim", "1",
         " PROBLEM NO.:         1     SUBPROBLEM NO.:      1", " ",
         " SIMULATION STEP PERFORMED", "1", " ", " #METH: ", "Stop Time: ",
         "Thu 02/25/2010 ", "09:21 AM")
  root <- tempfile(); make_run2(root, "sim71", L)
  x <- row_of2(suppressWarnings(nm_extract_runs(root, verbose = FALSE)), "sim71")
  expect_equal(x$run_state, "simulation")
  expect_equal(x$min_status, "SIMULATION")
  expect_true(is.na(x$min_ok))
})

test_that("NONMEM V/VI: tables printed in the .lst are not read as estimates; OFV is read", {
  L <- ex_lst2()
  L <- L[!grepl("#OBJV|OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT|#METH:|#OBJT", L)]
  i_min <- grep("MINIMUM VALUE OF OBJECTIVE FUNCTION", L)[1]
  if (is.na(i_min)) {
    i_fin <- grep("FINAL PARAMETER ESTIMATE", L)[1]
    L <- c(L[seq_len(i_fin - 3L)],
           " ********************               MINIMUM VALUE OF OBJECTIVE FUNCTION            ********************",
           L[(i_fin - 2L):length(L)])
    i_min <- grep("MINIMUM VALUE OF OBJECTIVE FUNCTION", L)[1]
  }
  # NONMEM V/VI print the value about ten lines below the banner
  j <- seq(i_min + 1L, length.out = 9L)
  L[j] <- ""
  L[i_min + 9L] <- " **************************************************      19955.455     *************************"
  tab <- c("1TABLE NO.  1      DATA RECS.     1 THROUGH   3", "",
           "           ID        TIME        DV        PRED",
           "+    1.00E+00    0.00E+00    1.20E+01    1.10E+01",
           "+    1.00E+00    1.00E+00    8.00E+00    7.90E+00",
           "+    1.00E+00    2.00E+00    5.00E+00    5.10E+00")
  L2 <- c(L, tab)
  e <- .nm_lst_estimates(L2)
  n_ref <- nrow(.nm_lst_estimates(L)$raw)
  expect_equal(nrow(e$raw), n_ref)
  root <- tempfile(); make_run2(root, "nm6", L2)
  x <- row_of2(suppressWarnings(nm_extract_runs(root, verbose = FALSE)), "nm6")
  expect_equal(x$ofv, 19955.455)
  expect_equal(x$n_theta, sum(grepl("^THETA", e$raw$raw_name)))
})

test_that("nm_ctl_check reads only the code records", {
  ctl <- c("$PROBLEM ORAL, 1 CMT, MDV in data",
           "$INPUT ID TIME AMT DV EVID",
           "$DATA data.csv IGNORE=@",
           "$SUBROUTINE ADVAN2 TRANS2",
           "$PK", "CL = THETA(1)*EXP(ETA(1))", "V = THETA(2)", "KA = THETA(3)", "S2 = V",
           "IF (MDV.EQ.0) DUM = 1",
           "$ERROR", "Y = F*(1 + EPS(1))",
           "$THETA (18.3 FIXED) ; CL", "$THETA (125 FIXED) ; V", "$THETA (0, 1.2)",
           "$OMEGA 0.1", "$SIGMA 0.04",
           "$EST METHOD=1 INTER",
           "$TABLE ID TIME CMT NOPRINT FILE=sdtab1")
  z <- nm_ctl_check(ctl)
  expect_false(any(grepl("THETA\\((18|125)\\)", z$message)))
  expect_false(any(grepl("^CMT is used", z$message)))
  # with PREDPP, NM-TRAN supplies EVID and MDV when $INPUT has none
  expect_false(any(grepl("^MDV is used", z$message)))
  ev <- sub("\\$INPUT ID TIME AMT DV EVID", "$INPUT ID TIME AMT DV", ctl)
  ev <- sub("IF \\(MDV.EQ.0\\) DUM = 1", "IF (EVID.EQ.0) DUM = 1", ev)
  expect_false(any(grepl("^EVID is used", nm_ctl_check(ev)$message)))

  # RATE with KA: no warning when a zero-order input (here D2) is modelled
  rt <- sub("\\$INPUT ID TIME AMT DV EVID", "$INPUT ID TIME AMT RATE DV EVID", ctl)
  expect_true(any(grepl("RATE is active", nm_ctl_check(rt)$message)))
  rt2 <- sub("S2 = V", "S2 = V\nD2 = 2", rt)
  rt2 <- unlist(strsplit(rt2, "\n"))
  expect_false(any(grepl("RATE is active", nm_ctl_check(rt2)$message)))
  expect_false(any(z$level == "ERROR"))

  # $PRED: MDV is not supplied by NM-TRAN
  pred <- c("$PROBLEM p", "$INPUT ID TIME DV", "$DATA d.csv",
            "$PRED", "IF (MDV.EQ.0) DUM = 1", "Y = THETA(1) + ETA(1) + EPS(1)",
            "$THETA 1", "$OMEGA 0.1", "$SIGMA 0.1", "$EST METHOD=1")
  expect_true(any(grepl("^MDV is used", nm_ctl_check(pred)$message)))
  expect_true(any(grepl("^EVID is used", nm_ctl_check(sub("MDV.EQ.0", "EVID.EQ.0", pred))$message)))

  # a THETA used in the code but not declared is still found
  bad <- sub("KA = THETA\\(3\\)", "KA = THETA(4)", ctl)
  expect_true(any(grepl("THETA\\(4\\) is used but not declared", nm_ctl_check(bad)$message)))

  # M3 checks still see $DATA and $EST
  m3 <- c("$PROBLEM m3", "$INPUT ID TIME DV BQL", "$DATA d.csv IGNORE=(BQL.EQ.1)",
          "$PRED", "LLOQ = 0.1", "IPRED = THETA(1)*EXP(ETA(1))", "W = 0.1",
          "IF (BQL.EQ.0) THEN", "F_FLAG = 0", "Y = IPRED + W*EPS(1)", "ELSE",
          "F_FLAG = 1", "Y = PHI((LLOQ - IPRED)/W)", "ENDIF",
          "$THETA 1", "$OMEGA 0.1", "$SIGMA 1 FIX", "$EST METHOD=1")
  zm <- nm_ctl_check(m3)
  expect_true(any(grepl("IGNORE refers to BQL", zm$message)))
  expect_true(any(grepl("without LAPLACE", zm$message)))
})

test_that("n_obs is read in the .lst when there is no .smr", {
  L <- ex_lst2()
  n_lst <- as.integer(sub(".*RECS:\\s*", "", grep("TOT. NO. OF OBS RECS:", L, value = TRUE)[1]))
  expect_equal(.nm_lst_n_obs(L), n_lst)
  expect_true(is.na(.nm_lst_n_obs(c("no header here"))))

  root <- tempfile(); make_run2(root, "obs", L)
  x <- row_of2(suppressWarnings(nm_extract_runs(root, verbose = FALSE)), "obs")
  expect_equal(x$n_obs, n_lst)
  expect_equal(x$bic, x$ofv + x$n_par * log(n_lst), tolerance = 1e-5)
})

test_that("nm_extract_run_table uses the .lst when the .ext is missing or empty", {
  skip_if_not_installed("flextable")
  skip_if_not_installed("officer")
  ex_dir <- system.file("extdata", "example", "run2.reg", package = "nmpilot")
  ref_root <- tempfile(); dir.create(file.path(ref_root, "run2.reg"), recursive = TRUE)
  file.copy(list.files(ex_dir, full.names = TRUE), file.path(ref_root, "run2.reg"))
  ref <- suppressMessages(nm_extract_run_table(file.path(ref_root, "run2")))
  expect_identical(ref$ext$source %||% "ext", "ext")

  # no .ext: only the .lst
  root <- tempfile(); d <- make_run2(root, "run2", ex_lst2())
  r <- suppressMessages(nm_extract_run_table(file.path(root, "run2")))
  expect_true(file.exists(file.path(d, "nonmem_results.docx")))
  expect_identical(r$ext$source, "lst")
  expect_true(is.na(r$ext_file))
  expect_equal(r$ext$ofv, ref$ext$ofv, tolerance = 1e-6)
  m <- merge(ref$params[, c("Parameter", "Estimate")], r$params[, c("Parameter", "Estimate")], by = "Parameter")
  expect_equal(nrow(m), nrow(ref$params))
  # the .lst prints 3 significant digits
  expect_equal(m$Estimate.y, m$Estimate.x, tolerance = 5e-3)

  # empty .ext: the .lst is used as well
  file.create(file.path(d, "run2.ext"))
  r2 <- suppressMessages(nm_extract_run_table(file.path(root, "run2")))
  expect_identical(r2$ext$source, "lst")

  # neither .ext nor .lst
  expect_error(nm_extract_run_table(file.path(tempfile(), "none")), "No .ext or .lst")
})
