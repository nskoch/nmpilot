# Action A.4: control streams and outputs written by others than nmpilot.
# Source: package bbr (Metrum Research Group, MIT licence), see
# fixtures/public/bbr/README.md. The expected values were read by hand in the
# .lst files (#OBJV, TOT. NO. OF OBS RECS, ETASHRINKSD, termination message).

bbr_dir <- function() test_path("fixtures", "public", "bbr")
bbr_ctl <- function(run) {
  f <- list.files(file.path(bbr_dir(), run), "\\.(ctl|mod)$", full.names = TRUE)
  f[1]
}
bbr_runs <- c("1", "1001", "200", "acop-fake-bayes", "acop-onlysim", "acop-iov", "iovmm")

test_that("the public control streams are all present", {
  for (r in bbr_runs) expect_true(file.exists(bbr_ctl(r)), info = r)
})

test_that("nm_ctl_check reports no ERROR on control streams that NM-TRAN accepted", {
  for (r in bbr_runs) {
    ck <- nm_ctl_check(bbr_ctl(r))
    expect_false(any(ck$level == "ERROR"), info = r)
  }
})

test_that("$TABLE: items below the options are accepted, not after an option of the same line", {
  # anomaly A.4-1 (run iovmm, accepted by NM-TRAN 7.4.4)
  base <- c("$PROBLEM t", "$INPUT ID TIME DV AMT", "$DATA d.csv IGNORE=@",
            "$SUBROUTINES ADVAN1 TRANS2", "$PK", "CL = THETA(1)*EXP(ETA(1))", "V = THETA(2)",
            "S1 = V", "$ERROR", "Y = F + EPS(1)", "$THETA (0, 5) (0, 50)", "$OMEGA 0.1",
            "$SIGMA 0.1", "$ESTIMATION METHOD=1")
  ok <- nm_ctl_check(write_ctl(base, "$TABLE NOPRINT ONEHEADER FILE=sdtab", "ID TIME CL V"))
  expect_false(any(grepl("\\$TABLE", ok$message[ok$level == "ERROR"])))
  ko <- nm_ctl_check(write_ctl(base, "$TABLE ID TIME NOPRINT FILE=sdtab CL", "V"))
  expect_true(any(grepl("same line of \\$TABLE.*: CL\\.", ko$message)))
})

test_that("$ABBR REPLACE ranges are read (nm_ctl_check and nm_identify_parameters)", {
  # anomaly A.4-2 (run acop-iov)
  ab <- .nm_abbr_replace(c("$ABBR REPLACE ETA(OCC_CL) = ETA(,3:6)",
                           "$ABBR REPLACE ETA(OCC_V)=ETA(7,9,11) ; list",
                           "$ABBR REPLACE THETA(CL)=THETA(2)",
                           "$ABBR REPLACE ETA(X)=ETA(12 TO 16 BY 2)"))
  expect_identical(ab$name, c("OCC_CL", "OCC_V", "CL", "X"))
  expect_identical(ab$index, list(3:6, c(7L, 9L, 11L), 2L, c(12L, 14L, 16L)))

  ck <- nm_ctl_check(bbr_ctl("acop-iov"))
  expect_false(any(grepl("ETA\\(", ck$message)))
  ip <- nm_identify_parameters(bbr_ctl("acop-iov"))
  expect_identical(nrow(ip$omega), 62L)                     # BLOCK(1) SAME(29), twice
  expect_identical(nrow(ip$unmapped_omega), 0L)
  expect_length(ip$missing_omega_declaration, 0)
  expect_identical(pt_row(ip, "CL")$ETA_IIV, paste0("ETA(", c(1, 3:32), ")", collapse = ", "))
  expect_identical(pt_row(ip, "V")$ETA_IIV, paste0("ETA(", c(2, 33:62), ")", collapse = ", "))
})

test_that("nm_identify_parameters: covariate helpers on the log scale, KOUT, MTIME and $MIX", {
  # anomaly A.4-5 (run 200)
  ip <- nm_identify_parameters(bbr_ctl("200"))
  expect_identical(pt_row(ip, "CL")$THETA_typical, "THETA(3)")
  expect_identical(pt_row(ip, "CL")$THETA_covariate, "THETA(6), THETA(7), THETA(8)")
  # a plain product is still a typical value
  ip <- nm_identify_parameters(write_ctl("$INPUT ID TIME DV WT", "$PK",
                                         "CL = THETA(1)*(WT/70)**0.75*EXP(ETA(1))",
                                         "$THETA 1", "$OMEGA 0.1"))
  expect_identical(pt_row(ip, "CL")$THETA_typical, "THETA(1)")
  expect_true(is.na(pt_row(ip, "CL")$THETA_covariate))

  # anomaly A.4-6 (runs 1001 and iovmm)
  ip <- nm_identify_parameters(bbr_ctl("1001"))
  expect_identical(pt_row(ip, "KOUT")$THETA_typical, "THETA(1)")
  expect_true(is.na(pt_row(ip, "KOUT")$THETA_covariate))
  expect_identical(ip$special_theta$role, "model_event_time")
  expect_identical(ip$special_theta$variable, "MTIME(1)")
  expect_identical(nrow(ip$unmapped_theta), 0L)
  ip <- nm_identify_parameters(bbr_ctl("iovmm"))
  expect_identical(unique(ip$special_theta$theta), 6L)
  expect_identical(unique(ip$special_theta$role), "mixture_probability")
  expect_identical(nrow(ip$unmapped_theta), 0L)
})

test_that("nm_identify_parameters reads the public control streams", {
  for (r in bbr_runs) {
    ip <- nm_identify_parameters(bbr_ctl(r))
    expect_true(nrow(ip$parameter_table) > 0, info = r)
    expect_length(ip$missing_theta_declaration, 0)
    expect_length(ip$missing_omega_declaration, 0)
  }
  ip <- nm_identify_parameters(bbr_ctl("1"))
  expect_setequal(ip$parameter_table$PK_parameter, c("CL", "KA", "V"))
  expect_identical(pt_row(ip, "CL")$THETA_typical, "THETA(2)")
  expect_identical(pt_row(ip, "CL")$ETA_IIV, "ETA(1)")
  expect_identical(pt_row(ip, "KA")$THETA_typical, "THETA(1)")
  expect_identical(pt_row(ip, "V")$ETA_IIV, "ETA(2)")
  expect_setequal(ip$residual_error$theta, c(4L, 5L))

  ip <- nm_identify_parameters(bbr_ctl("iovmm"))
  expect_identical(pt_row(ip, "CL")$ETA_IIV, "ETA(1), ETA(3), ETA(4), ETA(5)")
  expect_identical(pt_row(ip, "V")$ETA_IIV, "ETA(2)")
})

test_that("the .lst of the public runs are read (OFV, status, observations)", {
  ref <- list(
    "1"               = list(ofv = 2583.311, n_obs = 741L, state = "completed", est = "FOCEI", min_ok = TRUE),
    "1001"            = list(ofv = 3842.571, n_obs = 599L, state = "completed", est = "BAYES", min_ok = TRUE),
    "acop-fake-bayes" = list(ofv = 2675.293, n_obs = 760L, state = "completed", est = "FOCEI > BAYES > IMP", min_ok = NA),
    "acop-iov"        = list(ofv = 44158.939, n_obs = 6396L, state = "completed", est = "FOCEI", min_ok = TRUE),
    "acop-onlysim"    = list(ofv = NA_real_, n_obs = 760L, state = "simulation", est = NA_character_, min_ok = NA),
    "iovmm"           = list(ofv = 14722.149, n_obs = 12600L, state = "completed", est = "FOCE", min_ok = TRUE)
  )
  for (r in names(ref)) {
    p <- .nm_parse_lst(list.files(file.path(bbr_dir(), r), "\\.lst$", full.names = TRUE))
    e <- ref[[r]]
    expect_equal(p$ofv, e$ofv, tolerance = 1e-6, info = r)
    expect_identical(p$n_obs, e$n_obs, info = r)
    expect_identical(p$run_state, e$state, info = r)
    expect_identical(p$est_methods, e$est, info = r)
    expect_identical(p$min_ok, e$min_ok, info = r)
  }
  p <- .nm_parse_lst(file.path(bbr_dir(), "1", "1.lst"))
  expect_identical(p$min_status, "MINIMIZATION SUCCESSFUL")
  expect_identical(p$cov_status, "COMPLETED")
  expect_equal(p$eta_shrink[1:2], c(17.470, 2.0728), tolerance = 1e-4)
  # last estimation step (IMP) of the chain FOCEI > BAYES > IMP
  p <- .nm_parse_lst(file.path(bbr_dir(), "acop-fake-bayes", "acop-fake-bayes.lst"))
  expect_equal(p$eta_shrink[1:2], c(21.253, 4.4410), tolerance = 1e-4)
  # $MIX: smallest shrinkage of each ETA over the two subpopulations
  p <- .nm_parse_lst(file.path(bbr_dir(), "iovmm", "iovmm.lst"))
  expect_equal(p$eta_shrink[1:5], c(13.972, 3.5034, 18.301, 16.893, 16.012), tolerance = 1e-4)
})

test_that("more than 10 ETAs: shrinkage and OMEGA read over the continuation lines", {
  # anomalies A.4-3 and A.4-4 (run acop-iov: 62 ETAs, no .ext)
  p <- .nm_parse_lst(file.path(bbr_dir(), "acop-iov", "acop-iov.lst"))
  expect_identical(length(p$eta_shrink), 124L)              # ETASHRINKSD then ETASHRINKVR
  expect_equal(p$eta_shrink[c(1, 11, 58, 62)], c(4.4133, 1e-10, 32.700, 2.5124), tolerance = 1e-4)
  x <- nm_extract_runs(bbr_dir(), verbose = FALSE)
  r <- as.data.frame(x$runs)[x$runs$run == "acop-iov", ]
  expect_identical(r$n_eta_same, 58L)
  expect_identical(r$n_omega, 4L)
  expect_equal(r$max_eta_shrink, 32.7, tolerance = 1e-6)
  om <- x$params[x$params$run == "acop-iov" & x$params$type == "OMEGA", ]
  expect_identical(nrow(om), 62L)
  expect_equal(om$value[om$index %in% c(1, 3, 33, 62)], c(0.0629, 0.223, 0.309, 0.309), tolerance = 1e-6)
})

test_that("nm_extract_runs reads the folder of public runs", {
  x <- nm_extract_runs(bbr_dir(), verbose = FALSE)
  r <- as.data.frame(x$runs)
  expect_setequal(r$run, c("1", "1001", "acop-fake-bayes", "acop-iov", "acop-onlysim", "iovmm"))
  ofv <- stats::setNames(r$ofv, r$run)
  expect_equal(unname(ofv[c("1", "1001", "acop-fake-bayes", "acop-iov", "iovmm")]),
               c(2583.311, 3842.571, 2675.293, 44158.939, 14722.149), tolerance = 1e-6)
  expect_true(is.na(ofv[["acop-onlysim"]]))
  expect_identical(r$n_subpop[r$run == "iovmm"], 2L)
  expect_identical(r$n_theta[r$run == "1"], 5L)
})

test_that("nm_convert_ctl keeps the public control streams free of NM-TRAN errors", {
  td <- tempfile("pub_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  for (r in c("1", "1001", "200", "acop-iov")) {
    f <- file.path(td, basename(bbr_ctl(r)))
    file.copy(bbr_ctl(r), f)
    out <- nm_convert_ctl(f, verbose = FALSE)
    expect_true(file.exists(out), info = r)
    expect_false(any(nm_ctl_check(out)$level == "ERROR"), info = r)
  }
})
