test_that("nm_convert_ctl replaces Wings labels by THETA/ETA/ERR", {
  f <- write_ctl(
    "$THETA (0,5) ;POPCL", " (0,50) ;POPV", " (0,0.1) ;RUV_P",
    "$OMEGA 0.1 ; eta_cl", "$OMEGA 0.1 ; eta_v",
    "$SIGMA 1 FIX ; EPS1",
    "$PK",
    "  CL = POPCL*EXP(eta_cl)",
    "  V  = POPV*EXP(eta_v)",
    "$ERROR",
    "  Y = F + F*RUV_P*EPS1")
  out <- tempfile(fileext = ".ctl")
  res <- nm_convert_ctl(f, out_file = out, verbose = FALSE)
  expect_equal(normalizePath(res), normalizePath(out))
  txt <- readLines(out)
  expect_true(any(grepl("CL = THETA(1)*EXP(ETA(1))", txt, fixed = TRUE)))
  expect_true(any(grepl("V  = THETA(2)*EXP(ETA(2))", txt, fixed = TRUE)))
  expect_true(any(grepl("THETA(3)*ERR(1)", txt, fixed = TRUE)))
})

test_that("nm_convert_ctl does not overwrite without permission", {
  f <- write_ctl("$THETA (0,5) ;POPCL", "$PK", "  CL = POPCL")
  out <- tempfile(fileext = ".ctl")
  nm_convert_ctl(f, out_file = out, verbose = FALSE)
  expect_error(nm_convert_ctl(f, out_file = out, verbose = FALSE))
  expect_no_error(nm_convert_ctl(f, out_file = out, overwrite = TRUE, verbose = FALSE))
})

test_that("a multi-line $OMEGA starting on the record line is fully mapped", {
  f <- write_ctl(
    "$THETA (0,5) ;POPCL", " (0,50) ;POPV",
    "$OMEGA 0.1 ; eta_cl", " 0.1 ; eta_v",
    "$PK", "  CL = POPCL*EXP(eta_cl)", "  V  = POPV*EXP(eta_v)")
  out <- tempfile(fileext = ".ctl")
  nm_convert_ctl(f, out_file = out, verbose = FALSE)
  expect_true(any(grepl("V  = THETA(2)*EXP(ETA(2))", readLines(out), fixed = TRUE)))
})
