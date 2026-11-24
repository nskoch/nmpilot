# NM-TRAN errors found by running the control streams written by nm_write_ctl,
# nm_ctl_mu and nm_ctl_m3 through NONMEM 7.6 (action 26), and the static
# checks added to nm_ctl_check() for them.

no_error <- function(x) !any(nm_ctl_check(x)$level == "ERROR")

test_that("ODE models are written with TOL (NM-TRAN error 145)", {
  for (x in list(nm_write_ctl(dynamics = "ode", verbose = FALSE),
                 nm_write_ctl(dynamics = "ode", ode_subroutine = "ADVAN6", verbose = FALSE),
                 nm_write_ctl(elimination = "michaelis-menten", verbose = FALSE),
                 nm_write_ctl(absorption = "transit", cpts = 2, verbose = FALSE))) {
    expect_match(grep("^\\$SUBROUTINES", x, value = TRUE), "ADVAN(13|6) TOL=9")
    expect_true(no_error(x))
  }
})

test_that("FO is written METHOD=ZERO (NM-TRAN error 287)", {
  x <- nm_write_ctl(method = "FO", verbose = FALSE)
  e <- grep("^\\$ESTIMATION", x, value = TRUE)
  expect_match(e, "METHOD=ZERO")
  expect_false(grepl("METHOD=FO\\b", e))
  expect_true(no_error(x))
})

test_that("the two peripheral compartments of ADVAN12 have different rate constants", {
  x <- nm_write_ctl(cpts = 3, verbose = FALSE)
  th <- function(p) as.numeric(sub(".*\\(0, ([0-9.]+)\\).*", "\\1", grep(paste0("; ", p, "$"), x, value = TRUE)))
  expect_false(isTRUE(all.equal(th("Q3") / th("V3"), th("Q4") / th("V4"))))
})

test_that("MU linear: log(1) = 0 is not written as an initial estimate (NM-TRAN error 260)", {
  x <- nm_ctl_mu(nm_write_ctl(verbose = FALSE), mode = "linear", verbose = FALSE)
  ka <- grep("; log\\(KA\\)", x, value = TRUE)
  expect_match(ka, "^\\s*0\\.0001 ")
  expect_true(no_error(x))
  # a fixed THETA keeps its exact value
  y <- nm_write_ctl(theta_fix = "KA", verbose = FALSE)
  y <- nm_ctl_mu(y, mode = "linear", verbose = FALSE)
  expect_match(grep("; log\\(KA\\)", y, value = TRUE), "^\\s*0 FIX")
})

test_that("M3 items are added before the options of $TABLE (NM-TRAN error 20)", {
  base <- nm_write_ctl(input = c("ID", "TIME", "AMT", "DV", "EVID", "MDV", "CMT", "BQL"),
                       error = "combined", verbose = FALSE)
  x <- nm_ctl_m3(base, bql = "BQL", lloq = 0.1, verbose = FALSE)
  t1 <- grep("^\\$TABLE", x, value = TRUE)[1]
  tok <- strsplit(t1, "\\s+")[[1]]
  expect_true(all(c("M3LLOQ", "M3SD", "M3PBQL") %in% tok))
  expect_lt(max(match(c("M3LLOQ", "M3SD", "M3PBQL"), tok)), min(grep("=", tok)))
  expect_true(no_error(x))
  # a comment at the end of the line is kept
  expect_identical(.nm_ctl_table_add_items("$TABLE ID DV NOPRINT FILE=sd ; tab", "X"),
                   "$TABLE ID DV X NOPRINT FILE=sd ; tab")
  expect_identical(.nm_ctl_table_add_items("$TABLE ID DV", "X"), "$TABLE ID DV X")
})

test_that("nm_ctl_check reports the NM-TRAN errors met in action 26", {
  ctl <- c("$PROBLEM p", "$INPUT ID TIME AMT DV CMT", "$DATA d.csv IGNORE=@",
           "$SUBROUTINES ADVAN13", "$MODEL COMP=(DEPOT) COMP=(CENTRAL)",
           "$PK", "CL = THETA(1)*EXP(ETA(1))", "V = THETA(2)", "KA = THETA(3)", "S2 = V",
           "$DES", "DADT(1) = -KA*A(1)", "DADT(2) = KA*A(1) - CL/V*A(2)",
           "$ERROR", "Y = F + EPS(1)",
           "$THETA (0, 5) 50", "$THETA 0", "$OMEGA 0.1", "$SIGMA 0.1",
           "$ESTIMATION METHOD=FO MAXEVAL=0",
           "$TABLE ID TIME NOPRINT FILE=sd CL")
  m <- nm_ctl_check(ctl)$message
  expect_true(any(grepl("ADVAN13 requires TOL", m)))
  expect_true(any(grepl("METHOD=FO is not a NONMEM method", m)))
  expect_true(any(grepl("THETA\\(3\\): the initial estimate is 0", m)))
  expect_true(any(grepl("after an option on the same line of \\$TABLE.*CL", m)))
  expect_false(any(grepl("THETA\\((1|2)\\)", m)))      # (0, 5): lower bound 0, estimate 5

  ok <- ctl
  ok[4] <- "$SUBROUTINES ADVAN13 TOL=9"
  ok[17] <- "$THETA 0 FIX"
  ok[20] <- "$ESTIMATION METHOD=0 MAXEVAL=0"
  ok[21] <- "$TABLE ID TIME CL NOPRINT FILE=sd"
  expect_true(no_error(ok))
  # a $TOL record, METHOD=COND / CONDITIONAL / SAEM, BY in $TABLE
  ok2 <- c(ok[-4], "$SUBROUTINES ADVAN13", "$TOL 9")
  expect_false(any(grepl("requires TOL", nm_ctl_check(ok2)$message)))
  for (meth in c("COND", "CONDITIONAL", "SAEM", "IMP", "1", "ZERO"))
    expect_false(any(grepl("not a NONMEM method", nm_ctl_check(sub("METHOD=0", paste0("METHOD=", meth), ok))$message)))
  expect_length(.nm_ctl_nmtran_checks(c("$TABLE ID TIME NOPRINT BY ID FILE=x")), 0)
})

test_that("M3: both Y definitions are in one IF/ELSE block (NM-TRAN error 122)", {
  base <- nm_write_ctl(input = c("ID", "TIME", "AMT", "DV", "EVID", "MDV", "CMT", "BQL"),
                       error = "combined", verbose = FALSE)
  for (npde in c(FALSE, TRUE)) {
    x <- trimws(nm_ctl_m3(base, bql = "BQL", lloq = 0.1, npde_mode = npde, verbose = FALSE))
    i_if <- which(x == "IF (M3CENS.EQ.1) THEN"); i_else <- which(x == "ELSE"); i_end <- which(x == "ENDIF")
    expect_length(i_if, 1L)
    y <- grep("^Y\\s*=", x)
    expect_length(y, 2L)
    expect_true(y[1] > i_if && y[1] < i_else && y[2] > i_else && y[2] < min(i_end[i_end > i_else]))
    expect_match(x[y[2]], "ERR\\(1\\)")
    # M3 quantities defined before the block; no test of ICALL with a data item
    expect_true(max(grep("^M3(SD|Z|PBQL|CENS) =", x)) < i_if)
    expect_false(any(grepl("ICALL\\.NE\\.", x)))
    expect_true(no_error(x))
  }
})

test_that("nm_ctl_check reports Y with EPS defined only under ICALL.EQ.4 (NM-TRAN error 122)", {
  head <- c("$PROBLEM p", "$INPUT ID TIME DV BQL", "$DATA d.csv IGNORE=@", "$PRED", "IPRED = THETA(1)*EXP(ETA(1))")
  tail <- c("$THETA 1", "$OMEGA 0.1", "$SIGMA 0.1", "$ESTIMATION METHOD=1 LAPLACE")
  has122 <- function(code) any(grepl("error 122", nm_ctl_check(c(head, code, tail))$message))
  # former M3 code of nm_ctl_m3: refused by NM-TRAN 7.6
  expect_true(has122(c("IF (BQL.EQ.0.OR.ICALL.EQ.4) THEN", "F_FLAG = 0", "Y = IPRED + EPS(1)", "ENDIF",
                       "IF (BQL.EQ.1.AND.ICALL.NE.4) THEN", "F_FLAG = 1", "Y = PHI(-IPRED)", "ENDIF")))
  expect_true(has122(c("IF (ICALL.EQ.4) Y = IPRED + EPS(1)")))
  # two IF blocks without ICALL (accepted by NONMEM 7.3 on a real run)
  expect_false(has122(c("IF (BQL.EQ.0) THEN", "F_FLAG = 0", "Y = IPRED + EPS(1)", "ENDIF",
                        "IF (BQL.EQ.1) THEN", "F_FLAG = 1", "Y = PHI(-IPRED)", "ENDIF")))
  # simulation block in addition to the estimation code
  expect_false(has122(c("Y = IPRED + EPS(1)", "IF (ICALL.EQ.4) THEN", "Y = IPRED + EPS(1)", "ENDIF")))
  # EPS in the ELSE part of an ICALL block
  expect_false(has122(c("IF (ICALL.EQ.4) THEN", "Y = IPRED", "ELSE", "Y = IPRED + EPS(1)", "ENDIF")))
})
