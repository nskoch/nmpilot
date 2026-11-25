has_line <- function(x, pattern) any(grepl(pattern, x, fixed = TRUE))
check_levels <- function(x) unique(nm_ctl_check(x)$level)

test_that("nm_write_ctl builds a one-compartment oral model by default", {
  x <- nm_write_ctl()
  expect_s3_class(x, "nm_ctl")
  expect_true(has_line(x, "$SUBROUTINES ADVAN2 TRANS2"))
  expect_true(has_line(x, "KA = TVKA*EXP(ETA(3))"))
  expect_true(has_line(x, "S2 = V"))
  expect_equal(check_levels(x), "OK")
})

test_that("every model combination passes the static check", {
  g <- expand.grid(route = c("IV", "EV"),
                   absorption = c("first-order", "zero-order", "zero-order-d1", "transit"),
                   elimination = c("linear", "michaelis-menten"), cpts = 1:3,
                   dynamics = c("analytic", "ode"), error = c("additive", "proportional", "combined"),
                   omega_structure = c("diagonal", "block"), sigma_structure = c("auto", "block"),
                   stringsAsFactors = FALSE)
  g <- g[!(g$route == "IV" & g$absorption != "first-order"), ]
  lev <- vapply(seq_len(nrow(g)), function(i) {
    a <- g[i, ]
    x <- nm_write_ctl(route = a$route, absorption = a$absorption, elimination = a$elimination,
                      cpts = a$cpts, dynamics = a$dynamics, error = a$error,
                      omega_structure = a$omega_structure, sigma_structure = a$sigma_structure)
    paste(check_levels(x), collapse = ",")
  }, character(1))
  expect_true(all(lev == "OK"))
})

test_that("zero-order absorption uses a modelled rate R1 and keeps RATE", {
  x <- nm_write_ctl(absorption = "zero-order", cpts = 2,
                    input = c("ID", "TIME", "DV", "AMT", "RATE", "EVID", "MDV", "CMT"))
  expect_true(has_line(x, "$SUBROUTINES ADVAN3 TRANS4"))
  expect_true(has_line(x, "R1 = K0"))
  expect_true(has_line(x, "$INPUT ID TIME DV AMT RATE EVID MDV CMT"))
  y <- nm_write_ctl(absorption = "zero-order", elimination = "michaelis-menten")
  expect_false(any(grepl("DADT.*K0", y)))
  expect_true(has_line(y, "DADT(1) = - VMEL"))
})

test_that("RATE is dropped automatically with first-order absorption", {
  x <- nm_write_ctl(input = c("#ID", "TIME", "DV", "AMT", "RATE", "EVID", "MDV"))
  expect_true(has_line(x, "$INPUT ID TIME DV AMT RATE=DROP EVID MDV"))
  expect_true(has_line(x, "$DATA data.csv IGNORE=#"))
})

test_that("THETA, IIV and SIGMA options are written as asked", {
  x <- nm_write_ctl(route = "IV", cpts = 2, theta_init = c(CL = 3), theta_lower = c(V1 = NA),
                    theta_upper = c(Q = 50), theta_fix = "V2",
                    iiv = c(CL = 0.09, V1 = 0.04), iiv_scale = "omega",
                    error = "proportional", sigma_fix = "prop")
  expect_true(has_line(x, "(0, 3)  ; CL"))
  expect_true(has_line(x, "  50  ; V1"))
  expect_true(has_line(x, "(0, 5, 50)  ; Q"))
  expect_true(has_line(x, "100 FIX  ; V2"))
  expect_true(has_line(x, "$OMEGA 0.09 ; ETA(CL)"))
  expect_false(any(grepl("ETA(Q)", x, fixed = TRUE)))
  expect_true(has_line(x, "$SIGMA 0.09 FIX ; ERR(1) proportional"))
  expect_true(has_line(x, "$SIGMA 0 FIX ; ERR(2) additive"))
  expect_equal(check_levels(x), "OK")
})

test_that("nm_write_ctl stops on unknown parameters and existing files", {
  expect_error(nm_write_ctl(theta_init = c(KE = 1)), "unknown parameter")
  expect_error(nm_write_ctl(iiv = "Q"), "unknown parameter")
  expect_error(nm_write_ctl(theta_init = c(CL = -1)), "lower bound")
  f <- tempfile(fileext = ".ctl")
  nm_write_ctl(f, verbose = FALSE)
  expect_true(file.exists(f))
  expect_error(nm_write_ctl(f, verbose = FALSE), "already exists")
  expect_no_error(nm_write_ctl(f, overwrite = TRUE, verbose = FALSE))
})

test_that("nm_ctl_check counts declarations over all records", {
  ctl <- c("$INPUT ID TIME DV AMT", "$PK",
           "CL = THETA(1)*EXP(ETA(1)+ETA(4))", "V = THETA(2)*EXP(ETA(2)+ETA(5))",
           "KA = THETA(3)*EXP(ETA(3))",
           "$ERROR", "Y = F*(1+EPS(1)) + EPS(2)",
           "$THETA (0,5) (0,50) ; two values on the record line", "$THETA 1 FIX",
           "$OMEGA 0.1 0.2", "$OMEGA 0.3",
           "$OMEGA BLOCK(1) 0.02", "$OMEGA BLOCK SAME",
           "$SIGMA 0.04 0.01")
  expect_equal(check_levels(ctl), "OK")
  bad <- c(ctl[1:7], "$THETA (0,5) (0,50)", "$OMEGA 0.1 0.2 0.3", "$SIGMA 0.04")
  res <- nm_ctl_check(bad)
  expect_true(any(grepl("THETA(3) is used but not declared", res$message, fixed = TRUE)))
  expect_true(any(grepl("ETA(4) is used without", res$message, fixed = TRUE)))
  expect_true(any(grepl("EPS/ERR(2) is used without", res$message, fixed = TRUE)))
})

test_that("nm_ctl_check reports the differences from a reference", {
  ctl <- nm_write_ctl()
  res <- nm_ctl_check(nm_ctl_mu(ctl), reference = ctl)
  d <- attr(res, "diff")
  expect_equal(sum(startsWith(d, "- ")), 6L)
  expect_equal(sum(startsWith(d, "+ ")), 9L)
  expect_output(print(res), "Differences from reference")
})

test_that("nm_ctl_mu rewrites log-normal parameters and keeps blank lines", {
  ctl <- nm_write_ctl()
  mu <- nm_ctl_mu(ctl)
  expect_true(has_line(mu, "MU_1 = LOG(THETA(1))"))
  expect_true(has_line(mu, "CL = EXP(MU_1 + ETA(1))"))
  expect_equal(sum(mu == ""), sum(ctl == ""))
  expect_equal(check_levels(mu), "OK")
  lin <- nm_ctl_mu(ctl, parameters = "CL", mode = "linear")
  expect_true(has_line(lin, "MU_1 = THETA(1)"))
  expect_true(any(grepl("^\\s*1.609438 ; log\\(CL\\)", lin)))
  expect_error(nm_ctl_mu(ctl, parameters = "Q"), "not in a supported form")
  expect_error(nm_ctl_mu(mu), "no parameter written")
})

test_that("nm_ctl_iov adds ETAs with $ABBR REPLACE and BLOCK SAME", {
  ctl <- nm_write_ctl(input = c("ID", "TIME", "DV", "AMT", "EVID", "MDV", "CMT", "OCC"))
  x <- nm_ctl_iov(ctl, parameters = c("CL", "V"), n_occ = 3)
  expect_true(has_line(x, "$ABBR REPLACE ETA(OCC_BOVCL)=ETA(4,6,8)"))
  expect_true(has_line(x, "$ABBR REPLACE ETA(OCC_BOVV)=ETA(5,7,9)"))
  expect_true(has_line(x, "CL = TVCL*EXP(ETA(1)+ETA(OCC_BOVCL))"))
  expect_true(has_line(x, "$OMEGA BLOCK(2) SAME(2)"))
  expect_equal(check_levels(x), "OK")
  expect_error(nm_ctl_iov(nm_write_ctl(), "CL", n_occ = 2), "active item of \\$INPUT")
  expect_error(nm_ctl_iov(ctl, "CL", n_occ = 1), "at least 2")
  dat <- data.frame(ID = 1, OCC = c(1, 2, 4))
  expect_error(nm_ctl_iov(ctl, "CL", n_occ = 3, data = dat), "outside 1..3")
})

test_that("nm_ctl_m3 implements the censored likelihood", {
  ctl <- nm_write_ctl(input = c("ID", "TIME", "DV", "AMT", "EVID", "MDV", "CMT", "BQL"),
                      error = "combined", method = "SAEM")
  x <- suppressMessages(nm_ctl_m3(ctl, lloq = 0.1))
  expect_true(has_line(x, "M3LLOQ = 0.1"))
  expect_true(has_line(x, "M3SD = SQRT(IPRED**2*SIGMA(1,1)+SIGMA(2,2))"))
  expect_true(has_line(x, "IF (BQL.EQ.1) M3CENS = 1"))
  expect_true(has_line(x, "IF (M3CENS.EQ.1) THEN"))
  expect_true(has_line(x, "$ESTIMATION METHOD=COND INTERACTION LAPLACE"))
  expect_equal(check_levels(x), "OK")
  y <- nm_ctl_m3(ctl, method_col = "METH", method_lloq = c("1" = 0.1, "2" = 0.05), verbose = FALSE)
  expect_true(has_line(y, "IF (METH.EQ.2) M3LLOQ = 0.05"))
  expect_true(has_line(y, "ACTION REQUIRED: METH is not in $INPUT"))
  expect_error(nm_ctl_m3(ctl), "exactly one")
  expect_error(nm_ctl_m3(x, lloq = 0.1), "F_FLAG already exists")
})

test_that("transformations accept a file and write a file", {
  f <- tempfile(fileext = ".ctl"); out <- tempfile(fileext = ".ctl")
  nm_write_ctl(f, verbose = FALSE)
  nm_ctl_mu(f, output_file = out, verbose = FALSE)
  expect_true(any(grepl("MU_1", readLines(out), fixed = TRUE)))
  expect_error(nm_ctl_mu(tempfile(fileext = ".ctl")), "not found")
})
