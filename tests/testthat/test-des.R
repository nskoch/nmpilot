# Models written with $DES (ADVAN6/13): the functions that read or rewrite a
# control stream must keep $MODEL and $DES intact. The $DES model below is the
# example run2 (1 compartment, first-order absorption) rewritten with ADVAN13.

des_ctl <- function() {
  s <- readLines(system.file("extdata", "example", "run2.ctl", package = "nmpilot"))
  s <- sub("^\\$SUBROUTINES ADVAN2 TRANS2", "$SUBROUTINES ADVAN13 TOL=9", s)
  i <- grep("^\\$SUBROUTINES", s)
  s <- append(s, "$MODEL COMP=(DEPOT,DEFDOSE) COMP=(CENTRAL,DEFOBS)", after = i)
  j <- grep("^\\s*S2 = V", s)
  append(s, c("  K  = CL/V", "", "$DES", "  DADT(1) = -KA*A(1)",
              "  DADT(2) =  KA*A(1) - K*A(2)"), after = j)
}
keeps_des <- function(lines) {
  any(grepl("^\\$DES", lines)) && any(grepl("^\\$MODEL", lines)) &&
    any(grepl("ADVAN13", lines)) && sum(grepl("DADT\\(", lines)) == 2L
}
no_check_error <- function(lines) !any(unclass(nm_ctl_check(lines))$level == "ERROR")

test_that("nm_identify_parameters finds the same roles with $DES as with ADVAN2", {
  f <- tempfile(fileext = ".ctl"); writeLines(des_ctl(), f)
  a <- nm_identify_parameters(f)
  b <- nm_identify_parameters(system.file("extdata", "example", "run2.ctl", package = "nmpilot"))
  expect_equal(a$convergence, b$convergence)
  expect_equal(a$residual_error$theta, b$residual_error$theta)
})

test_that("control streams written by nmpilot keep $MODEL and $DES", {
  base <- des_ctl()
  expect_true(keeps_des(base)); expect_true(no_check_error(base))
  ext <- system.file("extdata", "example", "run2.reg", "run2.ext", package = "nmpilot")
  e <- .nm_vpc_get_estimates(.nm_vpc_read_ext_final(ext))
  orig <- .nm_vpc_build_orig_ctl(base, e, "run3_des")
  sim  <- .nm_vpc_build_sim_ctl(base, e, "run3_des", 1, 100)
  for (x in list(orig, sim)) { expect_true(keeps_des(x)); expect_true(no_check_error(x)) }
  tg <- .nm_find_target_theta(.nm_parse_theta_tokens(base), "CL")
  pr <- .nm_write_profile_ctl(base, tg, 4, 1, "C:/data/example_data.csv", TRUE, TRUE)
  expect_true(keeps_des(pr)); expect_true(no_check_error(pr))
  mu <- tempfile(fileext = ".ctl")
  f <- tempfile(fileext = ".ctl"); writeLines(base, f)
  suppressMessages(nm_ctl_mu(f, file = mu, verbose = FALSE))
  expect_true(keeps_des(readLines(mu))); expect_true(no_check_error(readLines(mu)))
})

test_that("nm_derivatives_ctl keeps $DES and adds the derivative tables", {
  root <- tempfile("des"); reg <- file.path(root, "run3_des.reg"); dir.create(reg, recursive = TRUE)
  writeLines(sub("\\$DATA example_data.csv", "$DATA ../example_data.csv", des_ctl()),
             file.path(reg, "run3_des.ctl"))
  ex <- system.file("extdata", "example", "run2.reg", package = "nmpilot")
  file.copy(file.path(ex, "run2.ext"), file.path(reg, "run3_des.ext"))
  file.copy(file.path(ex, "run2.lst"), file.path(reg, "run3_des.lst"))
  file.copy(file.path(ex, "run2.phi"), file.path(reg, "run3_des.phi"))
  out <- file.path(root, "deriv"); dir.create(out)
  suppressMessages(nm_derivatives_ctl(reg, out_dir = out))
  d <- readLines(list.files(out, "\\.ctl$", full.names = TRUE)[1])
  expect_true(keeps_des(d))
  expect_true(any(grepl("G_ETA1 = G\\(1,1\\)", d)))
  expect_true(any(grepl("MAXEVAL=0", d)))
})
