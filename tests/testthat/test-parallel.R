# Parallel NONMEM runs: a fake nmfe script (sh, or .bat on Windows) records how
# many runs are going on at the same time, sleeps, then writes a .lst and a .ext.

fake_nmfe <- function(root) {
  shared <- file.path(root, "shared")
  dir.create(shared, showWarnings = FALSE)
  if (.Platform$OS.type == "windows") {
    sh <- normalizePath(shared, winslash = "\\")
    f <- file.path(root, "fake_nmfe.bat")
    writeLines(c("@echo off",
                 # %RANDOM% is seeded from the clock: two cmd started in the same second
                 # get the same value, so the working folder name is added to the marker
                 "for %%d in (\"%CD%\") do set D=%%~nxd",
                 "set R=%D%.%RANDOM%",
                 sprintf("type nul > \"%s\\%%1.%%R%%.running\"", sh),
                 "set /a C=0",
                 sprintf("for %%%%f in (\"%s\\*.running\") do set /a C+=1", sh),
                 sprintf(">> \"%s\\conc.txt\" echo %%C%%", sh),
                 "ping -n 3 127.0.0.1 > nul",
                 "> %2 echo 0MINIMIZATION SUCCESSFUL",
                 ">> %2 echo  #OBJV:****************     123.456       ****************",
                 "for %%a in (%1) do set NAME=%%~na",
                 "> %NAME%.ext echo TABLE NO.     1: First Order Conditional Estimation",
                 ">> %NAME%.ext echo  ITERATION    THETA1       SIGMA(1,1)   OBJ",
                 ">> %NAME%.ext echo  -1000000000  2.5          1.0          123.456",
                 sprintf("del \"%s\\%%1.%%R%%.running\"", sh)), f)
  } else {
    f <- file.path(root, "fake_nmfe.sh")
    writeLines(c("#!/bin/sh",
                 sprintf("S='%s'", shared),
                 "touch \"$S/$1.$$.running\"",
                 "ls \"$S\"/*.running | wc -l >> \"$S/conc.txt\"",
                 "sleep 2",
                 "printf '0MINIMIZATION SUCCESSFUL\\n #OBJV:****************     123.456       ****************\\n' > \"$2\"",
                 "N=$(basename \"$1\" .ctl)",
                 "printf 'TABLE NO.     1: First Order Conditional Estimation\\n ITERATION    THETA1       SIGMA(1,1)   OBJ\\n  -1000000000  2.5   1.0   123.456\\n' > \"$N.ext\"",
                 "rm \"$S/$1.$$.running\""), f)
    Sys.chmod(f, "755")
  }
  list(nmfe = f, conc = file.path(shared, "conc.txt"))
}
max_conc <- function(fk) {
  x <- readLines(fk$conc, warn = FALSE, encoding = "latin1")
  max(suppressWarnings(as.integer(gsub("[^0-9]", "", x))), na.rm = TRUE)
}

test_that("n_cores is resolved from a number, \"auto\" or \"ask\"", {
  expect_equal(nmtools:::.nm_resolve_cores(1), 1L)
  expect_equal(suppressWarnings(nmtools:::.nm_resolve_cores(3.0)), 3L)
  expect_gte(nmtools:::.nm_resolve_cores("auto"), 1L)
  expect_error(nmtools:::.nm_resolve_cores(0), "n_cores")
  expect_error(nmtools:::.nm_resolve_cores("many"))
  old <- options(nmtools.n_cores = NULL)
  on.exit(options(old), add = TRUE)
  if (!interactive()) expect_error(nmtools:::.nm_resolve_cores("ask"), "interactive")
  # an answer given earlier in the session is reused without asking
  options(nmtools.n_cores = 3L)
  expect_message(k <- nmtools:::.nm_resolve_cores("ask"), "kept for this session")
  expect_equal(k, 3L)
  options(nmtools.n_cores = NULL)
  big <- nm_cores(verbose = FALSE)$logical + 5L
  expect_warning(nmtools:::.nm_resolve_cores(big), "logical cores")
  expect_named(nm_cores(verbose = FALSE), c("logical", "physical", "auto"))
})

test_that("nm_run runs several control streams in parallel, at most n_cores at a time", {
  skip_on_cran()
  root <- tempfile("par"); dir.create(root)
  fk <- fake_nmfe(root)
  writeLines(c("ID,TIME,DV", "1,1,5"), file.path(root, "data.csv"))
  ctls <- sprintf("r%d.ctl", 1:6)
  for (f in ctls) writeLines(c("$PROB x", "$INPUT ID TIME DV", "$DATA data.csv",
                               "$PRED", " Y = THETA(1) + EPS(1)", "$THETA 1", "$SIGMA 1",
                               "$EST METHOD=1"), file.path(root, f))
  t0 <- Sys.time()
  log <- nm_run(fk$nmfe, ctls, root, n_cores = 2, verbose = FALSE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_equal(nrow(log), 6L)
  expect_true(all(log$action == "launched"))
  expect_true(all(log$exit_status == 0L))
  expect_equal(log$ofv, rep(123.456, 6))
  expect_true(all(log$min_ok))
  expect_equal(max_conc(fk), 2L)
  expect_lt(el, 6 * 2)                          # faster than one after the other
  expect_true(file.exists(file.path(root, "r1.reg", "r1_nmtools.log")))
  # finished runs are skipped on a second call
  log2 <- nm_run(fk$nmfe, ctls, root, n_cores = 2, verbose = FALSE)
  expect_true(all(log2$action == "skipped_done"))
  # a flat layout falls back to one run after the other
  expect_warning(nm_run(fk$nmfe, ctls[1:2], root, run_subdir = "", n_cores = 2,
                        skip_terminated = FALSE, verbose = FALSE), "one after the other")
})

test_that("nm_profile, nm_jackknife and nm_bootstrap run in parallel", {
  skip_on_cran()
  root <- tempfile("par2"); dir.create(root)
  fk <- fake_nmfe(root)
  # profile
  writeLines(c("ID,TIME,DV", "1,1,5"), file.path(root, "data.csv"))
  writeLines(c("$PROB x", "$INPUT ID TIME DV", "$DATA data.csv IGNORE=@", "$PRED",
               "  CL = THETA(1)", "  Y = CL + EPS(1)", "$THETA (0, 2) ; CL", "$SIGMA 1",
               "$EST METHOD=1"), file.path(root, "base.ctl"))
  pr <- nm_profile(fk$nmfe, file.path(root, "base.ctl"), param = "CL",
                   from = 1, to = 3, by = 0.5, quiet = TRUE, n_cores = 2)
  expect_equal(pr$results$ofv, rep(123.456, 5))
  expect_equal(max_conc(fk), 2L)
  # jackknife: one subdirectory per run
  unlink(fk$conc)
  writeLines(c("#ID,TIME,AMT,DV", "1,0,100,.", "1,1,.,5.10", "2,0,100,.", "2,1,.,4.00",
               "3,0,100,.", "3,1,.,6.25", "4,0,100,.", "4,1,.,5.50"), file.path(root, "jdata.csv"))
  writeLines(c("$PROB x", "$INPUT ID TIME AMT DV", "$DATA jdata.csv IGNORE=#",
               "$PRED", "  Y = THETA(1) + EPS(1)", "$THETA 1", "$SIGMA 1", "$EST METHOD=1"),
             file.path(root, "run1.ctl"))
  writeLines(c("0MINIMIZATION SUCCESSFUL", " #OBJV:************     120.000   ********"),
             file.path(root, "run1.lst"))
  jk <- suppressWarnings(nm_jackknife(file.path(root, "run1.ctl"), nmfe_path = fk$nmfe,
                                      save_figures = FALSE, verbose = FALSE, n_cores = 2))
  expect_true(all(file.exists(file.path(jk$output_dir, "jk_run", paste0("jk_", 1:4),
                                        paste0("jk_", 1:4, ".lst")))))
  expect_equal(nrow(jk$summary), 4L)
  expect_equal(max_conc(fk), 2L)
  # bootstrap: working subdirectories slot1, slot2
  unlink(fk$conc)
  writeLines(c("ID,TIME,AMT,DV", paste(rep(1:4, each = 2), rep(c(0, 1), 4),
                                        rep(c(100, "."), 4), rep(c(".", 5), 4), sep = ",")),
             file.path(root, "bdata.csv"))
  writeLines(c("$PROB x", "$INPUT ID TIME AMT DV", "$DATA bdata.csv IGNORE=@", "$PRED",
               "  Y = THETA(1) + EPS(1)", "$THETA 1", "$SIGMA 1", "$EST METHOD=1"),
             file.path(root, "b1.ctl"))
  b <- nm_bootstrap(file.path(root, "b1.ctl"), nmfe_path = fk$nmfe, n = 5, seed = 1,
                    verbose = FALSE, n_cores = 2)
  expect_equal(b$n_run, 5L)
  expect_equal(b$estimates$replicate, 1:5)
  expect_true(all(b$estimates$status == "successful"))
  expect_equal(b$estimates$THETA1, rep(2.5, 5))
  expect_equal(max_conc(fk), 2L)
  expect_true(dir.exists(file.path(b$run_dir, "slot2")))
  expect_false(file.exists(file.path(b$run_dir, "slot1", "b1_boot.csv")))
})
