run_setup <- function() {
  root <- tempfile("an")
  dir.create(root)
  writeLines(c("ID,TIME,AMT,DV", "1,0,100,.", "1,1,.,5"), file.path(root, "data.csv"))
  # $DATA is relative to the .reg run directory (Wings convention).
  writeLines(c("$PROB x", "$INPUT ID TIME AMT DV", "$DATA ../data.csv IGNORE=@",
               "$PRED", "  Y = THETA(1) + EPS(1)", "$THETA 1", "$SIGMA 1"),
             file.path(root, "r1.ctl"))
  root
}

# Fake nmfe launcher: writes a .lst reporting a successful minimisation.
fake_nmfe <- function(dir) {
  if (.Platform$OS.type == "windows") {
    f <- file.path(dir, "fake_nmfe.bat")
    writeLines(c("@echo off", "echo MINIMIZATION SUCCESSFUL> %2"), f)
  } else {
    f <- file.path(dir, "fake_nmfe.sh")
    writeLines(c("#!/bin/sh", "echo \"MINIMIZATION SUCCESSFUL\" > \"$2\""), f)
    Sys.chmod(f, "755")
  }
  f
}

test_that("nm_run checks its inputs", {
  root <- run_setup()
  nmfe <- fake_nmfe(root)
  expect_error(nm_run(file.path(root, "none.bat"), "r1.ctl", root), "nmfe executable not found")
  expect_error(nm_run(nmfe, "r1.ctl", file.path(root, "nodir")), "does not exist")
  expect_error(nm_run(nmfe, character(0), root), "No .ctl")
  log <- nm_run(nmfe, "absent.ctl", root, verbose = FALSE)
  expect_equal(log$action, "skipped_no_ctl")
})

test_that("nm_run copies the control stream, fixes $DATA and launches nmfe", {
  root <- run_setup()
  nmfe <- fake_nmfe(root)
  expect_silent(log <- nm_run(nmfe, "r1.ctl", root, verbose = FALSE))
  expect_equal(log$action, "launched")
  expect_equal(log$exit_status, 0L)
  reg <- file.path(root, "r1.reg")
  expect_true(file.exists(file.path(reg, "r1.ctl")))
  expect_true(file.exists(file.path(reg, "r1.lst")))
  data_line <- grep("^\\$DATA", readLines(file.path(reg, "r1.ctl")), value = TRUE)
  data_path <- strsplit(trimws(data_line), "\\s+")[[1]][2]
  expect_true(file.exists(data_path))
  expect_true(log$min_ok)

  # A terminated run is skipped the second time.
  log2 <- nm_run(nmfe, "r1.ctl", root, verbose = FALSE)
  expect_equal(log2$action, "skipped_done")
})

test_that("Wings label syntax is detected only when labels replace THETA/ETA", {
  std <- write_ctl("$INPUT ID TIME DV WT", "$THETA (0,5) ; CL", " (0,50) ; V",
                   "$OMEGA 0.1 ; CL", "$PK", "  CL = THETA(1)*EXP(ETA(1))", "  V = THETA(2)",
                   "  S2 = V")
  expect_false(nmtools:::.nm_is_wings_ctl(std))
  wings <- write_ctl("$INPUT ID TIME DV", "$THETA (0,5) ;POPCL", "$OMEGA 0.1 ; eta_cl",
                     "$PK", "  CL = POPCL*EXP(eta_cl)")
  expect_true(nmtools:::.nm_is_wings_ctl(wings))
})

test_that("nm_run accepts $DATA relative to the control stream folder", {
  root <- run_setup()
  ctl <- file.path(root, "r1.ctl")
  l <- readLines(ctl); l[grepl("^\\$DATA", l)] <- "$DATA data.csv IGNORE=@"
  writeLines(l, ctl)
  log <- nm_run(fake_nmfe(root), "r1.ctl", root, verbose = FALSE)
  data_line <- grep("^\\$DATA", readLines(file.path(root, "r1.reg", "r1.ctl")), value = TRUE)
  data_path <- strsplit(trimws(data_line), "\\s+")[[1]][2]
  expect_true(file.exists(data_path))
  expect_equal(normalizePath(data_path), normalizePath(file.path(root, "data.csv")))
})

test_that("nm_run can run in the control stream folder (run_subdir = \"\")", {
  root <- run_setup()
  ctl <- file.path(root, "r1.ctl")
  l <- readLines(ctl); l[grepl("^\\$DATA", l)] <- "$DATA data.csv IGNORE=@"
  writeLines(l, ctl)
  before <- readLines(ctl)
  log <- nm_run(fake_nmfe(root), "r1.ctl", root, run_subdir = "", verbose = FALSE)
  expect_equal(log$action, "launched")
  expect_false(dir.exists(file.path(root, "r1.reg")))
  expect_true(file.exists(file.path(root, "r1.lst")))
  expect_identical(readLines(ctl), before)       # control stream left untouched
  expect_equal(nm_run(fake_nmfe(root), "r1.ctl", root, run_subdir = "",
                      verbose = FALSE)$action, "skipped_done")
  # Other suffix
  log3 <- nm_run(fake_nmfe(root), "r1.ctl", root, run_subdir = ".run", verbose = FALSE)
  expect_true(file.exists(file.path(root, "r1.run", "r1.lst")))
})

test_that("nm_run refuses a Wings label control stream in the flat layout", {
  root <- run_setup()
  writeLines(c("$PROB x", "$INPUT ID TIME AMT DV", "$DATA data.csv IGNORE=@",
               "$PRED", "  Y = POPA + EPS(1)", "$THETA 1 ;POPA", "$SIGMA 1"),
             file.path(root, "w1.ctl"))
  expect_error(nm_run(fake_nmfe(root), "w1.ctl", root, run_subdir = "", verbose = FALSE),
               "Wings label syntax")
})
