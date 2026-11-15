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
  log <- suppressMessages(nm_run(nmfe, "r1.ctl", root, verbose = FALSE))
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
