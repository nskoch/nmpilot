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
  expect_error(nm_run("r1.ctl", file.path(root, "none.bat"), root), "nmfe executable not found")
  expect_error(nm_run("r1.ctl", nmfe, file.path(root, "nodir")), "does not exist")
  expect_error(nm_run(character(0), nmfe, root), "No .ctl")
  log <- nm_run("absent.ctl", nmfe, root, verbose = FALSE)
  expect_equal(log$action, "skipped_no_ctl")
})

test_that("nm_run copies the control stream, fixes $DATA and launches nmfe", {
  root <- run_setup()
  nmfe <- fake_nmfe(root)
  expect_silent(log <- nm_run("r1.ctl", nmfe, root, verbose = FALSE))
  expect_s3_class(log, c("nm_run_log", "data.frame"))
  expect_output(print(log), "NONMEM runs: 1 \\(launched 1")
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
  log2 <- nm_run("r1.ctl", nmfe, root, verbose = FALSE)
  expect_equal(log2$action, "skipped_done")
})

test_that("Wings label syntax is detected only when labels replace THETA/ETA", {
  std <- write_ctl("$INPUT ID TIME DV WT", "$THETA (0,5) ; CL", " (0,50) ; V",
                   "$OMEGA 0.1 ; CL", "$PK", "  CL = THETA(1)*EXP(ETA(1))", "  V = THETA(2)",
                   "  S2 = V")
  expect_false(nmpilot:::.nm_is_wings_ctl(std))
  wings <- write_ctl("$INPUT ID TIME DV", "$THETA (0,5) ;POPCL", "$OMEGA 0.1 ; eta_cl",
                     "$PK", "  CL = POPCL*EXP(eta_cl)")
  expect_true(nmpilot:::.nm_is_wings_ctl(wings))
})

test_that("nm_run accepts $DATA relative to the control stream folder", {
  root <- run_setup()
  ctl <- file.path(root, "r1.ctl")
  l <- readLines(ctl); l[grepl("^\\$DATA", l)] <- "$DATA data.csv IGNORE=@"
  writeLines(l, ctl)
  log <- nm_run("r1.ctl", fake_nmfe(root), root, verbose = FALSE)
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
  log <- nm_run("r1.ctl", fake_nmfe(root), root, run_subdir = "", verbose = FALSE)
  expect_equal(log$action, "launched")
  expect_false(dir.exists(file.path(root, "r1.reg")))
  expect_true(file.exists(file.path(root, "r1.lst")))
  expect_identical(readLines(ctl), before)       # control stream left untouched
  expect_equal(nm_run("r1.ctl", fake_nmfe(root), root, run_subdir = "",
                      verbose = FALSE)$action, "skipped_done")
  # Other suffix
  log3 <- nm_run("r1.ctl", fake_nmfe(root), root, run_subdir = ".run", verbose = FALSE)
  expect_true(file.exists(file.path(root, "r1.run", "r1.lst")))
})

test_that("nm_run refuses a Wings label control stream in the flat layout", {
  root <- run_setup()
  writeLines(c("$PROB x", "$INPUT ID TIME AMT DV", "$DATA data.csv IGNORE=@",
               "$PRED", "  Y = POPA + EPS(1)", "$THETA 1 ;POPA", "$SIGMA 1"),
             file.path(root, "w1.ctl"))
  expect_error(nm_run("w1.ctl", fake_nmfe(root), root, run_subdir = "", verbose = FALSE),
               "Wings label syntax")
})

test_that("OFV and status of an evaluation run (MAXEVAL=0) are read", {
  lst <- tempfile(fileext = ".lst")
  writeLines(c(
    " ********************         FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION (EVALUATION)       *******************",
    " #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      *******************",
    "",
    " #OBJV:********************************************         -1147.341       ********************************************",
    " THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********"), lst)
  expect_equal(nmpilot:::.nm_lst_ofv(lst), -1147.341)
  st <- nmpilot:::.nm_lst_status(lst)
  expect_equal(st$minimization_status, "evaluation")
  expect_true(st$run_ok)
  expect_true(nmpilot:::.nm_is_run_terminal(lst))
})

test_that("status of a simulation-only run ($SIMULATION ONLYSIM) is read", {
  sim_lines <- c(
    " PROBLEM NO.:           1      SUBPROBLEM NO.:         500",
    "",
    " SIMULATION STEP PERFORMED",
    " SOURCE  1:",
    "    SEED1:      77690045   SEED2:             0",
    " Elapsed simulation  time in seconds:     0.01",
    " ESTIMATION STEP OMITTED:                 YES")
  done <- tempfile(fileext = ".lst")
  writeLines(c(sim_lines, " #CPUT: Total CPU Time in Seconds,        9.953",
               "Stop Time: ", "21/11/2026 ", "18:33"), done)
  st <- nmpilot:::.nm_lst_status(done)
  expect_equal(st$minimization_status, "simulation")
  expect_true(st$run_ok)
  expect_true(is.na(st$min_ok))
  expect_true(is.na(nmpilot:::.nm_lst_ofv(done)))
  expect_true(nmpilot:::.nm_is_run_terminal(done))
  expect_true(nmpilot:::.nm_is_run_terminal(done, accept_failed = FALSE))

  # output cut before NONMEM's final lines: not finished, not terminal
  cut <- tempfile(fileext = ".lst")
  writeLines(sim_lines, cut)
  st2 <- nmpilot:::.nm_lst_status(cut)
  expect_true(is.na(st2$minimization_status))
  expect_false(nmpilot:::.nm_is_run_terminal(cut))
})

test_that("a synchronised analysis folder triggers a warning", {
  expect_warning(nmpilot:::.nm_warn_synced_dir("D:/Dropbox/nm76g64/run/projet"), "synchronised")
  expect_warning(nmpilot:::.nm_warn_synced_dir("C:/Users/n/OneDrive - AP-HM/projet"), "synchronised")
  expect_silent(nmpilot:::.nm_warn_synced_dir("C:/nm76g64/run/projet"))
  old <- options(nmpilot.warn_synced = FALSE); on.exit(options(old))
  expect_silent(nmpilot:::.nm_warn_synced_dir("D:/Dropbox/projet"))
})
