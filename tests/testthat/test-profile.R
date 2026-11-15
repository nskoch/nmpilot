profile_setup <- function() {
  root <- tempfile("prof")
  dir.create(root)
  writeLines(c("ID,TIME,DV", "1,1,5"), file.path(root, "data.csv"))
  writeLines(c("$PROB x", "$INPUT ID TIME DV", "$DATA data.csv IGNORE=@", "$PRED",
               "  CL = THETA(1)", "  Y = CL + EPS(1)", "$THETA (0, 2) ; CL", " (0, 1) ; V",
               "$SIGMA 1", "$EST METHOD=1", "$COV", "$TABLE ID FILE=sdtab1"),
             file.path(root, "base.ctl"))
  writeLines("rem", file.path(root, "nmfe.bat"))
  root
}

# Fake launcher: OFV = 100 + 10 * (value - 2)^2, minimum at 2.
mock_launch <- function(nmfe_path, run_dir, run_name) {
  ctl <- readLines(file.path(run_dir, paste0(run_name, ".ctl")))
  v <- as.numeric(sub("^.*?([0-9.]+)\\s+FIX.*$", "\\1", ctl[grep("FIX", ctl)[1]], perl = TRUE))
  writeLines(c("MINIMIZATION SUCCESSFUL",
               sprintf(" #OBJV:**************     %.3f       **************", 100 + 10 * (v - 2)^2)),
             file.path(run_dir, paste0(run_name, ".lst")))
  0L
}

test_that("nm_profile runs the grid and derives the likelihood-profile CI", {
  root <- profile_setup()
  local_mocked_bindings(.nm_launch = mock_launch)
  pr <- nm_profile(file.path(root, "nmfe.bat"), file.path(root, "base.ctl"),
                   param = "CL", from = 1, to = 3, by = 0.5, quiet = TRUE)
  expect_s3_class(pr, "nm_profile")
  expect_equal(pr$target$theta_index, 1L)
  expect_equal(pr$results$value, c(1, 1.5, 2, 2.5, 3))
  expect_equal(pr$results$dOFV, 10 * (pr$results$value - 2)^2, tolerance = 1e-6)
  expect_equal(pr$min_value, 2)
  # Linear interpolation of dOFV = 3.84 between grid points.
  expect_equal(pr$ci, c(1.5 - 0.5 * (3.84 - 2.5) / 7.5, 2.5 + 0.5 * (3.84 - 2.5) / 7.5),
               tolerance = 1e-6)
  ctl1 <- readLines(file.path(pr$results$dir[1], paste0(pr$results$run[1], ".ctl")))
  expect_true(any(grepl("^\\$THETA 1.0 FIX", ctl1)))
  expect_false(any(grepl("^\\$COV|^\\$TABLE", ctl1)))
  expect_s3_class(plot(pr), "ggplot")
  expect_s3_class(plot(pr, type = "dofv"), "ggplot")
})

test_that("nm_profile (run = FALSE) only writes the control streams", {
  root <- profile_setup()
  pr <- nm_profile(file.path(root, "nmfe.bat"), file.path(root, "base.ctl"),
                   param = 2, from = 0.5, to = 1, by = 0.25, run = FALSE, quiet = TRUE)
  expect_equal(pr$target$theta_index, 2L)
  expect_equal(nrow(pr$results), 3L)
  expect_true(all(is.na(pr$results$ofv)))
  expect_true(all(file.exists(file.path(pr$results$dir, paste0(pr$results$run, ".ctl")))))
})

test_that("nm_profile checks its inputs", {
  root <- profile_setup()
  ctl <- file.path(root, "base.ctl")
  nmfe <- file.path(root, "nmfe.bat")
  expect_error(nm_profile(file.path(root, "none.bat"), ctl, "CL", 1, 2, 0.5), "launcher not found")
  expect_error(nm_profile(nmfe, file.path(root, "none.ctl"), "CL", 1, 2, 0.5), "Control file not found")
  expect_error(nm_profile(nmfe, ctl, "CL", 3, 1, 0.5, run = FALSE), "must be <=")
  expect_error(nm_profile(nmfe, ctl, "CL", 1, 2, -1, run = FALSE), "positive")
})
