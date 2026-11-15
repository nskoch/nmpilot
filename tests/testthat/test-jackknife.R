jk_setup <- function() {
  root <- tempfile("jk")
  dir.create(root)
  writeLines(c("#ID,TIME,AMT,DV", "1,0,100,.", "1,1,.,5.10", "2,0,100,.", "2,1,.,4.00",
               "3,0,100,.", "3,1,.,6.25"), file.path(root, "data.csv"))
  writeLines(c("$PROB x", "$INPUT ID TIME AMT DV", "$DATA data.csv IGNORE=#",
               "$PRED", "  Y = THETA(1) + THETA(2)*TIME + EPS(1)", "$THETA 1", " 2",
               "$SIGMA 1", "$EST METHOD=1", "$COV", "$TABLE ID FILE=sdtab1"),
             file.path(root, "run1.ctl"))
  writeLines(jk_lst(100, c(5, 2), se = c(0.5, 0.2)), file.path(root, "run1.lst"))
  root
}

# Minimal .lst with #OBJV, the THETA vector and optionally their SE.
jk_lst <- function(ofv, theta, se = NULL) {
  th_block <- function(v) c(" THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********", "",
                            paste(sprintf("        TH %d", seq_along(v)), collapse = ""), "",
                            paste(sprintf("%12.4E", v), collapse = ""), "",
                            " OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********")
  c("0MINIMIZATION SUCCESSFUL",
    sprintf(" #OBJV:********************************************     %.3f       ****", ofv),
    th_block(theta),
    if (!is.null(se)) c(" STANDARD ERROR OF ESTIMATE", th_block(se)))
}

test_that("nm_jackknife (run = FALSE) writes one data set and control stream per subject", {
  root <- jk_setup()
  jk <- nm_jackknife(file.path(root, "run1.ctl"), run = FALSE, verbose = FALSE)
  expect_s3_class(jk, "nm_jackknife")
  expect_equal(jk$n_subjects, 3L)
  d2 <- readLines(file.path(jk$output_dir, "jk_data", "jk_2.csv"))
  expect_equal(d2, c("1,0,100,.", "1,1,.,5.10", "3,0,100,.", "3,1,.,6.25"))  # verbatim lines
  c2 <- readLines(file.path(jk$output_dir, "jk_ctl", "jk_2.ctl"))
  expect_true(any(grepl("^\\$DATA jk_2.csv IGNORE=#", c2)))
  expect_false(any(grepl("^\\$TABLE", c2)))
  expect_true(any(grepl("^\\$COV", c2)))
  expect_output(print(jk), "dry run")
  jk2 <- nm_jackknife(file.path(root, "run1.ctl"), run = FALSE, strip_cov = TRUE,
                      output_dir = file.path(root, "jk2"), verbose = FALSE)
  expect_false(any(grepl("^\\$COV", readLines(file.path(root, "jk2", "jk_ctl", "jk_1.ctl")))))
})

test_that("nm_jackknife computes and ranks the influence of each subject", {
  root <- jk_setup()
  nmfe <- file.path(root, "nmfe.bat"); writeLines("", nmfe)
  fake <- list(jk_1 = list(102, c(5.5, 2)), jk_2 = list(90, c(5, 2.2)),
               jk_3 = list(100.5, c(5, 2)))
  launched <- character(0)
  local_mocked_bindings(.nm_launch = function(nmfe_path, run_dir, run_name) {
    launched <<- c(launched, run_name)
    f <- fake[[run_name]]
    writeLines(jk_lst(f[[1]], f[[2]]), file.path(run_dir, paste0(run_name, ".lst")))
    0L
  })
  jk <- nm_jackknife(file.path(root, "run1.ctl"), nmfe_path = nmfe, verbose = FALSE)
  expect_equal(jk$ref_ofv, 100)
  expect_equal(jk$ref_theta, c(5, 2))
  expect_equal(jk$ref_se, c(0.5, 0.2))
  s <- jk$summary
  expect_equal(s$id_removed, c("2", "1", "3"))          # ranked by |dOFV|
  expect_equal(s$abs_dOFV, c(10, 2, 0.5))
  expect_equal(s$dtheta_norm, c(0.1, 0.1, 0))
  expect_equal(s$dtheta_znorm, c(1, 1, 0))
  expect_equal(s$dtheta_rel_1, c(0, 0.1, 0))
  expect_true(all(c("rank_dtheta", "rank_dtheta_z") %in% names(s)))
  expect_true(file.exists(file.path(jk$output_dir, "jackknife_summary.csv")))
  expect_true(file.exists(file.path(jk$output_dir, "figures", "abs_dOFV_bar.png")))
  expect_output(print(jk), "most influential")
  for (t in c("bar", "rank", "scatter")) expect_s3_class(plot(jk, t), "ggplot")

  # Second call: the minimised runs are skipped (resume)
  launched <- character(0)
  nm_jackknife(file.path(root, "run1.ctl"), nmfe_path = nmfe, save_figures = FALSE,
               verbose = FALSE)
  expect_length(launched, 0)
})

test_that("nm_jackknife checks its inputs", {
  root <- jk_setup()
  expect_error(nm_jackknife(file.path(root, "none.ctl")), "not found")
  expect_error(nm_jackknife(file.path(root, "run1.ctl")), "nmfe_path")
  expect_error(nm_jackknife(file.path(root, "run1.ctl"), id_col = "SUBJ", run = FALSE,
                            verbose = FALSE), "not found in \\$INPUT")
})
