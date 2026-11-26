ex <- system.file("extdata", "example", package = "nmpilot")

test_that("nm_runs_table gives one row per run with parent, dOFV and status", {
  tab <- nm_runs_table(ex, parents = c(run2 = "run1"))
  expect_s3_class(tab, "nm_runs_table")
  expect_s3_class(tab, "data.frame")
  expect_identical(tab$run, c("run1", "run2"))             # tree order: parent first
  expect_identical(tab$parent, c(NA, "run1"))
  expect_equal(tab$dofv_parent[2], tab$ofv[2] - tab$ofv[1])
  expect_true(is.na(tab$dofv_parent[1]))
  expect_equal(min(tab$dofv_best), 0)
  expect_equal(tab$dofv_best, tab$ofv - min(tab$ofv))
  expect_true(all(c("description", "n_par", "aic", "bic", "minimization", "covariance",
                    "max_eta_shrink", "pct_rse_gt50", "status") %in% names(tab)))
  expect_match(tab$description[1], "BASE")

  # same values as nm_extract_runs() and nm_model_diagram()
  r <- nm_extract_runs(ex, verbose = FALSE)$runs
  i <- match(tab$run, r$run)
  expect_equal(tab$ofv, r$ofv[i])
  expect_equal(tab$bic, r$bic[i])
  expect_equal(tab$n_par, as.integer(r$n_par[i]))
  d <- nm_model_diagram(ex, parents = c(run2 = "run1"))$tree
  expect_identical(tab$status, as.character(d$status[match(tab$run, d$run)]))

  # the result of nm_extract_runs() is accepted as well
  tab2 <- nm_runs_table(nm_extract_runs(ex, verbose = FALSE), parents = c(run2 = "run1"))
  expect_equal(as.data.frame(tab2), as.data.frame(tab))

  expect_output(print(tab), "Summary of 2 run")
})

test_that("nm_runs_table sorts the rows and checks its arguments", {
  expect_identical(nm_runs_table(ex, sort = "ofv")$run, c("run2", "run1"))
  expect_identical(nm_runs_table(ex, sort = "run")$run, c("run1", "run2"))
  expect_true(all(is.na(nm_runs_table(ex)$parent)))         # no link written in the example
  expect_error(nm_runs_table(file.path(tempdir(), "no_such_dir")), "not found")
  expect_error(nm_runs_table(data.frame(a = 1)), "nm_extract_runs")
  expect_error(nm_runs_table(ex, output_file = file.path(tempdir(), "x.txt")), ".docx")
  expect_error(nm_runs_table(ex, sort = "dofv"))
})

test_that("nm_runs_table writes the Word file and refuses to overwrite it", {
  skip_if_not_installed("flextable")
  skip_if_not_installed("officer")
  d <- tempfile("rt"); on.exit(unlink(d, recursive = TRUE))
  f <- file.path(d, "sub", "runs.docx")
  expect_message(tab <- nm_runs_table(ex, parents = c(run2 = "run1"), output_file = f),
                 "Word file written")
  expect_true(file.exists(f))
  expect_identical(attr(tab, "output_file"), f)
  expect_error(nm_runs_table(ex, output_file = f), "overwrite")
  expect_silent(nm_runs_table(ex, output_file = f, overwrite = TRUE, verbose = FALSE))
  doc <- officer::read_docx(f)
  txt <- officer::docx_summary(doc)$text
  expect_true(any(grepl("Summary of the runs", txt)))
  cells <- officer::docx_summary(doc)
  expect_true(all(c("run1", "run2") %in% cells$text))
})
