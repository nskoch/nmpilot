test_that("nm_gof builds the standard diagnostic plots", {
  g <- nm_gof(write_sdtab())
  expect_s3_class(g, "nm_gof")
  expect_true(all(c("dv_pred", "dv_ipred", "cwres_time", "cwres_pred") %in% names(g$plots)))
  expect_true(all(vapply(g$plots, inherits, logical(1), "ggplot")))
  expect_equal(nrow(g$data), 6L)
  # plot() returns the figure (action 11b); `which` is the former name of `type`
  expect_s3_class(plot(g, type = "dv_pred"), "ggplot")
  expect_warning(p <- plot(g, which = "dv_pred"), "`which` is deprecated")
  expect_s3_class(p, "ggplot")
  skip_if_not_installed("patchwork")
  expect_s3_class(plot(g), "patchwork")
})

test_that("a run directory is accepted: estimation outputs are not taken for tables", {
  ex <- system.file("extdata", "example", package = "nmpilot")
  run <- tempfile("gofdir"); dir.create(run); on.exit(unlink(run, recursive = TRUE))
  file.copy(file.path(ex, "run2.reg"), run, recursive = TRUE)
  d <- file.path(run, "run2.reg")           # holds .ext, .phi, .cov, .cor, .shk and sdtab2
  expect_equal(basename(nmpilot:::.nm_locate_table(d)), "sdtab2")
  g <- nm_gof(d)
  expect_s3_class(g, "nm_gof")
  expect_equal(basename(g$file), "sdtab2")
  # a parameter table (no DV) next to sdtab: sdtab is still chosen
  writeLines(c("TABLE NO.  1", " ID CL V", " 1.0000E+00 5.0000E+00 2.0000E+01"),
             file.path(d, "patab2"))
  expect_equal(basename(nmpilot:::.nm_locate_table(d)), "sdtab2")
  # two tables with DV, none named sdtab: the user must choose
  d2 <- file.path(run, "two"); dir.create(d2)
  for (f in c("tabA", "tabB"))
    writeLines(c("TABLE NO.  1", " ID TIME DV PRED", " 1 0 1 1"), file.path(d2, f))
  expect_error(nmpilot:::.nm_locate_table(d2), "Several NONMEM tables")
  expect_equal(basename(nmpilot:::.nm_locate_table(d2, stem = "tabB")), "tabB")
})
