test_that("nm_plot_dv returns a ggplot (and a patchwork for both scales)", {
  dat <- suppressMessages(nm_read_csv(write_data_csv()))$data
  expect_s3_class(nm_plot_dv(dat, verbose = FALSE), "ggplot")
  skip_if_not_installed("patchwork")
  expect_s3_class(nm_plot_dv(dat, scales = "both", verbose = FALSE), "patchwork")
})

test_that("nm_data_previsual reads, summarises and plots a data set", {
  pv <- nm_data_previsual(write_data_csv(), verbose = FALSE)
  expect_true(all(c("data", "plot_dv", "summary") %in% names(pv)))
  expect_equal(nrow(pv$data), 8L)
  expect_s3_class(pv$plot_dv, "ggplot")
  expect_s3_class(pv, "nm_data_previsual")
  expect_output(print(pv), "Data set preview: 8 rows")
})

test_that("nm_data_previsual exports files only when asked", {
  out <- tempfile("pv")
  pv <- nm_data_previsual(write_data_csv(), export_jpg = TRUE, export_csv = TRUE,
                          output_dir = out, verbose = FALSE)
  expect_gt(length(list.files(out)), 0L)
})

test_that("nm_export_plot writes the image file", {
  f <- file.path(tempfile("plots"), "p.png")
  nm_export_plot(ggplot2::ggplot(), f, width = 3, height = 2, dpi = 72, verbose = FALSE)
  expect_true(file.exists(f))
})
