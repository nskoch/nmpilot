# nm_export_plot is the single entry point for writing figures (point B9)

test_that("nm_export_plot writes png, jpeg and pdf files and rejects other types", {
  skip_if_not_installed("ggplot2")
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  d <- tempfile("exp"); on.exit(unlink(d, recursive = TRUE))
  for (ext in c("png", "jpg", "pdf")) {
    f <- file.path(d, "sub", paste0("fig.", ext))
    out <- nm_export_plot(p, f, width = 3, height = 2, dpi = 72, verbose = FALSE)
    expect_true(file.exists(f))
    expect_identical(out, f)
  }
  expect_error(nm_export_plot(p, file.path(d, "fig.svg"), verbose = FALSE),
               "Unsupported file extension")
  expect_error(nm_export_plot(data.frame(), file.path(d, "fig.png")), "ggplot")
  expect_warning(expect_null(nm_export_plot(NULL, file.path(d, "fig.png"))), "NULL")
})

test_that("nm_gof(save_docx = TRUE) writes its figures through nm_export_plot", {
  skip_if_not_installed("officer")
  skip_if_not_installed("ggplot2")
  ex <- system.file("extdata", "example", package = "nmpilot")
  skip_if(!nzchar(ex))
  run <- tempfile("gof"); dir.create(run); on.exit(unlink(run, recursive = TRUE))
  file.copy(file.path(ex, "run2.reg"), run, recursive = TRUE)
  calls <- 0L
  real <- nm_export_plot
  testthat::local_mocked_bindings(nm_export_plot = function(...) {
    calls <<- calls + 1L
    real(...)
  })
  g <- suppressMessages(nm_gof(file.path(run, "run2.reg", "sdtab2"), save_docx = TRUE))
  expect_gt(calls, 0L)
  expect_equal(calls, length(g$plots))
  expect_true(length(list.files(file.path(run, "run2.reg"), pattern = "_GOF\\.docx$")) == 1L)
})
