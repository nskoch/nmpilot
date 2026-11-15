test_that("nm_gof builds the standard diagnostic plots", {
  g <- nm_gof(write_sdtab())
  expect_s3_class(g, "nm_gof")
  expect_true(all(c("dv_pred", "dv_ipred", "cwres_time", "cwres_pred") %in% names(g$plots)))
  expect_true(all(vapply(g$plots, inherits, logical(1), "ggplot")))
  expect_equal(nrow(g$data), 6L)
})
