compare_dirs <- function() {
  root <- tempfile("cmp")
  a <- file.path(root, "runA"); b <- file.path(root, "runB")
  dir.create(a, recursive = TRUE); dir.create(b)
  write_sdtab(file.path(a, "sdtab1"))
  write_sdtab(file.path(b, "sdtab1"), shift = 0.1)
  c(a, b)
}

test_that("nm_compare pairs observations and computes deltas", {
  d <- compare_dirs()
  cmp <- nm_compare(d[1], d[2], labels = c("A", "B"))
  expect_s3_class(cmp, "nm_compare")
  expect_named(cmp, c("data_obs", "data_subj", "metrics", "meta"), ignore.order = TRUE)
  expect_equal(nrow(cmp$data_obs), 6L)
  expect_equal(cmp$data_obs$dCWRES, rep(0.1, 6), tolerance = 1e-6)
  expect_equal(nrow(cmp$data_subj), 2L)
  m <- cmp$metrics
  expect_equal(m$run1[m$metric == "Mean |CWRES|"], 0.8 / 3, tolerance = 1e-6)
})

test_that("plot.nm_compare returns plots", {
  d <- compare_dirs()
  cmp <- nm_compare(d[1], d[2])
  p <- plot(cmp)
  expect_true(inherits(p, "ggplot") || is.list(p))
  expect_output(print(cmp))
})
