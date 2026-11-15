test_that("nm_read_csv reads a NONMEM data set", {
  r <- suppressMessages(nm_read_csv(write_data_csv()))
  expect_named(r, c("data", "mapping", "lloq"), ignore.order = TRUE)
  expect_equal(nrow(r$data), 8L)
  expect_true(all(c("ID", "TIME", "DV", "AMT", "BLQ") %in% names(r$data)))
  expect_true(all(is.na(r$data$DV[r$data$EVID == 1])))   # "." read as NA
  expect_equal(r$data$DV_norm[2], 5.2 / 100)
})

test_that("nm_summary counts subjects and observations", {
  r <- suppressMessages(nm_read_csv(write_data_csv()))
  s <- nm_summary(r$data, verbose = FALSE)
  expect_equal(s$n_subjects, 2L)
  expect_equal(s$n_records, 8L)
  expect_equal(s$n_observations, 6L)
})

test_that("nm_add_rate adds RATE after AMT", {
  out <- tempfile(fileext = ".csv")
  suppressMessages(nm_add_rate(write_data_csv(), out_file = out))
  d <- utils::read.csv(out)
  expect_equal(names(d)[which(names(d) == "AMT") + 1L], "RATE")
  expect_equal(d$RATE[d$EVID == 1], c(-2, -2))
  expect_equal(unique(d$RATE[d$EVID == 0]), 0)
})

test_that("nm_export_data writes '.' for missing values", {
  out <- file.path(tempdir(), "nmtools-export", "d.csv")
  invisible(capture.output(suppressMessages(
    nm_export_data(data.frame(ID = 1:2, DV = c(1.5, NA)), out))))
  expect_true(file.exists(out))
  expect_equal(readLines(out)[3], "2,.")
})
