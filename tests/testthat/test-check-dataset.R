# nm_check_dataset: NONMEM compatibility of a dataset (point B5)

write_csv_lines <- function(lines) {
  f <- tempfile(fileext = ".csv"); writeLines(lines, f); f
}

test_that("the example dataset is readable by NONMEM (no error)", {
  ex <- system.file("extdata", "example", package = "nmpilot")
  chk <- nm_check_dataset(file.path(ex, "example_data.csv"), verbose = FALSE)
  expect_s3_class(chk, "nm_dataset_check")
  expect_true(chk$ok)
  expect_equal(chk$n_rows, 600)
  expect_equal(chk$n_subjects, 60)
  expect_equal(unname(chk$columns[c("ID", "TIME", "DV")]), c("ID", "TIME", "DV"))   # "#ID" header
  expect_false(any(chk$findings$severity == "error"))
  # the data frame of nm_read_csv() gives the same result
  d <- suppressMessages(nm_read_csv(file.path(ex, "example_data.csv")))$data
  chk2 <- nm_check_dataset(d, verbose = FALSE)
  expect_true(chk2$ok)
  expect_equal(nrow(chk2$findings), nrow(chk$findings))
})

test_that("every problem of the checklist is reported with its rows", {
  f <- write_csv_lines(c(
    "ID,TIME,AMT,DV,EVID,MDV,CMT,SEX",
    "1,0,100,.,1,1,1,0",        # 1 dose
    "1,1,.,5,0,0,2,0",          # 2
    "1,0.5,.,4,0,0,2,0",        # 3 time decreases
    "1,2,.,.,0,0,2,0",          # 4 MDV 0 without DV
    "2,0,.,.,1,1,1,F",          # 5 dose without AMT; non-numeric SEX
    "2,1,50,3,0,0,2,0",         # 6 observation with AMT
    "2,2,.,3,0,1,2,0",          # 7 MDV 1 with DV (warning)
    "2,3,.,2,7,0,2,0",          # 8 EVID 7
    "3,0,100,.,1,0,1,0",        # 9 dose with MDV 0
    "3,1,.,6,0,0,2,0",          # 10 duplicate of 11
    "3,1,.,6,0,0,2,0",          # 11
    "3,2,.,7,0,0,2,0",          # 12 same time/CMT as 13, other DV
    "3,2,.,8,0,0,2,0",          # 13
    ".,3,.,1,0,0,2,0",          # 14 missing ID
    "3,.,.,1,0,0,2,0",          # 15 missing TIME
    "3,4,.,-1,0,0,2,0"))        # 16 negative DV
  chk <- nm_check_dataset(f, verbose = FALSE)
  fd <- chk$findings
  rows_of <- function(check, pattern) {
    k <- fd$check == check & grepl(pattern, fd$message)
    expect_true(any(k), info = paste(check, pattern))
    sort(unlist(fd$rows[k]))
  }
  expect_false(chk$ok)
  expect_equal(rows_of("time_order", "TIME decreases"), 3L)
  expect_equal(rows_of("missing_values", "Missing ID"), 14L)
  expect_equal(rows_of("missing_values", "Missing TIME"), 15L)
  expect_equal(rows_of("non_numeric", "SEX"), 5L)
  expect_equal(rows_of("evid_mdv_amt", "EVID must be"), 8L)
  expect_equal(rows_of("evid_mdv_amt", "without a positive AMT"), 5L)
  expect_equal(rows_of("evid_mdv_amt", "EVID 0\\) with a positive AMT"), 6L)
  expect_true(4L %in% rows_of("evid_mdv_amt", "MDV 0 with a missing DV"))
  expect_true(9L %in% rows_of("evid_mdv_amt", "not an observation"))   # row 8 too (EVID 7)
  expect_equal(rows_of("evid_mdv_amt", "value is ignored"), 7L)
  expect_equal(rows_of("duplicates", "Identical"), c(10L, 11L))
  expect_true(all(c(12L, 13L) %in% rows_of("duplicates", "different DV")))
  expect_equal(rows_of("outliers", "Negative"), 16L)
  # severities
  expect_equal(unique(fd$severity[fd$check %in% c("missing_values", "non_numeric", "time_order")]), "error")
  expect_equal(unique(fd$severity[fd$check %in% c("duplicates", "outliers")]), "warning")
  # errors are listed before warnings
  expect_true(all(diff(match(fd$severity, c("error", "warning"))) >= 0))
  # print and summary
  out <- capture.output(print(chk))
  expect_true(any(grepl("to fix before NONMEM", out)))
  s <- summary(chk)
  expect_true(all(c("check", "severity", "n_findings", "n_rows") %in% names(s)))
  expect_equal(sum(s$n_findings), nrow(fd))
})

test_that("missing columns, reset records, columns given with `cols`, data never modified", {
  chk <- nm_check_dataset(data.frame(TIME = 0:1, DV = c(1, 2)), verbose = FALSE)
  expect_false(chk$ok)
  expect_true(any(chk$findings$check == "columns" & grepl("No ID", chk$findings$message)))
  expect_true(any(grepl("No AMT", chk$findings$message)))
  # TIME restarts at a reset record (EVID 4): no error
  d <- data.frame(ID = 1, TIME = c(0, 1, 0, 1), AMT = c(10, NA, 10, NA), DV = c(NA, 2, NA, 2),
                  EVID = c(1, 0, 4, 0), MDV = c(1, 0, 1, 0))
  expect_true(nm_check_dataset(d, verbose = FALSE)$ok)
  # non-standard names given with `cols`
  d2 <- data.frame(SUBJ = c(1, 1), T = c(0, 1), CONC = c(NA, 2), AMT = c(10, NA))
  chk2 <- nm_check_dataset(d2, cols = c(ID = "SUBJ", TIME = "T", DV = "CONC"), verbose = FALSE)
  expect_true(chk2$ok)
  expect_error(nm_check_dataset(d2, cols = c(ID = "NOPE")), "not found")
  expect_error(nm_check_dataset(d2, cols = c(FOO = "SUBJ")), "standard items")
  # report only: the input file is unchanged
  f <- write_csv_lines(c("ID,TIME,AMT,DV", "1,1,.,5", "1,0,10,."))
  before <- readLines(f)
  invisible(nm_check_dataset(f, verbose = FALSE))
  expect_identical(readLines(f), before)
  expect_error(nm_check_dataset(file.path(tempdir(), "none.csv")), "not found")
})

test_that("an outlier among observations taken at the same time is reported", {
  set.seed(3)
  d <- do.call(rbind, lapply(1:20, function(i) data.frame(
    ID = i, TIME = c(0, 1, 4, 12), AMT = c(100, NA, NA, NA),
    DV = c(NA, exp(rnorm(3, log(c(5, 3, 1)), 0.2))), EVID = c(1, 0, 0, 0), MDV = c(1, 0, 0, 0))))
  chk <- nm_check_dataset(d, verbose = FALSE)
  expect_false(any(chk$findings$check == "outliers"))      # falling profile: no outlier
  d$DV[d$ID == 5 & d$TIME == 4] <- 300
  chk <- nm_check_dataset(d, verbose = FALSE)
  r <- unlist(chk$findings$rows[chk$findings$check == "outliers"])
  expect_equal(r, which(d$ID == 5 & d$TIME == 4))
})
