boot_setup <- function() {
  root <- tempfile("boot")
  dir.create(root)
  writeLines(c("ID,TIME,AMT,DV,SEX",
               paste(rep(1:6, each = 2), rep(c(0, 1), 6), rep(c(100, "."), 6),
                     rep(c(".", 5), 6), rep(c(0, 1), each = 6), sep = ",")),
             file.path(root, "boot.csv"))
  writeLines(c("$PROB x", "$INPUT ID TIME AMT DV SEX", "$DATA boot.csv IGNORE=@",
               "$PRED", "  Y = THETA(1) + EPS(1)", "$THETA 1", "$SIGMA 1",
               "$EST METHOD=1", "$COV", "$TABLE ID TIME FILE=sdtab1"),
             file.path(root, "b1.ctl"))
  file.path(root, "b1.ctl")
}

test_that("nm_bootstrap (run = FALSE) writes the bootstrap control stream", {
  b <- nm_bootstrap(boot_setup(), n = 5, seed = 1, run = FALSE, verbose = FALSE)
  expect_s3_class(b, "nm_bootstrap")
  expect_equal(b$n_run, 0L)
  ctl <- readLines(b$boot_ctl)
  expect_false(any(grepl("^\\$COV", ctl)))
  expect_false(any(grepl("^\\$TABLE", ctl)))
  expect_true(any(grepl("^\\$DATA b1_boot.csv IGNORE=@", ctl)))
  expect_true(any(grepl("^\\$WARN", ctl)))
  expect_true(file.exists(file.path(b$run_dir, "resample_index.rds")))
})

test_that("the resample plan is reproducible and respects strata", {
  f <- boot_setup()
  b1 <- nm_bootstrap(f, n = 20, seed = 42, strata = "SEX", run = FALSE, verbose = FALSE)
  b2 <- nm_bootstrap(f, n = 20, seed = 42, strata = "SEX", run = FALSE, verbose = FALSE)
  expect_identical(b1$plan, b2$plan)
  expect_length(b1$plan, 20L)
  for (p in b1$plan) {
    expect_length(p, 6L)
    expect_equal(sum(p %in% c("1", "2", "3")), 3L)   # SEX = 0
    expect_equal(sum(p %in% c("4", "5", "6")), 3L)   # SEX = 1
  }
})

test_that("nm_bootstrap checks its inputs", {
  f <- boot_setup()
  expect_error(nm_bootstrap(f, ci = 1.2, run = FALSE, verbose = FALSE), "ci")
  expect_error(nm_bootstrap(f, id_col = "SUBJ", run = FALSE, verbose = FALSE), "ID column")
  expect_error(nm_bootstrap(f, strata = "GRP", run = FALSE, verbose = FALSE), "strata column")
  expect_error(nm_bootstrap(f, n = 2, seed = 1, run = TRUE, verbose = FALSE), "nmfe_path")
})
