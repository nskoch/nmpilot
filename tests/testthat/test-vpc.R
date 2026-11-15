vpc_root <- function() {
  root <- tempfile("vpc")
  dir.create(root)
  set.seed(1)
  tab <- function(nsim) {
    tm <- rep(c(0.5, 1, 2, 4, 8, 12), 10)
    id <- rep(1:10, each = 6)
    do.call(rbind, lapply(seq_len(nsim), function(s)
      data.frame(NSIM = s, ID = id, TIME = tm, AMT = 0, MDV = 0,
                 DV = 10 * exp(-0.2 * tm) * exp(stats::rnorm(60, 0, 0.2)),
                 PRED = 10 * exp(-0.2 * tm), IPRED = 10 * exp(-0.2 * tm))))
  }
  write_fit <- function(df, path) {
    dir.create(dirname(path), showWarnings = FALSE)
    writeLines(c("TABLE NO.  1",
                 apply(df, 1, function(r) paste(sprintf("%12.4E", r), collapse = ""))), path)
  }
  write_fit(tab(1), file.path(root, "run1_orig.reg", "run1_orig.fit"))
  write_fit(tab(20), file.path(root, "run1_sim.reg", "run1_sim.fit"))
  root
}
vpc_cols <- c("NSIM", "ID", "TIME", "AMT", "MDV", "DV", "PRED", "IPRED")

test_that("nm_vpc builds VPC and pc-VPC plots", {
  root <- vpc_root()
  g1 <- suppressMessages(nm_vpc("run1", root, vpc_cols, show_hist = FALSE, bin_method = 3))
  expect_s3_class(g1, "ggplot")
  expect_match(g1$labels$title, "^VPC")
  g2 <- suppressMessages(nm_vpc("run1", root, vpc_cols, type = 2, show_hist = FALSE,
                                bin_method = 2, n_bins_quantile = 4))
  expect_match(g2$labels$title, "^pc-VPC")
})

test_that("nm_vpc checks its arguments and inputs", {
  root <- vpc_root()
  expect_error(nm_vpc("run1", root, vpc_cols, type = 3), "type")
  expect_error(nm_vpc("run1", root, vpc_cols, bin_method = 5), "bin_method")
  expect_error(nm_vpc("nope", root, vpc_cols, show_hist = FALSE), "_orig.fit not found")
})

test_that("nm_vpc_ctl writes _orig and _sim control streams and calls nm_run", {
  root <- tempfile("vpcctl")
  reg <- file.path(root, "est.reg")
  dir.create(reg, recursive = TRUE)
  writeLines(c("$PROB x", "$INPUT ID TIME AMT DV MDV", "$DATA ../data.csv IGNORE=@",
               "$SUBR ADVAN1 TRANS2", "$PK", "  NSIM = 1", "  CL = THETA(1)*EXP(ETA(1))",
               "  V = THETA(2)", "$ERROR", "  Y = F + F*EPS(1)", "$THETA (0,5)", " (0,50)",
               "$OMEGA 0.1", "$SIGMA 0.04", "$EST METHOD=1 INTER", "$COV",
               "$TABLE ID TIME FILE=sdtab1"), file.path(reg, "est.ctl"))
  writeLines(c("TABLE NO.     1: First Order",
               " ITERATION    THETA1       THETA2       SIGMA(1,1)   OMEGA(1,1)   OBJ",
               ext_row(-1000000000, c(4.2, 45, 0.04, 0.09, 200))), file.path(reg, "est.ext"))
  writeLines("MINIMIZATION SUCCESSFUL", file.path(reg, "est.lst"))

  called <- NULL
  local_mocked_bindings(nm_run = function(...) {
    called <<- list(...)
    data.frame(run = "mock")
  })
  nm_vpc_ctl("est", root, nmfe_path = "nmfe.bat", n_sim = 100, seed = 7, verbose = FALSE)

  expect_equal(basename(called$ctl_files), c("est_orig.ctl", "est_sim.ctl"))
  sim <- readLines(file.path(root, "est_sim.ctl"))
  orig <- readLines(file.path(root, "est_orig.ctl"))
  expect_true(any(grepl("$SIMULATION (7) ONLYSIM SUBPROBLEMS=100", sim, fixed = TRUE)))
  expect_true(any(grepl("NSIM = IREP", sim, fixed = TRUE)))
  expect_true(any(grepl("$ESTIMATION MAXEVAL=0", orig, fixed = TRUE)))
  expect_true(any(grepl("NSIM = 1", orig, fixed = TRUE)))
  expect_true(any(grepl("4.2 FIX", sim, fixed = TRUE)))
  expect_false(any(grepl("^\\$COV", sim)))
  expect_true(any(grepl("FILE=est_sim.fit", sim, fixed = TRUE)))
})

test_that("nm_vpc_ctl stops when the estimation run is missing", {
  root <- tempfile("vpcctl2")
  dir.create(root)
  expect_error(nm_vpc_ctl("est", root, nmfe_path = "nmfe.bat", verbose = FALSE),
               ".reg directory not found")
})
