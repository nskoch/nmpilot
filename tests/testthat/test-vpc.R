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
  expect_error(nm_vpc("nope", root, vpc_cols, show_hist = FALSE), "nope_orig.fit not found")
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
  # NM-TRAN refuses FIX repeated inside a BLOCK: once, on the record line.
  for (x in list(sim, orig)) {
    expect_true(any(grepl("^\\$OMEGA BLOCK\\(1\\) FIX$", x)))
    expect_true(any(grepl("^\\$SIGMA BLOCK\\(1\\) FIX$", x)))
    expect_equal(sum(grepl("FIX", x[grep("^\\$(OMEGA|SIGMA)", x)[1]:length(x)])), 2L)
  }
  expect_false(any(grepl("^\\$COV", sim)))
  expect_true(any(grepl("FILE=est_sim.fit", sim, fixed = TRUE)))
  # Estimation run in a .reg directory: the VPC runs use .reg directories too.
  expect_equal(called$run_subdir, ".reg")
})

test_that("nm_vpc_ctl stops when the estimation run is missing", {
  root <- tempfile("vpcctl2")
  dir.create(root)
  expect_error(nm_vpc_ctl("est", root, nmfe_path = "nmfe.bat", verbose = FALSE),
               "not found for run 'est'")
})

test_that("nm_vpc and nm_vpc_ctl work without .reg directories (flat layout)", {
  root <- vpc_root()
  for (s in c("orig", "sim")) {
    f <- file.path(root, paste0("run1_", s, ".reg"), paste0("run1_", s, ".fit"))
    file.rename(f, file.path(root, basename(f)))
  }
  unlink(file.path(root, c("run1_orig.reg", "run1_sim.reg")), recursive = TRUE)
  g <- suppressMessages(nm_vpc("run1", root, vpc_cols, show_hist = FALSE, bin_method = 3))
  expect_s3_class(g, "ggplot")

  # Estimation run made directly in the analysis directory
  root2 <- tempfile("vpcflat")
  dir.create(root2)
  writeLines(c("$PROB x", "$INPUT ID TIME AMT DV MDV", "$DATA data.csv IGNORE=@",
               "$SUBR ADVAN1 TRANS2", "$PK", "  NSIM = 1", "  CL = THETA(1)*EXP(ETA(1))",
               "  V = THETA(2)", "$ERROR", "  Y = F + F*EPS(1)", "$THETA (0,5)", " (0,50)",
               "$OMEGA 0.1", "$SIGMA 0.04", "$EST METHOD=1 INTER",
               "$TABLE ID TIME FILE=sdtab1"), file.path(root2, "est.ctl"))
  writeLines(c("TABLE NO.     1: First Order",
               " ITERATION    THETA1       THETA2       SIGMA(1,1)   OMEGA(1,1)   OBJ",
               ext_row(-1000000000, c(4.2, 45, 0.04, 0.09, 200))), file.path(root2, "est.ext"))
  writeLines("MINIMIZATION SUCCESSFUL", file.path(root2, "est.lst"))
  called <- NULL
  local_mocked_bindings(nm_run = function(...) {
    called <<- list(...)
    data.frame(run = "mock")
  })
  nm_vpc_ctl("est", root2, nmfe_path = "nmfe.bat", n_sim = 10, verbose = FALSE)
  expect_equal(called$run_subdir, "")
  expect_true(file.exists(file.path(root2, "est_sim.ctl")))
})

test_that(".nm_vpc_write_matrix writes FIX once per BLOCK", {
  m <- .nm_vpc_write_matrix("OMEGA", c(0.068259, 0, 0.0423479, 0, 0, 0.172495))
  expect_equal(m[1], "$OMEGA BLOCK(3) FIX")
  expect_equal(sum(grepl("FIX", m)), 1L)
  expect_equal(trimws(m[2:4]), c("0.068259", "0 0.0423479", "0 0 0.172495"))
})

test_that("bin midpoints stay right when a bin between two others is empty", {
  # breaks 0-2-4-...-24: bins 5 and 7 to 11 are empty, as in the example project
  df <- data.frame(TIME = c(0.5, 1, 2, 3, 4, 6, 8, 12, 24))
  b <- .nm_vpc_assign_bins(df, seq(0, 24, by = 2))
  expect_equal(b$TTIME, c(1, 1, 1, 2, 2, 3, 4, 6, 12))
  expect_equal(b$MED_T, c(1, 1, 1, 3.5, 3.5, 6, 8, 12, 24))
  expect_false(anyNA(b$MED_T))
})

test_that("nm_vpc draws no line with missing points when a bin is empty", {
  root <- vpc_root()
  # breaks 0-2-...-12: the bins 4-6 and 8-10 are empty
  g <- suppressMessages(nm_vpc("run1", root, vpc_cols, show_hist = FALSE,
                               bin_method = 1, TTIME_seq_manual = seq(0, 12, by = 2)))
  expect_no_warning(ggplot2::ggplot_build(g))
  obs_lines <- g$layers[vapply(g$layers, function(l) inherits(l$geom, "GeomLine"), TRUE)]
  expect_equal(sort(unique(obs_lines[[1]]$data$MED_T)), c(1, 4, 8, 12))
})
