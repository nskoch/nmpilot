make_run <- function(root, name, theta, ofv) {
  d <- file.path(root, name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    "$PROB test", "$INPUT ID TIME AMT DV MDV WT", "$DATA data.csv IGNORE=@",
    "$SUBR ADVAN1 TRANS2",
    "$PK",
    "  CL = THETA(1)*(WT/70)**0.75*EXP(ETA(1))",
    "  V  = THETA(2)*EXP(ETA(2))",
    "  KE = CL/V",
    "  S1 = V",
    "$ERROR", "  Y = F + F*THETA(3)*EPS(1)",
    "$THETA (0,5) ; CL", " (0,50) ; V", " (0,0.2) ; PROP",
    "$OMEGA 0.1 ; CL", " 0.1 ; V", "$SIGMA 1 FIX", "$EST METHOD=1 INTER"),
    file.path(d, paste0(name, ".ctl")))
  fmt <- function(v) paste(sprintf("%13.5E", v), collapse = "")
  write_ext(file.path(d, paste0(name, ".ext")),
    " ITERATION    THETA1       THETA2       THETA3       SIGMA(1,1)   OMEGA(1,1)   OMEGA(2,1)   OMEGA(2,2)   OBJ",
    c(paste0("  -1000000000", fmt(c(theta, 1, 0.09, 0, 0.04, ofv))),
      paste0("  -1000000001", fmt(c(theta / 10, 1e10, 0.018, 1e10, 0.01, 0))),
      paste0("  -1000000006", fmt(c(0, 0, 0, 1, 0, 1, 0, 0)))))
  d
}

runs_root <- function() {
  root <- tempfile("runs")
  make_run(root, "run1", c(4.2, 45, 0.15), 200)
  make_run(root, "run2", c(4.0, 40, 0.14), 190)
  root
}

test_that("nm_extract_runs reads .ext values and names parameters from the .ctl", {
  x <- nm_extract_runs(runs_root(), verbose = FALSE)
  expect_s3_class(x, "nm_runs")
  expect_equal(nrow(x$runs), 2L)
  p <- x$params[x$params$run == "run1", ]
  th <- p[p$type == "THETA", ]
  expect_equal(th$name[th$index == 1], "CL")    # not KE (derived)
  expect_equal(th$name[th$index == 2], "V")     # not S1 / KE
  expect_equal(th$name[th$index == 3], "PROP")  # residual error: $THETA label
  expect_equal(th$value[th$index == 1], 4.2)
  expect_equal(th$rse_pct[th$index == 2], 10)
  om <- p[p$type == "OMEGA", ]
  expect_equal(om$name[om$index == 1], "CL")
  expect_true(p$fixed[p$type == "SIGMA"])
})

test_that("OFV differences are computed against the best run", {
  x <- nm_extract_runs(runs_root(), verbose = FALSE)
  r <- x$runs[order(x$runs$run), ]
  expect_equal(r$ofv, c(200, 190))
  expect_equal(r$dOFV, c(10, 0))
  expect_equal(r$rank_ofv, c(2L, 1L))
})

test_that("nm_compare_runs aligns parameters across runs", {
  x <- nm_extract_runs(runs_root(), verbose = FALSE)
  cmp <- nm_compare_runs(x)
  expect_true(is.data.frame(cmp))
  expect_gte(ncol(cmp), 3L)
})

test_that("nm_extract_runs fails clearly when no run is found", {
  d <- tempfile("empty")
  dir.create(d)
  expect_error(nm_extract_runs(d, verbose = FALSE), "No NONMEM run")
})

test_that("residual-error THETAs are named from short $THETA labels only", {
  root <- tempfile("runs")
  d <- make_run(root, "run1", c(4.2, 45, 0.15), 200)
  ctl <- file.path(d, "run1.ctl")
  set_label <- function(lab) {
    l <- readLines(ctl); l[grepl("(0,0.2)", l, fixed = TRUE)] <- paste0(" (0,0.2) ; ", lab); writeLines(l, ctl)
    p <- nm_extract_runs(root, verbose = FALSE)$params
    p$name[p$type == "THETA" & p$index == 3]
  }
  expect_equal(set_label("3 Prop (sigma_prop)"), "Prop")
  expect_equal(set_label("SD1 = proportional SD"), "SD1")
  expect_equal(set_label("erreur proportionnelle"), "THETA3")   # free text
  expect_equal(set_label("CL"), "THETA3")                       # already a PK name
})

test_that("covariate-effect THETAs are named from their $THETA label (example project)", {
  ex <- system.file("extdata", "example", package = "nmtools")
  skip_if(!nzchar(ex))
  p <- nm_extract_runs(ex, verbose = FALSE)$params
  th <- p[p$run == "run2" & p$type == "THETA", ]
  expect_equal(th$name[order(th$index)], c("CL", "V", "KA", "PROP", "ADD", "CL_WT", "CL_GENO"))
})
