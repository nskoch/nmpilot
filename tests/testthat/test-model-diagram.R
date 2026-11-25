# nm_model_diagram (action 7): tree of the runs with their parent read in the
# control streams (";; Based on:") or given with `parents`.

md_project <- function() {
  td <- tempfile("md_"); dir.create(td)
  mk <- function(run, based, ofv, npar, min_ok = TRUE, cov_ok = TRUE, state = "completed") {
    writeLines(c(paste0("$PROBLEM ", run, " test model"),
                 if (!is.na(based)) based, "$PK", "CL = THETA(1)"),
               file.path(td, paste0(run, ".ctl")))
    data.frame(run = run, ofv = ofv, n_par = npar, min_ok = min_ok, cov_ok = cov_ok,
               run_state = state, dir = td, stringsAsFactors = FALSE)
  }
  r <- rbind(mk("run1", NA, 1000, 6),
             mk("run2", ";; 1. Based on: 1", 950, 8),            # PsN: number only
             mk("run3", ";; Based on: run1", 990, 7, cov_ok = FALSE),
             mk("run4", "; based on run2.ctl", 940, 9),            # Pirana-like, extension
             mk("run5", ";; 1. Based on: run2", 960, 9, min_ok = FALSE),
             mk("run10", ";; Based on: run4", NA, 10, state = "interrupted"),
             mk("run6", NA, 1100, 5))
  list(runs = r)
}

test_that("parents are read in the control streams and dOFV is relative to the parent", {
  d <- nm_model_diagram(md_project())
  t <- d$tree
  expect_s3_class(d, "nm_model_diagram")
  par <- stats::setNames(t$parent, t$run)
  expect_identical(unname(par[c("run1", "run2", "run3", "run4", "run5", "run10", "run6")]),
                   c(NA, "run1", "run1", "run2", "run2", "run4", NA))
  expect_true(all(t$parent_source[!is.na(t$parent)] == "control stream"))
  dofv <- stats::setNames(t$dofv_parent, t$run)
  expect_equal(unname(dofv[c("run2", "run3", "run4", "run5")]), c(-50, -10, -10, 10))
  expect_true(is.na(dofv[["run10"]]) && is.na(dofv[["run1"]]))
  expect_identical(t$dpar_parent[t$run == "run2"], 2L)
  # depth-first order: each run right after its parent
  expect_identical(t$run, c("run1", "run2", "run4", "run10", "run5", "run3", "run6"))
  expect_identical(t$x[match(c("run1", "run2", "run4", "run10"), t$run)], c(0, 1, 2, 3))
  st <- stats::setNames(as.character(t$status), t$run)
  expect_identical(unname(st[c("run1", "run3", "run5", "run10")]),
                   c("converged, covariance OK", "converged, no covariance",
                     "not converged", "not completed"))
  expect_identical(t$description[t$run == "run1"], "run1 test model")
  expect_s3_class(d$plot, "ggplot")
  expect_s3_class(plot(d), "ggplot")          # plot() returns the figure (action 11b)
  expect_output(print(d), "7 run\\(s\\), 2 root\\(s\\)")
})

test_that("`parents` completes, replaces or removes links, and is checked", {
  x <- md_project()
  d <- nm_model_diagram(x, parents = c(run6 = "run3", run4 = "run1"))
  t <- d$tree
  expect_identical(t$parent[t$run == "run6"], "run3")
  expect_identical(t$parent[t$run == "run4"], "run1")
  expect_identical(t$parent_source[t$run == "run4"], "argument")
  expect_equal(t$dofv_parent[t$run == "run4"], -60)
  d <- nm_model_diagram(x, parents = data.frame(run = "run2", parent = NA))
  expect_true(is.na(d$tree$parent[d$tree$run == "run2"]))
  expect_error(nm_model_diagram(x, parents = c(runX = "run1")), "Unknown run")
  expect_error(nm_model_diagram(x, parents = c(run2 = "runX")), "Unknown parent")
  expect_error(nm_model_diagram(x, parents = "run1"), "named vector")
})

test_that("circular links and unknown parents are reported", {
  x <- md_project()
  expect_warning(d <- nm_model_diagram(x, parents = c(run1 = "run4")), "Circular")
  expect_true(all(c("run1", "run2", "run4") %in% d$tree$run))
  writeLines(c("$PROBLEM x", ";; Based on: run99"), file.path(x$runs$dir[1], "run6.ctl"))
  expect_warning(d <- nm_model_diagram(x), "not found among the runs")
  expect_true(is.na(d$tree$parent[d$tree$run == "run6"]))
})

test_that("nm_model_diagram reads an analysis directory (example project)", {
  ex <- system.file("extdata", "example", package = "nmpilot")
  d <- nm_model_diagram(ex, parents = c(run2 = "run1"))
  t <- d$tree
  expect_identical(t$run, c("run1", "run2"))
  expect_equal(t$dofv_parent[2], t$ofv[2] - t$ofv[1])
  expect_true(grepl("BASE model", t$description[1]))
  f <- tempfile(fileext = ".png")
  nm_export_plot(d$plot, f, width = 5, height = 3, verbose = FALSE)
  expect_true(file.exists(f))
  expect_error(nm_model_diagram(tempfile()), "Directory not found")
})
