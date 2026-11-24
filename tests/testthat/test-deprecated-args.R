# Action 11 (interface harmonised before 1.0): the former argument names and
# argument orders are still accepted, with a warning.

test_that(".nm_deprecated_args maps old names and rejects unknown ones", {
  expect_identical(.nm_deprecated_args(list(), c(a = "b"), "f"), list())
  expect_warning(z <- .nm_deprecated_args(list(a = 1), c(a = "b"), "f"),
                 "In f\\(\\), `a` is deprecated; use `b` instead")
  expect_identical(z, list(b = 1))
  expect_error(.nm_deprecated_args(list(x = 1), c(a = "b"), "f"), "unused argument")
  expect_error(.nm_deprecated_args(list(1), c(a = "b"), "f"), "Unnamed")
})

test_that(".nm_run_arg reads a run given by its control stream, directory or stem", {
  expect_identical(.nm_run_arg("C:/p/run2.ctl"), list(analysis_dir = "C:/p", run_name = "run2"))
  expect_identical(.nm_run_arg("C:/p/run2.reg/run2.ctl"), list(analysis_dir = "C:/p", run_name = "run2"))
  expect_identical(.nm_run_arg("C:\\p\\run2.nm7\\"), list(analysis_dir = "C:/p", run_name = "run2"))
  expect_identical(.nm_run_arg("C:/p/run2"), list(analysis_dir = "C:/p", run_name = "run2"))
  expect_identical(.nm_run_arg("C:/p/run2.reg"), list(analysis_dir = "C:/p", run_name = "run2"))
  expect_error(.nm_run_arg(NA_character_), "path of a run")
})

test_that("nm_identify_parameters, nm_convert_ctl: `file` and `ctl_file` still work", {
  f <- write_ctl("$PK", "CL = THETA(1)*EXP(ETA(1))", "$THETA 1", "$OMEGA 0.1")
  a <- nm_identify_parameters(f)
  expect_warning(b <- nm_identify_parameters(file = f), "`file` is deprecated")
  expect_identical(a$parameter_table, b$parameter_table)
  expect_error(nm_identify_parameters(f, foo = 1), "unused argument")
  out <- tempfile(fileext = ".ctl")
  expect_warning(nm_convert_ctl(ctl_file = f, out_file = out, verbose = FALSE), "`ctl_file` is deprecated")
  expect_true(file.exists(out))
})

test_that("nm_run and nm_profile: former names and former order (nmfe_path first)", {
  root <- tempfile("dep_"); dir.create(root)
  writeLines(c("$PROBLEM t", "$INPUT ID TIME DV", "$DATA d.csv", "$PRED", "Y = THETA(1) + ETA(1) + EPS(1)",
               "$THETA (0, 5)", "$OMEGA 0.1", "$SIGMA 0.1", "$ESTIMATION METHOD=1"),
             file.path(root, "base.ctl"))
  writeLines("ID,TIME,DV", file.path(root, "d.csv"))
  nmfe <- file.path(root, "nmfe76.bat"); writeLines("rem", nmfe)
  # nm_run: missing control stream -> skipped, no launch needed
  new <- nm_run("absent.ctl", nmfe, root, verbose = FALSE)
  expect_warning(old <- nm_run(nmfe, "absent.ctl", root, verbose = FALSE), "former order")
  expect_identical(old$action, new$action)
  expect_warning(old2 <- nm_run(nmfe_path = nmfe, ctl_files = "absent.ctl", analysis_dir = root,
                                verbose = FALSE), "`ctl_files` is deprecated")
  expect_identical(old2$action, new$action)
  # nm_profile with run = FALSE: control streams only
  p_new <- nm_profile(file.path(root, "base.ctl"), nmfe, param = 1, from = 4, to = 5, by = 0.5,
                      run = FALSE, verbose = FALSE)
  expect_warning(expect_warning(
    p_old <- nm_profile(nmfe, ctl_path = file.path(root, "base.ctl"), param = 1, from = 4, to = 5,
                        by = 0.5, run = FALSE, quiet = TRUE, overwrite = TRUE),
    "`ctl_path` is deprecated"), "`quiet` is deprecated")
  expect_identical(p_old$results$value, p_new$results$value)
  expect_warning(nm_profile(nmfe, file.path(root, "base.ctl"), param = 1, from = 4, to = 5, by = 0.5,
                            run = FALSE, verbose = FALSE, overwrite = TRUE), "former order")
})

test_that("nm_bootstrap, nm_derivatives_ctl, nm_data_previsual: former output names", {
  root <- tempfile("depb_"); dir.create(root)
  writeLines(c("$PROBLEM t", "$INPUT ID TIME DV", "$DATA d.csv IGNORE=@", "$PRED",
               "Y = THETA(1) + ETA(1) + EPS(1)", "$THETA (0, 5)", "$OMEGA 0.1", "$SIGMA 0.1",
               "$ESTIMATION METHOD=1", "$TABLE ID FILE=sdtab"), file.path(root, "base.ctl"))
  writeLines(c("ID,TIME,DV", "1,0,1", "1,1,2", "2,0,1", "2,1,3"), file.path(root, "d.csv"))
  out <- file.path(root, "bt")
  # `strip_table` is matched to `strip_tables` by R's partial matching
  expect_warning(
    b <- nm_bootstrap(file.path(root, "base.ctl"), n = 2, seed = 1, run = FALSE, verbose = FALSE,
                      run_dir = out, strip_table = FALSE),
    "`run_dir` is deprecated")
  expect_true(dir.exists(out))
  expect_true(any(grepl("^\\$TABLE", readLines(b$boot_ctl))))
})

test_that("nm_vpc_ctl, nm_vpc, nm_extract_run_table: former run_name + analysis_dir forms", {
  expect_warning(expect_error(nm_vpc("nope", tempdir(), col_names = "ID", show_hist = FALSE)),
                 "`analysis_dir` is deprecated")
  expect_warning(expect_error(nm_vpc(run_name = file.path(tempdir(), "nope"), col_names = "ID",
                                     show_hist = FALSE)), "`run_name` is deprecated")
  skip_if_not_installed("flextable"); skip_if_not_installed("officer")
  expect_warning(expect_error(nm_extract_run_table(tempdir(), "none")), "`run_rep` \\+ `run_name` is deprecated")
  expect_warning(expect_error(nm_extract_run_table(run_path = file.path(tempdir(), "none"))),
                 "`run_path` is deprecated")
})
