# nmpilot (development version)

* New `nm_runs_table()` (action 8): summary table of all the runs of a
  project, one row per run with its description (`$PROBLEM`), its parent run
  (read in the control stream as in `nm_model_diagram()`, or given with
  `parents`), OFV, dOFV against the parent and against the best run, number
  of estimated parameters, AIC, BIC, minimisation and covariance status,
  largest ETA shrinkage and share of RSE above 50%. Written to a Word file in
  landscape format when `output_file` is given (`overwrite` to replace it);
  rows in tree order, by OFV or by run name (`sort`). In the task guide
  (step 5, no longer "planned") and on the site.

* `nm_suggest_model()`, sparse designs (no rich individual profile): the
  residual error model cannot be assessed before NONMEM, and the proportional
  model was always ranked first (the additive model never). The three error
  models now get the same support: the best structure is proposed with
  proportional, combined and additive error as candidates 1 to 3, which
  `nm_suggest_run()` fits by default, so that NONMEM decides. A test on the
  residuals within subjects around the pooled fit was tried and left out: it
  found an additive component in half of the additive-error simulations, but
  also in 60 to 90% of proportional-error simulations with a two-compartment
  or saturable structure.

* `nm_gof()` given a run directory: the estimation outputs of NONMEM (`.ext`,
  `.phi`, `.cov`, `.cor`, `.shk`...), which also start with "TABLE NO.", are
  no longer taken for `$TABLE` files; among several tables, the ones holding
  DV are kept, then `sdtab*`. `nm_gof()` failed on the run directories of the
  example project ("Several NONMEM tables found"), and so did the examples of
  the task guide for diagnostics and the Word export.

* `nm_suggest_model()`: the combined residual error model is now retained
  only when its BIC is better by more than 15 (10 before; point E28). The
  threshold was chosen with the simulation validation
  (`inst/validation/validate_suggest_model.R`, 8 scenarios, rich designs):
  correct error model 87% on average over the three true error models, as
  with 10, with the proportional model found more often (86% vs 82%), the
  additive model as often (85%) and the combined model a little less often
  (90% vs 95%); 20 replicates per scenario, one combined-error scenario. In
  sparse designs no BIC is computed and the threshold plays no role.

* New function `nm_check_dataset()` (point B5): checks that a dataset can be
  read by NONMEM (ID, TIME and DV columns; missing ID or TIME; non-numeric
  values; decreasing times within a subject, except at reset records;
  inconsistent EVID, MDV and AMT; duplicated records or observations;
  negative DV and outliers compared with the observations taken at the same
  time after dose). It reads the file as NONMEM does and only reports: the
  data are never modified. Listed in the task guide (stage 1).

* Output files are named `output_file` everywhere (point B11): in
  `nm_write_ctl()`, `nm_ctl_mu()`, `nm_ctl_iov()`, `nm_ctl_m3()`,
  `nm_export_data()`, `nm_export_plot()` and `nm_guide()` (formerly `file`),
  `nm_add_rate()` and `nm_convert_ctl()` (formerly `out_file`), `nm_gof()`
  (formerly `docx_file`), `nm_extract_run_table()` (formerly `output_name`)
  and `nm_vpc()` (formerly `plot_file`). The former names are still
  accepted, with a warning, until version 1.0. `overwrite` was already used
  everywhere.

* Every figure written by the package goes through `nm_export_plot()`
  (point B9): `nm_gof(save_docx = TRUE)` no longer calls `ggsave()`
  directly. Formats: `.jpg`, `.png`, `.pdf` and `.tiff`, chosen from the
  extension.

# nmpilot 0.9.0

* Task guide (`nm_guide()`): every argument of the code examples is
  explained by a comment at the end of its line, and `nm_model_diagram()` is
  listed as available, with an example. New test `test-guide.R`: every
  function and every argument cited by the examples must exist (a former
  argument name, accepted with a warning, is reported), and every example
  must parse.

* Automatic checks on GitHub (point A.14): `R CMD check` runs at each push on
  Windows, macOS and Linux (current R), and on Linux with R-devel and the
  previous R; test coverage is measured with covr and sent to Codecov. Badges
  in the README.

* Documentation website (point A.13): <https://nskoch.github.io/nmpilot/>,
  built with pkgdown and published at each push: the README and the task
  guide on the home page, the functions grouped by stage of the analysis,
  the vignettes and this file. Its address is in `DESCRIPTION`.

* One reader for the records `$INPUT`, `$THETA`, `$OMEGA` and `$SIGMA`
  (point A.11, third part), shared by `nm_ctl_check()`, the transformations,
  `nm_identify_parameters()`, `nm_convert_ctl()` and `nm_profile()`, which
  had each their own (new internal file `R/nm_ctl_read.R`). Results are
  unchanged on the 186 control streams of the comparison corpus, except for
  the following corrections:
  * `$THETA`: all the `$THETA` records are read by the check of the
    transformations (only the first one before); `(0.5 FIX)`, with FIX
    inside the parentheses, and `(0 1 5)`, with blanks as separators, are
    read; `(2.0)X8` gives 8 THETAs; `$THETAP` and `$THETAPV` (priors) are no
    longer counted as THETAs by `nm_profile()`;
  * `$OMEGA`/`$SIGMA`: FIX anywhere in a `BLOCK(n)` record fixes the whole
    block, as in NONMEM (only when FIX was on every line before); the label
    of an ETA of a BLOCK is the comment of the line of its variance (one
    value per line gave shifted labels); a comment-only line is no longer a
    label; `BLOCK(1) SAME` no longer swallows the values of the next record
    (one ETA was lost); new column `omega_comment` (full comment);
  * `$INPUT`: `DROP` written alone is not a data item;
  * `nm_convert_ctl()`: a name assigned after `IF(...)`, assigned as an
    array element (`MTIME(1) = ...`) or assigned in an earlier model record
    (`$PK` before `$DES` or `$ERROR`) is no longer replaced by
    `THETA(n)`/`ETA(n)`; control streams already in standard NONMEM syntax
    are left unchanged.
  Tests in `test-ctl-read.R`.

* Print and plot methods harmonised (point A.11, second part):
  * every `plot()` method returns the figure (a ggplot, or a patchwork for
    the combined panels of `nm_gof`), shown when not assigned and ready for
    `nm_export_plot()`. `plot.nm_gof()` and `plot.nm_model_diagram()`
    printed the figure themselves and returned the object;
  * the argument that chooses the figure is `type` in every method:
    `plot.nm_compare()` and `plot.nm_gof()` (formerly `which`),
    `plot.nm_optimize_base()` (formerly `area`); the former names still work
    until 1.0, with a warning;
  * the run log of `nm_run()` and `nm_vpc_ctl()` has the class `nm_run_log`,
    the result of `nm_extract_run_table()` the class `nm_run_table`, and the
    result of `nm_data_previsual()` the class `nm_data_previsual`, each with
    a print method giving a short summary. They remain a data frame or a
    list: existing code that uses their columns or elements is unchanged.

* Argument names harmonised before version 1.0 (point A.11, first part). The
  former names and forms still work until 1.0, with a warning that gives
  the new one:
  * the control stream is `ctl` everywhere: `nm_run(ctl, nmfe_path, ...)`
    (formerly `ctl_files`, after `nmfe_path`), `nm_profile(ctl, nmfe_path,
    ...)` (formerly `ctl_path`, after `nmfe_path`), `nm_identify_parameters(ctl)`
    (formerly `file`), `nm_convert_ctl(ctl)` (formerly `ctl_file`). The
    former order, with `nmfe_path` first, is recognised when the arguments
    are given by position;
  * a finished run is given by its path, `run` (control stream, run
    directory or stem): `nm_vpc(run, ...)` and `nm_vpc_ctl(run, ...)`
    (formerly `run_name` + `analysis_dir`), `nm_extract_run_table(run)`
    (formerly `run_rep` + `run_name`, or `run_path`);
  * the output directory is `output_dir`: `nm_bootstrap()` (formerly
    `run_dir`), `nm_derivatives_ctl()` and `nm_data_previsual()` (formerly
    `out_dir`);
  * `nm_profile()` takes `verbose` (formerly `quiet`, the opposite), and
    `nm_bootstrap()` `strip_tables` (formerly `strip_table`), as the other
    functions.
  The README, the vignettes, the guide and the tests use the new names.
  Tests in `test-deprecated-args.R`.

* New function `nm_model_diagram()` (point A.7): tree of the models tested in
  a project. Each run is linked to the run it was derived from, read in its
  control stream as a comment with the PsN and Pirana convention
  (`;; 1. Based on: run3`, or `;; Based on: 3` for `run3`) or given with
  `parents =`, which also corrects or removes links. Each run shows its OFV,
  its dOFV and its change in the number of estimated parameters relative to
  its parent, coloured by convergence state (covariance done, no
  covariance, not converged, not completed, simulation). Returns the table
  of the tree and a ggplot2 figure, which `nm_export_plot()` saves. Circular
  links and parents absent from the project are reported. Tests in
  `test-model-diagram.R`.

Reading of real outputs of NONMEM 7.1 to 7.6 (point A.3): a sample of 43 runs
(7.1, 7.2, 7.3, 7.4.1, 7.6, and a few NONMEM V/VI outputs; models with and
without `$DES`, ADVAN1 to ADVAN13) was read with `nm_extract_runs()`; the
following problems were found and fixed. NONMEM 7.5 was not checked (no
output of this version was available).

* `nm_extract_runs()` failed for a whole directory ("subscript out of
  bounds") when a run without `.ext` had 10 ETAs or more: from the 10th
  element on, NONMEM prints the rows of the matrices as `ET10`, `EP10`, which
  were not recognised. These labels are now read, and a correlation whose
  variance is missing is reported as `NA` instead of stopping
  `nm_extract_runs()` and `nm_extract_run_table()`.
* A simulation-only run of NONMEM 7.1, which prints an empty `#METH:` line,
  was reported as interrupted; it is now reported as a simulation.
* `nm_ctl_check()` reads the use of `THETA(n)`, `ETA(n)`, `EPS(n)`, `CMT`,
  `EVID`, `MDV` and `RATE` in the code records only: an initial estimate such
  as `$THETA (18.3 FIXED)` was taken for a use of `THETA(18)` ("used but not
  declared"), and a title such as `$PROBLEM ... 1 CMT` for a use of `CMT`.
  With PREDPP (`$SUBROUTINES`), `EVID` and `MDV` used in the code without
  being in `$INPUT` are no longer reported, as NM-TRAN supplies these items.
  `RATE` with `KA` is no longer reported when a zero-order input into any
  compartment (`D2`, `R3`...) is modelled.
* NONMEM V/VI outputs: the tables printed in the `.lst` were read as
  parameter estimates, and the OFV, printed about ten lines below its banner,
  was not found. Both are now read.
* `nm_extract_runs()` reads the number of observations (`n_obs`) in the
  `.lst` ("TOT. NO. OF OBS RECS") when there is no `.smr`, which is the case
  of the runs launched with `nm_run()` or nmfe without Wings for NONMEM; the
  BIC is then computed as well. On the 43 runs of the sample, the value read
  in the `.lst` equals the one of the `.smr` whenever both exist.
* `nm_extract_run_table()` no longer requires a `.ext` file: when it is
  missing or holds no estimates (empty file), the estimates printed in the
  `.lst` are used, as in `nm_extract_runs()`, and the table footer says so.
  This covers the runs of older NONMEM versions or of PsN without `.ext`
  (16 of the 43 runs of the sample have none). The `.lst` prints three
  significant digits: on three runs read both ways, the estimates differ by
  less than 0.4 % and the OFV is the same.
* A model written with `$DES` (the example run2 rewritten with ADVAN13) was
  run with NONMEM 7.6 through the whole chain: estimation (same OFV and
  estimates as run2), `nm_extract_runs()`, `nm_extract_run_table()`,
  `nm_vpc_ctl()` and `nm_vpc()`, `nm_profile()`, `nm_bootstrap()` and the
  derivatives run of `nm_derivatives_ctl()`.
* Control streams written by `nm_write_ctl()`, `nm_ctl_mu()`, `nm_ctl_iov()`
  and `nm_ctl_m3()` were run through NONMEM 7.6 (75 runs: every structure,
  error model, estimation method and transformation; action 26). Six
  errors were found and fixed:
  - `$DES` models (`ADVAN13`, `ADVAN6`) are written with `TOL=9`, required by
    NM-TRAN (error 145);
  - FO is written `METHOD=ZERO` instead of `METHOD=FO` (error 287);
  - `nm_ctl_m3()` adds `M3LLOQ`, `M3SD`, `M3PBQL` before the options of
    `$TABLE`, not after `FILE=` (error 20);
  - `nm_ctl_m3()` writes both definitions of `Y` in one `IF/THEN/ELSE`
    block, the M3 quantities being computed before it. The former code
    defined `Y` with the residual error in a block written
    `IF (BQL.EQ.0.OR.ICALL.EQ.4) THEN`, which NM-TRAN 7.6 refused (error 122,
    "EPS VARIABLES REQUIRED WITH POPULATION DATA", and warning 53): a block
    testing `ICALL.EQ.4` is simulation code for NM-TRAN. The censoring flag
    `M3CENS` keeps the test of `ICALL` out of the blocks defining `Y`. The
    four M3 control streams then ran with NONMEM 7.6 (estimation with
    LAPLACE converged, covariance step completed);
  - `nm_ctl_mu(mode = "linear")` writes 0.0001 instead of log(1) = 0 for a
    THETA whose initial value is 1 (error 260: initial estimate of 0);
  - three-compartment models with first-order absorption (ADVAN12): the
    default initial values of the second peripheral compartment are now
    Q4 = 2 and V4 = 400 (they were equal to Q3 and V3, which PREDPP refuses:
    "K32, OR K42 IS TOO CLOSE TO AN EIGENVALUE").
* `nm_ctl_check()` reports these NM-TRAN errors: `TOL` missing with a
  general ADVAN, unknown `METHOD`, `THETA` initial estimate of 0 not fixed,
  `$TABLE` item written after the options, `Y` with `EPS`/`ERR` defined only
  in `IF` blocks testing `ICALL.EQ.4`. No false alarm on the 38 real
  control streams of the sample nor on the example project.
* New tests for these cases and for models written with `$DES`
  (`test-lst-versions.R`, `test-des.R`, `test-nmtran-errors.R`), built from
  the example run or written by hand.

* Tests on control streams and outputs not written by nmpilot (point A.4):
  seven runs of the package bbr (Metrum Research Group, MIT licence; NONMEM
  7.4.4) are copied in `tests/testthat/fixtures/public/bbr`, with their
  licence and source. They cover `$DES`, `$PRIOR`, `MTIME`, BAYES, a
  FOCEI > BAYES > IMP chain, `$SIMULATION ONLYSIM`, `$ABBR REPLACE` with 62
  ETAs and `$MIX`. The new tests check `nm_ctl_check()`,
  `nm_identify_parameters()`, the reading of the `.lst`, `nm_extract_runs()`
  and `nm_convert_ctl()` on these files. They showed six problems, fixed:
  * `nm_ctl_check()` reported as an error `$TABLE` items written on a line
    below the options, which NM-TRAN accepts; only an item written after an
    option of the same line is reported now.
  * `$ABBR REPLACE` with a range (`ETA(OCC_CL)=ETA(,3:32)`, also
    `ETA(3 TO 32 BY 2)`): the ETAs were reported as unused by
    `nm_ctl_check()` and not linked to their parameter by
    `nm_identify_parameters()`. The names given by `$ABBR REPLACE` (ranges
    and lists) are now replaced by their elements, so that IOV ETAs are
    linked to their parameter. `nm_identify_parameters()` also counts the
    ETAs of `$OMEGA BLOCK(n) SAME(k)` (k blocks, not one).
  * Shrinkage read in the `.lst`: beyond 10 ETAs, the values continued on the
    following lines were ignored (20 values read out of 62).
  * Without `.ext`, the OMEGA matrix of the `.lst` was misread beyond 12
    ETAs (rows wrapped by NONMEM): 12 variances read out of 62.
  * `nm_identify_parameters()`: a covariate effect written in a helper
    variable added on the log scale (`CLEGFR = LOG(EGFR/90)*THETA(6)`,
    `CL = EXP(THETA(3) + CLEGFR + ETA(3))`) was taken for a typical value.
  * `nm_identify_parameters()`: a THETA present in every branch of a
    conditional definition (`IF(PATIENT.EQ.0) KOUT = EXP(THETA(1) + ETA(1))`,
    `IF(PATIENT.EQ.1) KOUT = EXP(THETA(1) + ETA(2))`) was taken for a
    covariate effect; it is the typical value. The THETAs of `MTIME(n)` and
    of the mixture probabilities `P(n)` of `$MIX` are returned in the new
    element `special_theta` instead of `unmapped_theta`.

# nmpilot 0.8.1

* Linux and macOS: `nm_bootstrap()` no longer depends on `cmd.exe`. Its
  sequential loop (`n_cores = 1`) now uses the launcher shared by `nm_run()`,
  `nm_profile()` and `nm_jackknife()`, which calls `cmd.exe` on Windows and
  the nmfe script directly on Linux and macOS. On these systems, a launcher
  that is not executable is reported with a clear message ("chmod +x")
  instead of a failure of the system call. The documentation of `nmfe_path`
  gives a Linux example. Checked on Linux with a simulated nmfe for
  `nm_run()`, `nm_profile()`, `nm_jackknife()`, `nm_bootstrap()`,
  `nm_vpc_ctl()` and `nm_vpc()`, one run at a time and in parallel; not yet
  with a real NONMEM installation on Linux or macOS.

# nmpilot 0.8.0

* The package is renamed **nmpilot** (formerly **nmtools**), to avoid
  confusion with an unrelated GitHub repository of the same name in the same
  field. Install it with `remotes::install_github("nskoch/nmpilot")` and load
  it with `library(nmpilot)`; remove the old package with
  `remove.packages("nmtools")`. No function is renamed.
* The options are renamed `nmpilot.n_cores` and `nmpilot.warn_synced`
  (formerly `nmtools.n_cores` and `nmtools.warn_synced`).
* The files written by parallel runs are renamed `<run>_nmpilot.bat` (or
  `.sh`), `<run>_nmpilot.log`, `<run>_nmpilot.err` and `<run>.nmpilot_done`.
  Files left by earlier versions (`*_nmtools.*`, `*.nmtools_done`) are not
  used and can be deleted.
* `nm_vpc_ctl()` wrote `FIX` on every row of the `$OMEGA` and `$SIGMA`
  blocks of the `_orig` and `_sim` control streams, which NM-TRAN refuses
  (error 52, "THIS OPTION HAS ALREADY BEEN SPECIFIED"): `FIX` is now written
  once, on the `BLOCK(n)` line. Found with NONMEM 7.6 on the example project.
* `nm_vpc()`: when a time bin between two others was empty, the bins after
  it were drawn at the wrong time (the midpoint of a bin was looked up by
  position instead of by bin number) and the last bin was dropped, with the
  warning "Removed 1 row containing missing values (geom_line)". On the
  example project with automatic bins, the statistics of the 12 h bin were
  drawn at 24 h and the 24 h bin was missing. Fixed.
* `nm_vpc()`: the first line of the `.fit` tables was always skipped, as if
  it were a header, whereas `nm_vpc_ctl()` writes them without header
  (`NOHEADER`): when a data set starts with an observation, that
  observation was lost. Header lines ("TABLE NO.", column names) are now
  recognised and left out, and no record is skipped.
* `nm_run()` (and the functions that launch runs through it, such as
  `nm_vpc_ctl()`): a simulation-only run (`$SIMULATION ONLYSIM`, e.g. the
  `_sim` run of a VPC) that finished is now reported with
  `minimization_status = "simulation"` and `run_ok = TRUE`, and is skipped
  as done by `skip_terminated = TRUE`. It was reported with `NA` and
  relaunched at each call.
* `nm_extract_runs()`: a simulation-only run (`$SIMULATION ONLYSIM`) that
  finished is now reported with `run_state = "simulation"`,
  `min_status = "SIMULATION"` and `min_ok = NA`; it was reported as
  `"interrupted"` with `min_ok = FALSE`. The subproblems of a run
  (`SUBPROBLEMS=n`) are no longer counted as problems in `n_problems` (the
  `_sim` run of a 500-replicate VPC had `n_problems = 501`).
* The entries below, up to 0.7.1, were released under the name nmtools.

# nmpilot 0.7.1

* No change in the code. This release is the first to contain the work of
  versions 0.6.0 (parallel NONMEM runs) and 0.7.0 (reading NONMEM outputs),
  described below. By mistake, the GitHub releases v0.6.0 and v0.7.0 (and
  their Zenodo archives) were made on the commit of version 0.5.0: they hold
  the code of 0.5.0. Use 0.7.1 or later.

# nmpilot 0.7.0

## Reading NONMEM outputs: several $EST and $PROB, IOV, $MIX, failed runs

* With several `$EST` records (for example ITS > SAEM > IMP > BAYES > FOCEI),
  the OFV, the minimisation status, the shrinkage and the covariance step are
  now those of the last estimation step; they came from the first steps.
  New column `est_methods` in `nm_extract_runs()` (e.g.
  `"ITS > SAEM > IMP(eval) > BAYES > FOCEI"`) and `n_est`. EM and Bayesian
  status lines ("STOCHASTIC PORTION WAS COMPLETED", "EXPECTATION ONLY PROCESS
  COMPLETED", "OPTIMIZATION WAS NOT TESTED") are recognised.
* The same rule is applied in `nm_run()`, `nm_profile()`, `nm_bootstrap()`,
  `nm_jackknife()` (OFV, status and THETA of the last step) and `nm_compare()`
  (OFV).
* `nm_run()`: the covariance status is read from the status lines only; the
  sentence "AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE
  OUTPUT", printed after "HOWEVER, PROBLEMS OCCURRED", made runs whose
  covariance step had succeeded appear as `"min_success_cov_failed"`.
* With several `$PROB`, the last problem that has an estimation is read (a
  simulation problem after the estimation is skipped); new column
  `n_problems`.
* `nm_extract_runs()`: when the `.ext` is empty, has no final row, or gives an
  OFV different from the `.lst` (an `.ext` left by another run with the same
  name), the estimates and standard errors are read in the `.lst`, with a
  warning and "EXT/LST MISMATCH" in `warnings`; new column `params_source`
  (`"ext"`, `"lst"` or `"smr"`).
* New column `run_state`: `"completed"`, `"interrupted"` (the output stops
  before the end of the last step), `"terminated"` (PROGRAM TERMINATED) or
  `"nmtran_error"`. A run that is not completed has `min_ok = FALSE`.
* IOV: ETAs declared by `$OMEGA ... SAME` (also `BLOC(n)`, `SAME(k)`) are
  marked in `$params` (new column `same`) and counted in `n_eta_same`.
* `n_omega` and `n_sigma` now count the estimated elements only (not FIX, not
  SAME, not the zeros outside the blocks): on the example project, 3 and 0
  instead of 6 and 1. ETAs fixed to 0 are left out of `max_eta_shrink`.
* `$MIX`: new columns `mix`, `n_subpop` and `mix_assign` (subjects per
  subpopulation, from the `.phm`); new element `$mixture` (probability of each
  subpopulation for each subject); mixing-proportion THETAs without another
  name are named `MIX_P`; the shrinkage of an ETA is its smallest value over
  the subpopulations (NONMEM prints 100 % in the subpopulations that do not
  use it). `nm_filter_runs()` filters `$mixture` too.
* Run files are matched without regard to case (NONMEM 7.1 wrote `RUN.EXT`
  next to `run.lst`), and `.mod` control streams are recognised.
* Checked on 11 real runs (NONMEM 7.1 to 7.4): chains of five `$EST`, SAEM >
  IMP, ITS > FOCE, IOV, `$MIX` with `METHOD=1` and `CONDITIONAL`, empty
  `.ext`, `.ext` from another run.

# nmpilot 0.6.0

## Parallel NONMEM runs

* New argument `n_cores` in `nm_run()`, `nm_profile()`, `nm_jackknife()`,
  `nm_bootstrap()`, `nm_vpc_ctl()` and `nm_suggest_run()` (default 1: runs
  one after the other, as before). With `n_cores > 1`, each run is started in
  the background from its own directory (NM-TRAN writes its work files in the
  current directory) and at most `n_cores` runs go on at the same time.
  `n_cores = "auto"` uses the physical cores minus one; `n_cores = "ask"`
  shows the cores and asks the number once per R session (the answer is kept
  in `options(nmtools.n_cores)`; set it to `NULL` to ask again).
* New `nm_cores()`: logical and physical cores of the computer, and the
  number of runs `"auto"` would use.
* Each parallel run writes the console output of nmfe to `<run>_nmtools.log`
  (and `<run>_nmtools.err`). If R is interrupted, the runs already started go
  on in the background and are listed.
* `nm_run()`: with `n_cores > 1`, several control streams in a flat layout
  (shared folder) are run one after the other, with a warning. `nm_jackknife()` now always runs
  each jackknife data set in its own subdirectory `jk_run/<name>/`.
  `nm_bootstrap()` uses working subdirectories `slot1`, `slot2`, ...
* Checked with NONMEM 7.6 on Windows with 5 runs at a time: `nm_run()` (4
  runs), `nm_profile()` (12 grid values), `nm_jackknife()` (60 runs) and
  `nm_bootstrap()` (20 replicates), all successful. On Windows
  the launcher calls nmfe without any redirection inside the batch file: a
  `2>&1` there keeps the first `2>trash.tmp` of nmfe in force and the
  compilation fails.

## Display fixes

* `nm_profile()`: the start message gives the initial value of the profiled
  THETA and, when the `.ext` of the control stream exists, its final estimate
  (it showed only the initial value, as "current value").
* `print.nm_jackknife()`: the reference OFV is printed with three decimals
  (it was rounded to four significant digits, e.g. -1224).

# nmpilot 0.5.0

## Direct comparison of two models (`nm_compare()`)

* Vuong (1989) and Clarke (2007) tests on the individual contributions to
  the objective function (iOFV of the `.phi` files), for non-nested models
  fitted to the same data, without correction and with AIC-like and
  BIC-like corrections for the number of parameters (read from the `.ext`
  files, or given with `n_par`). Likelihood ratio test when `nested = TRUE`.
* Subject bootstrap of the delta iOFV (no re-estimation): confidence
  interval of dOFV and probability that run 2 has a lower OFV (and AIC).
  Influence: number of subjects, largest contributions first, whose removal
  reverses the sign of dOFV.
* Delta iOFV versus covariates (`covariates =`, Spearman or Kruskal-Wallis
  with Bonferroni adjustment); parameters of the two runs matched by name
  (estimate, SE, RSE, shrinkage, ratio) and paired EBEs.
* New plots: `"gain_map"` (|CWRES_2| - |CWRES_1| versus TIME and PRED,
  loess with a 95% subject-bootstrap band), `"overlay_gof"` (diagnostics of
  both runs superimposed), `"bootstrap"`, `"iofv_covariate"`,
  `"parameters"`, `"eta_scatter"`.
* New `summary()` method: decision criteria (OFV, AIC, BIC, tests,
  bootstrap, influence, residual metrics, maximum RSE, condition number from
  the `.cor` file, maximum shrinkage, run status) with the run favoured by
  each.

## Improving a base model without iterative runs

* New `nm_derivatives_ctl()`: writes an evaluation run of a base model
  (final estimates fixed, `MAXEVAL=0`) that tabulates the derivatives of Y
  with respect to the ETAs, the EPSs and their cross derivatives (verbatim
  code at the end of `$ERROR`/`$PRED`, as PsN `linearize`). Checked with
  NONMEM 7.6.
* New `nm_optimize_base()`: level 1 from the outputs of the base run
  (resmod-like models on CWRES, EBE shrinkage, skewness and correlations,
  ETA-covariate screening); level 2 on the linearized model built from the
  derivatives table, each candidate re-estimated in R: residual error models,
  OMEGA blocks, removal of an ETA, inter-occasion variability, covariates on
  the ETAs (univariate tests, forward inclusion and backward elimination),
  then the combined proposal. The default objective (`objective = "focei"`)
  is NONMEM's FOCE with interaction (modes of the full individual
  objective, residual-variance term in the determinant) and reproduces the
  OFV of the base run at its estimates; `objective = "fast"` is a quicker
  approximation. `plot()` shows the dOFV of each candidate as better, no
  significant change or worse.

## Other changes

* `nm_run()` reads the OFV and status of an evaluation run (`MAXEVAL=0`,
  status `"evaluation"`), and warns when the analysis directory is inside a
  folder synchronised by a cloud client (Dropbox, OneDrive, Google Drive,
  iCloud, Box), whose client can lock files during compilation (silence with
  `options(nmtools.warn_synced = FALSE)`).

# nmpilot 0.4.1

* `nm_suggest_model()`: the residual error model is chosen on a common
  structural fit of each interval (iteratively reweighted least squares
  under a combined variance); the additive, proportional and combined
  variance models are then compared by BIC on the same leverage-corrected
  (PRESS) residuals, so that they are nested. Points of leverage above 0.7
  and predictions below twice the LLOQ are left out, and the combined model
  is kept only when its BIC is better by more than 10. On simulated rich
  data (20 replicates per scenario), the error model is recovered in 60 to
  100% of the datasets; better than 0.4.0 for Michaelis-Menten, additive and
  combined error, lower for some proportional-error scenarios
  (`inst/validation/validate_suggest_model.R`).

# nmpilot 0.4.0

## Data-driven suggestion of structural models

* New `nm_suggest_model()`: analyses a NONMEM dataset, rich or sparse,
  without running NONMEM and ranks candidate models: number of compartments
  (sums of exponentials compared by BIC), input (IV bolus or infusion,
  first-order, first-order with lag time, zero-order duration, transit
  compartments; full-profile fits with a censoring penalty for early BLQ
  samples), elimination (dose proportionality and curvature of the
  log-concentration decline) and residual error model (additive,
  proportional or combined; extended least squares on leverage-corrected
  residuals). Each candidate comes with initial values for THETA, IIV
  (variance of individual estimates) and SIGMA, the parameters likely to be
  hard to estimate and their reasons, and alerts on features outside the
  tested structures (double peaks, clearance changing between first dose
  and steady state, unexpected accumulation, less than dose-proportional
  exposure). Sparse data are analysed as naive pooled data. `print()` and
  `plot()` methods (`"profiles"`, `"support"`, `"pooled"`, `"design"`).
* New `nm_suggest_design()`: expected RSE and eta-shrinkage of a candidate
  for the sampling design of the dataset or another one (population Fisher
  information matrix, first-order approximation, analytic solution with
  superposition of doses).
* New `nm_suggest_run()`: writes the control streams of the best
  candidates, runs them with `nm_run()` and compares the BIC ranking with
  the data-driven score.
* New `nm_suggest_covariates()`: Wasserstein-distance screening of
  covariates on the individual estimates (CL, central volume, absorption)
  before any modelling.
* `nm_suggest_run()` gives each candidate a dataset that matches its
  compartment numbering (depot 1 and central 2 for first-order and transit
  absorption; central 1 for IV and zero-order input, with `RATE = -2` on the
  doses for a modelled duration). A copy is written only when the original
  dataset does not match; the original file is never modified. With
  `overwrite = TRUE`, the candidates are run again even when previous runs
  of the same name have terminated.
* `nm_write_ctl()`: new arguments `from` and `candidate` to write the
  control stream of a candidate of `nm_suggest_model()`; explicit arguments
  keep priority. New argument `conc_factor`: the scale is written
  `S = V/conc_factor` (e.g. 1000 for doses in mg, volumes in L and DV in
  ng/mL); with `from`, the `dose_conc_factor` of the suggestion is used.
* New validation script `inst/validation/validate_suggest_model.R`
  (recognition rate of the true structure on simulated datasets).

# nmpilot 0.3.0

## Writing control streams

Functions ported from the Shiny application "Générateur NONMEM .ctl"
(<https://nicolas-simon.shinyapps.io/nonmem_ctl_publi/>). All return an
`nm_ctl` object (control stream lines, printed as text) and can write a file.

* New `nm_write_ctl()`: writes a control stream for a 1- to 3-compartment PK
  model, IV or extravascular, with first-order, zero-order (modelled rate
  `R1 = K0`, `RATE = -1`, or modelled duration `D1`, `RATE = -2`) or
  transit-compartment absorption, linear or Michaelis-Menten elimination,
  analytic (`ADVAN`/`TRANS`) or ODE solution. `$INPUT` follows the column
  order of the dataset; THETA bounds and FIX, IIV as CV% or variance,
  diagonal or block `$OMEGA`, additive, proportional or combined error,
  nine estimation methods, xpose-style `sdtab`/`patab`/`cotab` tables.
* New `nm_ctl_mu()`: MU referencing of log-normal parameters (`"readable"`
  or `"linear"` mode).
* New `nm_ctl_iov()`: inter-occasion variability with `$ABBR REPLACE` and
  `$OMEGA BLOCK(p) SAME(n-1)`; optional check of the occasion column.
* New `nm_ctl_m3()`: M3 method for BQL data, with the LLOQ given as a
  constant, a dataset column, or by analytical method code.
* New `nm_ctl_check()`: static check of a control stream (THETA/ETA/EPS used
  but not declared, over all records including `BLOCK`, `SAME` and
  `$ABBR REPLACE`; `$INPUT` items; M3 requirements), with optional line
  differences from a reference control stream.

# nmpilot 0.2.3

* New `nm_jackknife()`: leave-one-subject-out jackknife of a NONMEM model
  (method of the standalone script nonmem-jackknife-influence,
  doi:10.5281/zenodo.22257994). Ranks the subjects by their likelihood
  influence |dOFV| and reports the relative THETA shifts (L2 norm, maximum,
  SE-standardised norm); `print()` and `plot()` methods (`"bar"`, `"rank"`,
  `"scatter"`). Resumes an interrupted analysis. Data lines are copied
  verbatim, as in `nm_bootstrap()`.
* `nm_identify_parameters()`: new column `derived_from` in
  `parameter_table`, naming the parameters a derived parameter is computed
  from (`"V"` for `S2 = V`, `"CL, V"` for `KE = CL/V`; `NA` otherwise).
* `nm_extract_runs()`: covariate-effect and residual-error THETAs are named
  from their `$THETA` comment label when it is a short name (`CL_WT`, `PROP`,
  `RUV_A`, `4 Prop (...)` -> `Prop`) instead of `THETA4`, `THETA5`...

# nmpilot 0.2.2

## Breaking changes

* Output column names are now in English. `nm_extract_run_table()`,
  `params` table: `Parameter`, `Estimate`, `RSE_pct`, `CI95`, `CI95_lower`,
  `CI95_upper` (were `Parametre`, `Estimation`, `RSE_pourcent`, `IC95`,
  `IC95_lower`, `IC95_upper`). `nm_summary()` and `nm_data_previsual()`:
  `unique_doses` (was `doses_uniques`).

## Wings for NONMEM is optional

* `nm_run()`: new argument `run_subdir`. `".reg"` (default) keeps the
  Wings-compatible layout `run1.reg/`; any other suffix can be used; `""`
  runs NONMEM in the folder of the control stream, without copying it.
  A control stream in Wings label syntax cannot be run with `run_subdir = ""`
  (clear error message).
* `nm_run()`: a relative `$DATA` path is looked up from the folder of the
  control stream (`$DATA data.csv`), then from the run directory
  (`$DATA ../data.csv`, Wings style).
* `nm_vpc_ctl()`, `nm_vpc()` and `nm_extract_run_table()` find a run in
  `<run>.reg/`, `<run>.nm7/`, a folder named `<run>`, or directly in the
  analysis directory. `nm_vpc_ctl()` gains `run_subdir`; by default the VPC
  runs use the layout of the estimation run.

## Bug fixes

* `nm_run()`: NONMEM is now launched from the run directory (NM-TRAN failed
  when started from another working directory).
* `nm_run()`: standard control streams with comments after `$THETA` or
  `$OMEGA` values are no longer taken for Wings control streams.
* `nm_run()`: `verbose = FALSE` now silences all messages.

## Documentation

* Example project in `inst/extdata/example`: simulated data (60 fictitious
  subjects, one-compartment oral model) and the NONMEM outputs of a BASE and
  a FINAL model.
* Two vignettes: `vignette("nmtools")` (getting started) and
  `vignette("wasserstein")`.
* `citation("nmtools")` cites the package and, for `nm_wasserstein()`, the
  original article.

# nmpilot 0.2.1

* `nm_wasserstein()`: new argument `continuous = "pairwise"` to compare every
  pair of quartiles of a continuous covariate (6 comparisons), as in the
  pairwise analysis of the original article; Bonferroni-adjusted p-values
  `p_W_adj` and `p_t_adj` are added to the results. The default
  (`"extremes"`, Q4 vs Q1) is unchanged.
* `nm_identify_parameters()`: a THETA assigned under an IF condition whose
  other branch is a plain constant (e.g. `IF (CT01.EQ.1) FALL = THETA(7)`,
  `ELSE FALL = 1`) is now a covariate effect; a covariate term added inside a
  sum written with extra parentheses (e.g.
  `((THETA(1)*FALL) + (THETA(8)*(CLCR/89)))`) is now a covariate effect.

# nmpilot 0.2.0

* New `nm_wasserstein()`: distribution-based covariate assessment with the
  Wasserstein distance (Simon, Hulot and von Fabeck, Front Pharmacol 2026,
  doi:10.3389/fphar.2026.1804989). For one or several runs, compares the
  distribution of an ETA between groups of subjects defined by covariates
  (permutation test on the 1-Wasserstein distance, with a t-test for
  reference) and reports a Wasserstein-based R2 with the eta-shrinkage.
  `print()` and `plot()` methods (ECDF by group, R2W bar chart).

# nmpilot 0.1.0

First public release.

* Data sets: `nm_read_csv()`, `nm_data_previsual()`, `nm_plot_dv()`,
  `nm_summary()`, `nm_add_rate()`, `nm_export_data()`.
* Running NONMEM: `nm_run()`, `nm_convert_ctl()` (Wings extended control
  streams).
* Results: `nm_extract_runs()`, `nm_filter_runs()`, `nm_compare_runs()`,
  `nm_extract_run_table()`.
* Parameter identification from the model code: `nm_identify_parameters()`,
  `nm_identify_directory()`.
* Diagnostics: `nm_gof()`, `nm_compare()`, `nm_vpc_ctl()`, `nm_vpc()`.
* Uncertainty: `nm_profile()`, `nm_bootstrap()`.
* Figures: `nm_export_plot()`.
