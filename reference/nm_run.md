# Run one or more NONMEM models

Sequentially launches a list of NONMEM control streams via the `nmfe`
executable (direct call, no Wings for NONMEM, no PsN required). By
default each run gets its own subdirectory `<run>.reg` inside
`analysis_dir`, holding a copy of the control stream and all NONMEM
output files; with `run_subdir = ""` the run is made directly in the
folder of its control stream, as when nmfe is called by hand. The data
file is not copied. A relative `$DATA` path is accepted: it is looked up
from the folder of the control stream (`$DATA data.csv`), then from the
run directory (`$DATA ../data.csv`, Wings style), and written as an
absolute path in the copied control stream.

## Usage

``` r
nm_run(
  ctl,
  nmfe_path,
  analysis_dir,
  skip_terminated = TRUE,
  accept_failed = TRUE,
  overwrite_ctl = TRUE,
  run_subdir = ".reg",
  verbose = TRUE,
  n_cores = 1L,
  ...
)
```

## Arguments

- ctl:

  Character vector of paths to `.ctl` files. May be absolute paths or
  file names relative to `analysis_dir`. (Formerly `ctl_files`, still
  accepted with a warning; the former order, `nmfe_path` first, is also
  recognised.)

- nmfe_path:

  Full path to the nmfe executable (e.g. `"C:/nm76g64/run/nmfe76.bat"`
  on Windows, `"/opt/nm760/run/nmfe76"` on Linux or macOS, where it is
  run directly and must be executable).

- analysis_dir:

  Root directory of the analysis. One run subdirectory is created here
  for each run (see `run_subdir`). A warning is issued when it lies in a
  folder synchronised by a cloud client (Dropbox, OneDrive, Google
  Drive, iCloud, Box): the client can lock the files written by NM-TRAN
  and make a run fail (silence it with
  `options(nmpilot.warn_synced = FALSE)`).

- skip_terminated:

  Logical. If `TRUE` (default), runs whose `.lst` file already exists
  and reports a terminal state are skipped.

- accept_failed:

  Logical. Only used when `skip_terminated = TRUE`. If `TRUE` (default),
  runs that terminated with an error are also skipped (not re-launched).
  If `FALSE`, only runs with `run_ok = TRUE` are skipped.

- overwrite_ctl:

  Logical. If `TRUE` (default), the `.ctl` file is overwritten in the
  `.reg` directory if already present.

- run_subdir:

  Suffix of the run subdirectories: `".reg"` (default, Wings-compatible
  layout `run001.reg/`), any other suffix (e.g. `".run"`), or `""` to
  run in the folder of the control stream itself, without copying it
  (outputs next to the `.ctl`). Control streams in Wings label syntax
  need a subdirectory, because the converted copy is written there.

- verbose:

  Logical. If `TRUE` (default), progress messages are printed to the
  console.

- n_cores:

  Number of runs launched at the same time (default 1: one after the
  other, NONMEM output shown in the console). A whole number, `"auto"`
  (physical cores minus one, see
  [`nm_cores()`](https://nskoch.github.io/nmpilot/reference/nm_cores.md))
  or `"ask"` (interactive sessions: the cores are shown and the number
  is asked once per R session, then reused;
  `options(nmpilot.n_cores = NULL)` asks again). In parallel, each run
  is started in the background from its run directory; the console
  output of nmfe goes to `<run>_nmpilot.log`. Runs written in the same
  folder (`run_subdir = ""`) are always run one after the other. If R is
  interrupted, the runs already started go on in the background.

- ...:

  Former argument names, accepted with a warning until version 1.0:
  `ctl_files` (now `ctl`).

## Value

A run log of class `nm_run_log` (a `data.frame`, printed as a short
summary), with one row per run and the following columns:

- run:

  Run name (without extension).

- reg_dir:

  Full path to the run directory (`<run>.reg`, or the folder of the
  control stream when `run_subdir = ""`).

- action:

  One of `"launched"`, `"skipped_done"`, or `"skipped_no_ctl"`.

- exit_status:

  Shell return code from `nmfe` (0 = success).

- lst_exists:

  Logical: does the `.lst` file exist after the run?

- ofv:

  Objective function value (OFV).

- minimization_status:

  Minimization status string: `"successful"`, `"terminated"`,
  `"rounding_errors"`, `"min_success_cov_failed"`, `"evaluation"` (run
  without minimization, `MAXEVAL=0`; `min_ok` is then `NA` and `run_ok`
  `TRUE`), `"simulation"` (simulation-only run, `$SIMULATION ONLYSIM`,
  finished; `ofv` and `min_ok` are then `NA` and `run_ok` `TRUE`), or
  `"unknown"`.

- min_ok:

  Logical: did minimization succeed?

- cov_ok:

  Logical: did the covariance step succeed?

- run_ok:

  Logical: is the run globally acceptable (minimization and covariance
  both successful)?

- cov_status:

  Detailed covariance step status: `"successful"`, `"aborted"`,
  `"not_requested"`, or `"unknown"`.

## Details

Wings for NONMEM is not needed. The default `.reg` layout mirrors the
Wings output structure (a run named `run001.ctl` produces a directory
`run001.reg/` containing `run001.ctl`, `run001.lst`, `run001.ext`,
etc.), so that projects started with Wings can be continued; the other
functions of the package read runs in any of these layouts.

Run status is parsed from the `.lst` file after completion. The
`skip_terminated` mechanism reads any existing `.lst` before launching,
allowing safe re-use of `nm_run()` on a partially completed batch.

## See also

[`nm_extract_runs()`](https://nskoch.github.io/nmpilot/reference/nm_extract_runs.md)
for parsing NONMEM output files after a run.

## Examples

``` r
if (FALSE) { # \dontrun{
log <- nm_run(
  ctl             = c("run001.ctl", "run002.ctl", "run003.ctl"),
  nmfe_path       = "C:/nm76g64/run/nmfe76.bat",
  analysis_dir    = "C:/analyses/my_pk_project",
  skip_terminated = TRUE,
  accept_failed   = TRUE,
  verbose         = TRUE
)

# Show only failed runs
log[!log$run_ok, ]
} # }
```
