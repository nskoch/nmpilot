# Profile the objective function over one model parameter

Runs a series of NONMEM estimations in which a single `$THETA` is fixed
to each value of a grid while every other parameter is re-estimated,
yielding the profile of the objective function value (OFV, *without* the
additive constant) for that parameter. The minimum of the curve locates
the conditional estimate; the width at `min(OFV) + ci_threshold` gives
an approximate likelihood-based confidence interval.

## Usage

``` r
nm_profile(
  ctl,
  nmfe_path,
  param,
  from,
  to,
  by,
  output_dir = NULL,
  strip_tables = TRUE,
  strip_cov = TRUE,
  overwrite = FALSE,
  ci_threshold = 3.84,
  digits = NULL,
  run = TRUE,
  verbose = TRUE,
  n_cores = 1L,
  ...
)

# S3 method for class 'nm_profile'
print(x, ...)

# S3 method for class 'nm_profile'
plot(x, type = c("ofv", "dofv"), ...)
```

## Arguments

- ctl:

  Path to the reference control file. (Formerly `ctl_path`, still
  accepted with a warning; the former order, `nmfe_path` first, is also
  recognised.)

- nmfe_path:

  Path to the NONMEM launcher (e.g. `nmfe76.bat` on Windows, the
  `nmfe76` script on Linux or macOS).

- param:

  The parameter to profile: a `$THETA` label (matched whole-word,
  case-insensitive, e.g. `"CL"`) or an index (`3`, `"3"`, or
  `"THETA(3)"`).

- from, to, by:

  Numeric grid: first value, last value, and increment.

- output_dir:

  Output directory. Defaults to `<ctl_dir>/<ctl>_<param>_profile`.

- strip_tables, strip_cov:

  Logical; drop `$TABLE` / `$COV` from the generated controls. Both
  `TRUE` by default (recommended for a profile).

- overwrite:

  Logical; if `FALSE` (default) a run whose `.lst` already exists and
  minimised successfully is skipped (resume after interruption).

- ci_threshold:

  Numeric; delta-OFV cut-off for the interval. Default `3.84`
  (chi-square, 1 df, two-sided 95%). Set `NA` to disable. Note: for a
  variance component tested at its boundary (\\\omega^2=0\\) the
  appropriate one-sided cut-off is `2.71`, not `3.84`.

- digits:

  Decimals used to format fixed values and run labels. Auto- derived
  from `by` when `NULL`.

- run:

  Logical; if `FALSE`, generate the control files without launching
  NONMEM.

- verbose:

  Logical; print progress messages (default `TRUE`). Formerly `quiet`
  (the opposite), still accepted with a warning.

- n_cores:

  Number of grid values run at the same time (default 1). A whole
  number, `"auto"` (physical cores minus one, see
  [`nm_cores()`](https://nskoch.github.io/nmpilot/reference/nm_cores.md))
  or `"ask"` (interactive sessions: the cores are shown and the number
  is asked once per R session, then reused;
  `options(nmpilot.n_cores = NULL)` asks again). Each run is then
  started in the background in its own directory; if R is interrupted,
  the runs already started go on in the background and a later call with
  the same arguments picks up their results.

- ...:

  Ignored.

- x:

  An `nm_profile` object.

- type:

  Which figure to draw: `"ofv"` (OFV vs value) or `"dofv"` (delta-OFV vs
  value, with the `ci_threshold` line).

## Value

An object of class `nm_profile`: a list with `results` (a data frame of
`run`, `value`, `ofv`, `dOFV`, `min_ok`, `cov_ok`, `cov_warning`, `dir`,
`lst`), plus `ci`, `min_ofv`, `min_value`, `param`, `param_safe`,
`target`, `ci_threshold`, and `output_dir`. Use
[`print()`](https://rdrr.io/r/base/print.html) for a summary and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) (with
`type = "ofv"` or `"dofv"`) for figures.

## Details

Only `$DATA` is remapped (to an absolute path) and `$TABLE`/`$COV` are
stripped by default; every other record is preserved byte-for-byte. Each
run is written to its own sub-directory under `output_dir/runs/`.

## See also

[`nm_run()`](https://nskoch.github.io/nmpilot/reference/nm_run.md),
[`nm_extract_runs()`](https://nskoch.github.io/nmpilot/reference/nm_extract_runs.md),
[`nm_compare()`](https://nskoch.github.io/nmpilot/reference/nm_compare.md)

## Examples

``` r
if (FALSE) { # \dontrun{
prof <- nm_profile(
  ctl       = "C:/nm76g64/run/test/run_test.reg/run_test.ctl",
  nmfe_path = "C:/nm76g64/run/nmfe76.bat",
  param     = "CL", from = 0.1, to = 3, by = 0.05
)
prof                       # summary (print.nm_profile)
plot(prof)                 # OFV vs value
plot(prof, type = "dofv")  # dOFV vs value, with the 3.84 threshold
} # }
```
