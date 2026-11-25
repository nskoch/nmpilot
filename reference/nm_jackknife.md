# Leave-one-subject-out jackknife of a NONMEM model

Removes one subject at a time from the data set, re-estimates the model
with NONMEM and measures how much each subject influences the fit. The
primary metric is the likelihood influence \\\|\Delta OFV\| =
\|OFV\_{(-i)} - OFV\_{ref}\|\\; the secondary metrics are the relative
THETA shifts \\(\theta\_{(-i)} - \theta\_{ref}) / \theta\_{ref}\\, their
L2 norm, their maximum, and an L2 norm standardised by the reference
standard errors when the reference covariance step is usable. Subjects
are ranked by \\\|\Delta OFV\|\\. The analysis is descriptive (an
influence diagnostic, not a formal outlier test).

## Usage

``` r
nm_jackknife(
  ctl,
  nmfe_path = NULL,
  lst = NULL,
  reference_data = NULL,
  id_col = "ID",
  comment_char = NULL,
  sep = ",",
  output_dir = NULL,
  overwrite_existing_runs = FALSE,
  strip_tables = TRUE,
  strip_cov = FALSE,
  dtheta_threshold = 0.01,
  save_figures = TRUE,
  run = TRUE,
  verbose = TRUE,
  n_cores = 1L
)

# S3 method for class 'nm_jackknife'
print(x, digits = 4, ...)
```

## Arguments

- ctl:

  Path to the reference control stream (`.ctl`), already run.

- nmfe_path:

  Path to the nmfe launcher (e.g. `"C:/nm76g64/run/nmfe76.bat"` on
  Windows, `"/opt/nm760/run/nmfe76"` on Linux or macOS). Needed when
  `run = TRUE`.

- lst:

  Reference output. Default: the `.lst` of the same name next to `ctl`.

- reference_data:

  Data file. `NULL` (default): read from `$DATA`.

- id_col:

  Label of the subject column in `$INPUT` (default `"ID"`).

- comment_char:

  Data comment/ignore character. `NULL` (default): the single-character
  `IGNORE` of `$DATA` (`"#"`, `"@"`, a letter...); `""`: none.

- sep:

  Field separator of the data file (default `","`).

- output_dir:

  Working directory (default `"<ctl_dir>/<run>_jackknife"`). Holds
  `jk_data/`, `jk_ctl/`, `jk_run/` (one subdirectory per run),
  `figures/` and the two CSV files.

- overwrite_existing_runs:

  `FALSE` (default): skip the runs already minimised (resume).

- strip_tables:

  Remove `$TABLE` from the jackknife control streams (default `TRUE`).

- strip_cov:

  Remove `$COVARIANCE` from the jackknife runs (default `FALSE`). `TRUE`
  is faster and does not change the influence metrics (the reference SE
  come from the reference run).

- dtheta_threshold:

  Relative THETA shift flagged as large (default 0.01).

- save_figures:

  Write the three figures to `output_dir/figures` (default `TRUE`).

- run:

  `FALSE`: only write the data sets and control streams (no NONMEM).

- verbose:

  Progress messages (default `TRUE`).

- n_cores:

  Number of jackknife runs at the same time (default 1). A whole number,
  `"auto"` (physical cores minus one, see
  [`nm_cores()`](https://nskoch.github.io/nmpilot/reference/nm_cores.md))
  or `"ask"` (interactive sessions: the cores are shown and the number
  is asked once per R session, then reused;
  `options(nmpilot.n_cores = NULL)` asks again). Each run is then
  started in the background in its own directory; if R is interrupted,
  the runs already started go on in the background and a later call with
  the same arguments picks up their results.

- x:

  An `nm_jackknife` object.

- digits:

  Significant digits in the printed table (default 4).

- ...:

  Ignored.

## Value

An object of class `nm_jackknife`: the reference OFV, THETA and SE,
`all_runs` (every replicate, including the failed ones), `summary`
(successful runs ranked by \\\|\Delta OFV\|\\, with `rank`,
`rank_dtheta` and, when available, `rank_dtheta_z`) and the settings.
`jackknife_all_runs.csv` and `jackknife_summary.csv` are written to
`output_dir`.

## Details

The method is that of the standalone script nonmem-jackknife-influence
(Simon N, doi:10.5281/zenodo.22257994). The data lines are copied
verbatim (positional reading through `$INPUT`), so number formatting and
rows removed by a conditional `IGNORE=(...)` are preserved. Runs already
minimised are skipped on a second call (resume), unless
`overwrite_existing_runs = TRUE`.

## See also

[`nm_bootstrap()`](https://nskoch.github.io/nmpilot/reference/nm_bootstrap.md),
[`nm_profile()`](https://nskoch.github.io/nmpilot/reference/nm_profile.md)

## Examples

``` r
if (FALSE) { # \dontrun{
jk <- nm_jackknife("C:/study/run2.reg/run2.ctl",
                   nmfe_path = "C:/nm76g64/run/nmfe76.bat")
jk
plot(jk, "bar")
plot(jk, "scatter")
} # }
```
