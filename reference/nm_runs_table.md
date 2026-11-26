# Summary table of all the runs of a project

Builds one table with one row per run of a project and, when
`output_file` is given, writes it to a Word document. Each run is shown
with its description (text of `$PROBLEM`), the run it was derived from
(its parent, read in the control stream as in
[`nm_model_diagram()`](https://nskoch.github.io/nmpilot/reference/nm_model_diagram.md)),
its OFV, the change of OFV from its parent and from the best run, its
number of estimated parameters, its AIC and BIC, its minimisation and
covariance status, its largest ETA shrinkage and the share of parameters
estimated with a relative standard error above 50%.

## Usage

``` r
nm_runs_table(
  x,
  parents = NULL,
  sort = c("tree", "ofv", "run"),
  output_file = NULL,
  overwrite = FALSE,
  title = "Summary of the runs",
  verbose = TRUE
)
```

## Arguments

- x:

  Result of
  [`nm_extract_runs()`](https://nskoch.github.io/nmpilot/reference/nm_extract_runs.md),
  or the path of an analysis directory, which is then read with
  `nm_extract_runs(x)`.

- parents:

  Optional links that complete or replace those read in the control
  streams, as in
  [`nm_model_diagram()`](https://nskoch.github.io/nmpilot/reference/nm_model_diagram.md):
  a named character vector (`c(run2 = "run1")`) or a data frame with
  columns `run` and `parent`.

- sort:

  Order of the rows: `"tree"` (default; each run right after its parent,
  as in
  [`nm_model_diagram()`](https://nskoch.github.io/nmpilot/reference/nm_model_diagram.md)),
  `"ofv"` (best OFV first) or `"run"` (run names, numbers sorted as
  numbers).

- output_file:

  Path of the Word file (.docx). `NULL` (default): nothing is written.

- overwrite:

  Logical. Replace an existing `output_file` (default `FALSE`).

- title:

  Title written above the table in the Word file.

- verbose:

  Print a message when the file is written.

## Value

An object of class `nm_runs_table`: a data frame with one row per run
and the columns `run`, `description`, `parent`, `ofv`, `dofv_parent`
(OFV of the run minus OFV of its parent), `dofv_best` (minus the lowest
OFV of the project), `n_par`, `aic`, `bic`, `minimization`,
`covariance`, `max_eta_shrink` (%), `pct_rse_gt50` (% of the parameters
with an RSE above 50%) and `status` (as in
[`nm_model_diagram()`](https://nskoch.github.io/nmpilot/reference/nm_model_diagram.md)).
The Word file, when written, has the same rows in landscape format; its
path is kept in the attribute `output_file`.

## See also

[`nm_extract_runs()`](https://nskoch.github.io/nmpilot/reference/nm_extract_runs.md),
[`nm_model_diagram()`](https://nskoch.github.io/nmpilot/reference/nm_model_diagram.md),
[`nm_compare_runs()`](https://nskoch.github.io/nmpilot/reference/nm_compare_runs.md),
[`nm_extract_run_table()`](https://nskoch.github.io/nmpilot/reference/nm_extract_run_table.md)
(table of the estimates of one run)

## Examples

``` r
ex <- system.file("extdata", "example", package = "nmpilot")
tab <- nm_runs_table(ex, parents = c(run2 = "run1"))
tab
#> Summary of 2 run(s)
#>  Run  Parent OFV      dOFV (parent) dOFV (best) N par AIC     BIC    
#>  run1        -1147.34               76.69       8     -1131.3 -1097.0
#>  run2 run1   -1224.04 -76.69        0.00        10    -1204.0 -1161.1
#>  Covariance Max shrinkage (%) RSE > 50% (%)
#>  COMPLETED  14.7              0.0          
#>  COMPLETED  14.4              0.0          
# \donttest{
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  f <- file.path(tempdir(), "runs_summary.docx")
  nm_runs_table(ex, parents = c(run2 = "run1"), output_file = f, overwrite = TRUE)
}
#> Word file written: /tmp/Rtmp4BwTni/runs_summary.docx
#> Summary of 2 run(s)
#>  Run  Parent OFV      dOFV (parent) dOFV (best) N par AIC     BIC    
#>  run1        -1147.34               76.69       8     -1131.3 -1097.0
#>  run2 run1   -1224.04 -76.69        0.00        10    -1204.0 -1161.1
#>  Covariance Max shrinkage (%) RSE > 50% (%)
#>  COMPLETED  14.7              0.0          
#>  COMPLETED  14.4              0.0          
#> Word file: /tmp/Rtmp4BwTni/runs_summary.docx 
# }
```
