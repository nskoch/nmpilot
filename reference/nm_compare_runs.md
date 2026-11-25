# Compare the parameters of several runs side by side

Compare the parameters of several runs side by side

## Usage

``` r
nm_compare_runs(
  x,
  what = c("value", "rse_pct", "se", "shrinkage"),
  types = NULL,
  format = FALSE,
  digits = 3,
  align_by = c("name", "index"),
  run_filter = c("all", "successful", "successful_no_warning", "successful_no_comment",
    "clean")
)
```

## Arguments

- x:

  An `nm_runs` object or a params tibble.

- what:

  `"value"`, `"rse_pct"`, `"se"` or `"shrinkage"`.

- types:

  Parameter types to keep. `NULL` = all.

- format:

  If `TRUE`, cells are formatted as "value (RSE%)".

- digits:

  Significant digits when `format = TRUE`.

- align_by:

  `"name"` (default) for semantic alignment; `"index"` for a strict
  comparison by NONMEM index.

- run_filter:

  Run filter: `"all"` (default), `"successful"`,
  `"successful_no_warning"`, `"successful_no_comment"` or `"clean"`.

## Value

A wide tibble with one row per parameter and one column per run.

## Examples

``` r
if (FALSE) { # \dontrun{
nm_compare_runs(res, types = "THETA", format = TRUE)
nm_compare_runs(res, types = "THETA", format = TRUE, run_filter = "successful")
nm_compare_runs(res, types = "THETA", format = TRUE, run_filter = "clean")
res_clean <- nm_filter_runs(res, "clean")
res_clean$runs
} # }
```
