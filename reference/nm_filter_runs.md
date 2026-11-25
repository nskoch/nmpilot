# Filter NONMEM runs and recompute relative criteria

Filter NONMEM runs and recompute relative criteria

## Usage

``` r
nm_filter_runs(
  x,
  run_filter = c("all", "successful", "successful_no_warning", "successful_no_comment",
    "clean")
)
```

## Arguments

- x:

  An object of class `nm_runs`.

- run_filter:

  Filter to apply: `"all"` (default), `"successful"`,
  `"successful_no_warning"`, `"successful_no_comment"` or `"clean"`.

## Value

A filtered `nm_runs` object. The `$runs` and `$params` tables are
filtered together; dOFV, rank_ofv, dAIC and dBIC are recomputed on the
remaining runs.
