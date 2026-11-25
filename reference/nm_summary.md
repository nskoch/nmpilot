# Diagnostic summary of a NONMEM-like dataset

Diagnostic summary of a NONMEM-like dataset

## Usage

``` r
nm_summary(data, verbose = TRUE)
```

## Arguments

- data:

  A data.frame from
  [`nm_read_csv()`](https://nskoch.github.io/nmpilot/reference/nm_read_csv.md).

- verbose:

  Print the summary to the console (default `TRUE`).

## Value

A one-row data.frame with: `n_subjects`, `n_records`, `n_observations`,
`n_blq`, `pct_blq`, `time_min`, `time_max`, `time_unit`, `unique_doses`
(distinct positive doses), `n_cmt_obs`.
