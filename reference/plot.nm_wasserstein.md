# Plots for `nm_wasserstein` objects

Plots for `nm_wasserstein` objects

## Usage

``` r
# S3 method for class 'nm_wasserstein'
plot(x, type = c("ecdf", "r2w"), covariate = NULL, ...)
```

## Arguments

- x:

  An `nm_wasserstein` object.

- type:

  `"ecdf"` (default): empirical cumulative distribution of the ETA in
  the two groups, one panel per run and comparison, annotated with W,
  its p-value, t and its p-value. `"r2w"`: bar chart of the
  Wasserstein-based R2 per run.

- covariate:

  Optional covariate(s) to show (ECDF plot only).

- ...:

  Ignored.

## Value

A ggplot object.
