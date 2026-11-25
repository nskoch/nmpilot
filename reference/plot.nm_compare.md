# Diagnostic plots for `nm_compare` objects

Diagnostic plots for `nm_compare` objects

## Usage

``` r
# S3 method for class 'nm_compare'
plot(x, type = "all", n_band = 200L, ...)
```

## Arguments

- x:

  An `nm_compare` object.

- type:

  Character. One or several of `"cwres_scatter"`, `"npde_scatter"`,
  `"iofv_caterpillar"`, `"delta_density"`, `"qq_npde"`, `"gain_map"`
  (Delta\|CWRES\| versus TIME and PRED, loess with a 95 %
  subject-bootstrap band), `"overlay_gof"` (diagnostics of both runs
  superimposed), `"bootstrap"` (subject bootstrap of dOFV),
  `"iofv_covariate"` (delta iOFV versus the covariates), `"parameters"`
  (estimates with 95 % Wald intervals), `"eta_scatter"` (EBEs of the two
  runs), or `"all"` (default). Formerly `which`, still accepted with a
  warning.

- n_band:

  Number of subject resamples for the bands of `"gain_map"` and
  `"overlay_gof"`.

- ...:

  Unused (except the former name `which`).

## Value

A `ggplot` object (single plot) or a list of `ggplot` objects (multiple
plots). Requires the `ggplot2` package.
