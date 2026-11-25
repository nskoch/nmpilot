# Proposals to improve a base model without iterative NONMEM runs

Examines a base run and proposes changes to the residual error model,
the random effects and the covariate model.

## Usage

``` r
nm_optimize_base(
  run,
  covariates = NULL,
  deriv = NULL,
  nmfe_path = NULL,
  data = NULL,
  occ = NULL,
  level = 2L,
  cov_form = c("power", "exponential"),
  p_forward = 0.05,
  p_backward = 0.01,
  max_levels = 5L,
  objective = c("focei", "fast"),
  verbose = TRUE
)

# S3 method for class 'nm_optimize_base'
print(x, ...)

# S3 method for class 'nm_optimize_base'
plot(x, type = c("all", "residual", "random", "covariate"), ...)
```

## Arguments

- run:

  Run directory (e.g. `"run1.reg"`) or path prefix of the base run.

- covariates:

  Covariates to test (columns of the data set). `NULL`: no covariate
  analysis.

- deriv:

  Table written by the run of
  [`nm_derivatives_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_derivatives_ctl.md)
  (`<run>_deriv.dtab`). `NULL`: searched for next to the run, otherwise
  the function runs level 1 only (unless `nmfe_path` is given).

- nmfe_path:

  Path of nmfe: if given and no derivatives table exists, the
  derivatives control stream is written and run with
  [`nm_run()`](https://nskoch.github.io/nmpilot/reference/nm_run.md).

- data:

  Data set (path or data.frame) holding the covariates. Default: the
  `$DATA` file of the run.

- occ:

  Occasion column for inter-occasion variability (level 2; it must be in
  the derivatives table: `nm_derivatives_ctl(extra = "OCC")`).

- level:

  `1` or `2` (2 runs both levels when the table is available).

- cov_form:

  Form of continuous covariate effects: `"power"` (`(COV/median)^theta`)
  or `"exponential"` (`exp(theta*(COV-median))`).

- p_forward, p_backward:

  Significance levels of the forward inclusion and backward elimination
  of covariates (linearized model).

- max_levels:

  A numeric covariate with at most this number of distinct values is
  treated as categorical.

- objective:

  Objective of the linearized model. `"focei"` (default): the FOCE-I
  objective of NONMEM with interaction (individual modes of the full
  individual objective, residual-variance term in the determinant); at
  the NONMEM estimates it reproduces the OFV of the base run. `"fast"`:
  an approximation (modes given by the linear predictor with the
  residual variance held fixed), about three times faster, whose OFV
  differs from NONMEM's by a few units when the residual error depends
  on the ETAs; the dOFV are close.

- verbose:

  Print progress.

- x:

  An `nm_optimize_base` object.

- ...:

  Ignored.

- type:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): `"all"`,
  `"residual"`, `"random"` or `"covariate"`. Formerly `area`, still
  accepted with a warning.

## Value

An object of class `nm_optimize_base`: lists `level1` (`residual`,
`random`, `pairs`, `covariates`) and `level2` (`check`, `residual`,
`random`, `covariates`, `scm`, `proposal`), and `proposal_text`.

## Details

**Level 1** uses only the outputs of the base run (`.phi`, `.ext`, a
table with `CWRES`, `IPRED` and `TIME`, and the data set): resmod-like
models on the CWRES (power, time-varying, autocorrelated residuals,
individual residual variance), EBE correlations, distribution and
shrinkage, and an ETA-covariate screening.

**Level 2** needs one extra NONMEM evaluation run written by
[`nm_derivatives_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_derivatives_ctl.md).
Its table of derivatives gives a first-order linearization of the model
around the EBEs; every candidate is re-estimated in R under the FOCE-I
objective of the linearized model: residual models (additive,
proportional, combined, power, time-varying), OMEGA structure (full
block, pairwise blocks, removal of an ETA, inter-occasion variability),
covariates on the ETAs (univariate tests, then forward inclusion and
backward elimination), and finally the combined proposal.

The changes in objective function (`dOFV`) are those of the linearized
model: they approximate the changes of the real model, and the proposal
must be confirmed by a NONMEM run of the modified model. The quality of
the approximation is reported in `check` (linearized vs NONMEM OFV of
the base model).

## References

Khandelwal A, Harling K, Jonsson EN, Hooker AC, Karlsson MO. A fast
method for testing covariates in population PK/PD models. AAPS J
2011;13:464-72. Ibrahim MMA, Nordgren R, Kjellsson MC, Karlsson MO.
Model-based residual post-processing for residual model identification.
AAPS J 2018;20:81.

## See also

[`nm_derivatives_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_derivatives_ctl.md),
[`nm_wasserstein()`](https://nskoch.github.io/nmpilot/reference/nm_wasserstein.md)

## Examples

``` r
ex <- system.file("extdata", "example", "run1.reg", package = "nmpilot")
o <- nm_optimize_base(ex, covariates = c("WT", "GENO"), level = 1, verbose = FALSE)
o
#> Optimization of base run run1 (NONMEM OFV -1147.34)
#> 
#> == Level 1: outputs of the base run ==
#> 
#> Residual error (models on CWRES, dOFV vs CWRES ~ N(mu, sigma) with random intercept):
#>            model  dOFV estimate                                   meaning
#>            power -4.28   0.0479        delta (variance ~ IPRED^(2*delta))
#>     time_varying -4.49   0.9100                 SD ratio TIME > 4 vs <= 4
#>  autocorrelation -2.24   0.0371 half-life of the correlation (time units)
#>       iiv_on_eps -2.26   0.0778          SD of the ETA on EPS (log scale)
#>  suggestion
#>           -
#>           -
#>           -
#>           -
#> 
#> Random effects:
#>        eta  omega shrinkage_pct skewness p_shapiro suggestion
#>  ETA1 (CL) 0.1430           1.7     0.54      0.17          -
#>   ETA2 (V) 0.0843           4.9    -0.16      0.60          -
#>  ETA3 (KA) 0.1780          14.0    -0.18      0.96          -
#> 
#> EBE correlations:
#>      eta_1     eta_2     r     p max_shrinkage_pct omega_cov suggestion
#>  ETA1 (CL)  ETA2 (V)  0.18 0.180               4.9         0          -
#>  ETA1 (CL) ETA3 (KA) -0.16 0.220              14.0         0          -
#>   ETA2 (V) ETA3 (KA) -0.28 0.029              14.0         0          -
#> 
#> ETA-covariate screening (Bonferroni-adjusted):
#>        eta covariate statistic       p   p_adj shrinkage_pct
#>   ETA2 (V)        WT     0.638 4.1e-08 2.5e-07           4.9
#>  ETA1 (CL)      GENO     0.549 1.1e-06 6.6e-06           1.7
#>  ETA3 (KA)        WT    -0.260 4.4e-02 2.7e-01          14.0
#>  ETA1 (CL)        WT     0.173 1.9e-01 1.0e+00           1.7
#>   ETA2 (V)      GENO     0.075 6.8e-01 1.0e+00           4.9
#>  ETA3 (KA)      GENO     0.083 9.3e-01 1.0e+00          14.0
#>           suggestion
#>  candidate covariate
#>  candidate covariate
#>                    -
#>                    -
#>                    -
#>                    -
#> 
#> Summary:
#>   - [level 1, covariate] WT on ETA2 (V): candidate covariate (adjusted p 2.5e-07)
#>   - [level 1, covariate] GENO on ETA1 (CL): candidate covariate (adjusted p 6.6e-06)
```
