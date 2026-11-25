# Expected precision of a suggested model for a sampling design

Computes, before any NONMEM run, the expected relative standard errors
(RSE) of the parameters of a candidate model from
[`nm_suggest_model()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_model.md)
and the expected eta-shrinkage, for the sampling design of the dataset
or for another design. The population Fisher information matrix is
computed with the first-order (FO) approximation, block-diagonal between
fixed effects and variance parameters (as in PopED's reduced FIM), with
numerical derivatives; predictions use the analytic solution of the
linear model, with superposition of all doses (`ADDL`/`II` expanded).
The expected shrinkage of each ETA is
`1 - sqrt(1 - mean(C_kk / omega_k))`, where `C` is the posterior
covariance of the ETAs of each subject in the linearised model.

## Usage

``` r
nm_suggest_design(x, candidate = 1L, times = NULL, n_subjects = NULL)

# S3 method for class 'nm_suggest_design'
print(x, ...)

# S3 method for class 'nm_suggest_design'
plot(x, ...)
```

## Arguments

- x:

  An `nm_model_suggestion` object.

- candidate:

  Candidate number (default 1).

- times:

  Optional sampling times after the first dose, applied to every subject
  in place of the observed times (to compare designs).

- n_subjects:

  Optional number of subjects: the information is scaled from the
  subjects of the dataset to this number.

- ...:

  Ignored.

## Value

An object of class `nm_suggest_design`: `parameters` (parameter, type,
value, se, rse_pct), `shrinkage` (parameter, shrinkage_pct),
`n_subjects`, `n_obs`, `spec`, `times`.

## Details

The values depend on the initial values of the candidate (THETA, IIV,
SIGMA), which are approximate: they rank designs and parameters rather
than predict the RSE reported by NONMEM. Michaelis-Menten candidates are
not supported (no analytic solution).

## See also

[`nm_suggest_model()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_model.md)

## Examples

``` r
ex <- system.file("extdata", "example", "example_data.csv", package = "nmpilot")
s <- nm_suggest_model(ex, verbose = FALSE, design_eval = FALSE, max_profiles = 15)
nm_suggest_design(s)
#> Expected precision: 1 cpt | first-order | proportional error | 60 subjects | design of the dataset
#>    parameter  type   value      se rse_pct
#>           CL theta  5.7200 0.23700     4.1
#>            V theta 55.6000 2.57000     4.6
#>           KA theta  1.1100 0.12500    11.3
#>    OMEGA(CL) omega  0.0975 0.01880    19.3
#>     OMEGA(V) omega  0.1160 0.02340    20.2
#>    OMEGA(KA) omega  0.6930 0.13900    20.1
#>  SIGMA(prop) sigma  0.0377 0.00281     7.4
#> 
#> Expected eta-shrinkage (%)
#>  parameter shrinkage_pct
#>         CL           2.7
#>          V           5.0
#>         KA           4.6
# a sparse design with 3 samples per subject
nm_suggest_design(s, times = c(1, 4, 12))
#> Expected precision: 1 cpt | first-order | proportional error | 60 subjects | design: times given
#>    parameter  type   value     se rse_pct
#>           CL theta  5.7200 0.2650     4.6
#>            V theta 55.6000 3.1800     5.7
#>           KA theta  1.1100 0.1480    13.3
#>    OMEGA(CL) omega  0.0975 0.0259    26.6
#>     OMEGA(V) omega  0.1160 0.0377    32.6
#>    OMEGA(KA) omega  0.6930 0.1930    27.8
#>  SIGMA(prop) sigma  0.0377 0.0159    42.3
#> 
#> Expected eta-shrinkage (%)
#>  parameter shrinkage_pct
#>         CL          12.2
#>          V          19.2
#>         KA          15.1
```
