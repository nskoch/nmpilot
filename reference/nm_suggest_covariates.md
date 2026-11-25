# Covariate screening on individual PK estimates before modelling

Applies the Wasserstein-distance comparison of
[`nm_wasserstein()`](https://nskoch.github.io/nmpilot/reference/nm_wasserstein.md)
to the individual estimates obtained by
[`nm_suggest_model()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_model.md)
(NCA clearance, central volume and absorption parameter of the rich
intervals, log scale), before any NONMEM run. For each covariate,
subjects are split into groups (each level against the first for a
categorical covariate; Q4 against Q1, or all pairs of quartiles, for a
continuous one) and the distributions of the log estimates are compared
with the 1-Wasserstein distance, a permutation test and a Welch t-test.
It points to covariates worth testing in the model; it is no substitute
for the covariate model (individual NCA estimates carry their own
estimation noise, and only subjects with a rich interval are used).

## Usage

``` r
nm_suggest_covariates(
  x,
  covariates,
  data = NULL,
  parameters = c("CL", "Vc"),
  n_perm = 2000L,
  seed = 2025L,
  max_levels = 5L,
  continuous = c("extremes", "pairwise")
)

# S3 method for class 'nm_suggest_covariates'
print(x, ...)
```

## Arguments

- x:

  An `nm_model_suggestion` object (individual estimates are needed, i.e.
  at least 3 rich intervals).

- covariates:

  Character vector of covariate columns.

- data:

  Covariate data (data.frame or CSV path with an `ID` column). Default:
  the dataset analysed by
  [`nm_suggest_model()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_model.md).

- parameters:

  Individual estimates to analyse among `"CL"`, `"Vc"`, `"KA"`, `"DUR"`,
  `"MTT"` (default `c("CL", "Vc")`).

- n_perm:

  Number of permutations (default 2000).

- seed:

  Random seed of the permutation tests; the user's random state is
  restored.

- max_levels:

  A numeric covariate with at most this many distinct values is treated
  as categorical (default 5).

- continuous:

  `"extremes"` (Q4 vs Q1, default) or `"pairwise"`.

- ...:

  Ignored.

## Value

A data.frame of class `nm_suggest_covariates`: `parameter`, `covariate`,
`type`, `comparison`, `n_1`, `n_0`, `ratio_median` (median of group 1 /
median of group 0), `W`, `p_W`, `p_t`, `p_W_adj` (Bonferroni over the
comparisons of a covariate and parameter).

## References

Simon N, Hulot JS, von Fabeck K. Distribution-based covariate assessment
using Wasserstein distance in population pharmacokinetic models. *Front
Pharmacol.* 2026;17:1804989.
[doi:10.3389/fphar.2026.1804989](https://doi.org/10.3389/fphar.2026.1804989)

## See also

[`nm_suggest_model()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_model.md),
[`nm_wasserstein()`](https://nskoch.github.io/nmpilot/reference/nm_wasserstein.md)

## Examples

``` r
ex <- system.file("extdata", "example", "example_data.csv", package = "nmpilot")
s <- nm_suggest_model(ex, verbose = FALSE, design_eval = FALSE, max_profiles = 15)
nm_suggest_covariates(s, covariates = c("WT", "GENO"), n_perm = 500)
#> Covariate screening on individual estimates (log scale, Wasserstein distance)
#>  parameter covariate        type comparison n_1 n_0 ratio_median     W   p_W
#>         CL        WT  continuous   Q4 vs Q1   3   4        0.834 0.255 0.902
#>         Vc        WT  continuous   Q4 vs Q1   3   4        1.309 0.404 0.084
#>         CL      GENO categorical     1 vs 0   2  13        1.577 0.428 0.070
#>         Vc      GENO categorical     1 vs 0   2  13        0.949 0.185 0.932
#>    p_t p_W_adj
#>  0.596   0.902
#>  0.111   0.084
#>  0.438   0.070
#>  0.952   0.932
```
