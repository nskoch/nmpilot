# Distribution-based covariate assessment with the Wasserstein distance

Implements the method of Simon, Hulot and von Fabeck (2026). For each
run and each covariate, the empirical Bayes estimates of one ETA (read
from the `.phi` file) are split into two groups of subjects and compared
with the 1-Wasserstein distance `W` and a permutation test, together
with a Welch t-test for reference. A covariate that the model accounts
for leaves no residual distributional signal in the ETA: comparing a
base and a final model shows whether the signal has been removed. A
Wasserstein-based R2 (`R2W`) compares the distance of each run's ETA
distribution to 0 with the one of a reference run.

## Usage

``` r
nm_wasserstein(
  runs,
  covariates,
  eta = "CL",
  data = NULL,
  ref = 1,
  n_perm = 10000L,
  seed = 2025L,
  max_levels = 5L,
  continuous = c("extremes", "pairwise")
)

# S3 method for class 'nm_wasserstein'
print(x, ...)
```

## Arguments

- runs:

  Character vector of run directories (`.reg`, `.nm7` or any directory
  holding the `.phi` file) or run path prefixes. Names, if any, are used
  as labels; otherwise the run names are used.

- covariates:

  Character vector of covariate columns of the data set.

- eta:

  The ETA to analyse: a PK parameter name (e.g. `"CL"`, identified from
  the control stream with
  [`nm_identify_parameters()`](https://nskoch.github.io/nmpilot/reference/nm_identify_parameters.md)),
  `"ETA(n)"` or an integer index.

- data:

  Covariate data: a data.frame or the path of a CSV file with an `ID`
  column (one row per subject is kept, the first one). If `NULL`
  (default), the data set named in `$DATA` of the first run's control
  stream is read.

- ref:

  Reference run for `R2W`: index or label (default 1).

- n_perm:

  Number of permutations (default 10000).

- seed:

  Random seed for the permutation tests (default 2025). Each test uses
  `seed + k`, so results are reproducible and independent of the order
  of the other tests. The user's random state is restored.

- max_levels:

  A numeric covariate with at most this many distinct values is treated
  as categorical (default 5).

- continuous:

  How continuous covariates are split into quartiles (defined once on
  the subjects of the covariate data): `"extremes"` (default) compares
  the highest quartile Q4 with the lowest Q1; `"pairwise"` compares
  every pair of quartiles (6 comparisons), as in the pairwise analysis
  of the original article. Adjusted p-values (`p_W_adj`, `p_t_adj`) use
  the Bonferroni correction over the comparisons of the same covariate
  within a run.

- x:

  An `nm_wasserstein` object.

- ...:

  Ignored.

## Value

An object of class `nm_wasserstein`, a list with:

- `tests`:

  one row per run x covariate comparison: `run`, `covariate`, `type`,
  `comparison`, `n_1`, `n_0`, `W`, `p_W`, `t_stat`, `p_t`, and the
  Bonferroni-adjusted `p_W_adj` and `p_t_adj` (equal to `p_W` and `p_t`
  when a covariate gives a single comparison).

- `r2w`:

  one row per run: `run`, `eta`, `n`, `W_to_zero`, `R2W`
  (`1 - W_to_zero / W_to_zero[ref]`), `omega` and `shrinkage_pct`
  (eta-shrinkage `1 - sd(EBE)/sqrt(omega)`, when the `.ext` is found). A
  large difference in shrinkage between runs can make the ETA
  distribution narrower without any real gain.

- `data`:

  the EBE merged with the comparison groups (long format), used by
  [`plot.nm_wasserstein()`](https://nskoch.github.io/nmpilot/reference/plot.nm_wasserstein.md).

- `settings`:

  the arguments used.

## References

Simon N, Hulot JS, von Fabeck K. Distribution-based covariate assessment
using Wasserstein distance in population pharmacokinetic models. *Front
Pharmacol.* 2026;17:1804989.
[doi:10.3389/fphar.2026.1804989](https://doi.org/10.3389/fphar.2026.1804989)

## See also

[`nm_identify_parameters()`](https://nskoch.github.io/nmpilot/reference/nm_identify_parameters.md),
[`nm_compare()`](https://nskoch.github.io/nmpilot/reference/nm_compare.md)

## Examples

``` r
if (FALSE) { # \dontrun{
w <- nm_wasserstein(
  runs       = c(BASE = "C:/nm76g64/run/mtx/base.reg",
                 FINAL = "C:/nm76g64/run/mtx/final.reg"),
  covariates = c("CLCR", "CT01"),
  eta        = "CL"
)
w                      # tests and R2W
plot(w)                # ECDF of the ETA by group, per run and covariate
plot(w, type = "r2w")  # Wasserstein-based R2

# All pairs of CLCR quartiles, with Bonferroni-adjusted p-values
wp <- nm_wasserstein(c(BASE = "base.reg", FINAL = "final.reg"),
                     covariates = "CLCR", eta = "CL", continuous = "pairwise")
wp$tests
} # }
```
