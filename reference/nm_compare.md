# Compare two NONMEM runs

Paired comparison of two NONMEM runs at the observation level (CWRES,
NPDE) and at the subject level (individual OFV from `.phi`). Provides a
summary of residual diagnostics, per-subject iOFV deltas, and identifies
which subjects drive the OFV difference.

## Usage

``` r
nm_compare(
  run1,
  run2,
  labels = c("Run 1", "Run 2"),
  join_keys = c("ID", "TIME"),
  n_par = NULL,
  nested = FALSE,
  covariates = NULL,
  data = NULL,
  n_boot = 2000L,
  seed = 2026L
)

# S3 method for class 'nm_compare'
print(x, ...)

# S3 method for class 'nm_compare'
summary(object, ...)
```

## Arguments

- run1, run2:

  Path to a run directory or run prefix (without extension). Each
  location must contain at least a `sdtab*` table. `.phi` and `.lst`
  files are used if available for iOFV and global OFV respectively.

- labels:

  Character vector of length 2. Display labels for the two runs.

- join_keys:

  Character vector. Columns used to pair observations between runs.
  Default `c("ID", "TIME")`. Add `"OCC"` for BOV models, `"CMT"` for
  multi-DV models.

- n_par:

  Numbers of estimated parameters of the two runs (length 2). Default:
  read from the `.ext` files (parameters not fixed).

- nested:

  `TRUE` if one model is nested in the other: a likelihood ratio test is
  added to the Vuong and Clarke tests.

- covariates:

  Covariates (columns of the data set) against which the per-subject
  delta iOFV is examined.

- data:

  Data set (path or data.frame) holding the covariates. Default: the
  `$DATA` file of run 1.

- n_boot:

  Number of subject resamples for the bootstrap of the delta iOFV.

- seed:

  Seed of the bootstrap.

- x:

  An `nm_compare` object.

- ...:

  Unused.

- object:

  An `nm_compare` object.

## Value

An object of class `nm_compare`, a list with:

- `data_obs`:

  Observation-level paired data.frame with `_1`/`_2` suffixed columns
  and delta columns (`dCWRES`, `dAbsCWRES`, `dNPDE`, `dAbsNPDE`).

- `data_subj`:

  Per-subject summary including `iOFV_1`, `iOFV_2`, `diOFV`, and
  residual aggregates.

- `metrics`:

  Global metrics table (OFV, RMSE(CWRES), Mean \|CWRES\|, proportions of
  outliers, NPDE moments). Wilcoxon paired p-value stored as attribute
  `wilcox_absCWRES_p`.

- `meta`:

  Labels, global OFV, numbers of parameters, sample sizes, dataset
  compatibility check.

- `tests`:

  Vuong, Clarke (and LRT) tests, subject bootstrap and influence (`NULL`
  without `.phi` files).

- `parameters`:

  Estimates, SE, RSE and shrinkage of the two runs matched by parameter
  name, with their ratio.

- `etas`:

  Paired EBEs of the ETAs present in both runs.

- `covariate_tests`:

  Delta iOFV versus each covariate.

- `run_info`:

  Run-level criteria (RSE, condition number, shrinkage, status).

[`summary()`](https://rdrr.io/r/base/summary.html): a data.frame of
decision criteria (class `nm_compare_summary`), one row per criterion
with the value for each run and the run it favours.

## Details

Observations are merged on `join_keys` (default `c("ID", "TIME")`). If
the two datasets differ (e.g. different BLQ handling), only the
intersection is retained and a warning is issued.

**Tests on the individual contributions.** With
`d_i = (iOFV_1,i - iOFV_2,i)/2`, the difference of the individual
log-likelihoods, the Vuong (1989) statistic is
`z = (sum(d_i) - c) / (sqrt(N) * sd(d_i))` and the Clarke (2007) test is
a sign test on `d_i - c/N`, where `N` is the number of subjects and `c`
the correction for the number of parameters: none, `k_2 - k_1`
(AIC-like) or `(k_2 - k_1) log(N) / 2` (BIC-like). Both apply to
non-nested models fitted to the same data; the subjects are the
independent units. The individual contributions are those of the
estimation method (FOCE: the iOFV of the `.phi` file).

**Bootstrap.** The subjects are resampled with their delta iOFV, without
re-estimation: the distribution measures how much the difference depends
on a few subjects, conditionally on the estimates (optimistic for the
larger model).

**Influence.** `n_reverse` is the number of subjects with the largest
contributions whose removal reverses the sign of dOFV (as `dOFV.vs.id`
in xpose4).

## References

Vuong QH. Likelihood ratio tests for model selection and non-nested
hypotheses. Econometrica 1989;57:307-33.

Clarke KA. A simple distribution-free test for nonnested model
selection. Political Analysis 2007;15:347-63.

## See also

[`nm_extract_runs()`](https://nskoch.github.io/nmpilot/reference/nm_extract_runs.md),
[`plot.nm_compare()`](https://nskoch.github.io/nmpilot/reference/plot.nm_compare.md),
`print.nm_compare()`

## Examples

``` r
if (FALSE) { # \dontrun{
  cmp <- nm_compare("runs/run42", "runs/run43",
                    labels = c("Base", "+ CLcr"))
  print(cmp)
  plot(cmp)
  plot(cmp, type = "iofv_caterpillar")
  summary(cmp)
  plot(cmp, type = c("gain_map", "overlay_gof"))
} # }
```
