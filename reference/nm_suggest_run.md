# Write, run and rank the candidates of a model suggestion

Writes the control streams of the chosen candidates of
[`nm_suggest_model()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_model.md)
with
[`nm_write_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_write_ctl.md)
(`from = x`), launches them with
[`nm_run()`](https://nskoch.github.io/nmpilot/reference/nm_run.md),
reads their objective function value and ranks them by BIC (OFV + number
of estimated parameters x log(number of observations)). The rank by BIC
is compared with the rank of the data-driven score (Spearman
correlation), which shows whether the suggestion anticipated the NONMEM
result. A lag time of a candidate is added after `$PK` as
`ALAG1 = THETA(n)` with its initial value, since
[`nm_write_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_write_ctl.md)
writes none.

## Usage

``` r
nm_suggest_run(
  x,
  analysis_dir,
  nmfe_path = NULL,
  candidates = 1:3,
  prefix = "sugg",
  data = NULL,
  run = TRUE,
  overwrite = FALSE,
  verbose = TRUE,
  n_cores = 1L,
  ...
)

# S3 method for class 'nm_suggest_run'
print(x, ...)
```

## Arguments

- x:

  An `nm_model_suggestion` object.

- analysis_dir:

  Analysis directory where the control streams are written and run. It
  must contain the dataset, or `data` must point to it.

- nmfe_path:

  Path to the nmfe launcher (see
  [`nm_run()`](https://nskoch.github.io/nmpilot/reference/nm_run.md)).

- candidates:

  Candidate numbers (default: the first 3).

- prefix:

  Prefix of the control stream names (`<prefix><n>.ctl`).

- data:

  Dataset file name written in `$DATA` (default: the file analysed by
  [`nm_suggest_model()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_model.md)),
  relative to `analysis_dir` or absolute.

- run:

  If `FALSE`, only the control streams are written.

- overwrite:

  Replace existing control streams and adapted datasets, and run the
  candidates again even when a previous run of the same name has
  terminated (otherwise
  [`nm_run()`](https://nskoch.github.io/nmpilot/reference/nm_run.md)
  skips it and its old results are read).

- verbose:

  Print progress messages.

- n_cores:

  Number of candidates run at the same time, passed to
  [`nm_run()`](https://nskoch.github.io/nmpilot/reference/nm_run.md)
  (default 1).

- ...:

  Further arguments passed to
  [`nm_write_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_write_ctl.md)
  (e.g. `method`).

## Value

An object of class `nm_suggest_run`: a data.frame with one row per
candidate (`candidate`, `model`, `score`, `rank_score`, `ctl`, `data`
(dataset used), `ofv`, `n_par`, `n_obs`, `bic`, `rank_bic`, `min_ok`),
with the Spearman correlation between score and BIC ranks in attribute
`spearman` and the
[`nm_run()`](https://nskoch.github.io/nmpilot/reference/nm_run.md) log
in attribute `log`.

## Details

Each candidate gets a dataset that matches its compartment numbering:
doses in the depot (CMT 1) and observations in the central compartment
(CMT 2) for first-order and transit absorption (ADVAN2/4/12 or ODE);
doses and observations in the central compartment (CMT 1) for IV and
zero-order input (ADVAN1/3/11), with `RATE = -2` on the doses for a
modelled zero-order duration (`D1`). When the original dataset does not
match, a copy `<data>_<prefix><n>.csv` is written in `analysis_dir`
(records copied verbatim, only CMT and RATE changed); the original file
is never modified. Datasets with several observed compartments (parent
and metabolite, PK and PD) or doses by several routes are not handled:
the function stops.

## See also

[`nm_suggest_model()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_model.md),
[`nm_write_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_write_ctl.md),
[`nm_run()`](https://nskoch.github.io/nmpilot/reference/nm_run.md)

## Examples

``` r
if (FALSE) { # \dontrun{
s <- nm_suggest_model("C:/nm76g64/run/etude/data.csv")
r <- nm_suggest_run(s, analysis_dir = "C:/nm76g64/run/etude",
                    nmfe_path = "C:/nm76g64/run/nmfe76.bat", candidates = 1:4)
r
} # }
```
