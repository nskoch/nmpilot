# Nonparametric (case) bootstrap of a NONMEM model

Nonparametric (case) bootstrap of a NONMEM model

## Usage

``` r
nm_bootstrap(
  ctl,
  nmfe_path = NULL,
  n = 500L,
  ci = 0.95,
  seed = NULL,
  id_col = "ID",
  strata = NULL,
  comment_char = NULL,
  sep = ",",
  output_dir = NULL,
  keep_datasets = FALSE,
  keep_outputs = FALSE,
  strip_cov = TRUE,
  strip_tables = TRUE,
  add_warnings = TRUE,
  warnings_line = "$WARN WARN=NONE",
  include_rounding = FALSE,
  resume = TRUE,
  run = TRUE,
  verbose = TRUE,
  n_cores = 1L,
  ...
)
```

## Arguments

- ctl:

  Path to the base control stream (.ctl).

- nmfe_path:

  Path to the NONMEM launcher: nmfe\*.bat on Windows (e.g.
  "C:/nm74g64/run/nmfe74.bat"), the nmfe script on Linux or macOS (e.g.
  "/opt/nm760/run/nmfe76").

- n:

  Number of bootstrap replicates (default 500).

- ci:

  Confidence level(s) for the percentile CI (default 0.95; e.g.
  c(0.95,0.90)).

- seed:

  RNG seed for reproducibility (default: a random seed, stored in the
  result).

- id_col:

  Label of the subject ID column, as it appears in \$INPUT (default
  "ID").

- strata:

  Optional \$INPUT column to stratify the resampling on (constant within
  ID).

- comment_char:

  Data-file comment/ignore character. NULL = auto-derive from the
  single-char IGNORE in \$DATA ("@", "#", a letter...); "" = none.

- sep:

  Data-file field separator (default ",").

- output_dir:

  Working directory (default "\<ctl_dir\>/\_bootstrap"). Formerly
  `run_dir`, still accepted with a warning.

- keep_datasets:

  Keep each replicate data file? (default FALSE = space-optimal).

- keep_outputs:

  Keep each replicate's NONMEM outputs? (default FALSE; .ext is read
  then removed).

- strip_cov:

  Remove \$COVARIANCE (default TRUE; useless for a bootstrap).

- strip_tables:

  Remove \$TABLE records (default TRUE; fixed FILE= names would
  collide). Formerly `strip_table`, which still works.

- add_warnings:

  Insert a warning-suppression record (default TRUE).

- warnings_line:

  Exact record inserted (default "\$WARN WARN=NONE"; canonical form
  "\$WARNINGS").

- include_rounding:

  Count "minimisation successful, but rounding errors" runs as
  successful (default FALSE).

- resume:

  Skip replicates already in bootstrap_estimates.csv (default TRUE).

- run:

  FALSE = only build the optimised ctl + resample plan (dry run, no
  NONMEM).

- verbose:

  Progress messages (default TRUE).

- n_cores:

  Number of replicates run at the same time (default 1). A whole number,
  `"auto"` (physical cores minus one, see
  [`nm_cores()`](https://nskoch.github.io/nmpilot/reference/nm_cores.md))
  or `"ask"` (interactive sessions: the cores are shown and the number
  is asked once per R session, then reused;
  `options(nmpilot.n_cores = NULL)` asks again). Each run is then
  started in the background in its own directory; if R is interrupted,
  the runs already started go on in the background and a later call with
  the same arguments picks up their results. In parallel, each run uses
  a working subdirectory `slot<k>/` of `output_dir`; the estimates are
  appended to `bootstrap_estimates.csv` as the replicates end (the file
  is sorted by replicate when the result is assembled).

- ...:

  Former argument names, accepted with a warning until version 1.0:
  `run_dir` (now `output_dir`). The former name `strip_table` still
  works (R matches it to `strip_tables`).

## Value

An object of class "nm_bootstrap" (see print()/plot()).
