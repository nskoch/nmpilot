# Visual Predictive Check (VPC or pc-VPC) from NONMEM simulation output

Reads the `_orig.fit` and `_sim.fit` files produced by
[`nm_vpc_ctl`](https://nskoch.github.io/nmpilot/reference/nm_vpc_ctl.md),
computes observed and simulated percentiles in time bins, and returns a
`ggplot2` VPC or pc-VPC plot.

## Usage

``` r
nm_vpc(
  run,
  analysis_dir = NULL,
  col_names,
  type = 1,
  pi_low = 0.05,
  pi_up = 0.95,
  ci_low = 0.025,
  ci_up = 0.975,
  bin_method = 4,
  TTIME_seq_manual = NULL,
  n_bins_quantile = 8,
  bin_rule = "FD",
  plot_title = NULL,
  plot_file = NULL,
  show_hist = TRUE,
  ...
)
```

## Arguments

- run:

  Path of the estimation run, as given to
  [`nm_vpc_ctl`](https://nskoch.github.io/nmpilot/reference/nm_vpc_ctl.md):
  its control stream (`"C:/project/run001.ctl"`), its run directory or
  its stem (`"C:/project/run001"`). The function looks for
  `run001_orig.fit` and `run001_sim.fit`, looked for in
  `run001_orig.reg/` (resp. `_sim.reg/`), `.nm7/`, a folder of the same
  name, then directly in the folder of the run.

- analysis_dir:

  Deprecated (former form `run_name` + `analysis_dir`): when given,
  `run` is read as the run name in this directory, with a warning.

- col_names:

  Character vector of column names matching the `$TABLE` output columns
  in the `.fit` files (e.g.
  `c("NSIM","ID","TIME","AMT","MDV","DV","PRED","IPRED")`).

- type:

  Integer. `1` = standard VPC (default), `2` = pc-VPC
  (prediction-corrected, Bergstrand 2011).

- pi_low:

  Lower bound of the prediction interval (default: 0.05).

- pi_up:

  Upper bound of the prediction interval (default: 0.95).

- ci_low:

  Lower bound of the confidence interval around simulated percentiles
  (default: 0.025).

- ci_up:

  Upper bound of the confidence interval around simulated percentiles
  (default: 0.975).

- bin_method:

  Integer. Binning method: `1` = manual, `2` = equiprobable quantiles,
  `3` = nominal time points, `4` = automatic (default: `4`).

- TTIME_seq_manual:

  Numeric vector of bin breaks. Required when `bin_method = 1`.

- n_bins_quantile:

  Number of bins for quantile method (`bin_method = 2`). Default: 8.

- bin_rule:

  Rule for automatic binning (`bin_method = 4`): `"Sturges"`, `"Scott"`,
  or `"FD"` (default).

- plot_title:

  Character. Title of the VPC plot. If `NULL` (default), a title is
  generated from the run name and `type`.

- plot_file:

  Character. If provided, the plot is saved to this path (e.g.
  `"vpc_run001.jpeg"`) through
  [`nm_export_plot`](https://nskoch.github.io/nmpilot/reference/nm_export_plot.md)
  (`.jpg`, `.png`, `.pdf` or `.tiff`).

- show_hist:

  Logical. If `TRUE` (default), prints a diagnostic histogram of
  observed TIME with bin boundaries before the VPC plot.

- ...:

  Former argument name, accepted with a warning until version 1.0:
  `run_name` (now `run`).

## Value

A `ggplot2` object, displayed when not assigned.

## Details

**pc-VPC correction (Bergstrand 2011, Eq. 1, lb = 0):** \$\$DV\_{pc} =
DV \times \tilde{PRED}\_{bin} / PRED\$\$ where \\\tilde{PRED}\_{bin}\\
(`MEDpred_ref`) is the median of all simulated PRED values pooled across
replicates within each bin. The same reference is applied to both
observed and simulated data.

**Required packages:** `ggplot2` and `dplyr` must be installed (listed
in `Imports`).

## See also

[`nm_vpc_ctl`](https://nskoch.github.io/nmpilot/reference/nm_vpc_ctl.md)
for generating the simulation runs.

## Examples

``` r
if (FALSE) { # \dontrun{
col_names <- c("NSIM", "ID", "TIME", "AMT", "MDV", "DV", "PRED", "IPRED")

# Standard VPC
p <- nm_vpc(
  run          = "C:/analyses/my_pk_project/run001.ctl",
  col_names    = col_names,
  type         = 1,
  bin_method   = 4
)

# pc-VPC with manual bins, saved to file
p <- nm_vpc(
  run              = "C:/analyses/my_pk_project/run001.ctl",
  col_names        = col_names,
  type             = 2,
  bin_method       = 1,
  TTIME_seq_manual = c(0, 1, 2, 4, 8, 12, 24),
  plot_file        = "pcvpc_run001.jpeg"
)
} # }
```
