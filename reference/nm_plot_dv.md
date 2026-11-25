# Plot DV vs TIME with median and interval band

Plot DV vs TIME with median and interval band

## Usage

``` r
nm_plot_dv(
  data,
  y = c("DV", "DV_norm"),
  time_var = c("TIME", "TAD"),
  strat = NULL,
  ci_probs = c(0.05, 0.95),
  cmt = NULL,
  xlim = NULL,
  scales = c("linear", "log", "both"),
  blq_handling = c("show", "exclude", "keep"),
  verbose = TRUE
)
```

## Arguments

- data:

  A data.frame from
  [`nm_read_csv()`](https://nskoch.github.io/nmpilot/reference/nm_read_csv.md).

- y:

  Y variable: `"DV"` (default) or `"DV_norm"`.

- time_var:

  X variable: `"TIME"` (default) or `"TAD"` (time after dose, computed
  by
  [`nm_read_csv()`](https://nskoch.github.io/nmpilot/reference/nm_read_csv.md)
  when EVID is present).

- strat:

  Name of a stratification column (e.g. `"DOSE"`).

- ci_probs:

  Quantiles of the interval band (default 90% interval,
  `c(0.05, 0.95)`).

- cmt:

  CMT value(s) to keep (`NULL` = all).

- xlim:

  X-axis limits `c(min, max)` (`NULL` = full range).

- scales:

  `"linear"`, `"log"` or `"both"` (side-by-side panels).

- blq_handling:

  BLQ handling: `"show"` (grey crosses, default), `"exclude"` or
  `"keep"`.

- verbose:

  Print diagnostics (default `TRUE`).

## Value

A ggplot (or a patchwork object if `scales = "both"`).

## Examples

``` r
if (FALSE) { # \dontrun{
parsed <- nm_read_csv("data.csv", lloq = 0.5)
p <- nm_plot_dv(parsed$data, strat = "DOSE", cmt = 2, scales = "both")
print(p)
} # }
```
