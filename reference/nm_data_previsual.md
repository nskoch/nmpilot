# Pre-visualise a NONMEM-like dataset

Pre-visualise a NONMEM-like dataset

## Usage

``` r
nm_data_previsual(
  file,
  sep = ",",
  na = c("", "NA", "."),
  time_unit = c("h", "d"),
  reconstruct_time = FALSE,
  clock_col = NULL,
  date_type = NULL,
  dv_col = NULL,
  strat = NULL,
  ci_probs = c(0.05, 0.95),
  dose_norm = FALSE,
  scales = c("linear", "log", "both"),
  cmt = NULL,
  xlim = NULL,
  blq_handling = c("show", "exclude", "half_lloq", "keep"),
  lloq = NULL,
  blq_col = NULL,
  export_jpg = FALSE,
  export_csv = FALSE,
  output_dir = NULL,
  jpg_width = 9,
  jpg_height = 5,
  jpg_dpi = 300,
  verbose = TRUE,
  ...
)
```

## Arguments

- file:

  Path to the CSV file (or an already loaded data.frame).

- sep:

  Column separator (default `","`).

- na:

  Strings to be treated as missing values.

- time_unit:

  Output unit for TIME: `"h"` (default) or `"d"`.

- reconstruct_time:

  If `TRUE`, rebuild a relative TIME (t0 = first dose per subject) from
  a clock-time column and, optionally, a DATE column. If `FALSE`
  (default), the TIME column is used as is (`as.numeric`). Detecting a
  DATE/DAT1/DAT2/DAT3 column NEVER triggers the reconstruction
  automatically.

- clock_col:

  Name of the column holding the clock time (`hh:mm:ss`, `hh:mm` or
  decimal hours), used when `reconstruct_time = TRUE`. `NULL` (default)
  = column `"TIME"`.

- date_type:

  How the DATE column is interpreted during the reconstruction (only if
  `reconstruct_time = TRUE`):

  - `NULL` (default): clock time only, anchored on a single day

  - `"continuous"`: DATE = integer day counter (1, 2, 3...)

  - `"DATE"`: calendar date month-day-year (mm-dd-yy)

  - `"DAT1"`: calendar date day-month-year (dd-mm-yy)

  - `"DAT2"`: calendar date year-month-day (yy-mm-dd)

  - `"DAT3"`: calendar date year-day-month (yy-dd-mm)

- dv_col:

  Name of the DV column to use (manual override).

- strat:

  Name of a stratification column for the plot (e.g. `"DOSE"`, `"OCC"`,
  `"STUDY"`). `NULL` = none.

- ci_probs:

  Two quantiles defining the interval band around the median (default
  `c(0.05, 0.95)`, i.e. a 90% interval).

- dose_norm:

  If `TRUE`, also compute DV_norm = DV / AMT (dose given to the subject)
  and produce the matching plot. Requires an identified AMT column.

- scales:

  `"linear"`, `"log"` or `"both"`: Y axes produced. `"both"` gives
  side-by-side linear and log panels (patchwork).

- cmt:

  CMT value(s) to keep for the plot. Requires an identified CMT column.
  `NULL` (default) = all observed compartments (CMT != 0 when EVID is
  available). E.g. `cmt = 2` for the central compartment.

- xlim:

  X-axis limits `c(min, max)`. `NULL` (default) = full data range.

- blq_handling:

  How BLQ values are handled:

  - `"show"` (default): BLQ drawn as grey crosses below a dashed line at
    LLOQ; excluded from the median and quantiles

  - `"exclude"`: BLQ removed from the plot and statistics

  - `"half_lloq"`: DV replaced by LLOQ/2 (requires `lloq`)

  - `"keep"`: nothing done, plotted as is

- lloq:

  Numeric LLOQ. Needed to identify BLQ records when no flag column is
  found, and for `"half_lloq"`. `NULL` = no threshold.

- blq_col:

  Name of the censoring flag column (1 = BLQ, 0 = quantified). If
  `NULL`, auto-detected among `"BLQ"`, `"BQL"`, `"LLQ"`, `"CENS"`.

- export_jpg:

  If `TRUE`, export the plots as JPG files in `output_dir`.

- export_csv:

  If `TRUE`, export the final data.frame (relative TIME, DV_norm) as CSV
  in `output_dir`.

- output_dir:

  Output directory (default: directory of the source file). Formerly
  `out_dir`, still accepted with a warning.

- jpg_width:

  Width in inches of the JPG export (default 9).

- jpg_height:

  Height in inches of the JPG export (default 5).

- jpg_dpi:

  Resolution in dpi (default 300).

- verbose:

  Print the identification report (default `TRUE`).

- ...:

  Former argument name, accepted with a warning until version 1.0:
  `out_dir` (now `output_dir`).

## Value

An object of class `nm_data_previsual` (a list, printed as a short
summary):

- `data`: data.frame with TIME (relative if rebuilt) and DV_norm

- `mapping`: reserved labels -\> original column names

- `plot_dv`: ggplot of DV vs TIME (according to `scales`)

- `plot_norm`: ggplot of DV/dose vs TIME (`NULL` if `dose_norm = FALSE`)

- `summary`: diagnostic data.frame (N subjects, N obs, time range, %
  BLQ...)

- `files`: paths of the written files (`NULL` if nothing was exported)

## Examples

``` r
if (FALSE) { # \dontrun{
# TIME already numeric/relative: used as is (default behaviour)
res <- nm_data_previsual("data.csv")

# Rebuild from a clock-time column TIME_R + DATE = day counter
res <- nm_data_previsual("data.csv", reconstruct_time = TRUE,
                         clock_col = "TIME_R", date_type = "continuous")

# Rebuild from a calendar DATE in day-month-year order (DAT1)
res <- nm_data_previsual("data.csv", reconstruct_time = TRUE,
                         clock_col = "TIME", date_type = "DAT1")

res <- nm_data_previsual("data.csv", strat = "DOSE", dose_norm = TRUE,
                         scales = "both", lloq = 0.5, export_jpg = TRUE)
print(res$plot_dv); print(res$summary)
} # }
```
