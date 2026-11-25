# Goodness-of-fit diagnostic plots from a NONMEM \$TABLE file

Reads a NONMEM output table (space/tab-delimited, possibly with repeated
headers and multiple "TABLE NO." blocks), keeps the observation rows,
and builds the standard goodness-of-fit diagnostics: observed vs
predictions and residual plots (CWRES / NPDE / IWRES). Plots are only
built when the required columns are present. Optionally writes all plots
to a Word (.docx) file in the analysis directory.

## Usage

``` r
nm_gof(
  file,
  pred = c("PRED", "CPRED"),
  time = c("TIME", "TAD"),
  by = NULL,
  loess = TRUE,
  reg = TRUE,
  se = FALSE,
  alpha = 0.5,
  log_pred = FALSE,
  stem = NULL,
  save_docx = FALSE,
  output_file = NULL,
  width = 6.5,
  height = 5,
  dpi = 150,
  ...
)
```

## Arguments

- file:

  Path to a NONMEM table file, OR a run directory (.reg / .nm7 / any
  folder) inside which a NONMEM table is auto-detected. Since the table
  is not necessarily named "sdtab", detection is by content ("TABLE NO."
  header), optionally filtered by `stem`.

- pred:

  Preference for the population prediction column used on x-axes. The
  first available name is used. Default `c("PRED","CPRED")`. Pass e.g.
  `pred = "CPRED"` to force one.

- time:

  Preference for the time column used on residual-vs-time plots. Default
  `c("TIME","TAD")`.

- by:

  Optional character vector of column names for faceting (`facet_wrap`),
  e.g. `c("SEX","CMT")`.

- loess:

  Add a loess trend (default TRUE).

- reg:

  Add a linear regression (default TRUE).

- se:

  Show confidence band around the linear regression (default FALSE).

- alpha:

  Point transparency (default 0.5).

- log_pred:

  Log10 axes for the DV-vs-prediction plots (default FALSE).

- stem:

  Optional filename stem to disambiguate auto-detection in a run
  directory (e.g. "sdtab", "patab"). Not required to be "sdtab".

- save_docx:

  If TRUE, write all plots to a .docx in the analysis directory (default
  FALSE).

- output_file:

  Output .docx path. If NULL and `save_docx=TRUE`, a name is derived
  from the table file, written next to it. (Formerly `docx_file`, still
  accepted with a warning.)

- width, height, dpi:

  Image size (inches) and resolution for the docx export.

- ...:

  Former argument name, accepted with a warning until version 1.0.

## Value

An object of class `"nm_gof"`: a list with `plots` (named list of ggplot
objects), `data` (the filtered observation rows), `file`, `columns`
(resolved column mapping), `skipped` (plots that could not be built and
why), and `by`.

## Examples

``` r
if (FALSE) { # \dontrun{
g <- nm_gof("C:/nm74g64/run/lopi_04.reg/lopi_04_03_04.fit")
g                       # summary
plot(g)                 # all panels
plot(g, type = "dv_pred")
nm_gof(".../run.reg", pred = "CPRED", by = "SEX", save_docx = TRUE)
} # }
```
