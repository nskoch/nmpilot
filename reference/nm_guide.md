# Interactive task guide

Builds and opens an HTML page that lists the steps of a population PK
analysis with nmpilot. Buttons group the tasks by stage (data
preparation, exploration, control streams, running, results,
diagnostics, model selection, validation, reporting); selecting a task
shows the matching nmpilot function(s), a short description and an
example of R code that can be copied. Tasks whose function is not
available yet are shown as "planned".

## Usage

``` r
nm_guide(output_file = NULL, open = interactive(), ...)
```

## Arguments

- output_file:

  Output path of the HTML page. Default: `nmpilot_guide.html` in the
  temporary directory of the R session. (Formerly `file`, still accepted
  with a warning.)

- open:

  Open the page in the web browser (default: in interactive sessions
  only).

- ...:

  Former argument name, accepted with a warning until version 1.0.

## Value

The path of the HTML file, invisibly.

## Details

The page is generated from the file `guide/guide_tasks.csv` installed
with the package (one row per task: stage, task, functions, description,
code, status), so that it follows the package. It is a single
self-contained file: no internet connection is needed. The examples use
the example project of the package, copied to a temporary folder (see
the "Setup" block of the page); `nmfe` must be set to the path of your
NONMEM launcher.

## Examples

``` r
p <- nm_guide(open = FALSE)
file.exists(p)
#> [1] TRUE
if (FALSE) { # \dontrun{
nm_guide()     # opens the guide in the browser
} # }
```
