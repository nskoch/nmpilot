# Parameter table of one NONMEM run, exported to Word

Reads the `.ext` file of a run and writes a Word document with the
parameter estimates, 95% confidence intervals and RSE. When there is no
`.ext` file, or when it holds no estimates (empty file), the estimates
printed in the `.lst` file are used, as in
[`nm_extract_runs()`](https://nskoch.github.io/nmpilot/reference/nm_extract_runs.md);
the table footer then says so. ETA and EPS standard deviations are
derived from the variances (delta method for their SE); ETA correlations
are derived from the OMEGA matrix.

## Usage

``` r
nm_extract_run_table(
  run = NULL,
  run_name = NULL,
  run_path = NULL,
  output_name = "nonmem_results.docx",
  extensions = c("nm7", "reg"),
  ...
)
```

## Arguments

- run:

  Path of the run: run directory, control stream or stem.

- run_name:

  Deprecated: with `run_name`, `run` is read as the parent directory
  (former argument `run_rep`) and `run_name` as the run name. A full
  path given here by mistake is converted automatically.

- run_path:

  Deprecated: former name of `run`.

- output_name:

  Name of the Word file, written in the run directory (default
  `"nonmem_results.docx"`).

- extensions:

  Run directory extensions to try first, in order (default
  `c("nm7", "reg")`), before the folder `<run>` and the flat layout.

- ...:

  Former argument name, accepted with a warning until version 1.0:
  `run_rep` (with `run_name`).

## Value

Invisibly, an object of class `nm_run_table` (a list, printed as a short
summary) with `flextable`, `params` (parameter table), `ext` (parsed
.ext, or estimates read in the .lst, with `source = "lst"`),
`output_file`, `ext_file` (`NA` without .ext), `lst_file`, `run_dir` and
`extension`.

## Details

The run is given by its path, `run`: its run directory (`.nm7`, `.reg`),
its control stream, or its stem (e.g. `C:/project/run1` for
`C:/project/run1.ext`). The `.ext` file is looked for in `<run>.nm7/`,
`<run>.reg/`, a folder named `<run>`, then directly in the folder of the
run (flat layout, as produced by a plain nmfe call). Wings for NONMEM is
not required.

The former call forms, `run_rep` (parent directory) + `run_name`, and
`run_path`, are still accepted with a warning until version 1.0.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- nm_extract_run_table("C:/nm76g64/run/para/run_1_prop_kaBbloc.reg")
res <- nm_extract_run_table("C:/nm76g64/run/para/run_1_prop_kaBbloc.ctl")
res$params
} # }
```
