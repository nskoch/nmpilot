# Convert a Wings for NONMEM extended control stream to standard NONMEM

Wings for NONMEM (WFN) supports an extended control stream format where
parameter names (e.g. `POP_CL`, `PPV_CL`, `RUV_PROP`) are used directly
in model records (`$PK`, `$PRED`, `$ERROR`, `$DES`) instead of
`THETA(n)`, `ETA(n)`, `ERR(n)`. This function parses the parameter names
from `$THETA`, `$OMEGA` and `$SIGMA` comment fields and substitutes them
back into standard NONMEM syntax, producing a control stream compatible
with direct `nmfe` execution.

## Usage

``` r
nm_convert_ctl(ctl, output_file = NULL, overwrite = FALSE, verbose = TRUE, ...)
```

## Arguments

- ctl:

  Path to the Wings extended control stream (`.ctl`). (Formerly
  `ctl_file`, still accepted with a warning.)

- output_file:

  Output path for the converted control stream. If `NULL` (default), the
  output file is named `<run>_converted.ctl` in the same directory as
  `ctl`. (Formerly `out_file`, still accepted with a warning.)

- overwrite:

  Logical. If `TRUE`, overwrite `output_file` if it already exists
  (default `FALSE`).

- verbose:

  Logical. If `TRUE` (default), print a summary of substitutions made.

- ...:

  Former argument name, accepted with a warning until version 1.0:
  `ctl_file` (now `ctl`).

## Value

Invisibly returns the path to the converted control stream.

## Details

The conversion follows the WFN specification:

- Parameter names are extracted from the first non-numeric word after
  the `;` comment delimiter in each `$THETA`, `$OMEGA` and `$SIGMA`
  line.

- Names are substituted by `THETA(n)`, `ETA(n)` and `ERR(n)` in model
  records (`$PK`, `$PRED`, `$ERROR`, `$DES`, `$AES`).

- Once a name appears on the left-hand side of an assignment, it is no
  longer substituted on subsequent right-hand sides (WFN rule).

- `$OMEGA BLOCK` structures are handled correctly.

## See also

[`nm_run`](https://nskoch.github.io/nmpilot/reference/nm_run.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert a WFN extended control stream
out <- nm_convert_ctl("C:/analyses/run001_wfn.ctl")

# Then run with nmpilot
log <- nm_run(
  ctl          = out,
  nmfe_path    = "C:/nm76g64/run/nmfe76.bat",
  analysis_dir = "C:/analyses"
)
} # }
```
