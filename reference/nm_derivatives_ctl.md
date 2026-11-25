# Control stream of derivatives for the linearization of a base model

Writes a copy of a base control stream that evaluates the model once at
the final estimates (`MAXEVAL=0`, all parameters `FIX`) and tabulates,
for every record, the individual prediction and the derivatives of `Y`
with respect to the ETAs (`G_ETAk`), the EPSs (`H_EPSj`) and the cross
derivatives (`D_EPSj_ETAk`, interaction). The derivatives are read
through verbatim code (`G(k,1)`, `HH(j,1)`, `HH(j,k+1)` in `$ERROR`; `G`
and `H` in `$PRED`) placed at the end of the block (`"LAST`), as in the
PsN `linearize` tool. The table is the input of level 2 of
[`nm_optimize_base()`](https://nskoch.github.io/nmpilot/reference/nm_optimize_base.md).

## Usage

``` r
nm_derivatives_ctl(run, output_dir = NULL, name = NULL, extra = NULL, ...)
```

## Arguments

- run:

  Run directory of the base model (e.g. `"run1.reg"`) or path prefix
  (`"dir/run1"`) whose `.ctl` and `.ext` files are used.

- output_dir:

  Directory of the new control stream (formerly `out_dir`, still
  accepted with a warning). Default: the analysis directory (parent of a
  `.reg`/`.nm7` run directory), so that the `$DATA` path keeps the
  convention of
  [`nm_run()`](https://nskoch.github.io/nmpilot/reference/nm_run.md).

- name:

  Name of the new run. Default `"<run>_deriv"`.

- extra:

  Additional data items to tabulate (e.g. an occasion column `"OCC"` for
  inter-occasion variability). Items absent from `$INPUT` are dropped
  with a warning.

- ...:

  Former argument name, accepted with a warning until version 1.0:
  `out_dir` (now `output_dir`).

## Value

Invisibly, the path of the control stream written. The table file is
`<name>.dtab`.

## Details

Not testable without NONMEM: the verbatim code follows the NONMEM help
(ERROR: `G(i,1)` is the partial derivative of F with respect to ETA(i);
`HH(j,1)` with respect to EPS(j); `HH(j,i+1)` the second derivative with
respect to EPS(j) and ETA(i)). Models with `F_FLAG` (M3 likelihood) or
simulation-only records are not supported.

## See also

[`nm_optimize_base()`](https://nskoch.github.io/nmpilot/reference/nm_optimize_base.md)
