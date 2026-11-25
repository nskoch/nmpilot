# Identify the semantic role of THETA and ETA in a NONMEM control stream

Parses `$THETA`, `$OMEGA`, `$PK` and `$ERROR` of a control stream (or
`$PRED` when there is no `$PK`) and links each THETA and ETA to the PK
parameter(s) it acts on, with a semantic role: `typical_value`,
`covariate_effect`, `iiv_scale`, `residual_error`, `iiv` or
`structural_other`. Executable equations are the primary source;
comments are only used as labels and consistency checks.

## Usage

``` r
nm_identify_parameters(ctl, likely_pk_only = TRUE, ...)
```

## Arguments

- ctl:

  Path to a NONMEM control stream (`.ctl`). (Formerly `file`, still
  accepted with a warning.)

- likely_pk_only:

  If `TRUE` (default), `mapping_detail` keeps only variables that look
  like PK parameters (CL, V, Q, KA, F, ALAG...). `mapping_detail_all`
  always keeps every variable.

- ...:

  Former argument name, accepted with a warning until version 1.0:
  `file` (now `ctl`).

## Value

An object of class `nm_parameter_map`, a list with, among others:
`parameter_table` (one row per PK parameter with its THETA/ETA by role;
`derived_from` names the parameters of the table that a derived
parameter is computed from, e.g. `"V"` for `S2 = V` or `"CL, V"` for
`KE = CL/V`, and is `NA` for the others), `mapping_detail` and
`mapping_detail_all` (one row per PK parameter / THETA or ETA
association; `depth` is the number of assignment steps between the
THETA/ETA and the parameter, 0 when it appears directly in its
definition), `residual_error`, `special_theta` (THETAs of model event
times `MTIME(n)` in `$PK` and of mixture probabilities `P(n)` in `$MIX`,
with their role `model_event_time` or `mixture_probability`; they are
not counted in `unmapped_theta`), `theta`, `omega`, `omega_elements`
(variances and BLOCK covariances), `unmapped_theta`, `unmapped_omega`,
`missing_theta_declaration`, `missing_omega_declaration` and
`label_substitutions` (the labels replaced and their
`THETA(n)`/`ETA(n)`/`ERR(n)` reference; empty for a standard control
stream).

## Details

Covariate effects are recognised in power, exponential, linear,
saturable and conditional forms (e.g. `(WT/70)**THETA(4)`,
`EXP(THETA(4)*WT)`, `THETA(1)*(1 + THETA(4)*(WT-70))`,
`C_RES/(THETA(4) + C_RES)`, `IF(SEX.EQ.1) TVCL = TVCL*THETA(4)`,
including `IF ... THEN ... ELSE ... ENDIF` blocks). A variable defined
several times keeps the THETA/ETA of all its definitions. When a
variable is only defined under conditions (e.g.
`IF(SEX.EQ.1) CL = THETA(1)` / `IF(SEX.EQ.0) CL = THETA(2)`), the THETA
of its first definition in the file is the reference (`typical_value`)
and the others are covariate effects. A covariate is any data item
declared in `$INPUT` that is not a NONMEM-defined item (`ID`, `TIME`,
`DV`, `AMT`, `MDV`, `EVID`, `CMT`...), whatever its name; items dropped
with `=DROP` are ignored. Without `$INPUT`, any symbol used in `$PK` but
not defined there is taken as a data item.

Wings for NONMEM extended control streams, where model records use the
labels of the `$THETA`, `$OMEGA` and `$SIGMA` comments (e.g.
`CL = POPCL*EXP(eta_cl)`), are handled: such labels are first replaced
by `THETA(n)`, `ETA(n)` or `ERR(n)`, as in
[`nm_convert_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_convert_ctl.md).
A label is replaced only if it is used but never assigned in a model
record and is not a `$INPUT` item.

## See also

[`nm_identify_directory()`](https://nskoch.github.io/nmpilot/reference/nm_identify_directory.md),
[`nm_extract_runs()`](https://nskoch.github.io/nmpilot/reference/nm_extract_runs.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- nm_identify_parameters("run.ctl")
x$parameter_table
x$mapping_detail
x$omega_elements
} # }
```
