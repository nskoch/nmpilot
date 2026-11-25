# Extract the results of several NONMEM runs

Extract the results of several NONMEM runs

## Usage

``` r
nm_extract_runs(path, runs = NULL, verbose = TRUE)
```

## Arguments

- path:

  Analysis directory, or the directory of a single run.

- runs:

  Optional names of the runs to keep.

- verbose:

  Print progress.

## Value

An object of class `nm_runs` with:

- `$runs`: one row per run

- `$params`: one row per parameter (`same = TRUE` for the OMEGA elements
  of an `$OMEGA ... SAME` record)

- `$mixture`: for `$MIX` runs with a `.phm` file, one row per subject
  and subpopulation (`PMIX`, `assigned` = most probable subpopulation)

- `$path`: the analysis directory

- `$created`: time of extraction

## Details

When the control stream of a run is found, parameters are named from the
model code with
[`nm_identify_parameters()`](https://nskoch.github.io/nmpilot/reference/nm_identify_parameters.md):
typical-value THETAs and ETAs after the PK parameter they act on (`CL`,
`V`...), covariate-effect and residual-error THETAs after their `$THETA`
comment label when it is a short name (`CL_WT`, `PROP`, `ADD`...). Other
THETAs keep their NONMEM name (`THETA6`...).

## Several estimation steps, several problems, failed runs

With several `$EST` records, the OFV, the minimisation status, the
shrinkage and the covariance step are those of the last estimation step
(`est_methods` lists the steps, e.g.
`"ITS > SAEM > IMP(eval) > FOCEI"`). With several `$PROB`, the last
problem that has an estimation is read (`n_problems`). The estimates are
read in the `.ext`; when it is empty, has no final row, or reports an
OFV different from the `.lst` (an `.ext` left by another run), they are
read in the `.lst` (`params_source`) and a warning is given. `run_state`
is `"completed"`, `"interrupted"` (the output stops before the end of
the last step), `"terminated"` (PROGRAM TERMINATED), `"nmtran_error"` or
`"simulation"` (simulation-only run, `$SIMULATION ONLYSIM`, finished: no
OFV, `min_ok` is `NA`; subproblems are not counted in `n_problems`).
`n_omega` and `n_sigma` count the estimated elements only (not FIX, not
SAME, not the zeros outside the blocks); `n_eta_same` counts the ETAs
declared `SAME` (inter-occasion variability); ETAs fixed to 0 are left
out of `max_eta_shrink`. For `$MIX` runs, `n_subpop` is the number of
subpopulations, the mixing-proportion THETAs are named `MIX_P` when they
have no other name, the shrinkage of an ETA is its smallest value over
the subpopulations, and `mix_assign` counts the subjects assigned to
each subpopulation.
