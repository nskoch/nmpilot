# Generate VPC control streams and launch NONMEM runs

From an estimation run, generates two NONMEM control streams required
for a VPC:

- `<run>_orig.ctl`: re-runs the model at `MAXEVAL=0` with final
  estimates fixed, producing observed predictions (`_orig.fit`).

- `<run>_sim.ctl`: simulates `n_sim` replicates using
  `$SIMULATION ONLYSIM SUBPROBLEMS=n_sim`, producing `_sim.fit`.

Both control streams are then launched via
[`nm_run`](https://nskoch.github.io/nmpilot/reference/nm_run.md).

## Usage

``` r
nm_vpc_ctl(
  run,
  analysis_dir = NULL,
  nmfe_path,
  n_sim = 500,
  seed = 12345,
  skip_terminated = TRUE,
  run_subdir = NULL,
  verbose = TRUE,
  n_cores = 1L,
  ...
)
```

## Arguments

- run:

  Path of the estimation run: its control stream
  (`"C:/project/run001.ctl"`), its run directory
  (`"C:/project/run001.reg"`) or its stem (`"C:/project/run001"`). Its
  `.ctl`, `.ext` and `.lst` files are looked for in `run001.reg/`,
  `run001.nm7/`, `run001/`, then directly in the folder of the run (flat
  layout); this folder is the analysis directory, where the VPC control
  streams are written.

- analysis_dir:

  Deprecated (former form `run_name` + `analysis_dir`): when given,
  `run` is read as the run name in this directory, with a warning.

- nmfe_path:

  Full path to the nmfe executable (e.g. `"C:/nm76g64/run/nmfe76.bat"`
  on Windows, `"/opt/nm760/run/nmfe76"` on Linux or macOS).

- n_sim:

  Number of simulation replicates (default: 500).

- seed:

  Random seed for `$SIMULATION` (default: 12345).

- skip_terminated:

  Passed to
  [`nm_run`](https://nskoch.github.io/nmpilot/reference/nm_run.md). If
  `TRUE` (default), already completed VPC runs are skipped.

- run_subdir:

  Layout of the two VPC runs, passed to
  [`nm_run`](https://nskoch.github.io/nmpilot/reference/nm_run.md). By
  default (`NULL`) the layout of the estimation run is reused: `""`
  (flat) if its files were found directly in the analysis directory,
  `".reg"` otherwise.

- verbose:

  Passed to
  [`nm_run`](https://nskoch.github.io/nmpilot/reference/nm_run.md) and
  used locally for messages.

- n_cores:

  Passed to
  [`nm_run`](https://nskoch.github.io/nmpilot/reference/nm_run.md): 2
  runs the `_orig` and `_sim` control streams at the same time (default
  1).

- ...:

  Former argument name, accepted with a warning until version 1.0:
  `run_name` (now `run`).

## Value

A `data.frame` run log (from
[`nm_run`](https://nskoch.github.io/nmpilot/reference/nm_run.md)) with
one row per VPC run (`_orig` and `_sim`).

## Details

The `$PK` block of the source control stream must contain a line of the
form `NSIM = <value>` which will be replaced by `NSIM = 1` in
`_orig.ctl` and `NSIM = IREP` in `_sim.ctl`.

Final parameter estimates are read from the `.ext` file of the
estimation run (row `-1000000000`) and written as `FIX` values in both
generated control streams.

The generated `.ctl` files are written to the analysis directory (not
inside a run subdirectory). They can be inspected or edited before
running if needed. Wings for NONMEM is not required.

## See also

[`nm_run`](https://nskoch.github.io/nmpilot/reference/nm_run.md)

## Examples

``` r
if (FALSE) { # \dontrun{
log <- nm_vpc_ctl(
  run          = "C:/analyses/my_pk_project/run001.ctl",
  nmfe_path    = "C:/nm76g64/run/nmfe76.bat",
  n_sim        = 500,
  seed         = 12345
)
print(log)
# Then run the VPC:
# fig <- nm_vpc("C:/analyses/my_pk_project/run001.ctl")
} # }
```
