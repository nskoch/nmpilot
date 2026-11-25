# Add a RATE column to a NONMEM dataset

Reads the CSV verbatim (all columns as character) so that nothing is
reformatted, adds a RATE column based on the dosing records, and writes
a new file. The original columns (values, format, '#ID' header) are kept
unchanged.

## Usage

``` r
nm_add_rate(
  file,
  output_file = NULL,
  sep = ",",
  dose_value = -2,
  base_value = 0,
  after = "AMT",
  dose_cmt = NULL,
  cmt_col = NULL,
  evid_col = NULL,
  amt_col = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- file:

  Path to the original CSV.

- output_file:

  Name (or path) of the output CSV. `NULL` (default) -\>
  `"<base>_rate.csv"` in the same directory as `file`. If only a file
  name is given (no directory), it is placed next to the source file.
  (Formerly `out_file`, still accepted with a warning.)

- sep:

  Column separator (default `","`).

- dose_value:

  RATE value on the targeted dosing records: -2 (default, estimated
  duration) or -1 (estimated rate).

- base_value:

  RATE value on all other records (default 0).

- after:

  Name of the column after which RATE is inserted (default `"AMT"`).
  `NULL` = added as the last column.

- dose_cmt:

  CMT value(s) identifying the route(s) that receive `dose_value`.
  `NULL` (default): all dosing records receive `dose_value`
  (single-route dataset). E.g. `dose_cmt = 1` to affect only the oral
  depot, leaving IV boluses (other CMT) at `base_value`.

- cmt_col:

  Name of the CMT column (default: auto-detected). Required only when
  `dose_cmt` is supplied.

- evid_col:

  Name of the EVID column (default: auto-detected).

- amt_col:

  Name of the AMT column (default: auto-detected; used only when no EVID
  column is found).

- overwrite:

  If `TRUE`, overwrite an existing RATE column (default `FALSE`: the
  function stops if RATE already exists).

- verbose:

  Print a diagnostic report (default `TRUE`).

- ...:

  Former argument name, accepted with a warning until version 1.0.

## Value

The path of the written file (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
# Single route: RATE=-2 on all doses -> data_nmol_rate.csv
nm_add_rate("C:/nm74g64/run/Bienvenu/data_nmol.csv")

# Estimated rate rather than duration
nm_add_rate("data.csv", dose_value = -1)

# Multiple routes: -2 only on the oral depot (CMT=1),
# IV boluses (CMT=2) keep RATE=0
nm_add_rate("data.csv", dose_cmt = 1)

# Chosen output name (placed next to the source file)
nm_add_rate("data.csv", output_file = "data_dur.csv")
} # }
```
