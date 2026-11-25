# Check that a dataset is compatible with NONMEM

Reads a dataset as NONMEM would (every value as text, `.` or an empty
cell as missing) and reports the problems that make NM-TRAN refuse the
data or NONMEM misread them. The data are never modified.

## Usage

``` r
nm_check_dataset(data, cols = NULL, sep = ",", outlier_k = 4, verbose = TRUE)
```

## Arguments

- data:

  Path of a CSV file, or a data frame (for example the `data` element of
  [`nm_read_csv()`](https://nskoch.github.io/nmpilot/reference/nm_read_csv.md)).

- cols:

  Optional named character vector giving the column of a standard item
  when its name differs, e.g. `c(ID = "SUBJ", DV = "CONC")`. Names are
  otherwise matched without regard to case.

- sep:

  Column separator of the CSV file (default `","`).

- outlier_k:

  Number of median absolute deviations (log DV) beyond which an
  observation is reported as an outlier (default 4).

- verbose:

  Print the report (default `TRUE`).

## Value

An object of class `nm_dataset_check`: a list with `findings` (one row
per problem found: `check`, `severity` "error" or "warning", `column`,
`n_rows`, `rows` (row numbers in the data, header excluded), `ids` and
`message`), `n_rows`, `n_subjects`, `columns` (standard item -\> column
used) and `ok` (`TRUE` when no error was found). Printed (invisibly
returned) when `verbose = TRUE`.

## Details

Checks:

- **columns**: ID, TIME and DV present (error); AMT absent (warning: no
  dose record);

- **missing ID or TIME** values (error);

- **non-numeric values** in a column (error; NONMEM reads numbers only,
  unless the column is dropped in `$INPUT`);

- **times**: TIME decreasing within a subject (error), except at a reset
  record (EVID 3 or 4);

- **EVID / MDV / AMT**: EVID outside 0 to 4 (error); dose record (EVID 1
  or 4) without a positive AMT (error); observation (EVID 0) with a
  positive AMT (error); MDV 0 with a missing DV (error); MDV 0 on a
  record that is not an observation (error); observation with MDV 1 but
  a DV (warning: the value is ignored); without EVID, a record with AMT
  \> 0 is a dose and MDV must be 1;

- **duplicates**: identical records (warning); two observations of the
  same subject, time and compartment with different DV (warning);

- **outliers** in the observed DV: negative values (warning) and values
  far from the observations taken at the same time after the last dose
  (log scale, more than `outlier_k` median absolute deviations from
  their median; per compartment when CMT exists; groups of at least 10
  observations) (warning).

The header line may start with `#` or `@` (`#ID,TIME,...`), as read with
`IGNORE=#` or `IGNORE=@` in `$DATA`.

## Examples

``` r
ex <- system.file("extdata", "example", package = "nmpilot")
chk <- nm_check_dataset(file.path(ex, "example_data.csv"))
#> <nm_dataset_check> example_data.csv 
#>   600 records, 60 subjects: 0 error(s), 1 warning(s) -> readable by NONMEM
#>   [WARNING] Observed DV far from the others taken at the same time after dose (more than 4 median absolute deviations from their median, log scale, per CMT). (4 record(s): rows 59, 98-100; ID 6, 10)
chk$ok
#> [1] TRUE
chk$findings
#>      check severity column n_rows         rows   ids
#> 1 outliers  warning     DV      4 59, 98, .... 6, 10
#>                                                                                                                                             message
#> 1 Observed DV far from the others taken at the same time after dose (more than 4 median absolute deviations from their median, log scale, per CMT).
```
