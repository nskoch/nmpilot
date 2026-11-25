# Read and normalise a NONMEM-like dataset

Read and normalise a NONMEM-like dataset

## Usage

``` r
nm_read_csv(
  file,
  sep = ",",
  na = c("", "NA", "."),
  time_unit = c("h", "d"),
  lloq = NULL,
  blq_col = NULL,
  date_format = NULL,
  dv_col = NULL,
  verbose = TRUE
)
```

## Arguments

- file:

  Path to the CSV file.

- sep:

  Column separator (default `","`).

- na:

  Strings to be treated as missing values.

- time_unit:

  Output unit for TIME: `"h"` (default) or `"d"`.

- lloq:

  LLOQ value. If supplied and no BLQ/CENS column exists, records with DV
  \< lloq are flagged as BLQ.

- blq_col:

  Name of the BLQ flag column (1 = BLQ, 0 = quantified). By default,
  `"BLQ"` or `"CENS"` is auto-detected.

- date_format:

  Date format forced in
  [`lubridate::parse_date_time()`](https://lubridate.tidyverse.org/reference/parse_date_time.html)
  (e.g. `"dmy"` for DD/MM/YYYY). By default, chosen from the name of the
  DATE/DAT1-3 column.

- dv_col:

  Name of the DV column if non-standard (e.g. `"DV_nmol/L"`). If `NULL`
  (default), a column whose name starts with "DV" is auto-detected.
  Required when several DV columns exist.

- verbose:

  Print the identification report (default `TRUE`).

## Value

A list:

- `data`: data.frame with relative TIME, BLQ (logical), and DV_norm when
  AMT is available. Reserved columns are renamed in upper case.

- `mapping`: reserved labels -\> original column names.

- `lloq`: LLOQ value used (`NULL` if not supplied).

## Examples

``` r
if (FALSE) { # \dontrun{
res <- nm_read_csv("data.csv", lloq = 0.5)
head(res$data)
} # }
```
