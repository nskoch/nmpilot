# Identify THETA/ETA roles for all control streams of a directory

Applies
[`nm_identify_parameters()`](https://nskoch.github.io/nmpilot/reference/nm_identify_parameters.md)
to every `.ctl` file found in `path` and stacks the parameter tables. A
file that cannot be parsed gives one row with its error message.

## Usage

``` r
nm_identify_directory(path, recursive = TRUE, likely_pk_only = TRUE)
```

## Arguments

- path:

  Directory to search.

- recursive:

  Search subdirectories (default `TRUE`).

- likely_pk_only:

  Passed to
  [`nm_identify_parameters()`](https://nskoch.github.io/nmpilot/reference/nm_identify_parameters.md).

## Value

A data.frame with a `file` column followed by the columns of
`parameter_table` (or an `error` column for files that failed).
