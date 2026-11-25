# Export a data.frame to CSV (NONMEM convention: na = ".")

Export a data.frame to CSV (NONMEM convention: na = ".")

## Usage

``` r
nm_export_data(data, output_file, ...)
```

## Arguments

- data:

  A data.frame, typically from
  [`nm_read_csv()`](https://nskoch.github.io/nmpilot/reference/nm_read_csv.md).

- output_file:

  Output path (.csv). The directory is created if needed. (Formerly
  `file`, still accepted with a warning.)

- ...:

  Former argument name, accepted with a warning until version 1.0.

## Value

The path of the written file (invisibly).
