# Export a data.frame to CSV (NONMEM convention: na = ".")

Export a data.frame to CSV (NONMEM convention: na = ".")

## Usage

``` r
nm_export_data(data, file)
```

## Arguments

- data:

  A data.frame, typically from
  [`nm_read_csv()`](https://nskoch.github.io/nmpilot/reference/nm_read_csv.md).

- file:

  Output path (.csv). The directory is created if needed.

## Value

The path of the written file (invisibly).
