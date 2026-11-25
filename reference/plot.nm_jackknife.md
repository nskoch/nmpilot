# Plots for `nm_jackknife` objects

Plots for `nm_jackknife` objects

## Usage

``` r
# S3 method for class 'nm_jackknife'
plot(x, type = c("bar", "rank", "scatter"), ...)
```

## Arguments

- x:

  An `nm_jackknife` object.

- type:

  `"bar"` (default): \\\|\Delta OFV\|\\ per subject; `"rank"`: ranked
  \\\|\Delta OFV\|\\; `"scatter"`: \\\|\Delta OFV\|\\ against the
  relative THETA shift norm.

- ...:

  Ignored.

## Value

A ggplot object.
