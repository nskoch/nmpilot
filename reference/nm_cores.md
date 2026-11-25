# Number of processor cores available for parallel NONMEM runs

Reports the logical and physical cores of the computer, and the number
of parallel runs that `n_cores = "auto"` would use (physical cores minus
one, at least 1). Several runs share the memory and the disk: on a
laptop, or with large data sets, fewer runs than cores can be faster.

## Usage

``` r
nm_cores(verbose = TRUE)
```

## Arguments

- verbose:

  Print the result.

## Value

Invisibly, a list with `logical`, `physical` and `auto`.

## See also

[`nm_run()`](https://nskoch.github.io/nmpilot/reference/nm_run.md),
[`nm_bootstrap()`](https://nskoch.github.io/nmpilot/reference/nm_bootstrap.md),
[`nm_jackknife()`](https://nskoch.github.io/nmpilot/reference/nm_jackknife.md),
[`nm_profile()`](https://nskoch.github.io/nmpilot/reference/nm_profile.md)

## Examples

``` r
nm_cores()
#> Cores: 4 logical, 4 physical. n_cores = "auto" would use 3 parallel runs.
```
