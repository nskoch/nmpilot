# Plot a model suggestion

Plot a model suggestion

## Usage

``` r
# S3 method for class 'nm_model_suggestion'
plot(
  x,
  type = c("profiles", "support", "pooled", "design"),
  candidate = 1L,
  n_max = 12L,
  ...
)
```

## Arguments

- x:

  An `nm_model_suggestion` object.

- type:

  `"profiles"` (observed concentrations of the analysed intervals with
  the best fit of each, log scale), `"support"` (relative support of
  each component), `"pooled"` (binned pooled data with their fit) or
  `"design"` (expected RSE and eta-shrinkage of a candidate).

- candidate:

  Candidate shown by `type = "design"`.

- n_max:

  Maximum number of intervals shown by `type = "profiles"`.

- ...:

  Ignored.

## Value

A ggplot object.
