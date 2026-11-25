# Export a ggplot to an image file

Single export entry point for all nmpilot plots: plot methods return the
ggplot object and leave saving to this function.

## Usage

``` r
nm_export_plot(plot, file, width = 9, height = 5, dpi = 300, verbose = TRUE)
```

## Arguments

- plot:

  A ggplot (or patchwork) object.

- file:

  Output path. The file type is taken from the extension:
  `.jpg`/`.jpeg`, `.png`, `.pdf` or `.tif`/`.tiff`. The directory is
  created if needed.

- width:

  Width in inches (default 9).

- height:

  Height in inches (default 5).

- dpi:

  Resolution in dots per inch (default 300; ignored for PDF).

- verbose:

  Print the path of the written file (default `TRUE`).

## Value

The path of the written file (invisibly), or `NULL` if `plot` is `NULL`.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- nm_plot_dv(nm_read_csv("data.csv")$data)
nm_export_plot(p, "figures/dv_vs_time.png")
} # }
```
