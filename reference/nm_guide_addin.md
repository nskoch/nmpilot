# RStudio addin: open the task guide

Called by the RStudio addin "nmpilot task guide" (menu Addins). Writes
the guide with
[`nm_guide()`](https://nskoch.github.io/nmpilot/reference/nm_guide.md)
in the temporary directory and opens it in the RStudio Viewer pane, or
in the web browser outside RStudio.

## Usage

``` r
nm_guide_addin()
```

## Value

The path of the HTML file, invisibly.
