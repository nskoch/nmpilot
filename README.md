# nmtools

<!-- badges: start -->
<!-- badges: end -->

**nmtools** is an R package that drives a population pharmacokinetic (popPK)
analysis with NONMEM from R, from the data set to the final model diagnostics.
It calls the NONMEM `nmfe` launcher directly and does not need PsN or
Wings for NONMEM, while still reading the Wings conventions (`.reg` run
directories, `.smr` files, extended control streams with parameter labels).

> **Status:** version 0.1, under active development. The interface may still
> change. Bug reports and suggestions are welcome in the
> [issues](https://github.com/nskoch/nmtools/issues).

## Installation

```r
# install.packages("remotes")
remotes::install_github("nskoch/nmtools")
```

Running models requires a working NONMEM installation. The package is
developed and used on Windows with NONMEM 7.4 and 7.6 (`nmfe74.bat`,
`nmfe76.bat`). Reading outputs, diagnostics and plots work without NONMEM.

## What the package covers

| Step | Functions |
|---|---|
| Read and explore the data set | `nm_read_csv()`, `nm_data_previsual()`, `nm_plot_dv()`, `nm_summary()`, `nm_add_rate()`, `nm_export_data()` |
| Run NONMEM | `nm_run()`, `nm_convert_ctl()` (Wings extended control stream to standard NONMEM) |
| Read the results of one or several runs | `nm_extract_runs()`, `nm_filter_runs()`, `nm_compare_runs()`, `nm_extract_run_table()` (Word table) |
| Identify what each THETA / ETA means | `nm_identify_parameters()`, `nm_identify_directory()` |
| Goodness-of-fit diagnostics | `nm_gof()` |
| Compare two runs | `nm_compare()` (paired CWRES / NPDE, individual OFV) |
| VPC and pc-VPC | `nm_vpc_ctl()`, `nm_vpc()` |
| Parameter uncertainty | `nm_profile()` (likelihood profile), `nm_bootstrap()` (case bootstrap) |
| Export figures | `nm_export_plot()` |

## A typical workflow

```r
library(nmtools)

# 1. Data set: read, summarise, plot
dat <- nm_read_csv("C:/nm76g64/run/study/data.csv")
nm_summary(dat$data)
nm_plot_dv(dat$data, scales = "both")

# 2. Run one or several models. Each run goes to <run>.reg/ in the analysis
#    directory; $DATA paths are read relative to that run directory
#    (e.g. $DATA ../data.csv).
log <- nm_run(
  nmfe_path    = "C:/nm76g64/run/nmfe76.bat",
  ctl_files    = c("run1.ctl", "run2.ctl"),
  analysis_dir = "C:/nm76g64/run/study"
)

# 3. Results of all the runs of the project
runs <- nm_extract_runs("C:/nm76g64/run/study")
runs                          # OFV, dOFV, minimisation and covariance status
nm_compare_runs(runs)         # parameters side by side
nm_extract_run_table(run_path = "C:/nm76g64/run/study/run2.reg")   # Word table

# 4. Diagnostics, and comparison of two models
gof <- nm_gof("C:/nm76g64/run/study/run2.reg/sdtab2")
plot(gof)
cmp <- nm_compare("C:/nm76g64/run/study/run1.reg",
                  "C:/nm76g64/run/study/run2.reg",
                  labels = c("Base", "+ WT on CL"))
plot(cmp)

# 5. VPC / pc-VPC
nm_vpc_ctl("run2", "C:/nm76g64/run/study",
           nmfe_path = "C:/nm76g64/run/nmfe76.bat", n_sim = 500)
nm_vpc("run2", "C:/nm76g64/run/study",
       col_names = c("NSIM", "ID", "TIME", "AMT", "MDV", "DV", "PRED", "IPRED"),
       type = 2)

# 6. Uncertainty
prof <- nm_profile("C:/nm76g64/run/nmfe76.bat",
                   "C:/nm76g64/run/study/run2.reg/run2.ctl",
                   param = "CL", from = 2, to = 6, by = 0.25)
plot(prof, type = "dofv")
boot <- nm_bootstrap("C:/nm76g64/run/study/run2.reg/run2.ctl",
                     nmfe_path = "C:/nm76g64/run/nmfe76.bat",
                     n = 500, seed = 1)
```

## Parameter names from the model code

`nm_identify_parameters()` reads `$PK` (or `$PRED`) and links every THETA and
ETA to the PK parameter it acts on, with a role: typical value, covariate
effect (power, exponential, linear, saturable or conditional), IIV scale or
residual error. `nm_extract_runs()` uses it so that results are labelled `CL`,
`V`, `KA`... rather than `THETA1`, `THETA2`. Wings control streams that use
labels (`CL = POPCL*EXP(eta_cl)`) are handled.

```r
x <- nm_identify_parameters("run2.ctl")
x$parameter_table
```

## License

MIT © Nicolas Simon
