# nmtools

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24244340.svg)](https://doi.org/10.5281/zenodo.24244340)
<!-- badges: end -->

**nmtools** is an R package that drives a population pharmacokinetic (popPK)
analysis with NONMEM from R, from the data set to the final model diagnostics.
It calls the NONMEM `nmfe` launcher directly and does not need PsN or
Wings for NONMEM. Wings is optional: runs can sit next to their control
streams or in run subdirectories, and the Wings conventions (`.reg` run
directories, `.smr` files, extended control streams with parameter labels)
are also read, so that projects started with Wings can be continued.

> **Status:** version 0.3, under active development. The interface may still
> change. Bug reports and suggestions are welcome in the
> [issues](https://github.com/nskoch/nmtools/issues).

## Installation

```r
# install.packages("remotes")
remotes::install_github("nskoch/nmtools")

# with the vignettes:
remotes::install_github("nskoch/nmtools", build_vignettes = TRUE)
```

Running models requires a working NONMEM installation. The package is
developed and used on Windows with NONMEM 7.4 and 7.6 (`nmfe74.bat`,
`nmfe76.bat`). Reading outputs, diagnostics and plots work without NONMEM.

## What the package covers

| Step | Functions |
|---|---|
| Read and explore the data set | `nm_read_csv()`, `nm_data_previsual()`, `nm_plot_dv()`, `nm_summary()`, `nm_add_rate()`, `nm_export_data()` |
| Write a control stream | `nm_write_ctl()` (compartmental PK model), `nm_ctl_mu()` (MU referencing), `nm_ctl_iov()` (inter-occasion variability), `nm_ctl_m3()` (M3 method for BQL data), `nm_ctl_check()` (static check) |
| Run NONMEM | `nm_run()`, `nm_convert_ctl()` (Wings extended control stream to standard NONMEM) |
| Read the results of one or several runs | `nm_extract_runs()`, `nm_filter_runs()`, `nm_compare_runs()`, `nm_extract_run_table()` (Word table) |
| Identify what each THETA / ETA means | `nm_identify_parameters()`, `nm_identify_directory()` |
| Goodness-of-fit diagnostics | `nm_gof()` |
| Compare two runs | `nm_compare()` (paired CWRES / NPDE, individual OFV) |
| VPC and pc-VPC | `nm_vpc_ctl()`, `nm_vpc()` |
| Parameter uncertainty | `nm_profile()` (likelihood profile), `nm_bootstrap()` (case bootstrap) |
| Influential subjects | `nm_jackknife()` (leave-one-subject-out jackknife: OFV change and THETA shifts) |
| Covariate assessment on the ETA distributions | `nm_wasserstein()` (Wasserstein distance, Simon et al. 2026) |
| Export figures | `nm_export_plot()` |

## Example project and vignettes

The package ships a small example project, with simulated data (60
fictitious subjects, one-compartment oral model) and the NONMEM outputs of a
BASE and a FINAL model, so that every function that reads results can be
tried without NONMEM:

```r
ex <- system.file("extdata", "example", package = "nmtools")
runs <- nm_extract_runs(ex)
```

Two vignettes use it: `vignette("nmtools")` (getting started, from the data
set to the comparison of two models) and `vignette("wasserstein")`.

## A typical workflow

```r
library(nmtools)

# 1. Data set: read, summarise, plot
dat <- nm_read_csv("C:/nm76g64/run/study/data.csv")
nm_summary(dat$data)
nm_plot_dv(dat$data, scales = "both")

# 2. Write a first control stream from the columns of the data set
nm_write_ctl("C:/nm76g64/run/study/run1.ctl",
             input = "C:/nm76g64/run/study/data.csv", error = "combined")

# 3. Run one or several models. By default each run goes to <run>.reg/ in the
#    analysis directory; run_subdir = "" runs it next to its control stream.
#    A relative $DATA path is read from the control stream folder
#    ($DATA data.csv) or from the run directory ($DATA ../data.csv).
log <- nm_run(
  nmfe_path    = "C:/nm76g64/run/nmfe76.bat",
  ctl_files    = c("run1.ctl", "run2.ctl"),
  analysis_dir = "C:/nm76g64/run/study"
)

# 4. Results of all the runs of the project
runs <- nm_extract_runs("C:/nm76g64/run/study")
runs                          # OFV, dOFV, minimisation and covariance status
nm_compare_runs(runs)         # parameters side by side
nm_extract_run_table(run_path = "C:/nm76g64/run/study/run2.reg")   # Word table

# 5. Diagnostics, and comparison of two models
gof <- nm_gof("C:/nm76g64/run/study/run2.reg/sdtab2")
plot(gof)
cmp <- nm_compare("C:/nm76g64/run/study/run1.reg",
                  "C:/nm76g64/run/study/run2.reg",
                  labels = c("Base", "+ WT on CL"))
plot(cmp)

# 6. VPC / pc-VPC
nm_vpc_ctl("run2", "C:/nm76g64/run/study",
           nmfe_path = "C:/nm76g64/run/nmfe76.bat", n_sim = 500)
nm_vpc("run2", "C:/nm76g64/run/study",
       col_names = c("NSIM", "ID", "TIME", "AMT", "MDV", "DV", "PRED", "IPRED"),
       type = 2)

# 7. Uncertainty
prof <- nm_profile("C:/nm76g64/run/nmfe76.bat",
                   "C:/nm76g64/run/study/run2.reg/run2.ctl",
                   param = "CL", from = 2, to = 6, by = 0.25)
plot(prof, type = "dofv")
boot <- nm_bootstrap("C:/nm76g64/run/study/run2.reg/run2.ctl",
                     nmfe_path = "C:/nm76g64/run/nmfe76.bat",
                     n = 500, seed = 1)
```

## Run layouts (Wings not required)

`nm_run()` writes each run either in its own subdirectory or directly in the
folder of its control stream:

| `run_subdir` | Outputs of `run1.ctl` |
|---|---|
| `".reg"` (default) | `run1.reg/run1.lst`, `run1.reg/run1.ext`... (same layout as Wings) |
| `".run"` or any suffix | `run1.run/run1.lst`... |
| `""` | `run1.lst`, `run1.ext`... next to `run1.ctl` (as a manual nmfe call) |

`nm_extract_run_table()`, `nm_vpc_ctl()` and `nm_vpc()` find a run in
`<run>.reg/`, `<run>.nm7/`, a folder named `<run>`, or directly in the
analysis directory. `nm_extract_runs()` reads every run of a folder, in
subdirectories or not. `nm_wasserstein()` and `nm_compare()` take a run
directory of any name or a run prefix (`"C:/study/run1"` for
`C:/study/run1.ext`, `run1.phi`...), and `nm_gof()` a table file or the
folder that holds it. Control streams
in Wings label syntax are converted automatically by `nm_run()` (this needs a
run subdirectory) or by `nm_convert_ctl()`.

## Writing control streams

`nm_write_ctl()` writes a complete control stream for a 1- to 3-compartment
PK model: IV or extravascular administration; first-order, zero-order
(modelled rate `R1 = K0` with `RATE = -1`, or modelled duration `D1` with
`RATE = -2`) or transit-compartment absorption; linear or Michaelis-Menten
elimination; analytic (`ADVAN`/`TRANS`) or ODE solution. `$INPUT` follows
the column order of the data set. Three functions then transform a control
stream, generated or written by hand, and `nm_ctl_check()` looks for
THETA/ETA/EPS used without a declaration and other static inconsistencies.

```r
ctl <- nm_write_ctl(input = c("ID", "TIME", "DV", "AMT", "EVID", "MDV", "CMT", "OCC", "BQL"),
                    cpts = 2, error = "combined",
                    theta_init = c(CL = 3, V2 = 30), iiv = c(CL = 30, V2 = 25, KA = 50))
ctl <- nm_ctl_iov(ctl, parameters = "CL", n_occ = 2)   # $ABBR REPLACE + BLOCK SAME
ctl <- nm_ctl_m3(ctl, lloq = 0.1)                       # M3 likelihood, LAPLACE
nm_ctl_check(ctl)
nm_ctl_mu(nm_write_ctl())                               # MU referencing for SAEM/IMP
```

The same functions are available as a Shiny application:
<https://nicolas-simon.shinyapps.io/nonmem_ctl_publi/>.

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

## Wasserstein-based covariate assessment

`nm_wasserstein()` implements the method of Simon N, Hulot JS, von Fabeck K.
Distribution-based covariate assessment using Wasserstein distance in
population pharmacokinetic models. *Front Pharmacol.* 2026;17:1804989
([doi:10.3389/fphar.2026.1804989](https://doi.org/10.3389/fphar.2026.1804989)).
It compares the distribution of an ETA between groups of subjects defined by
a covariate, for one or several runs, and reports a Wasserstein-based R2.

```r
w <- nm_wasserstein(c(BASE = "base.reg", FINAL = "final.reg"),
                    covariates = c("CLCR", "CT01"), eta = "CL")
w
plot(w)
plot(w, type = "r2w")
```

## Citation

```r
citation("nmtools")
```

gives the reference of the package and, for `nm_jackknife()` and
`nm_wasserstein()`, the reference of the original method.

## License

MIT © Nicolas Simon
