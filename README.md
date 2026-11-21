# nmpilot

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24244339.svg)](https://doi.org/10.5281/zenodo.24244339)
<!-- badges: end -->

**nmpilot** is an R package that drives a population pharmacokinetic (popPK)
analysis with NONMEM from R, from the data set to the final model diagnostics.
It calls the NONMEM `nmfe` launcher directly and does not need PsN or
Wings for NONMEM. Wings is optional: runs can sit next to their control
streams or in run subdirectories, and the Wings conventions (`.reg` run
directories, `.smr` files, extended control streams with parameter labels)
are also read, so that projects started with Wings can be continued.

> **Status:** version 0.7.1, under active development. The interface may still
> change. Bug reports and suggestions are welcome in the
> [issues](https://github.com/nskoch/nmpilot/issues).

> **New name:** the package was called **nmtools** up to version 0.7.1. It was
> renamed **nmpilot** to avoid confusion with an unrelated GitHub repository of
> the same name in the same field. If you installed nmtools, remove it with
> `remove.packages("nmtools")`, install nmpilot as below, and replace
> `library(nmtools)` by `library(nmpilot)` in your scripts.

## Installation

```r
# install.packages("remotes")
remotes::install_github("nskoch/nmpilot")

# with the vignettes:
remotes::install_github("nskoch/nmpilot", build_vignettes = TRUE)
```

Running models requires a working NONMEM installation. The package is
developed and used on Windows with NONMEM 7.4 and 7.6 (`nmfe74.bat`,
`nmfe76.bat`); Linux and macOS are not tested yet. Reading outputs,
diagnostics and plots work without NONMEM, on outputs of NONMEM 7.1 to 7.6.

## What the package covers

| Step | Functions |
|---|---|
| Read and explore the data set | `nm_read_csv()`, `nm_data_previsual()`, `nm_plot_dv()`, `nm_summary()`, `nm_add_rate()`, `nm_export_data()` |
| Suggest a structural model before any run | `nm_suggest_model()` (compartments, input, elimination, error model, initial values, hard-to-estimate parameters, alerts), `nm_suggest_design()` (expected RSE and shrinkage), `nm_suggest_run()` (run and rank the candidates), `nm_suggest_covariates()` (covariate screening on individual estimates) |
| Write a control stream | `nm_write_ctl()` (compartmental PK model), `nm_ctl_mu()` (MU referencing), `nm_ctl_iov()` (inter-occasion variability), `nm_ctl_m3()` (M3 method for BQL data), `nm_ctl_check()` (static check) |
| Run NONMEM | `nm_run()` (`n_cores`: several runs at the same time), `nm_cores()`, `nm_convert_ctl()` (Wings extended control stream to standard NONMEM) |
| Read the results of one or several runs | `nm_extract_runs()` (also several `$EST` or `$PROB`, IOV, `$MIX`, failed or interrupted runs), `nm_filter_runs()`, `nm_compare_runs()`, `nm_extract_run_table()` (Word table) |
| Identify what each THETA / ETA means | `nm_identify_parameters()`, `nm_identify_directory()` |
| Goodness-of-fit diagnostics | `nm_gof()` |
| Compare two runs | `nm_compare()` (paired CWRES / NPDE, individual OFV; Vuong and Clarke tests for non-nested models, bootstrap, gain map, overlaid diagnostics) |
| Improve a base model from one evaluation run | `nm_optimize_base()` (residual error, random effects and covariates ranked on a linearised model), `nm_derivatives_ctl()` |
| VPC and pc-VPC | `nm_vpc_ctl()`, `nm_vpc()` |
| Parameter uncertainty | `nm_profile()` (likelihood profile), `nm_bootstrap()` (case bootstrap); both accept `n_cores` |
| Influential subjects | `nm_jackknife()` (leave-one-subject-out jackknife: OFV change and THETA shifts) |
| Covariate assessment on the ETA distributions | `nm_wasserstein()` (Wasserstein distance, Simon et al. 2026) |
| Export figures | `nm_export_plot()` |

## Example project and vignettes

The package ships a small example project, with simulated data (60
fictitious subjects, one-compartment oral model) and the NONMEM outputs of a
BASE and a FINAL model, so that every function that reads results can be
tried without NONMEM:

```r
ex <- system.file("extdata", "example", package = "nmpilot")
runs <- nm_extract_runs(ex)
```

Two vignettes use it: `vignette("nmpilot")` (getting started, from the data
set to the comparison of two models) and `vignette("wasserstein")`.

## A typical workflow

```r
library(nmpilot)

# 1. Data set: read, summarise, plot
dat <- nm_read_csv("C:/nm76g64/run/study/data.csv")
nm_summary(dat$data)
nm_plot_dv(dat$data, scales = "both")

# 2. Suggest a structural model and write its control stream
s <- nm_suggest_model("C:/nm76g64/run/study/data.csv", dose_conc_factor = 1000)
s                                   # ranked candidates, initial values, alerts
nm_write_ctl("C:/nm76g64/run/study/run1.ctl", from = s, candidate = 1)
#    or write a control stream from the columns of the data set
nm_write_ctl("C:/nm76g64/run/study/run2.ctl",
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

## Several NONMEM runs at the same time

`nm_run()`, `nm_profile()`, `nm_jackknife()`, `nm_bootstrap()`,
`nm_vpc_ctl()` and `nm_suggest_run()` take an `n_cores` argument (default 1:
one run after the other). Each run is then started in the background from
its own folder, at most `n_cores` at a time. `n_cores = "auto"` uses the
physical cores minus one; `n_cores = "ask"` shows the cores and asks the
number once per R session.

```r
nm_cores()          # logical and physical cores
log <- nm_run("C:/nm76g64/run/nmfe76.bat", sprintf("run%d.ctl", 1:4),
              "C:/nm76g64/run/study", n_cores = 4)
jk  <- nm_jackknife("C:/nm76g64/run/study/run2.reg/run2.ctl",
                    nmfe_path = "C:/nm76g64/run/nmfe76.bat", n_cores = "auto")
```

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

## Suggesting a structural model

`nm_suggest_model()` looks at the data before any NONMEM run: NCA and sums of
exponentials on each rich dosing interval (or on naive pooled data when the
sampling is sparse) give the support for 1, 2 or 3 compartments, for the
input model (first-order, lag time, zero-order, transit), for linear or
saturable elimination and for an additive, proportional or combined error.
Candidates come with initial values in the names used by `nm_write_ctl()`,
the parameters likely to be hard to estimate, and, for linear models, the
expected RSE and eta-shrinkage for the sampling design.

```r
s <- nm_suggest_model("data.csv", dose_conc_factor = 1000)   # mg and ng/mL
s$candidates
plot(s, "profiles")
nm_suggest_design(s, candidate = 1, times = c(1, 4, 12))     # another design
nm_suggest_covariates(s, covariates = c("WT", "SEX"))
r <- nm_suggest_run(s, analysis_dir = "C:/nm76g64/run/study",
                    nmfe_path = "C:/nm76g64/run/nmfe76.bat", candidates = 1:3)
```

The scores rank candidates; they are no substitute for model comparison in
NONMEM. The validation script `inst/validation/validate_suggest_model.R`
measures how often the true structure is recovered on simulated datasets.

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
citation("nmpilot")
```

gives the reference of the package and, for `nm_jackknife()` and
`nm_wasserstein()`, the reference of the original method. Each release is
archived on Zenodo with its own DOI; the DOI badge above
(10.5281/zenodo.24244339) cites all versions and resolves to the latest one
(v0.7.1: 10.5281/zenodo.24406632). The archives up to v0.7.1 carry the
former name, nmtools.
The releases v0.6.0 and v0.7.0 were made by mistake on the commit of
version 0.5.0 and hold its code: use v0.7.1 or later.

## License

MIT © Nicolas Simon
