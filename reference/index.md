# Package index

## 

Data preparation & inspection

- [`nm_read_csv()`](https://nskoch.github.io/nmpilot/reference/nm_read_csv.md)
  : Read and normalise a NONMEM-like dataset
- [`nm_check_dataset()`](https://nskoch.github.io/nmpilot/reference/nm_check_dataset.md)
  : Check that a dataset is compatible with NONMEM
- [`nm_summary()`](https://nskoch.github.io/nmpilot/reference/nm_summary.md)
  : Diagnostic summary of a NONMEM-like dataset
- [`nm_plot_dv()`](https://nskoch.github.io/nmpilot/reference/nm_plot_dv.md)
  : Plot DV vs TIME with median and interval band
- [`nm_data_previsual()`](https://nskoch.github.io/nmpilot/reference/nm_data_previsual.md)
  : Pre-visualise a NONMEM-like dataset
- [`nm_add_rate()`](https://nskoch.github.io/nmpilot/reference/nm_add_rate.md)
  : Add a RATE column to a NONMEM dataset
- [`nm_export_data()`](https://nskoch.github.io/nmpilot/reference/nm_export_data.md)
  : Export a data.frame to CSV (NONMEM convention: na = ".")

## 

Exploration & base model

- [`nm_suggest_model()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_model.md)
  [`print(`*`<nm_model_suggestion>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_suggest_model.md)
  : Suggest structural PK models from a NONMEM dataset
- [`plot(`*`<nm_model_suggestion>`*`)`](https://nskoch.github.io/nmpilot/reference/plot.nm_model_suggestion.md)
  : Plot a model suggestion
- [`nm_suggest_design()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_design.md)
  [`print(`*`<nm_suggest_design>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_suggest_design.md)
  [`plot(`*`<nm_suggest_design>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_suggest_design.md)
  : Expected precision of a suggested model for a sampling design
- [`nm_suggest_run()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_run.md)
  [`print(`*`<nm_suggest_run>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_suggest_run.md)
  : Write, run and rank the candidates of a model suggestion
- [`nm_suggest_covariates()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_covariates.md)
  [`print(`*`<nm_suggest_covariates>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_suggest_covariates.md)
  : Covariate screening on individual PK estimates before modelling

## 

Control streams

- [`nm_write_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_write_ctl.md)
  : Write a NONMEM control stream for a population PK model
- [`nm_ctl_check()`](https://nskoch.github.io/nmpilot/reference/nm_ctl_check.md)
  : Static check of a NONMEM control stream
- [`nm_ctl_mu()`](https://nskoch.github.io/nmpilot/reference/nm_ctl_mu.md)
  : MU referencing of a NONMEM control stream
- [`nm_ctl_iov()`](https://nskoch.github.io/nmpilot/reference/nm_ctl_iov.md)
  : Add inter-occasion variability to a NONMEM control stream
- [`nm_ctl_m3()`](https://nskoch.github.io/nmpilot/reference/nm_ctl_m3.md)
  : Implement the M3 method for BQL data in a NONMEM control stream
- [`nm_convert_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_convert_ctl.md)
  : Convert a Wings for NONMEM extended control stream to standard
  NONMEM

## 

Running NONMEM

- [`nm_run()`](https://nskoch.github.io/nmpilot/reference/nm_run.md) :
  Run one or more NONMEM models
- [`nm_cores()`](https://nskoch.github.io/nmpilot/reference/nm_cores.md)
  : Number of processor cores available for parallel NONMEM runs

## 

Results

- [`nm_extract_runs()`](https://nskoch.github.io/nmpilot/reference/nm_extract_runs.md)
  : Extract the results of several NONMEM runs
- [`nm_filter_runs()`](https://nskoch.github.io/nmpilot/reference/nm_filter_runs.md)
  : Filter NONMEM runs and recompute relative criteria
- [`nm_compare_runs()`](https://nskoch.github.io/nmpilot/reference/nm_compare_runs.md)
  : Compare the parameters of several runs side by side
- [`nm_extract_run_table()`](https://nskoch.github.io/nmpilot/reference/nm_extract_run_table.md)
  : Parameter table of one NONMEM run, exported to Word
- [`nm_identify_parameters()`](https://nskoch.github.io/nmpilot/reference/nm_identify_parameters.md)
  : Identify the semantic role of THETA and ETA in a NONMEM control
  stream
- [`nm_identify_directory()`](https://nskoch.github.io/nmpilot/reference/nm_identify_directory.md)
  : Identify THETA/ETA roles for all control streams of a directory
- [`nm_model_diagram()`](https://nskoch.github.io/nmpilot/reference/nm_model_diagram.md)
  : Diagram of the models tested in a project
- [`nm_runs_table()`](https://nskoch.github.io/nmpilot/reference/nm_runs_table.md)
  : Summary table of all the runs of a project

## 

Diagnostics

- [`nm_gof()`](https://nskoch.github.io/nmpilot/reference/nm_gof.md) :
  Goodness-of-fit diagnostic plots from a NONMEM \$TABLE file

## 

Model selection & comparison

- [`nm_compare()`](https://nskoch.github.io/nmpilot/reference/nm_compare.md)
  [`print(`*`<nm_compare>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_compare.md)
  [`summary(`*`<nm_compare>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_compare.md)
  : Compare two NONMEM runs

- [`plot(`*`<nm_compare>`*`)`](https://nskoch.github.io/nmpilot/reference/plot.nm_compare.md)
  :

  Diagnostic plots for `nm_compare` objects

- [`nm_derivatives_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_derivatives_ctl.md)
  : Control stream of derivatives for the linearization of a base model

- [`nm_optimize_base()`](https://nskoch.github.io/nmpilot/reference/nm_optimize_base.md)
  [`print(`*`<nm_optimize_base>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_optimize_base.md)
  [`plot(`*`<nm_optimize_base>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_optimize_base.md)
  : Proposals to improve a base model without iterative NONMEM runs

- [`nm_wasserstein()`](https://nskoch.github.io/nmpilot/reference/nm_wasserstein.md)
  [`print(`*`<nm_wasserstein>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_wasserstein.md)
  : Distribution-based covariate assessment with the Wasserstein
  distance

- [`plot(`*`<nm_wasserstein>`*`)`](https://nskoch.github.io/nmpilot/reference/plot.nm_wasserstein.md)
  :

  Plots for `nm_wasserstein` objects

## 

Model validation

- [`nm_vpc_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_vpc_ctl.md)
  : Generate VPC control streams and launch NONMEM runs

- [`nm_vpc()`](https://nskoch.github.io/nmpilot/reference/nm_vpc.md) :
  Visual Predictive Check (VPC or pc-VPC) from NONMEM simulation output

- [`nm_bootstrap()`](https://nskoch.github.io/nmpilot/reference/nm_bootstrap.md)
  : Nonparametric (case) bootstrap of a NONMEM model

- [`nm_jackknife()`](https://nskoch.github.io/nmpilot/reference/nm_jackknife.md)
  [`print(`*`<nm_jackknife>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_jackknife.md)
  : Leave-one-subject-out jackknife of a NONMEM model

- [`plot(`*`<nm_jackknife>`*`)`](https://nskoch.github.io/nmpilot/reference/plot.nm_jackknife.md)
  :

  Plots for `nm_jackknife` objects

- [`nm_profile()`](https://nskoch.github.io/nmpilot/reference/nm_profile.md)
  [`print(`*`<nm_profile>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_profile.md)
  [`plot(`*`<nm_profile>`*`)`](https://nskoch.github.io/nmpilot/reference/nm_profile.md)
  : Profile the objective function over one model parameter

## 

Reporting

- [`nm_export_plot()`](https://nskoch.github.io/nmpilot/reference/nm_export_plot.md)
  : Export a ggplot to an image file

## Task guide

- [`nm_guide()`](https://nskoch.github.io/nmpilot/reference/nm_guide.md)
  : Interactive task guide
- [`nm_guide_addin()`](https://nskoch.github.io/nmpilot/reference/nm_guide_addin.md)
  : RStudio addin: open the task guide
