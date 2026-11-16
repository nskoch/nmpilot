# nmtools 0.4.1

* `nm_suggest_model()`: the residual error model is chosen on a common
  structural fit of each interval (iteratively reweighted least squares
  under a combined variance); the additive, proportional and combined
  variance models are then compared by BIC on the same leverage-corrected
  (PRESS) residuals, so that they are nested. Points of leverage above 0.7
  and predictions below twice the LLOQ are left out, and the combined model
  is kept only when its BIC is better by more than 10. On simulated rich
  data (20 replicates per scenario), the error model is recovered in 60 to
  100% of the datasets; better than 0.4.0 for Michaelis-Menten, additive and
  combined error, lower for some proportional-error scenarios
  (`inst/validation/validate_suggest_model.R`).

# nmtools 0.4.0

## Data-driven suggestion of structural models

* New `nm_suggest_model()`: analyses a NONMEM dataset, rich or sparse,
  without running NONMEM and ranks candidate models: number of compartments
  (sums of exponentials compared by BIC), input (IV bolus or infusion,
  first-order, first-order with lag time, zero-order duration, transit
  compartments; full-profile fits with a censoring penalty for early BLQ
  samples), elimination (dose proportionality and curvature of the
  log-concentration decline) and residual error model (additive,
  proportional or combined; extended least squares on leverage-corrected
  residuals). Each candidate comes with initial values for THETA, IIV
  (variance of individual estimates) and SIGMA, the parameters likely to be
  hard to estimate and their reasons, and alerts on features outside the
  tested structures (double peaks, clearance changing between first dose
  and steady state, unexpected accumulation, less than dose-proportional
  exposure). Sparse data are analysed as naive pooled data. `print()` and
  `plot()` methods (`"profiles"`, `"support"`, `"pooled"`, `"design"`).
* New `nm_suggest_design()`: expected RSE and eta-shrinkage of a candidate
  for the sampling design of the dataset or another one (population Fisher
  information matrix, first-order approximation, analytic solution with
  superposition of doses).
* New `nm_suggest_run()`: writes the control streams of the best
  candidates, runs them with `nm_run()` and compares the BIC ranking with
  the data-driven score.
* New `nm_suggest_covariates()`: Wasserstein-distance screening of
  covariates on the individual estimates (CL, central volume, absorption)
  before any modelling.
* `nm_suggest_run()` gives each candidate a dataset that matches its
  compartment numbering (depot 1 and central 2 for first-order and transit
  absorption; central 1 for IV and zero-order input, with `RATE = -2` on the
  doses for a modelled duration). A copy is written only when the original
  dataset does not match; the original file is never modified. With
  `overwrite = TRUE`, the candidates are run again even when previous runs
  of the same name have terminated.
* `nm_write_ctl()`: new arguments `from` and `candidate` to write the
  control stream of a candidate of `nm_suggest_model()`; explicit arguments
  keep priority. New argument `conc_factor`: the scale is written
  `S = V/conc_factor` (e.g. 1000 for doses in mg, volumes in L and DV in
  ng/mL); with `from`, the `dose_conc_factor` of the suggestion is used.
* New validation script `inst/validation/validate_suggest_model.R`
  (recognition rate of the true structure on simulated datasets).

# nmtools 0.3.0

## Writing control streams

Functions ported from the Shiny application "Générateur NONMEM .ctl"
(<https://nicolas-simon.shinyapps.io/nonmem_ctl_publi/>). All return an
`nm_ctl` object (control stream lines, printed as text) and can write a file.

* New `nm_write_ctl()`: writes a control stream for a 1- to 3-compartment PK
  model, IV or extravascular, with first-order, zero-order (modelled rate
  `R1 = K0`, `RATE = -1`, or modelled duration `D1`, `RATE = -2`) or
  transit-compartment absorption, linear or Michaelis-Menten elimination,
  analytic (`ADVAN`/`TRANS`) or ODE solution. `$INPUT` follows the column
  order of the dataset; THETA bounds and FIX, IIV as CV% or variance,
  diagonal or block `$OMEGA`, additive, proportional or combined error,
  nine estimation methods, xpose-style `sdtab`/`patab`/`cotab` tables.
* New `nm_ctl_mu()`: MU referencing of log-normal parameters (`"readable"`
  or `"linear"` mode).
* New `nm_ctl_iov()`: inter-occasion variability with `$ABBR REPLACE` and
  `$OMEGA BLOCK(p) SAME(n-1)`; optional check of the occasion column.
* New `nm_ctl_m3()`: M3 method for BQL data, with the LLOQ given as a
  constant, a dataset column, or by analytical method code.
* New `nm_ctl_check()`: static check of a control stream (THETA/ETA/EPS used
  but not declared, over all records including `BLOCK`, `SAME` and
  `$ABBR REPLACE`; `$INPUT` items; M3 requirements), with optional line
  differences from a reference control stream.

# nmtools 0.2.3

* New `nm_jackknife()`: leave-one-subject-out jackknife of a NONMEM model
  (method of the standalone script nonmem-jackknife-influence,
  doi:10.5281/zenodo.22257994). Ranks the subjects by their likelihood
  influence |dOFV| and reports the relative THETA shifts (L2 norm, maximum,
  SE-standardised norm); `print()` and `plot()` methods (`"bar"`, `"rank"`,
  `"scatter"`). Resumes an interrupted analysis. Data lines are copied
  verbatim, as in `nm_bootstrap()`.
* `nm_identify_parameters()`: new column `derived_from` in
  `parameter_table`, naming the parameters a derived parameter is computed
  from (`"V"` for `S2 = V`, `"CL, V"` for `KE = CL/V`; `NA` otherwise).
* `nm_extract_runs()`: covariate-effect and residual-error THETAs are named
  from their `$THETA` comment label when it is a short name (`CL_WT`, `PROP`,
  `RUV_A`, `4 Prop (...)` -> `Prop`) instead of `THETA4`, `THETA5`...

# nmtools 0.2.2

## Breaking changes

* Output column names are now in English. `nm_extract_run_table()`,
  `params` table: `Parameter`, `Estimate`, `RSE_pct`, `CI95`, `CI95_lower`,
  `CI95_upper` (were `Parametre`, `Estimation`, `RSE_pourcent`, `IC95`,
  `IC95_lower`, `IC95_upper`). `nm_summary()` and `nm_data_previsual()`:
  `unique_doses` (was `doses_uniques`).

## Wings for NONMEM is optional

* `nm_run()`: new argument `run_subdir`. `".reg"` (default) keeps the
  Wings-compatible layout `run1.reg/`; any other suffix can be used; `""`
  runs NONMEM in the folder of the control stream, without copying it.
  A control stream in Wings label syntax cannot be run with `run_subdir = ""`
  (clear error message).
* `nm_run()`: a relative `$DATA` path is looked up from the folder of the
  control stream (`$DATA data.csv`), then from the run directory
  (`$DATA ../data.csv`, Wings style).
* `nm_vpc_ctl()`, `nm_vpc()` and `nm_extract_run_table()` find a run in
  `<run>.reg/`, `<run>.nm7/`, a folder named `<run>`, or directly in the
  analysis directory. `nm_vpc_ctl()` gains `run_subdir`; by default the VPC
  runs use the layout of the estimation run.

## Bug fixes

* `nm_run()`: NONMEM is now launched from the run directory (NM-TRAN failed
  when started from another working directory).
* `nm_run()`: standard control streams with comments after `$THETA` or
  `$OMEGA` values are no longer taken for Wings control streams.
* `nm_run()`: `verbose = FALSE` now silences all messages.

## Documentation

* Example project in `inst/extdata/example`: simulated data (60 fictitious
  subjects, one-compartment oral model) and the NONMEM outputs of a BASE and
  a FINAL model.
* Two vignettes: `vignette("nmtools")` (getting started) and
  `vignette("wasserstein")`.
* `citation("nmtools")` cites the package and, for `nm_wasserstein()`, the
  original article.

# nmtools 0.2.1

* `nm_wasserstein()`: new argument `continuous = "pairwise"` to compare every
  pair of quartiles of a continuous covariate (6 comparisons), as in the
  pairwise analysis of the original article; Bonferroni-adjusted p-values
  `p_W_adj` and `p_t_adj` are added to the results. The default
  (`"extremes"`, Q4 vs Q1) is unchanged.
* `nm_identify_parameters()`: a THETA assigned under an IF condition whose
  other branch is a plain constant (e.g. `IF (CT01.EQ.1) FALL = THETA(7)`,
  `ELSE FALL = 1`) is now a covariate effect; a covariate term added inside a
  sum written with extra parentheses (e.g.
  `((THETA(1)*FALL) + (THETA(8)*(CLCR/89)))`) is now a covariate effect.

# nmtools 0.2.0

* New `nm_wasserstein()`: distribution-based covariate assessment with the
  Wasserstein distance (Simon, Hulot and von Fabeck, Front Pharmacol 2026,
  doi:10.3389/fphar.2026.1804989). For one or several runs, compares the
  distribution of an ETA between groups of subjects defined by covariates
  (permutation test on the 1-Wasserstein distance, with a t-test for
  reference) and reports a Wasserstein-based R2 with the eta-shrinkage.
  `print()` and `plot()` methods (ECDF by group, R2W bar chart).

# nmtools 0.1.0

First public release.

* Data sets: `nm_read_csv()`, `nm_data_previsual()`, `nm_plot_dv()`,
  `nm_summary()`, `nm_add_rate()`, `nm_export_data()`.
* Running NONMEM: `nm_run()`, `nm_convert_ctl()` (Wings extended control
  streams).
* Results: `nm_extract_runs()`, `nm_filter_runs()`, `nm_compare_runs()`,
  `nm_extract_run_table()`.
* Parameter identification from the model code: `nm_identify_parameters()`,
  `nm_identify_directory()`.
* Diagnostics: `nm_gof()`, `nm_compare()`, `nm_vpc_ctl()`, `nm_vpc()`.
* Uncertainty: `nm_profile()`, `nm_bootstrap()`.
* Figures: `nm_export_plot()`.
