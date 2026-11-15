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
