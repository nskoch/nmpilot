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
