Setup used by the examples

Copy

    library(nmpilot)
    # Copy the example project to a temporary folder (runs write files there)
    file.copy(system.file("extdata", "example", package = "nmpilot"),
              tempdir(), recursive = TRUE)
    ex   <- file.path(tempdir(), "example")
    nmfe <- "C:/nm76g64/run/nmfe76.bat"   # path to your nmfe launcher

Examples that only read results need `ex`; examples that run NONMEM also
need `nmfe`. Objects such as `dat`, `s`, `runs` or `cmp` are created by
earlier tasks of the guide.

1\. Data preparation & inspection 7

2\. Exploration & base model 4

3\. Control streams 6

4\. Running NONMEM 2

5\. Results 6

6\. Diagnostics 1

7\. Model selection & comparison 4

8\. Model validation 5

9\. Reporting 2

Read and standardise a dataset

Summarise the dataset

Plot concentrations vs time

One-step preview (read, summary, plots)

Add a RATE column

Export in NONMEM format

Check NONMEM compatibility

Suggest a structural model and initial values

Expected precision of the design

Run the top candidates

Pre-screen covariates

Write a control stream

Static check

MU referencing

Inter-occasion variability

M3 method for BLQ data

Convert a Wings label control stream

Run one or several models

Number of cores available

Summarise all runs of a project

Parameters side by side

Parameter table of one run (Word)

Identify the role of THETAs and ETAs

Project summary in Word planned

Diagram of tested models

Goodness-of-fit plots

Compare two models (LRT, Vuong, Clarke, bootstrap, influence)

Comparison plots (gain map, overlaid diagnostics)

Optimise a base model without iterative runs

Covariate effects on ETAs (Wasserstein)

VPC

pcVPC

Bootstrap

Jackknife

Likelihood profile

Export a figure

Diagnostics in Word

## Read and standardise a dataset

`nm_read_csv`

Reads a NONMEM CSV file, recognises the reserved columns (ID, TIME, DV,
AMT, EVID, MDV, CMT...), converts '.' to missing values, computes the
relative time, the time after dose (TAD) and the dose-normalised DV, and
flags BLQ records.

Help in R: `?nm_read_csv`

Copy

    dat <- nm_read_csv(
      file.path(ex, "example_data.csv")   # file: path of the NONMEM CSV data set
    )
    head(dat$data)     # standardised data (reserved NONMEM names)
    dat$mapping        # reserved labels -> original column names

## Summarise the dataset

`nm_summary`

Number of subjects, records and observations, proportion of BLQ, time
range, distinct doses and observed compartments.

Help in R: `?nm_summary`

Copy

    nm_summary(
      dat$data    # data: data frame returned by nm_read_csv() ($data)
    )

## Plot concentrations vs time

`nm_plot_dv`

Observed concentrations against TIME or TAD, with the median and an
interval band, on linear and/or log scale, optionally stratified by a
covariate.

Help in R: `?nm_plot_dv`

Copy

    nm_plot_dv(dat$data,          # data: data frame from nm_read_csv()
               scales = "both")   # y axis: "linear", "log" or "both" (side by side)
    nm_plot_dv(dat$data,
               time_var = "TAD",  # x axis: "TIME" (default) or "TAD" (time after dose)
               scales   = "log")
    nm_plot_dv(dat$data,
               y     = "DV_norm", # y: "DV" (default) or "DV_norm" (DV / dose)
               strat = "GENO")    # column used to colour the profiles

## One-step preview (read, summary, plots)

`nm_data_previsual`

Reads the dataset, summarises it and draws the raw and dose-normalised
profiles in one call, with configurable BLQ handling and optional
JPEG/CSV export.

Help in R: `?nm_data_previsual`

Copy

    pv <- nm_data_previsual(
      file.path(ex, "example_data.csv"),  # path of the CSV data set
      scales    = "both",                 # y axis: "linear", "log" or "both"
      dose_norm = TRUE                    # also plot DV_norm = DV / AMT
    )
    pv$summary     # summary table
    pv$plot_dv     # concentrations vs time
    pv$plot_norm   # dose-normalised concentrations

## Add a RATE column

`nm_add_rate`

Adds a RATE column (e.g. -2 for an estimated zero-order duration) on all
doses or only on some routes, keeping the exact format of the file.

Help in R: `?nm_add_rate`

Copy

    # writes example_data_rate.csv next to the source file
    nm_add_rate(
      file.path(ex, "example_data.csv"),  # file: path of the original CSV
      dose_value = -2                     # RATE on doses: -2 (duration estimated)
    )                                     #   or -1 (rate estimated)

## Export in NONMEM format

`nm_export_data`

Writes a data.frame as CSV with the NONMEM convention for missing values
('.').

Help in R: `?nm_export_data`

Copy

    nm_export_data(
      dat$data,                       # data: data frame to write
      file.path(ex, "data_nm.csv")    # file: output path (folder created if needed)
    )

## Check NONMEM compatibility

`nm_check_dataset`

Checks that a dataset can be read by NONMEM: missing ID or TIME,
non-numeric values, decreasing times, inconsistent EVID/MDV/AMT,
duplicates, outliers. Report only: the data are never modified.

Help in R: `?nm_check_dataset`

Copy

    chk <- nm_check_dataset(
      file.path(ex, "example_data.csv")   # data: CSV path (or a data frame)
    )
    chk$ok          # TRUE when no error was found
    chk$findings    # one row per problem: check, severity, rows, IDs, message
    summary(chk)    # number of findings per check
    # columns with other names: cols = c(ID = "SUBJ", DV = "CONC")

## Suggest a structural model and initial values

`nm_suggest_model` `plot.nm_model_suggestion`

Analyses the data before any NONMEM run and ranks candidate models
(compartments, input, elimination, residual error), each with initial
values for THETA, IIV and SIGMA and alerts on what no candidate
describes.

Help in R: `?nm_suggest_model`, `?plot.nm_model_suggestion`

Copy

    s <- nm_suggest_model(
      file.path(ex, "example_data.csv")   # data: path of the CSV (or a data frame)
    )
    s                  # ranked candidates, initial values, alerts
    plot(s, "profiles")   # type: "profiles", "support", "pooled" or "design"
    plot(s, "support")

## Expected precision of the design

`nm_suggest_design`

Expected RSE of the parameters and eta-shrinkage of a candidate for the
sampling design of the dataset or for another design, from the
population Fisher information matrix.

Help in R: `?nm_suggest_design`

Copy

    nm_suggest_design(s,              # s: result of nm_suggest_model()
                      candidate = 1)  # candidate number (default 1)
    # sparse design: 3 samples, 100 subjects
    nm_suggest_design(s,
                      times      = c(1, 4, 12),  # sampling times after the dose
                      n_subjects = 100)          # number of subjects

## Run the top candidates

`nm_suggest_run`

Writes the control streams of the chosen candidates, runs them with
nm_run(), reads their OFV and ranks them by BIC.

Help in R: `?nm_suggest_run`

Copy

    r <- nm_suggest_run(
      s,                     # result of nm_suggest_model()
      analysis_dir = ex,     # folder where the control streams are written and run
      nmfe_path    = nmfe,   # path of the nmfe launcher
      candidates   = 1:3,    # candidate numbers to run (default: the first 3)
      n_cores      = 3       # number of runs at the same time
    )
    r                        # OFV, BIC, ranks, Spearman correlation

## Pre-screen covariates

`nm_suggest_covariates`

Compares the individual estimates (CL, central volume, absorption)
between covariate groups with the Wasserstein distance, before any
modelling.

Help in R: `?nm_suggest_covariates`

Copy

    nm_suggest_covariates(
      s,                            # result of nm_suggest_model()
      covariates = c("WT", "GENO")  # covariate columns of the data set
    )

## Write a control stream

`nm_write_ctl`

Writes a complete control stream for a compartmental PK model (1 to 3
compartments, IV or extravascular, first-order, zero-order or transit
absorption, linear or Michaelis-Menten elimination), possibly from a
candidate of nm_suggest_model().

Help in R: `?nm_write_ctl`

Copy

    # from candidate 1 of a suggestion
    nm_write_ctl(file.path(ex, "run3.ctl"),  # file: control stream to write
                 from      = s,              # result of nm_suggest_model()
                 candidate = 1)              # candidate number in `from`
    # or from the dataset columns
    nm_write_ctl(input = file.path(ex, "example_data.csv"),  # data set (or its column names)
                 error = "combined")   # "additive", "proportional" or "combined"

## Static check

`nm_ctl_check`

Checks a control stream without NM-TRAN: THETA, ETA and EPS used without
declaration, items missing from \$INPUT, M3 requirements; optional line
differences from a reference.

Help in R: `?nm_ctl_check`

Copy

    nm_ctl_check(
      file.path(ex, "run2.ctl")   # ctl: control stream (path or lines)
    )
    ctl <- nm_write_ctl(input = file.path(ex, "example_data.csv"))
    nm_ctl_check(nm_ctl_mu(ctl),     # ctl: the transformed control stream
                 reference = ctl)    # original: line differences are added

## MU referencing

`nm_ctl_mu`

Rewrites log-normal parameters as MU_m = LOG(THETA(n)) and P =
EXP(MU_m + ETA(m)), as required by EM and Bayesian methods (SAEM, IMP,
BAYES).

Help in R: `?nm_ctl_mu`

Copy

    ctl <- nm_write_ctl(input = file.path(ex, "example_data.csv"),
                        error = "combined")
    nm_ctl_mu(ctl)       # ctl: control stream (path or lines); all parameters
    nm_ctl_mu(ctl,
              parameters = c("CL", "V"),  # parameters to transform
              mode       = "linear")      # "readable" (default) or "linear"

## Inter-occasion variability

`nm_ctl_iov`

Adds IOV to chosen parameters with \$ABBR REPLACE and \$OMEGA BLOCK
SAME; the occasion item must be in \$INPUT.

Help in R: `?nm_ctl_iov`

Copy

    cols <- c("ID", "TIME", "AMT", "DV", "EVID", "MDV", "CMT", "OCC", "BQL")
    ctl  <- nm_write_ctl(input = cols,         # column names of the data set
                         error = "combined")
    nm_ctl_iov(ctl,
               parameters = c("CL", "V"),  # parameters receiving IOV
               n_occ      = 3,             # number of occasions (at least 2)
               occ        = "OCC")         # occasion column, values 1 to n_occ

## M3 method for BLQ data

`nm_ctl_m3`

Applies the M3 method (Beal 2001): BLQ observations contribute P(DV \<
LLOQ); the LLOQ is a constant, a dataset column or depends on the
analytical method.

Help in R: `?nm_ctl_m3`

Copy

    cols <- c("ID", "TIME", "AMT", "DV", "EVID", "MDV", "CMT", "OCC", "BQL")
    ctl  <- nm_write_ctl(input = cols, error = "combined")
    nm_ctl_m3(ctl,
              bql  = "BQL",   # BQL flag column (0 = quantified, 1 = below LLOQ)
              lloq = 0.02)    # LLOQ, on the DV scale

## Convert a Wings label control stream

`nm_convert_ctl`

Converts a Wings for NONMEM 'extended' control stream, where the model
uses the parameter labels, into a standard one (THETA(n), ETA(n),
ERR(n)).

Help in R: `?nm_convert_ctl`

Copy

    # a Wings control stream of your own (not in the example project)
    nm_convert_ctl(
      "C:/myproject/run_wings.ctl"   # ctl: Wings control stream with labels
    )                                # -> run_wings_converted.ctl

## Run one or several models

`nm_run`

Runs control streams with nmfe, without Wings or PsN, each in its own
run directory; finished runs are skipped. n_cores runs several models at
the same time.

Help in R: `?nm_run`

Copy

    log <- nm_run(
      ctl          = c("run1.ctl", "run2.ctl"),  # control streams (relative to analysis_dir)
      nmfe_path    = nmfe,                       # path of the nmfe launcher
      analysis_dir = ex,                         # analysis folder (one subfolder per run)
      n_cores      = 2                           # runs at the same time ("auto", "ask")
    )
    log[, c("run", "action", "ofv", "min_ok", "cov_status")]

## Number of cores available

`nm_cores`

Logical and physical cores of the computer and the number of
simultaneous runs used by n_cores = "auto".

Help in R: `?nm_cores`

Copy

    nm_cores()   # no argument: logical and physical cores of this computer

## Summarise all runs of a project

`nm_extract_runs` `nm_filter_runs`

Reads every run of an analysis folder: OFV, dOFV, AIC/BIC, minimisation
and covariance status, warnings, and all parameters named from the model
code (CL, V, KA...).

Help in R: `?nm_extract_runs`, `?nm_filter_runs`

Copy

    runs <- nm_extract_runs(
      ex      # path: analysis folder (or the folder of a single run)
    )
    runs      # one line per run
    nm_filter_runs(runs,           # x: result of nm_extract_runs()
                   "successful")$runs   # filter: "all", "successful",
                                        #   "successful_no_warning", "clean"...

## Parameters side by side

`nm_compare_runs`

One row per parameter and one column per run, aligned by name or by
NONMEM index; values, RSE, SE or shrinkage.

Help in R: `?nm_compare_runs`

Copy

    nm_compare_runs(runs,           # result of nm_extract_runs()
                    format = TRUE)  # cells written "value (RSE%)"
    nm_compare_runs(runs,
                    what = "shrinkage")   # "value", "rse_pct", "se" or "shrinkage"

## Parameter table of one run (Word)

`nm_extract_run_table`

Parameter table of a run (estimates, SE, RSE, 95% CI, SD of ETAs and
EPS, ETA correlations) exported to a Word document in the run directory.

Help in R: `?nm_extract_run_table`

Copy

    res <- nm_extract_run_table(
      file.path(ex, "run2.reg")   # run: run folder, control stream or stem
    )
    res$params          # the table in R
    res$output_file     # the Word file

## Identify the role of THETAs and ETAs

`nm_identify_parameters` `nm_identify_directory`

Reads \$PK or \$PRED and links each THETA and ETA to the PK parameter it
acts on, with its role (typical value, covariate effect, IIV, residual
error).

Help in R: `?nm_identify_parameters`, `?nm_identify_directory`

Copy

    x <- nm_identify_parameters(
      file.path(ex, "run2.ctl")   # ctl: path of the control stream
    )
    x$parameter_table
    nm_identify_directory(ex)     # path: every control stream of a folder

## Project summary in Word planned

Summary table of all runs of a project exported to Word (today
nm_extract_runs() gives it in R only).

Not available yet: this task is planned for a future version of nmpilot.

## Diagram of tested models

`nm_model_diagram`

Tree of the runs with, for each, its parent run, OFV, dOFV against the
parent and convergence status. Parents are read from the ";; Based on:"
comments of the control streams, or given with parents.

Help in R: `?nm_model_diagram`

Copy

    d <- nm_model_diagram(
      ex,                         # x: analysis folder (or nm_extract_runs() result)
      parents = c(run2 = "run1")  # optional links run = parent, added to (or
    )                             #   replacing) the ";; Based on:" comments
    d                       # table of the runs
    plot(d)                 # the diagram

## Goodness-of-fit plots

`nm_gof`

Standard diagnostic plots from a NONMEM table: observations vs
population and individual predictions, CWRES, NPDE and IWRES vs time and
predictions.

Help in R: `?nm_gof`

Copy

    gof <- nm_gof(
      file.path(ex, "run2.reg")   # file: run folder, or the NONMEM table itself
    )
    plot(gof)
    nm_gof(file.path(ex, "run2.reg"),
           by = "GENO")           # columns used to split the plots into panels

## Compare two models (LRT, Vuong, Clarke, bootstrap, influence)

`nm_compare` `summary.nm_compare`

Paired comparison of two runs: CWRES and NPDE, individual OFVs, Vuong
and Clarke tests for non-nested models, LRT if nested, subject bootstrap
of dOFV, influence, paired parameters and ETAs; summary() gives a
decision grid.

Help in R: `?nm_compare`, `?summary.nm_compare`

Copy

    cmp <- nm_compare(
      file.path(ex, "run1.reg"),      # run1: first run (folder or path prefix)
      file.path(ex, "run2.reg"),      # run2: second run
      labels     = c("Base", "Final"),   # display names of the two runs
      covariates = c("WT", "GENO")       # covariates examined against the
    )                                    #   per-subject difference of iOFV
    summary(cmp)            # decision grid
    cmp$tests$table

## Comparison plots (gain map, overlaid diagnostics)

`plot.nm_compare`

Where the second model improves or degrades the fit (gain map with
bootstrap band), diagnostics of both runs overlaid, caterpillar of
individual dOFV, parameters and ETAs.

Help in R: `?plot.nm_compare`

Copy

    plot(cmp,                  # result of nm_compare()
         type   = "gain_map",  # figure(s) to draw (see ?plot.nm_compare)
         n_band = 200)         # subject resamples for the confidence band
    plot(cmp, type = c("overlay_gof", "iofv_caterpillar", "eta_scatter"))

## Optimise a base model without iterative runs

`nm_derivatives_ctl` `nm_optimize_base` `plot.nm_optimize_base`

Ranks the possible changes of a base model (residual error, random
effects, covariates) by their expected dOFV, from its outputs (level 1)
and one NONMEM evaluation run (level 2, linearised model).

Help in R: `?nm_derivatives_ctl`, `?nm_optimize_base`,
`?plot.nm_optimize_base`

Copy

    # level 1 only, no NONMEM run
    nm_optimize_base(file.path(ex, "run1.reg"),      # run: base run folder
                     covariates = c("WT", "GENO"),   # covariate columns to test
                     level      = 1)                 # 1 = residual analysis only
    # level 2: one evaluation run, then the linearised model
    p <- nm_derivatives_ctl(file.path(ex, "run1.reg"))  # writes run1_deriv.ctl
    nm_run(basename(p),   # ctl: the control stream written above
           nmfe,          # nmfe_path
           dirname(p))    # analysis_dir: the folder of that control stream
    o <- nm_optimize_base(file.path(ex, "run1.reg"), covariates = c("WT", "GENO"))
    o
    plot(o)

## Covariate effects on ETAs (Wasserstein)

`nm_wasserstein`

Compares the distribution of an ETA between covariate groups with the
Wasserstein distance and a permutation test, with the Wasserstein R2
between runs (Simon et al., Front Pharmacol 2026).

Help in R: `?nm_wasserstein`

Copy

    w <- nm_wasserstein(
      c(BASE  = file.path(ex, "run1.reg"),   # runs: run folders, named
        FINAL = file.path(ex, "run2.reg")),  #   (the names label the results)
      covariates = c("WT", "GENO"),          # covariate columns compared
      eta        = "CL"                      # ETA analysed: a parameter name,
    )                                        #   "ETA(n)" or an index
    w
    plot(w)                 # type = "ecdf" (default)
    plot(w, type = "r2w")

## VPC

`nm_vpc_ctl` `nm_vpc`

Writes and runs the \_orig and \_sim control streams of a VPC from an
estimated run, then draws observed percentiles against the confidence
intervals of simulated percentiles.

Help in R: `?nm_vpc_ctl`, `?nm_vpc`

Copy

    nm_vpc_ctl(file.path(ex, "run2.ctl"),  # run: estimation run (ctl, folder or stem)
               nmfe_path = nmfe,     # path of the nmfe launcher
               n_sim     = 500,      # number of simulated replicates
               n_cores   = 2)        # 2 = the _orig and _sim runs at the same time
    cols <- c("NSIM", "ID", "TIME", "AMT", "MDV", "DV", "PRED", "IPRED")
    nm_vpc(file.path(ex, "run2.ctl"),  # run: same as in nm_vpc_ctl()
           col_names = cols)           # names of the $TABLE columns of the .fit files

## pcVPC

`nm_vpc`

Prediction-corrected VPC (Bergstrand 2011) from the same simulation runs
as the VPC.

Help in R: `?nm_vpc`

Copy

    cols <- c("NSIM", "ID", "TIME", "AMT", "MDV", "DV", "PRED", "IPRED")
    nm_vpc(file.path(ex, "run2.ctl"),
           col_names  = cols,   # names of the $TABLE columns of the .fit files
           type       = 2,      # 1 = VPC (default), 2 = prediction-corrected VPC
           bin_method = 3)      # 3 = nominal times (1 manual, 2 quantiles, 4 auto)

## Bootstrap

`nm_bootstrap`

Non-parametric bootstrap by subject: reproducible and optionally
stratified resampling, resume after interruption, percentile confidence
intervals and success rate.

Help in R: `?nm_bootstrap`

Copy

    boot <- nm_bootstrap(
      file.path(ex, "run2.reg", "run2.ctl"),  # ctl: control stream of the final model
      nmfe_path = nmfe,     # path of the nmfe launcher
      n         = 500,      # number of replicates
      seed      = 1,        # seed, for a reproducible resampling
      strata    = "GENO",   # optional column to stratify the resampling on
      n_cores   = "auto"    # replicates at the same time (physical cores - 1)
    )
    boot
    plot(boot)

## Jackknife

`nm_jackknife`

Leave-one-subject-out runs ranked by \|dOFV\|, with relative and
standardised THETA shifts. A descriptive influence diagnostic, not an
outlier test.

Help in R: `?nm_jackknife`

Copy

    jk <- nm_jackknife(
      file.path(ex, "run2.reg", "run2.ctl"),  # ctl: control stream, already run
      nmfe_path = nmfe,     # path of the nmfe launcher
      n_cores   = "auto"    # runs at the same time (physical cores - 1)
    )
    jk
    plot(jk, "bar")         # type: "bar" (default), "rank" or "scatter"
    plot(jk, "scatter")

## Likelihood profile

`nm_profile`

Fixes a THETA at each value of a grid, re-runs the model and reads the
OFV; the confidence interval is where dOFV crosses 3.84.

Help in R: `?nm_profile`

Copy

    prof <- nm_profile(
      file.path(ex, "run2.reg", "run2.ctl"),  # ctl: reference control stream
      nmfe,                  # nmfe_path: path of the nmfe launcher
      param   = "CL",        # parameter profiled: $THETA label or index
      from    = 4.6,         # grid of fixed values: first value,
      to      = 5.7,         #   last value,
      by      = 0.1,         #   and step
      n_cores = "auto"       # grid values run at the same time
    )
    prof
    plot(prof, type = "dofv")   # "ofv" or "dofv" (with the CI threshold line)

## Export a figure

`nm_export_plot`

Saves a ggplot to JPEG, PNG, PDF or TIFF according to the file
extension, creating the folder if needed.

Help in R: `?nm_export_plot`

Copy

    g <- nm_plot_dv(dat$data)
    nm_export_plot(g,                                       # plot: the figure
                   file.path(ex, "figures", "dv_time.png"), # file: type from the extension
                   width  = 8,                              # width in inches
                   height = 5)                              # height in inches

## Diagnostics in Word

`nm_gof`

Writes all goodness-of-fit plots of a run to a Word document next to its
table.

Help in R: `?nm_gof`

Copy

    nm_gof(file.path(ex, "run2.reg"),  # file: run folder (or table)
           save_docx = TRUE)            # write all the plots to a Word file
                                        #   in the analysis folder
