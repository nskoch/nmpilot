# Suggest structural PK models from a NONMEM dataset

Analyses a NONMEM dataset (rich or sparse data) without running NONMEM
and ranks candidate structural models: number of compartments (1 to 3),
input model (IV bolus or infusion, first-order, first-order with lag
time, zero-order with estimated duration, transit compartments),
elimination (linear or Michaelis-Menten) and residual error model
(additive, proportional or combined). Each candidate comes with initial
values for THETA, IIV and SIGMA in the parameter names of
[`nm_write_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_write_ctl.md),
the parameters likely to be hard to estimate and, for linear models, the
precision expected from the sampling design
([`nm_suggest_design()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_design.md)).

## Usage

``` r
nm_suggest_model(
  data,
  cols = NULL,
  route = c("auto", "EV", "IV"),
  lloq = NULL,
  comment_char = c("#", "@"),
  sep = ",",
  dose_conc_factor = 1,
  min_obs_rich = 6L,
  n_bins = 10L,
  max_cpts = 3L,
  n_candidates = 6L,
  max_profiles = 60L,
  design_eval = TRUE,
  seed = 1L,
  verbose = TRUE
)

# S3 method for class 'nm_model_suggestion'
print(x, ...)
```

## Arguments

- data:

  Path to a NONMEM CSV file, or a data.frame.

- cols:

  Optional named list mapping roles to column names: `id`, `time`, `dv`,
  `amt`, `evid`, `mdv`, `cmt`, `rate`, `addl`, `ii`, `ss`, `blq`.
  Missing roles are guessed from standard names.

- route:

  `"auto"` (default: from the CMT of doses and observations, then from
  the shape of the profiles), `"EV"` or `"IV"`.

- lloq:

  Optional LLOQ (DV scale). Concentrations below it are BLQ. Without it,
  the lowest quantified concentration is used as LLOQ for the BLQ
  samples.

- comment_char:

  Characters marking ignored lines (first field).

- sep:

  Field separator of the CSV file.

- dose_conc_factor:

  Unit factor between amount/volume and the DV units, such that CL =
  factor \* Dose / AUC: 1 (default) when DV is in dose units per volume
  unit (mg and mg/L), 1000 for doses in mg and DV in ng/mL (volumes then
  in L). It is passed to
  [`nm_write_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_write_ctl.md)
  by `from`, which writes the scale `S = V/factor`, so that the initial
  values match the control stream.

- min_obs_rich:

  Minimum number of quantified concentrations in a dosing interval for
  it to be analysed individually (default 6).

- n_bins:

  Number of time-after-dose bins for the pooled analysis.

- max_cpts:

  Maximum number of compartments tested (1 to 3).

- n_candidates:

  Number of candidate models returned.

- max_profiles:

  Maximum number of rich intervals analysed individually (a reproducible
  random subset is drawn above it).

- design_eval:

  If `TRUE` (default), the expected precision (RSE) and eta-shrinkage of
  each linear candidate are computed for the sampling design of the
  dataset with
  [`nm_suggest_design()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_design.md).

- seed:

  Seed of the subset of intervals; the user's random state is restored.

- verbose:

  Print progress messages.

- x:

  An `nm_model_suggestion` object.

- ...:

  Ignored.

## Value

An object of class `nm_model_suggestion`, a list with:

- `candidates`:

  one row per candidate: `rank`, `cpts`, `route`, `input`,
  `elimination`, `error`, `score`, `parameters`, `theta_init`, `iiv_cv`,
  `sigma_init`, `max_rse_pct`, `max_shrinkage_pct`, `hard_to_estimate`,
  `notes`.

- `specs`:

  one list per candidate with the arguments used by
  `nm_write_ctl(from = , candidate = )`.

- `calls`:

  the matching
  [`nm_write_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_write_ctl.md)
  calls, as text.

- `support`:

  relative support of each level of the four components.

- `design`:

  expected RSE and shrinkage of each candidate (`nm_suggest_design`
  objects, `NULL` for Michaelis-Menten).

- `alerts`:

  data.frame `type`, `message`.

- `nca`:

  one row per analysed interval.

- `individual`:

  individual estimates (CL, central volume, absorption parameter) used
  for the IIV and by
  [`nm_suggest_covariates()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_covariates.md).

- `evidence`, `profiles`, `data_summary`, `flags`, `notes`, `settings`:

  details of the analysis.

## Details

Method (heuristic; it chooses a starting model and does not replace
model comparison in NONMEM):

1.  Records are split into dosing intervals; `ADDL`/`II` doses are
    expanded. An interval is "rich" when it holds at least
    `min_obs_rich` quantified concentrations, with at least 3 points
    after the peak for extravascular data.

2.  NCA on each rich interval (linear-up log-down AUC, `lambda_z` by
    best adjusted R2, CL = Dose/AUCinf after a first dose, Dose/AUCtau
    at steady state).

3.  Number of compartments: sums of 1 to 3 exponentials fitted to the
    post-peak phase (log scale); BIC weights averaged over intervals. A
    phase is kept only when it is identifiable (rates at least 1.5-fold
    apart, each amplitude at least 1% of the total, fast phase decaying
    within the sampling range).

4.  Input (extravascular data): first-order, first-order + lag,
    zero-order (duration) and gamma (transit) inputs, each convolved
    with the disposition retained for the interval and fitted to the
    whole profile; BLQ samples before the first quantified one enter as
    a censoring penalty (M3-like). BIC weights averaged.

5.  Elimination: dose proportionality (power model on AUC, or on
    dose-normalised concentrations by time bin for sparse data) and
    curvature of the log-concentration decline (concave = saturable).

6.  Residual error: the best structural fit of each interval is refitted
    by iteratively reweighted least squares under a combined variance
    `a^2 + b^2 f^2`; on these common, leverage-corrected (PRESS)
    residuals, the variance models `a^2` (additive), `b^2 f^2`
    (proportional) and `a^2 + b^2 f^2` (combined) are fitted by maximum
    likelihood and compared by BIC, the combined model being kept only
    when its BIC is better by more than 10. Points of leverage above 0.7
    and predictions below twice the LLOQ are left out. The estimates
    give the initial SIGMA values (variances).

7.  IIV: variance of the log of the individual estimates of CL, central
    volume and absorption parameter (median per subject), given as CV%.
    It includes the estimation noise of each individual value and tends
    to overestimate the true IIV, which is the safe side for initial
    values.

8.  With fewer than 3 rich intervals, all observations are pooled by
    time after dose (naive pooled, dose-normalised, binned) and steps 2
    to 5 are run on the binned geometric means (indicative only); the
    residual error model and the IIV cannot be separated in that case
    and default values are used.

The candidate score is the product of the relative supports of the four
components; it sums to 1 over all combinations and is a ranking aid, not
a probability.

Alerts (`$alerts`) flag features that none of the tested structures
describes: double peaks (enterohepatic recycling, several absorption
sites), a clearance that differs between first dose and steady state
(auto-induction, inhibition, time-varying CL), an accumulation different
from the one predicted by the first-dose terminal slope, and a
less-than-proportional exposure (saturable absorption or
bioavailability).

## See also

[`nm_write_ctl()`](https://nskoch.github.io/nmpilot/reference/nm_write_ctl.md),
[`nm_suggest_design()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_design.md),
[`nm_suggest_run()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_run.md),
[`nm_suggest_covariates()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_covariates.md),
[`plot.nm_model_suggestion()`](https://nskoch.github.io/nmpilot/reference/plot.nm_model_suggestion.md)

## Examples

``` r
ex <- system.file("extdata", "example", "example_data.csv", package = "nmpilot")
s <- nm_suggest_model(ex, verbose = FALSE, max_profiles = 15)
s
#> Structural model suggestion
#> ===========================
#> Subjects: 60 | quantified obs: 540 | BLQ: 0 | median obs/subject: 9
#> Design: mixed (15 rich intervals from 15 subjects) | evidence: individual
#> Route: EV | dose levels: 100
#> 
#> Relative support
#>   Compartments:  1 = 0.93, 2 = 0.07, 3 = 0.00 
#>   Input:         first-order = 0.38, first-order+lag = 0.16, zero-order-d1 = 0.31, transit = 0.16 
#>   Elimination:   linear = 0.95, michaelis-menten = 0.05 
#>   Error model:   additive = 0.00, proportional = 1.00, combined = 0.00 
#> 
#> Candidate models
#> 1. 1 cpt | first-order | linear | proportional error | score 0.334
#>    THETA: CL=5.72, V=55.6, KA=1.11
#>    IIV (CV%): CL=32, V=35, KA=100 | SIGMA: prop=0.0377
#>    expected precision: max RSE 20%, max eta-shrinkage 5%
#> 2. 1 cpt | zero-order-d1 | linear | proportional error | score 0.269
#>    THETA: CL=5.72, V=55.6, DUR=1.27
#>    IIV (CV%): CL=32, V=35, DUR=46 | SIGMA: prop=0.0377
#>    expected precision: max RSE 21%, max eta-shrinkage 6.1%
#>    hard to estimate: DUR (fewer than 3 samples before the peak)
#>    note: Zero-order duration (ADVAN1/3/11): doses and observations in CMT 1 and RATE = -2 on doses; nm_suggest_run() adapts the dataset.
#> 3. 1 cpt | transit | linear | proportional error | score 0.140
#>    THETA: CL=5.72, V=55.6, MTT=0.759
#>    IIV (CV%): CL=32, V=35, MTT=100 | SIGMA: prop=0.0377
#>    expected precision: max RSE 20%, max eta-shrinkage 4.2%
#>    hard to estimate: MTT (fewer than 3 samples before the peak); transit_n (number of transit compartments rarely identifiable: test 2 to 5, or estimate it (Savic 2007))
#> 4. 1 cpt | first-order+lag | linear | proportional error | score 0.137
#>    THETA: CL=5.72, V=55.6, KA=1.6
#>    IIV (CV%): CL=32, V=35, KA=100 | SIGMA: prop=0.0377
#>    expected precision: max RSE 57%, max eta-shrinkage 5.3%
#>    hard to estimate: ALAG1 (expected RSE 57%)
#>    note: ALAG1 = 0.0832 to add by hand (nm_write_ctl writes no lag time).
#> 5. 2 cpt | first-order | linear | proportional error | score 0.026
#>    THETA: CL=5.72, V2=29, Q=43, V3=12.6, KA=1.11
#>    IIV (CV%): CL=32, V2=35, KA=100 | SIGMA: prop=0.0377
#>    expected precision: max RSE 349%, max eta-shrinkage 9.1%
#>    hard to estimate: V2 (expected RSE 50%); Q (expected RSE 349%); V3 (expected RSE 120%)
#> 6. 2 cpt | zero-order-d1 | linear | proportional error | score 0.021
#>    THETA: CL=5.72, V1=29, Q=43, V2=12.6, DUR=1.27
#>    IIV (CV%): CL=32, V1=35, DUR=46 | SIGMA: prop=0.0377
#>    expected precision: max RSE 350%, max eta-shrinkage 7.9%
#>    hard to estimate: DUR (fewer than 3 samples before the peak); V1 (expected RSE 51%); Q (expected RSE 350%); V2 (expected RSE 120%)
#>    note: Zero-order duration (ADVAN1/3/11): doses and observations in CMT 1 and RATE = -2 on doses; nm_suggest_run() adapts the dataset.
#> 
#> Alerts
#>  ! [double peak] 4 of 15 intervals show a second peak after the maximum: enterohepatic recycling or several absorption sites are not described by the candidate models.
#> 
#> Notes
#>  - 15 rich intervals analysed (random subset of the 60 available, `max_profiles`). 
#>  - Single dose level (or range < 2-fold): saturable elimination cannot be assessed from dose proportionality. 
#>  - SIGMA initial values come from the residuals of individual fits (common reweighted fit, leverage-corrected residuals): they include part of the structural misfit. 
#> 
#> Control stream of candidate 1: nm_write_ctl(from = x, candidate = 1)
s$candidates[, c("rank", "cpts", "input", "elimination", "error", "score")]
#>   rank cpts           input elimination        error score
#> 1    1    1     first-order      linear proportional 0.334
#> 2    2    1   zero-order-d1      linear proportional 0.269
#> 3    3    1         transit      linear proportional 0.140
#> 4    4    1 first-order+lag      linear proportional 0.137
#> 5    5    2     first-order      linear proportional 0.026
#> 6    6    2   zero-order-d1      linear proportional 0.021
nm_write_ctl(from = s, candidate = 1, verbose = FALSE)
#> $PROBLEM 1cpt EV SUB/TRANS abs=first-order elim=linear
#> $INPUT ID TIME AMT DV EVID MDV CMT WT GENO
#> $DATA example_data.csv IGNORE=#
#> $SUBROUTINES ADVAN2 TRANS2
#> $PK
#> TVCL = THETA(1)
#> CL = TVCL*EXP(ETA(1))
#> TVV = THETA(2)
#> V = TVV*EXP(ETA(2))
#> TVKA = THETA(3)
#> KA = TVKA*EXP(ETA(3))
#> S2 = V   ; scale: concentration in the central compartment
#> 
#> $ERROR
#>   IPRED = F
#>   W     = SQRT(IPRED**2*SIGMA(1,1) + SIGMA(2,2))
#>   IF (W.EQ.0) W = 1
#>   IRES  = DV - IPRED
#>   IWRES = IRES/W
#>   Y     = IPRED + IPRED*ERR(1) + ERR(2)
#> 
#> $THETA ; initial estimates and bounds
#>   (0, 5.72)  ; CL
#>   (0, 55.6)  ; V
#>   (0, 1.11)  ; KA
#> 
#> $OMEGA 0.097 ; ETA(CL)
#> $OMEGA 0.116 ; ETA(V)
#> $OMEGA 0.693 ; ETA(KA)
#> 
#> $SIGMA 0.038 ; ERR(1) proportional
#> $SIGMA 0 FIX ; ERR(2) additive
#> 
#> $ESTIMATION METHOD=COND INTERACTION MAXEVAL=9999 SIGDIGITS=3 PRINT=5
#> $COVARIANCE PRINT=E MATRIX=S
#> 
#> $TABLE ID TIME AMT EVID MDV CMT DV PRED IPRED RES IRES WRES IWRES CWRES NPDE ESAMPLE=300 SEED=12334 NOAPPEND NOPRINT ONEHEADER FILE=sdtab
#> $TABLE ID ETA1 ETA2 ETA3 FIRSTONLY NOPRINT NOAPPEND ONEHEADER FILE=patab
```
