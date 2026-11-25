# Write a NONMEM control stream for a population PK model

Builds a complete control stream (`$PROBLEM` to `$TABLE`) for a standard
compartmental PK model: 1 to 3 compartments, IV or extravascular
administration, first-order, zero-order (modelled rate or duration) or
transit-compartment absorption, linear or Michaelis-Menten elimination,
analytic (`ADVAN`/`TRANS`) or ODE solution. The `$INPUT` record follows
the actual column order of the dataset when `input` is given.

## Usage

``` r
nm_write_ctl(
  file = NULL,
  data = "data.csv",
  input = NULL,
  cols = NULL,
  cpts = 1L,
  route = c("EV", "IV"),
  absorption = c("first-order", "zero-order", "zero-order-d1", "transit"),
  transit_n = 2L,
  elimination = c("linear", "michaelis-menten"),
  dynamics = c("analytic", "ode"),
  ode_subroutine = c("ADVAN13", "ADVAN6"),
  error = c("additive", "proportional", "combined"),
  method = c("FOCEI", "FOCE", "FO", "LAPLACE", "LAPLACEI", "IMP", "ITS", "SAEM", "BAYES"),
  theta_init = NULL,
  theta_lower = NULL,
  theta_upper = NULL,
  theta_fix = NULL,
  iiv = NULL,
  iiv_scale = c("cv", "omega"),
  omega_structure = c("diagonal", "block"),
  omega_rho = 0.1,
  sigma_init = c(prop = 0.09, add = 0.01),
  sigma_fix = NULL,
  sigma_structure = c("auto", "diagonal", "block"),
  covariates = NULL,
  drop = NULL,
  ignore = NULL,
  npde_esample = 300L,
  npde_seed = 12334L,
  overwrite = FALSE,
  verbose = TRUE,
  from = NULL,
  candidate = 1L,
  conc_factor = 1
)
```

## Arguments

- file:

  Output path. If `NULL` (default), nothing is written.

- data:

  File name written in `$DATA`.

- input:

  Dataset column names, in file order, or the path to the dataset (its
  first line is read). If `NULL`, `$INPUT` lists the mapped columns
  only.

- cols:

  Named list mapping NONMEM roles (`id`, `time`, `dv`, `amt`, `evid`,
  `mdv`, `cmt`, `rate`) to column names. Missing roles are guessed from
  `input`; without `input`, standard names are used.

- cpts:

  Number of compartments (1, 2 or 3).

- route:

  `"EV"` (extravascular) or `"IV"`.

- absorption:

  For `route = "EV"`: `"first-order"`, `"zero-order"` (rate K0
  estimated), `"zero-order-d1"` (duration D1 estimated) or `"transit"`.

- transit_n:

  Number of transit compartments.

- elimination:

  `"linear"` or `"michaelis-menten"`.

- dynamics:

  `"analytic"` or `"ode"`. Analytic is used only when an `ADVAN` exists
  for the model; otherwise ODE is used.

- ode_subroutine:

  `"ADVAN13"` or `"ADVAN6"`, written with `TOL=9`.

- error:

  `"additive"`, `"proportional"` or `"combined"`.

- method:

  Estimation method: `"FO"`, `"FOCE"`, `"FOCEI"`, `"LAPLACE"`,
  `"LAPLACEI"`, `"IMP"`, `"ITS"`, `"SAEM"` or `"BAYES"`. `$COVARIANCE`
  is added for FO, FOCE and Laplace methods. FO is written `METHOD=ZERO`
  (NONMEM does not accept `METHOD=FO`).

- theta_init:

  Named numeric vector of initial values; defaults are used for the
  other parameters.

- theta_lower, theta_upper:

  Named numeric vectors of bounds. The lower bound defaults to 0 and the
  upper bound to none; `NA` removes a bound.

- theta_fix:

  Names of the parameters to fix at their initial value.

- iiv:

  Interindividual variability (exponential ETA). `NULL` (default): an
  ETA on every structural parameter (CV 20%, 30% for absorption
  parameters). A character vector: ETAs on these parameters with default
  values. A named numeric vector: ETAs on these parameters with these
  values. `character(0)`: no ETA.

- iiv_scale:

  `"cv"` (values of `iiv` are CV%, converted exactly with OMEGA =
  log(1 + CV^2)) or `"omega"` (variances).

- omega_structure:

  `"diagonal"` or `"block"` (one `$OMEGA BLOCK`).

- omega_rho:

  Initial correlation of the `$OMEGA BLOCK`.

- sigma_init:

  Named numeric vector `c(prop = , add = )` of initial variances.

- sigma_fix:

  Components to fix: any of `"prop"`, `"add"`.

- sigma_structure:

  `"auto"`, `"diagonal"` or `"block"` (combined error only).

- covariates:

  Covariate names written in the `cotab` table.

- drop:

  Column names written `=DROP` in `$INPUT`. For extravascular
  first-order or transit absorption, a mapped `RATE` column is dropped
  automatically.

- ignore:

  Extra `IGNORE` conditions, e.g. `"BQL.EQN.1"`.

- npde_esample, npde_seed:

  `ESAMPLE` and `SEED` of the `sdtab` table (NPDE).

- overwrite:

  Logical. Replace an existing `file`.

- verbose:

  Logical. Print a message when the file is written.

- from:

  Optional `nm_model_suggestion` object from
  [`nm_suggest_model()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_model.md).
  The structure (`cpts`, `route`, `absorption`, `transit_n`,
  `elimination`, `error`), the initial values (`theta_init`, `iiv`,
  `sigma_init`) and the dataset (`data`, `input`) of candidate
  `candidate` are used for every argument not given explicitly. A lag
  time of the candidate is not written (no `ALAG1` in this function): a
  message gives its value. When the structure is changed explicitly,
  only the initial values whose names exist in the new model are kept.

- candidate:

  Candidate number in `from` (default 1).

- conc_factor:

  Unit factor between amount/volume and the DV units: the scale is
  written `S = V/conc_factor` (e.g. 1000 for doses in mg, volumes in L
  and concentrations in ng/mL). Default 1 (`S = V`). With `from`, the
  `dose_conc_factor` of the suggestion is used, so that its initial
  values (volumes in L, KM in amount/L) match the DV units.

## Value

An object of class `nm_ctl` (a character vector of control stream lines,
printed as text), invisibly when `file` is given.

## Details

Structural parameter names, used in `theta_init`, `theta_lower`,
`theta_upper`, `theta_fix` and `iiv`:

- elimination: `CL` (linear) or `VMAX`, `KM` (Michaelis-Menten);

- central volume: `V` (1 compartment), `V1` (IV or zero-order input, 2-3
  compartments), `V2` (first-order absorption, 2-3 compartments);

- peripheral compartments: `Q`, `V2` / `Q2`, `V3`, `Q3` / `V3`, `Q3`,
  `V4`, `Q4` (names of the matching `ADVAN`);

- absorption: `KA` (first-order), `K0` (zero-order rate, `RATE = -1` in
  the dataset), `DUR` (zero-order duration D1, `RATE = -2` in the
  dataset), `MTT` (mean transit time).

Combinations without an analytic solution (Michaelis-Menten elimination,
transit compartments) are written with `$DES` and `ode_subroutine`, with
`TOL=9` on the `$SUBROUTINES` line (NM-TRAN requires `TOL` for these
subroutines). Residual error uses `ERR(1)` (proportional) and `ERR(2)`
(additive); the component absent from the chosen model is fixed to 0.

## See also

[`nm_ctl_check()`](https://nskoch.github.io/nmpilot/reference/nm_ctl_check.md),
[`nm_ctl_mu()`](https://nskoch.github.io/nmpilot/reference/nm_ctl_mu.md),
[`nm_ctl_iov()`](https://nskoch.github.io/nmpilot/reference/nm_ctl_iov.md),
[`nm_ctl_m3()`](https://nskoch.github.io/nmpilot/reference/nm_ctl_m3.md),
[`nm_suggest_model()`](https://nskoch.github.io/nmpilot/reference/nm_suggest_model.md)

## Examples

``` r
# One-compartment oral model with combined error
nm_write_ctl(error = "combined")
#> $PROBLEM 1cpt EV SUB/TRANS abs=first-order elim=linear
#> $INPUT ID TIME DV AMT EVID MDV CMT
#> $DATA data.csv IGNORE=@
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
#>   (0, 5)  ; CL
#>   (0, 50)  ; V
#>   (0, 1)  ; KA
#> 
#> $OMEGA 0.039 ; ETA(CL)
#> $OMEGA 0.039 ; ETA(V)
#> $OMEGA 0.086 ; ETA(KA)
#> 
#> $SIGMA 0.09 ; ERR(1) proportional
#> $SIGMA 0.01 ; ERR(2) additive
#> 
#> $ESTIMATION METHOD=COND INTERACTION MAXEVAL=9999 SIGDIGITS=3 PRINT=5
#> $COVARIANCE PRINT=E MATRIX=S
#> 
#> $TABLE ID TIME AMT EVID MDV CMT DV PRED IPRED RES IRES WRES IWRES CWRES NPDE ESAMPLE=300 SEED=12334 NOAPPEND NOPRINT ONEHEADER FILE=sdtab
#> $TABLE ID ETA1 ETA2 ETA3 FIRSTONLY NOPRINT NOAPPEND ONEHEADER FILE=patab

# Two-compartment IV model, no IIV on Q, fixed V2, written to a file
f <- tempfile(fileext = ".ctl")
nm_write_ctl(f, route = "IV", cpts = 2,
             theta_init = c(CL = 3, V1 = 20), theta_fix = "V2",
             iiv = c(CL = 30, V1 = 25, V2 = 20), verbose = FALSE)
nm_ctl_check(f)
#> [OK] Check - No obvious inconsistency found by the static checks.

# From a data-driven suggestion
ex <- system.file("extdata", "example", "example_data.csv", package = "nmpilot")
s <- nm_suggest_model(ex, verbose = FALSE, design_eval = FALSE, max_profiles = 15)
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
