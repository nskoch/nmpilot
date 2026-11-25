# MU referencing of a NONMEM control stream

Rewrites log-normal individual parameters, `P = THETA(n)*EXP(ETA(m))` or
`P = TVP*EXP(ETA(m))` with `TVP = THETA(n)`, as `MU_m = LOG(THETA(n))`
and `P = EXP(MU_m + ETA(m))`, as required by the EM and Bayesian methods
(SAEM, IMP, BAYES). The former `TVP` line is kept as a comment.
Covariate effects on the typical value are not handled: only parameters
written exactly in one of the two forms above are transformed.

## Usage

``` r
nm_ctl_mu(
  ctl,
  parameters = NULL,
  mode = c("readable", "linear"),
  file = NULL,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- ctl:

  Path to a control stream, or its lines (character vector or `nm_ctl`
  object).

- parameters:

  Names of the parameters to transform. `NULL` (default): every
  parameter written in a supported form.

- mode:

  `"readable"` (`MU_m = LOG(THETA(n))`, THETAs unchanged) or `"linear"`
  (`MU_m = THETA(n)`, the THETA initial value and bounds are
  log-transformed; requires one value per line in `$THETA`; an initial
  value of 1, whose log is 0, is written 0.0001 because NM-TRAN refuses
  an initial estimate of 0 unless it is fixed).

- file:

  Output path. If `NULL` (default), nothing is written.

- overwrite:

  Logical. Replace an existing `file`.

- verbose:

  Logical. Print a message when the file is written.

## Value

An object of class `nm_ctl` (control stream lines), invisibly when
`file` is given.

## Examples

``` r
ctl <- nm_write_ctl()
nm_ctl_mu(ctl)
#> $PROBLEM 1cpt EV SUB/TRANS abs=first-order elim=linear
#> $INPUT ID TIME DV AMT EVID MDV CMT
#> $DATA data.csv IGNORE=@
#> $SUBROUTINES ADVAN2 TRANS2
#> $PK
#> ; MU referencing - former definition: TVCL = THETA(1)
#> MU_1 = LOG(THETA(1))
#> CL = EXP(MU_1 + ETA(1))
#> ; MU referencing - former definition: TVV = THETA(2)
#> MU_2 = LOG(THETA(2))
#> V = EXP(MU_2 + ETA(2))
#> ; MU referencing - former definition: TVKA = THETA(3)
#> MU_3 = LOG(THETA(3))
#> KA = EXP(MU_3 + ETA(3))
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
#> $SIGMA 0 FIX ; ERR(1) proportional
#> $SIGMA 0.01 ; ERR(2) additive
#> 
#> $ESTIMATION METHOD=COND INTERACTION MAXEVAL=9999 SIGDIGITS=3 PRINT=5
#> $COVARIANCE PRINT=E MATRIX=S
#> 
#> $TABLE ID TIME AMT EVID MDV CMT DV PRED IPRED RES IRES WRES IWRES CWRES NPDE ESAMPLE=300 SEED=12334 NOAPPEND NOPRINT ONEHEADER FILE=sdtab
#> $TABLE ID ETA1 ETA2 ETA3 FIRSTONLY NOPRINT NOAPPEND ONEHEADER FILE=patab
nm_ctl_mu(ctl, parameters = c("CL", "V"), mode = "linear")
#> $PROBLEM 1cpt EV SUB/TRANS abs=first-order elim=linear
#> $INPUT ID TIME DV AMT EVID MDV CMT
#> $DATA data.csv IGNORE=@
#> $SUBROUTINES ADVAN2 TRANS2
#> $PK
#> ; MU referencing - former definition: TVCL = THETA(1)
#> MU_1 = THETA(1)
#> CL = EXP(MU_1 + ETA(1))
#> ; MU referencing - former definition: TVV = THETA(2)
#> MU_2 = THETA(2)
#> V = EXP(MU_2 + ETA(2))
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
#>   1.609438 ; log(CL)
#>   3.912023 ; log(V)
#>   (0, 1)  ; KA
#> 
#> $OMEGA 0.039 ; ETA(CL)
#> $OMEGA 0.039 ; ETA(V)
#> $OMEGA 0.086 ; ETA(KA)
#> 
#> $SIGMA 0 FIX ; ERR(1) proportional
#> $SIGMA 0.01 ; ERR(2) additive
#> 
#> $ESTIMATION METHOD=COND INTERACTION MAXEVAL=9999 SIGDIGITS=3 PRINT=5
#> $COVARIANCE PRINT=E MATRIX=S
#> 
#> $TABLE ID TIME AMT EVID MDV CMT DV PRED IPRED RES IRES WRES IWRES CWRES NPDE ESAMPLE=300 SEED=12334 NOAPPEND NOPRINT ONEHEADER FILE=sdtab
#> $TABLE ID ETA1 ETA2 ETA3 FIRSTONLY NOPRINT NOAPPEND ONEHEADER FILE=patab
```
