# Add inter-occasion variability to a NONMEM control stream

Adds an inter-occasion (between-occasion) random effect to the selected
parameters with `$ABBR REPLACE` (NONMEM 7.3 or later): one ETA per
parameter and occasion, selected through the occasion data item, and an
`$OMEGA BLOCK(p)` followed by `$OMEGA BLOCK(p) SAME(n_occ - 1)`. A
parameter with IIV becomes `P = TVP*EXP(ETA(m)+ETA(OCC_BOVP))`; a
parameter without IIV (`P = TVP`) becomes `P = TVP*EXP(ETA(OCC_BOVP))`.

## Usage

``` r
nm_ctl_iov(
  ctl,
  parameters,
  n_occ,
  occ = "OCC",
  cv = 20,
  structure = c("independent", "correlated"),
  rho = 0.1,
  data = NULL,
  output_file = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- ctl:

  Path to a control stream, or its lines (character vector or `nm_ctl`
  object).

- parameters:

  Names of the parameters receiving IOV.

- n_occ:

  Number of occasions (at least 2).

- occ:

  Name of the occasion data item; it must be active in `$INPUT` and take
  the values 1 to `n_occ`. Letters and digits only.

- cv:

  Initial IOV, as CV% (one value, or one per parameter, possibly named).

- structure:

  `"independent"` (diagonal block) or `"correlated"` (common initial
  correlation `rho`).

- rho:

  Initial correlation for `structure = "correlated"`.

- data:

  Optional dataset (data frame or CSV path) used to check that `occ`
  holds the integers 1 to `n_occ`.

- output_file:

  Output path. If `NULL` (default), nothing is written. (Formerly
  `file`, still accepted with a warning.)

- overwrite:

  Logical. Replace an existing `file`.

- verbose:

  Logical. Print a message when the file is written.

- ...:

  Former argument name, accepted with a warning until version 1.0.

## Value

An object of class `nm_ctl` (control stream lines), invisibly when
`file` is given.

## Examples

``` r
ctl <- nm_write_ctl(input = c("ID", "TIME", "DV", "AMT", "EVID", "MDV", "CMT", "OCC"))
nm_ctl_iov(ctl, parameters = c("CL", "V"), n_occ = 3)
#> $PROBLEM 1cpt EV SUB/TRANS abs=first-order elim=linear
#> $INPUT ID TIME DV AMT EVID MDV CMT OCC
#> $DATA data.csv IGNORE=@
#> $SUBROUTINES ADVAN2 TRANS2
#> ; inter-occasion variability (IOV)
#> $ABBR REPLACE ETA(OCC_BOVCL)=ETA(4,6,8)
#> $ABBR REPLACE ETA(OCC_BOVV)=ETA(5,7,9)
#> $PK
#> TVCL = THETA(1)
#> CL = TVCL*EXP(ETA(1)+ETA(OCC_BOVCL))
#> TVV = THETA(2)
#> V = TVV*EXP(ETA(2)+ETA(OCC_BOVV))
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
#> ; IOV: 3 occasions, independent structure
#> $OMEGA BLOCK(2)
#> 0.039221
#> 0 0.039221
#> $OMEGA BLOCK(2) SAME(2)
#> $SIGMA 0 FIX ; ERR(1) proportional
#> $SIGMA 0.01 ; ERR(2) additive
#> 
#> $ESTIMATION METHOD=COND INTERACTION MAXEVAL=9999 SIGDIGITS=3 PRINT=5
#> $COVARIANCE PRINT=E MATRIX=S
#> 
#> $TABLE ID TIME AMT EVID MDV CMT DV PRED IPRED RES IRES WRES IWRES CWRES NPDE ESAMPLE=300 SEED=12334 NOAPPEND NOPRINT ONEHEADER FILE=sdtab
#> $TABLE ID ETA1 ETA2 ETA3 FIRSTONLY NOPRINT NOAPPEND ONEHEADER FILE=patab
```
