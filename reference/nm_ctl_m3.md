# Implement the M3 method for BQL data in a NONMEM control stream

Rewrites `$ERROR` so that quantifiable observations keep their usual
likelihood and BQL observations contribute the probability
`P(DV < LLOQ) = PHI((LLOQ - IPRED)/SD)` (Beal 2001, method M3), with
`F_FLAG`. The dataset is not modified: comments listing its requirements
are inserted after `$DATA` (BQL flag coded 0/1, BQL records kept, EVID=0
and MDV=0 on BQL records). `LAPLACE` is added to `$ESTIMATION` (EM,
Bayesian and FO methods are replaced by FOCE-I Laplace), and `M3LLOQ`,
`M3SD`, `M3PBQL` are added to the first `$TABLE`.

## Usage

``` r
nm_ctl_m3(
  ctl,
  bql = "BQL",
  lloq = NULL,
  lloq_col = NULL,
  method_col = NULL,
  method_lloq = NULL,
  residual_sd = NULL,
  npde_mode = FALSE,
  file = NULL,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- ctl:

  Path to a control stream, or its lines (character vector or `nm_ctl`
  object).

- bql:

  Name of the BQL flag data item (0 = quantifiable, 1 = BQL).

- lloq:

  Constant LLOQ, on the DV scale.

- lloq_col:

  Name of a dataset column holding the LLOQ.

- method_col:

  Name of a dataset column holding a numeric code of the analytical
  method.

- method_lloq:

  Named numeric vector of LLOQs by method code, e.g.
  `c("1" = 0.1, "2" = 0.05)`.

- residual_sd:

  NONMEM expression of the residual SD on the DV scale. `NULL`
  (default): inferred from the `Y =` line for additive, proportional,
  combined and `IPRED + W*EPS(1)` models.

- npde_mode:

  Logical. Add `NPDE_MODE` and `DV_LOQ` handling for the NONMEM NPDE
  diagnostics.

- file:

  Output path. If `NULL` (default), nothing is written.

- overwrite:

  Logical. Replace an existing `file`.

- verbose:

  Logical. Print a message when the file is written.

## Value

An object of class `nm_ctl` (control stream lines), invisibly when
`file` is given.

## Details

The LLOQ comes from exactly one of: a constant (`lloq`), a dataset
column (`lloq_col`), or a code of the analytical method (`method_col`
and `method_lloq`), for studies pooling assays with different LLOQs.

## References

Beal SL. Ways to fit a PK model with some data below the quantification
limit. J Pharmacokinet Pharmacodyn. 2001;28(5):481-504.

## Examples

``` r
ctl <- nm_write_ctl(input = c("ID", "TIME", "DV", "AMT", "EVID", "MDV", "CMT", "BQL"),
                    error = "combined")
nm_ctl_m3(ctl, lloq = 0.1)
#> $PROBLEM 1cpt EV SUB/TRANS abs=first-order elim=linear
#> $INPUT ID TIME DV AMT EVID MDV CMT BQL
#> $DATA data.csv IGNORE=@
#> 
#> ; =====================================================================
#> ; M3 BQL - DATASET REQUIREMENTS (the dataset is not modified)
#> ; - BQL must be coded 0 = quantifiable, 1 = BQL.
#> ; - Records with BQL = 1 stay in the analysis: do NOT IGNORE them.
#> ; - BQL records are observations: EVID = 0 and MDV = 0.
#> ; - DV must stay numeric on BQL records; its value is not used as an
#> ;   observed concentration in the censored M3 contribution.
#> ; - Constant LLOQ used in the control stream: 0.1.
#> ; - IMPORTANT: this value enters the likelihood of the BQL records.
#> ; - BQL is in $INPUT.
#> ; =====================================================================
#> 
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
#>   ; --- M3 BQL likelihood ---
#>   M3LLOQ = 0.1
#>   M3SD = SQRT(IPRED**2*SIGMA(1,1)+SIGMA(2,2))
#>   IF (M3SD.LE.0) M3SD = 1.0E-12
#>   M3Z = (M3LLOQ - IPRED)/M3SD
#>   M3PBQL = PHI(M3Z) + 1e-30
#>   M3CENS = 0
#>   IF (BQL.EQ.1) M3CENS = 1
#>   IF (ICALL.EQ.4) M3CENS = 0
#>   IF (M3CENS.EQ.1) THEN
#>     F_FLAG = 1
#>     Y = M3PBQL
#>     MDVRES = 1
#>   ELSE
#>     F_FLAG = 0
#>     Y     = IPRED + IPRED*ERR(1) + ERR(2)
#>   ENDIF
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
#> $ESTIMATION METHOD=COND INTERACTION MAXEVAL=9999 SIGDIGITS=3 PRINT=5 LAPLACE
#> $COVARIANCE PRINT=E MATRIX=S
#> 
#> $TABLE ID TIME AMT EVID MDV CMT DV PRED IPRED RES IRES WRES IWRES CWRES NPDE M3LLOQ M3SD M3PBQL ESAMPLE=300 SEED=12334 NOAPPEND NOPRINT ONEHEADER FILE=sdtab
#> $TABLE ID ETA1 ETA2 ETA3 FIRSTONLY NOPRINT NOAPPEND ONEHEADER FILE=patab
```
