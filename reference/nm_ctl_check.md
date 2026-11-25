# Static check of a NONMEM control stream

Looks for inconsistencies that can be detected without running NM-TRAN:
`THETA(n)`, `ETA(n)` and `EPS(n)`/`ERR(n)` used without a matching
`$THETA`, `$OMEGA` or `$SIGMA` declaration (and declared parameters that
are never used), data items used in the code but absent from `$INPUT`,
an active `RATE` with first-order absorption, and the requirements of
the M3 method (BQL records not ignored, LLOQ defined, Laplace
estimation). It also reports four errors refused by NM-TRAN: `TOL`
missing with a general ADVAN (6, 8, 9, 13-18), an unknown `METHOD` (FO
is written `METHOD=ZERO`), a `THETA` initial estimate of 0 that is not
fixed, a `$TABLE` item written after an option of the same line (items
written on a line below the options are accepted), and `Y` with
`EPS`/`ERR` defined only in `IF` blocks testing `ICALL.EQ.4` (simulation
code for NM-TRAN).

## Usage

``` r
nm_ctl_check(ctl, reference = NULL)
```

## Arguments

- ctl:

  Path to a control stream, or its lines (character vector or `nm_ctl`
  object).

- reference:

  Optional original control stream (same forms as `ctl`). When given,
  the line differences from `reference` to `ctl` are added.

## Value

An object of class `nm_ctl_check`: a data frame with columns `level`
(`"ERROR"`, `"WARNING"` or `"OK"`), `category` and `message`. With
`reference`, attribute `"diff"` holds the line differences (`" "`
unchanged, `"- "` removed, `"+ "` added).

## Details

Declarations are counted over all `$THETA`, `$OMEGA` and `$SIGMA`
records, including values written on the record line, `BLOCK(n)`,
`DIAGONAL(n)`, `SAME(k)` and `$ABBR REPLACE` lists and ranges
(`ETA(4,6,8)`, `ETA(,3:32)`).

This check does not replace NM-TRAN: it catches common mistakes early.

## Examples

``` r
ctl <- nm_write_ctl()
nm_ctl_check(ctl)
#> [OK] Check - No obvious inconsistency found by the static checks.
nm_ctl_check(nm_ctl_mu(ctl), reference = ctl)
#> [OK] Check - No obvious inconsistency found by the static checks.
#> 
#> Differences from reference: 6 line(s) removed, 9 added
#> - TVCL = THETA(1)
#> - CL = TVCL*EXP(ETA(1))
#> - TVV = THETA(2)
#> - V = TVV*EXP(ETA(2))
#> - TVKA = THETA(3)
#> - KA = TVKA*EXP(ETA(3))
#> + ; MU referencing - former definition: TVCL = THETA(1)
#> + MU_1 = LOG(THETA(1))
#> + CL = EXP(MU_1 + ETA(1))
#> + ; MU referencing - former definition: TVV = THETA(2)
#> + MU_2 = LOG(THETA(2))
#> + V = EXP(MU_2 + ETA(2))
#> + ; MU referencing - former definition: TVKA = THETA(3)
#> + MU_3 = LOG(THETA(3))
#> + KA = EXP(MU_3 + ETA(3))
```
