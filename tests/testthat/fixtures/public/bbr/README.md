# Public control streams and NONMEM outputs from bbr

Source: R package **bbr**, Metrum Research Group,
<https://github.com/metrumresearchgroup/bbr>, folder `inst/model/nonmem`,
commit `2fded1feff1a7c453c69031beeecb99f80205a6e` (2026-10-23).

Licence: MIT, copyright (c) 2020 Metrum Research Group. The licence text is
in `LICENSE.md` (this folder), as the MIT licence requires.

Files copied unchanged (control stream, `.lst`, and `.ext` when it exists and
is small); runs made with NONMEM 7.4.4:

| Folder            | Files                 | Features                                              |
|-------------------|-----------------------|-------------------------------------------------------|
| `1`               | ctl, lst, ext         | ADVAN2 TRANS2, FOCEI, `$COV`                          |
| `1001`            | ctl, lst              | ADVAN13, `$DES`, `$PRIOR NWPRI`, `MTIME`, BAYES       |
| `200`             | ctl                   | ADVAN13, `$DES`, log-scale covariate model            |
| `acop-fake-bayes` | ctl, lst, ext         | FOCEI > BAYES > IMP chain                             |
| `acop-onlysim`    | ctl, lst              | `$SIMULATION ONLYSIM`                                 |
| `acop-iov`        | mod, lst              | `$ABBR REPLACE ETA(OCC_CL)=ETA(,3:32)`, 62 ETAs, SAME(29) |
| `iovmm`           | mod, lst, ext         | `$MIX` (2 subpopulations), IOV, `$TABLE` on several lines |

Left out: `example2_saemimp` (NONMEM example distributed by ICON), the data
sets, `acop-iov.ext` (943 kB) and the other output files.
