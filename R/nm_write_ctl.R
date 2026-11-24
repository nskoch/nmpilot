# =============================================================================
# nm_write_ctl.R
# Write a NONMEM control stream for a standard population PK model.
#
# Ported from the Shiny application "Générateur NONMEM .ctl"
# (https://nicolas-simon.shinyapps.io/nonmem_ctl_publi/). The model-building
# helpers below are the application's helpers, renamed to the nmpilot
# conventions (.nm_ctl_*), with English comments in the generated code.
# Transformations of an existing control stream (MU referencing, IOV, M3) and
# the static check live in nm_ctl_transform.R.
# =============================================================================


# =============================================================================
# Internal helpers: model structure
# =============================================================================

# Single source of truth for the model structure: parameter names, central
# compartment, peripheral compartments and absorption parameter.
.nm_ctl_struct <- function(cpts, route, absorption, elimination, transit_n = 2L) {
  linear    <- (elimination == "linear")
  transit   <- (route == "EV" && absorption == "transit")
  first_ord <- (route == "EV" && absorption == "first-order")
  zero_d1   <- (route == "EV" && absorption == "zero-order-d1")
  zero_k0   <- (route == "EV" && absorption == "zero-order")
  # "central" family: dose given directly in the central compartment
  # (bolus or zero-order input), same ADVAN as IV (ADVAN1/3/11)
  central_family <- (route == "IV" || zero_d1 || zero_k0)
  central_cmt <- if (first_ord || transit) 2L else 1L

  if (central_family) {                 # ADVAN1/3/11 naming (TRANS2/TRANS4)
    if (cpts == 1L)      { vc <- "V";  peri <- list() }
    else if (cpts == 2L) { vc <- "V1"; peri <- list(list(vol = "V2", q = "Q",  cmt = 2L)) }
    else                 { vc <- "V1"; peri <- list(list(vol = "V2", q = "Q2", cmt = 2L),
                                                    list(vol = "V3", q = "Q3", cmt = 3L)) }
  } else if (first_ord) {               # ADVAN2/4/12 naming
    if (cpts == 1L)      { vc <- "V";  peri <- list() }
    else if (cpts == 2L) { vc <- "V2"; peri <- list(list(vol = "V3", q = "Q",  cmt = 3L)) }
    else                 { vc <- "V2"; peri <- list(list(vol = "V3", q = "Q3", cmt = 3L),
                                                    list(vol = "V4", q = "Q4", cmt = 4L)) }
  } else {                              # transit compartments (ODE)
    if (cpts == 1L)      { vc <- "V"; peri <- list() }
    else if (cpts == 2L) { vc <- "V"; peri <- list(list(vol = "V2", q = "Q",  cmt = central_cmt + 1L)) }
    else                 { vc <- "V"; peri <- list(list(vol = "V2", q = "Q",  cmt = central_cmt + 1L),
                                                   list(vol = "V3", q = "Q2", cmt = central_cmt + 2L)) }
  }
  peri <- lapply(peri, function(p) {
    p$kcp <- paste0("K", central_cmt, p$cmt); p$kpc <- paste0("K", p$cmt, central_cmt); p
  })

  params <- if (linear) "CL" else c("VMAX", "KM")
  params <- c(params, vc)
  for (p in peri) params <- c(params, p$q, p$vol)
  absorp <- NULL
  if (first_ord)    absorp <- "KA"
  else if (zero_k0) absorp <- "K0"
  else if (zero_d1) absorp <- "DUR"
  else if (transit) absorp <- "MTT"
  params <- c(params, absorp)

  list(cpts = cpts, route = route, absorption = absorption, elimination = elimination,
       linear = linear, transit = transit, zero_d1 = zero_d1, zero_k0 = zero_k0,
       first_ord = first_ord, central_family = central_family,
       central_cmt = central_cmt, vc = vc, peri = peri, params = params,
       absorp = absorp, transit_n = as.integer(transit_n))
}

# $SUBROUTINES line for analytic models, or NULL when no ADVAN fits.
.nm_ctl_advan <- function(st) {
  if (st$central_family) {
    return(switch(st$cpts, "$SUBROUTINES ADVAN1 TRANS2",
                  "$SUBROUTINES ADVAN3 TRANS4", "$SUBROUTINES ADVAN11 TRANS4"))
  }
  if (st$first_ord) {
    return(switch(st$cpts, "$SUBROUTINES ADVAN2 TRANS2",
                  "$SUBROUTINES ADVAN4 TRANS4", "$SUBROUTINES ADVAN12 TRANS4"))
  }
  NULL
}

# Default initial values of the structural THETAs.
.nm_ctl_theta_default <- function(p) switch(p,
  "CL" = 5, "V" = 50, "V1" = 50, "VMAX" = 5, "KM" = 2,
  # the two peripheral compartments of ADVAN12 (Q3/V3, Q4/V4) get different
  # rate constants: equal ones are refused by PREDPP ("K32, OR K42 IS TOO
  # CLOSE TO AN EIGENVALUE")
  "Q" = 5, "Q2" = 5, "Q3" = 5, "Q4" = 2,
  "V2" = 100, "V3" = 200, "V4" = 400,
  "KA" = 1, "K0" = 1, "MTT" = 2, "DUR" = 2, 1)

# Absorption parameters get a larger default IIV.
.nm_ctl_is_absorption <- function(p) p %in% c("KA", "MTT", "K0", "DUR")

.nm_ctl_fmt_num <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x)) return(as.character(x))
  format(x, trim = TRUE, scientific = FALSE)
}

# Compact format of initial variances: at most 3 decimals.
.nm_ctl_fmt_var <- function(x, digits = 3L) {
  x <- as.numeric(x)
  if (!is.finite(x)) return(as.character(x))
  format(round(x, digits), trim = TRUE, scientific = FALSE, nsmall = 0)
}

# Exact conversion of a CV% of a log-normal parameter P = TVP*EXP(ETA)
# into the variance OMEGA = log(1 + CV^2).
.nm_ctl_cv_to_omega <- function(cv_pct) log(1 + (as.numeric(cv_pct) / 100)^2)

# $THETA record. spec: list per parameter with init, lower, upper, fix.
.nm_ctl_theta_lines <- function(params, spec) {
  out <- "$THETA ; initial estimates and bounds"
  for (p in params) {
    s <- spec[[p]]
    if (isTRUE(s$fix)) {
      txt <- sprintf("  %s FIX  ; %s", .nm_ctl_fmt_num(s$init), p)
    } else if (is.na(s$lower) && is.na(s$upper)) {
      txt <- sprintf("  %s  ; %s", .nm_ctl_fmt_num(s$init), p)
    } else {
      lo <- if (is.na(s$lower)) "-INF" else .nm_ctl_fmt_num(s$lower)
      txt <- if (is.na(s$upper))
        sprintf("  (%s, %s)  ; %s", lo, .nm_ctl_fmt_num(s$init), p)
      else
        sprintf("  (%s, %s, %s)  ; %s", lo, .nm_ctl_fmt_num(s$init),
                .nm_ctl_fmt_num(s$upper), p)
    }
    out <- c(out, txt)
  }
  out
}

# $OMEGA records, diagonal or one BLOCK with a common initial correlation.
.nm_ctl_omega_lines <- function(eta_params, var_init, structure = "diagonal", rho = 0.1) {
  n <- length(eta_params)
  if (n == 0L) return(character(0))
  if (structure == "diagonal" || n == 1L)
    return(sprintf("$OMEGA %s ; ETA(%s)", vapply(var_init, .nm_ctl_fmt_var, character(1)), eta_params))
  L <- diag(var_init, n)
  for (j in 2:n) for (i in 1:(j - 1)) L[j, i] <- rho * sqrt(var_init[i] * var_init[j])
  c(sprintf("$OMEGA BLOCK(%d)", n),
    vapply(seq_len(n), function(i) paste(vapply(L[i, 1:i], .nm_ctl_fmt_var, character(1)),
                                         collapse = " "), character(1)),
    paste0("; ETAs: ", paste(eta_params, collapse = ", ")))
}

# $SIGMA records. ERR(1) proportional, ERR(2) additive; a component absent
# from the error model is fixed to 0. BLOCK only for a combined model.
.nm_ctl_sigma_lines <- function(error, structure, init_prop, init_add, fix_prop, fix_add) {
  prop_inactive <- (error == "additive")
  add_inactive  <- (error == "proportional")
  prop_txt <- if (prop_inactive) "0 FIX" else trimws(paste(.nm_ctl_fmt_var(init_prop), if (fix_prop) "FIX" else ""))
  add_txt  <- if (add_inactive)  "0 FIX" else trimws(paste(.nm_ctl_fmt_var(init_add),  if (fix_add)  "FIX" else ""))

  # FIX inside a SIGMA BLOCK would fix the whole block: use diagonal instead.
  fixed_active <- (!prop_inactive && fix_prop) || (!add_inactive && fix_add)
  if (structure == "block" && error == "combined" && !fixed_active) {
    return(c("$SIGMA BLOCK(2)", .nm_ctl_fmt_var(init_prop),
             paste("0", .nm_ctl_fmt_var(init_add)),
             "; [1,1] = proportional (ERR1), [2,2] = additive (ERR2)"))
  }
  note <- NULL
  if (structure == "block" && fixed_active)
    note <- "; SIGMA BLOCK replaced by diagonal: an active ERR component is FIX"
  else if (structure == "block" && error != "combined")
    note <- "; SIGMA BLOCK not applicable: a single active ERR component"
  c(sprintf("$SIGMA %s ; ERR(1) proportional", prop_txt),
    sprintf("$SIGMA %s ; ERR(2) additive", add_txt), note)
}

# $PK: individual parameters and scale; micro-constants only for ODE models
# (the TRANS subroutine computes them for analytic models).
.nm_ctl_pk_lines <- function(st, use_ode, eta_map, rate_mapped, conc_factor = 1) {
  tdx <- stats::setNames(seq_along(st$params), st$params)
  out <- "$PK"
  for (p in st$params) {
    out <- c(out, sprintf("TV%s = THETA(%d)", p, tdx[[p]]))
    out <- c(out, if (p %in% names(eta_map))
      sprintf("%s = TV%s*EXP(ETA(%d))", p, p, eta_map[[p]]) else sprintf("%s = TV%s", p, p))
  }
  out <- c(out, if (conc_factor == 1)
    sprintf("S%d = %s   ; scale: concentration in the central compartment", st$central_cmt, st$vc) else
    sprintf("S%d = %s/%s   ; scale: concentration in DV units (amount/volume x %s)",
            st$central_cmt, st$vc, format(conc_factor, scientific = FALSE), format(conc_factor, scientific = FALSE)))
  if (use_ode) {
    if (st$linear) out <- c(out, sprintf("K = CL/%s", st$vc))
    for (p in st$peri) out <- c(out, sprintf("%s = %s/%s", p$kcp, p$q, st$vc),
                                sprintf("%s = %s/%s", p$kpc, p$q, p$vol))
    if (st$transit) {
      N <- max(1L, st$transit_n)
      out <- c(out, sprintf("KTR = %d/MTT   ; (N+1)/MTT (Savic 2007): depot + %d transit compartments",
                            N + 1L, N))
    }
  }
  no_rate <- if (!rate_mapped)
    "; WARNING: no RATE column mapped -> add it to the dataset and to $INPUT" else NULL
  # Zero-order input with modelled duration D1 (analytic and ODE).
  if (st$zero_d1)
    out <- c(out, "D1 = DUR   ; modelled zero-order duration (DUR = THETA)",
             "; REQUIRED: RATE = -2 on dose records, IN THE DATASET (data item RATE)", no_rate)
  # Zero-order input with modelled rate R1 (analytic and ODE): the input
  # lasts AMT/K0 and stops when the dose is exhausted.
  if (st$zero_k0)
    out <- c(out, sprintf("R%d = K0   ; modelled zero-order rate (K0 = THETA)", st$central_cmt),
             "; REQUIRED: RATE = -1 on dose records, IN THE DATASET (data item RATE)", no_rate)
  out
}

.nm_ctl_model_lines <- function(st) {
  m <- "$MODEL"
  if (st$first_ord || st$transit) m <- c(m, "COMP=(DEPOT)", "COMP=(CENTRAL)")
  else                            m <- c(m, "COMP=(CENTRAL)")
  for (p in st$peri) m <- c(m, sprintf("COMP=(PERIPH%d)", p$cmt - st$central_cmt))
  if (st$transit) for (k in seq_len(max(1L, st$transit_n))) m <- c(m, sprintf("COMP=(TRANSIT%d)", k))
  m
}

.nm_ctl_des_lines <- function(st) {
  c_i <- st$central_cmt
  rhs <- list(); prelude <- character(0)
  add_rhs <- function(i, term) rhs[[as.character(i)]] <<- c(rhs[[as.character(i)]], term)

  # input into the central compartment (depot / transit). Zero-order inputs
  # (RATE=-1 with R1, RATE=-2 with D1) are handled by PREDPP: nothing in $DES.
  if (st$first_ord) {
    add_rhs(1L, "-KA*A(1)"); add_rhs(c_i, "KA*A(1)")
  } else if (st$transit) {
    N <- max(1L, st$transit_n); tstart <- c_i + length(st$peri) + 1L
    add_rhs(1L, "-KTR*A(1)")
    add_rhs(tstart, sprintf("KTR*(A(1) - A(%d))", tstart))
    if (N > 1) for (k in (tstart + 1):(tstart + N - 1)) add_rhs(k, sprintf("KTR*(A(%d) - A(%d))", k - 1, k))
    add_rhs(c_i, sprintf("KTR*A(%d)", tstart + N - 1))
  }
  if (st$linear) add_rhs(c_i, sprintf("- K*A(%d)", c_i))
  else {
    prelude <- c(sprintf("CC = A(%d)/%s", c_i, st$vc), "VMEL = VMAX*CC/(KM + CC)")
    add_rhs(c_i, "- VMEL")
  }
  for (p in st$peri) {
    add_rhs(p$cmt, sprintf("%s*A(%d) - %s*A(%d)", p$kcp, c_i, p$kpc, p$cmt))
    add_rhs(c_i,   sprintf("- %s*A(%d) + %s*A(%d)", p$kcp, c_i, p$kpc, p$cmt))
  }
  d <- c("$DES", prelude)
  for (i in sort(as.integer(names(rhs))))
    d <- c(d, sprintf("DADT(%d) = %s", i, paste(rhs[[as.character(i)]], collapse = " ")))
  d
}


# =============================================================================
# Internal helpers: $INPUT / $DATA
# =============================================================================

# IGNORE character from the first column name (#ID, @ID, cID...).
.nm_ctl_ignore_char <- function(headers, default = "@") {
  if (is.null(headers) || !length(headers)) return(default)
  first <- headers[1]
  if (is.na(first) || !nzchar(first) || toupper(first) == "ID") return(default)
  if (nchar(first) >= 2L && toupper(substring(first, 2L)) == "ID") return(substr(first, 1L, 1L))
  default
}

# Extra IGNORE conditions: "BQL.EQN.1" -> "IGNORE=(BQL.EQN.1)".
.nm_ctl_ignore_conditions <- function(x) {
  if (is.null(x) || !length(x)) return(character(0))
  z <- trimws(unlist(strsplit(paste(as.character(x), collapse = "\n"), "[\r\n;]+")))
  z <- z[nzchar(z)]
  vapply(z, function(a) {
    a <- trimws(sub("^IGNORE\\s*=\\s*", "", a, ignore.case = TRUE))
    if (!grepl("^\\(.*\\)$", a)) a <- paste0("(", a, ")")
    paste0("IGNORE=", a)
  }, character(1), USE.NAMES = FALSE)
}

# Guess the NONMEM role of dataset columns from their names.
.nm_ctl_guess_cols <- function(nms) {
  clean <- sub("^[#@;c](?=id$)", "", tolower(nms), perl = TRUE)
  pick <- function(patterns) {
    for (pt in patterns) {
      j <- which(grepl(pt, clean, perl = TRUE))
      if (length(j)) return(nms[j[1]])
    }
    NULL
  }
  list(
    id   = pick(c("^id$", "^subject$", "^subj$", "^usubjid$", "pat(id)?$")),
    time = pick(c("^time$", "^times?$", "^t$", "time_?(min|hr|h|hour|sec|s)$")),
    dv   = pick(c("^dv$", "^y$", "^obs$", "conc", "concentration", "pk_?obs")),
    amt  = pick(c("^amt$", "dose(amt)?$", "^amount$")),
    evid = pick(c("^evid$", "^event$")),
    mdv  = pick(c("^mdv$")),
    cmt  = pick(c("^cmt$", "comp(artment)?$")),
    rate = pick(c("^rate$", "inf(us(e|ion))?_?rate$", "k0$"))
  )
}

# Positional $INPUT following the actual column order of the dataset.
.nm_ctl_input_line <- function(headers, cols, ignore_char = "@", drop_cols = NULL) {
  roles <- c(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID",
             mdv = "MDV", cmt = "CMT", rate = "RATE")
  mapped <- function(r) !is.null(cols[[r]]) && length(cols[[r]]) && nzchar(cols[[r]][1])
  if (is.null(headers) || !length(headers)) {
    inp <- unname(roles[vapply(names(roles), mapped, logical(1))])
    return(if (length(inp)) paste("$INPUT", paste(inp, collapse = " ")) else NULL)
  }
  name2label <- list()
  for (r in names(roles)) if (mapped(r)) name2label[[cols[[r]][1]]] <- roles[[r]]
  labels <- vapply(headers, function(h) {
    lab <- if (!is.null(name2label[[h]])) name2label[[h]] else NULL
    if (!is.null(drop_cols) && h %in% drop_cols)
      return(sprintf("%s=DROP", if (is.null(lab)) h else lab))
    if (!is.null(lab)) return(lab)
    if (nzchar(ignore_char) && startsWith(h, ignore_char) && toupper(substring(h, 2L)) == "ID")
      return("ID")
    h
  }, character(1), USE.NAMES = FALSE)
  paste("$INPUT", paste(labels, collapse = " "))
}

# Column names from a character vector or from the first line of a file.
.nm_ctl_headers <- function(input) {
  if (is.null(input)) return(NULL)
  if (length(input) == 1L && file.exists(input)) {
    first <- sub("\r$", "", readLines(input, n = 1L, warn = FALSE))
    h <- trimws(strsplit(first, if (grepl(",", first)) "," else "\\s+")[[1]])
    h <- gsub('^"|"$', "", h)
    return(h[nzchar(h)])
  }
  as.character(input)
}


# =============================================================================
# Internal helpers: reading and writing control streams
# =============================================================================

# Control stream lines from a file path, a single string or a vector of lines.
.nm_ctl_lines <- function(ctl) {
  if (inherits(ctl, "nm_ctl")) return(as.character(unclass(ctl)))
  if (!is.character(ctl) || !length(ctl)) stop("`ctl` must be a file path or a character vector of lines.")
  if (length(ctl) == 1L && !grepl("\n", ctl, fixed = TRUE) && file.exists(ctl))
    return(sub("\r$", "", readLines(ctl, warn = FALSE)))
  if (length(ctl) == 1L && !grepl("\n", ctl, fixed = TRUE) && !grepl("^\\s*\\$", ctl))
    stop("Control stream file not found: ", ctl)
  x <- unlist(lapply(ctl, function(a) if (nzchar(a)) strsplit(a, "\n", fixed = TRUE)[[1]] else ""),
              use.names = FALSE)
  sub("\r$", "", x)
}

# Return an nm_ctl object, optionally written to `file`.
.nm_ctl_output <- function(lines, file = NULL, overwrite = FALSE, verbose = TRUE) {
  out <- structure(as.character(lines), class = c("nm_ctl", "character"))
  if (is.null(file)) return(out)
  if (file.exists(file) && !isTRUE(overwrite))
    stop("File already exists: ", file, " (use overwrite = TRUE to replace it).")
  writeLines(out, file)
  if (isTRUE(verbose)) message("Control stream written: ", file)
  invisible(out)
}

#' @export
print.nm_ctl <- function(x, ...) {
  cat(unclass(x), sep = "\n")
  invisible(x)
}


# =============================================================================
# Exported function
# =============================================================================

#' Write a NONMEM control stream for a population PK model
#'
#' Builds a complete control stream (`$PROBLEM` to `$TABLE`) for a standard
#' compartmental PK model: 1 to 3 compartments, IV or extravascular
#' administration, first-order, zero-order (modelled rate or duration) or
#' transit-compartment absorption, linear or Michaelis-Menten elimination,
#' analytic (`ADVAN`/`TRANS`) or ODE solution. The `$INPUT` record follows
#' the actual column order of the dataset when `input` is given.
#'
#' Structural parameter names, used in `theta_init`, `theta_lower`,
#' `theta_upper`, `theta_fix` and `iiv`:
#' * elimination: `CL` (linear) or `VMAX`, `KM` (Michaelis-Menten);
#' * central volume: `V` (1 compartment), `V1` (IV or zero-order input,
#'   2-3 compartments), `V2` (first-order absorption, 2-3 compartments);
#' * peripheral compartments: `Q`, `V2` / `Q2`, `V3`, `Q3` / `V3`, `Q3`, `V4`,
#'   `Q4` (names of the matching `ADVAN`);
#' * absorption: `KA` (first-order), `K0` (zero-order rate, `RATE = -1` in
#'   the dataset), `DUR` (zero-order duration D1, `RATE = -2` in the
#'   dataset), `MTT` (mean transit time).
#'
#' Combinations without an analytic solution (Michaelis-Menten elimination,
#' transit compartments) are written with `$DES` and `ode_subroutine`, with
#' `TOL=9` on the `$SUBROUTINES` line (NM-TRAN requires `TOL` for these
#' subroutines).
#' Residual error uses `ERR(1)` (proportional) and `ERR(2)` (additive); the
#' component absent from the chosen model is fixed to 0.
#'
#' @param file Output path. If `NULL` (default), nothing is written.
#' @param data File name written in `$DATA`.
#' @param input Dataset column names, in file order, or the path to the
#'   dataset (its first line is read). If `NULL`, `$INPUT` lists the mapped
#'   columns only.
#' @param cols Named list mapping NONMEM roles (`id`, `time`, `dv`, `amt`,
#'   `evid`, `mdv`, `cmt`, `rate`) to column names. Missing roles are guessed
#'   from `input`; without `input`, standard names are used.
#' @param cpts Number of compartments (1, 2 or 3).
#' @param route `"EV"` (extravascular) or `"IV"`.
#' @param absorption For `route = "EV"`: `"first-order"`, `"zero-order"`
#'   (rate K0 estimated), `"zero-order-d1"` (duration D1 estimated) or
#'   `"transit"`.
#' @param transit_n Number of transit compartments.
#' @param elimination `"linear"` or `"michaelis-menten"`.
#' @param dynamics `"analytic"` or `"ode"`. Analytic is used only when an
#'   `ADVAN` exists for the model; otherwise ODE is used.
#' @param ode_subroutine `"ADVAN13"` or `"ADVAN6"`, written with `TOL=9`.
#' @param error `"additive"`, `"proportional"` or `"combined"`.
#' @param method Estimation method: `"FO"`, `"FOCE"`, `"FOCEI"`,
#'   `"LAPLACE"`, `"LAPLACEI"`, `"IMP"`, `"ITS"`, `"SAEM"` or `"BAYES"`.
#'   `$COVARIANCE` is added for FO, FOCE and Laplace methods. FO is written
#'   `METHOD=ZERO` (NONMEM does not accept `METHOD=FO`).
#' @param theta_init Named numeric vector of initial values; defaults are
#'   used for the other parameters.
#' @param theta_lower,theta_upper Named numeric vectors of bounds. The lower
#'   bound defaults to 0 and the upper bound to none; `NA` removes a bound.
#' @param theta_fix Names of the parameters to fix at their initial value.
#' @param iiv Interindividual variability (exponential ETA). `NULL`
#'   (default): an ETA on every structural parameter (CV 20%, 30% for
#'   absorption parameters). A character vector: ETAs on these parameters
#'   with default values. A named numeric vector: ETAs on these parameters
#'   with these values. `character(0)`: no ETA.
#' @param iiv_scale `"cv"` (values of `iiv` are CV%, converted exactly with
#'   OMEGA = log(1 + CV^2)) or `"omega"` (variances).
#' @param omega_structure `"diagonal"` or `"block"` (one `$OMEGA BLOCK`).
#' @param omega_rho Initial correlation of the `$OMEGA BLOCK`.
#' @param sigma_init Named numeric vector `c(prop = , add = )` of initial
#'   variances.
#' @param sigma_fix Components to fix: any of `"prop"`, `"add"`.
#' @param sigma_structure `"auto"`, `"diagonal"` or `"block"` (combined
#'   error only).
#' @param covariates Covariate names written in the `cotab` table.
#' @param drop Column names written `=DROP` in `$INPUT`. For extravascular
#'   first-order or transit absorption, a mapped `RATE` column is dropped
#'   automatically.
#' @param ignore Extra `IGNORE` conditions, e.g. `"BQL.EQN.1"`.
#' @param npde_esample,npde_seed `ESAMPLE` and `SEED` of the `sdtab` table
#'   (NPDE).
#' @param overwrite Logical. Replace an existing `file`.
#' @param verbose Logical. Print a message when the file is written.
#' @param from Optional `nm_model_suggestion` object from
#'   [nm_suggest_model()]. The structure (`cpts`, `route`, `absorption`,
#'   `transit_n`, `elimination`, `error`), the initial values (`theta_init`,
#'   `iiv`, `sigma_init`) and the dataset (`data`, `input`) of candidate
#'   `candidate` are used for every argument not given explicitly. A lag
#'   time of the candidate is not written (no `ALAG1` in this function): a
#'   message gives its value. When the structure is changed explicitly, only
#'   the initial values whose names exist in the new model are kept.
#' @param candidate Candidate number in `from` (default 1).
#' @param conc_factor Unit factor between amount/volume and the DV units:
#'   the scale is written `S = V/conc_factor` (e.g. 1000 for doses in mg,
#'   volumes in L and concentrations in ng/mL). Default 1 (`S = V`). With
#'   `from`, the `dose_conc_factor` of the suggestion is used, so that its
#'   initial values (volumes in L, KM in amount/L) match the DV units.
#'
#' @return An object of class `nm_ctl` (a character vector of control
#'   stream lines, printed as text), invisibly when `file` is given.
#'
#' @seealso [nm_ctl_check()], [nm_ctl_mu()], [nm_ctl_iov()], [nm_ctl_m3()],
#'   [nm_suggest_model()]
#'
#' @examples
#' # One-compartment oral model with combined error
#' nm_write_ctl(error = "combined")
#'
#' # Two-compartment IV model, no IIV on Q, fixed V2, written to a file
#' f <- tempfile(fileext = ".ctl")
#' nm_write_ctl(f, route = "IV", cpts = 2,
#'              theta_init = c(CL = 3, V1 = 20), theta_fix = "V2",
#'              iiv = c(CL = 30, V1 = 25, V2 = 20), verbose = FALSE)
#' nm_ctl_check(f)
#'
#' # From a data-driven suggestion
#' ex <- system.file("extdata", "example", "example_data.csv", package = "nmpilot")
#' s <- nm_suggest_model(ex, verbose = FALSE, design_eval = FALSE, max_profiles = 15)
#' nm_write_ctl(from = s, candidate = 1, verbose = FALSE)
#' @export
nm_write_ctl <- function(file = NULL,
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
                         method = c("FOCEI", "FOCE", "FO", "LAPLACE", "LAPLACEI",
                                    "IMP", "ITS", "SAEM", "BAYES"),
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
                         conc_factor = 1) {

  # ---- candidate of nm_suggest_model() -------------------------------------------
  from_theta <- FALSE; from_iiv <- FALSE
  if (!is.null(from)) {
    if (!inherits(from, "nm_model_suggestion"))
      stop("`from` must be an nm_model_suggestion object (see nm_suggest_model()).")
    if (length(candidate) != 1L || !candidate %in% seq_along(from$specs))
      stop("`candidate` must be between 1 and ", length(from$specs), ".")
    sp <- from$specs[[candidate]]
    if (missing(cpts)) cpts <- sp$cpts
    if (missing(route)) route <- sp$route
    if (missing(absorption)) absorption <- sp$absorption
    if (missing(transit_n)) transit_n <- sp$transit_n
    if (missing(elimination)) elimination <- sp$elimination
    if (missing(error)) error <- sp$error
    if (missing(theta_init)) { theta_init <- sp$theta; from_theta <- TRUE }
    if (missing(iiv) && missing(iiv_scale)) { iiv <- sp$iiv; from_iiv <- TRUE }
    if (missing(sigma_init)) sigma_init <- sp$sigma
    if (missing(conc_factor)) conc_factor <- from$settings$dose_conc_factor
    df_file <- from$data_summary$data_file
    if (missing(input))
      input <- if (!is.null(df_file) && file.exists(df_file)) df_file else from$data_summary$columns
    if (missing(data) && !is.null(df_file)) data <- basename(df_file)
    if (!is.null(sp$alag) && verbose)
      message("Candidate ", candidate, " has a lag time: add ALAG1 (initial value ", sp$alag,
              ") by hand; nm_write_ctl() writes no lag time.")
  }

  if (!is.numeric(conc_factor) || length(conc_factor) != 1L || !is.finite(conc_factor) || conc_factor <= 0)
    stop("`conc_factor` must be a positive number.")
  route <- match.arg(route); absorption <- match.arg(absorption)
  elimination <- match.arg(elimination); dynamics <- match.arg(dynamics)
  ode_subroutine <- match.arg(ode_subroutine); error <- match.arg(error)
  method <- match.arg(toupper(method[1]), c("FOCEI", "FOCE", "FO", "LAPLACE", "LAPLACEI",
                                            "IMP", "ITS", "SAEM", "BAYES"))
  iiv_scale <- match.arg(iiv_scale); omega_structure <- match.arg(omega_structure)
  sigma_structure <- match.arg(sigma_structure)
  cpts <- as.integer(cpts)
  if (length(cpts) != 1L || !cpts %in% 1:3) stop("`cpts` must be 1, 2 or 3.")
  transit_n <- as.integer(transit_n)
  if (absorption == "transit" && (length(transit_n) != 1L || is.na(transit_n) || transit_n < 1L))
    stop("`transit_n` must be a positive integer.")
  if (route == "IV") absorption <- "first-order"   # ignored for IV

  st <- .nm_ctl_struct(cpts, route, absorption, elimination, transit_n)
  params <- st$params
  # values taken from a suggestion: keep only those of the (possibly changed) model
  if (from_theta) theta_init <- theta_init[names(theta_init) %in% params]
  if (from_iiv) iiv <- iiv[names(iiv) %in% params]

  # ---- checks of user-supplied parameter names -------------------------------
  check_names <- function(x, arg) {
    nm <- if (is.character(x)) x else names(x)
    if (length(x) && (is.null(nm) || any(!nzchar(nm))))
      stop("`", arg, "` must be named with structural parameter names.")
    bad <- setdiff(nm, params)
    if (length(bad))
      stop("`", arg, "`: unknown parameter(s) ", paste(bad, collapse = ", "),
           ". Parameters of this model: ", paste(params, collapse = ", "), ".")
  }
  check_names(theta_init, "theta_init"); check_names(theta_lower, "theta_lower")
  check_names(theta_upper, "theta_upper"); check_names(theta_fix, "theta_fix")
  if (!is.null(iiv) && length(iiv)) check_names(iiv, "iiv")
  if (!is.null(sigma_fix) && length(setdiff(sigma_fix, c("prop", "add"))))
    stop("`sigma_fix` must contain \"prop\" and/or \"add\".")
  if (is.null(names(sigma_init)) || !all(c("prop", "add") %in% names(sigma_init)))
    stop("`sigma_init` must be a named vector c(prop = , add = ).")

  # ---- THETA specification -------------------------------------------------------
  spec <- lapply(params, function(p) {
    init <- if (!is.null(theta_init) && p %in% names(theta_init)) as.numeric(theta_init[[p]])
            else .nm_ctl_theta_default(p)
    lower <- if (!is.null(theta_lower) && p %in% names(theta_lower)) as.numeric(theta_lower[[p]]) else 0
    upper <- if (!is.null(theta_upper) && p %in% names(theta_upper)) as.numeric(theta_upper[[p]]) else NA_real_
    if (!is.finite(init)) stop("Initial value of ", p, " must be finite.")
    if (!is.na(lower) && init <= lower) stop("Initial value of ", p, " must be above its lower bound.")
    if (!is.na(upper) && init >= upper) stop("Initial value of ", p, " must be below its upper bound.")
    list(init = init, lower = lower, upper = upper, fix = p %in% theta_fix)
  })
  names(spec) <- params

  # ---- IIV ---------------------------------------------------------------------
  default_cv <- function(p) if (.nm_ctl_is_absorption(p)) 30 else 20
  if (is.null(iiv)) {
    eta_params <- params
    omega <- vapply(eta_params, function(p) .nm_ctl_cv_to_omega(default_cv(p)), numeric(1))
  } else if (is.character(iiv)) {
    eta_params <- params[params %in% iiv]
    omega <- vapply(eta_params, function(p) .nm_ctl_cv_to_omega(default_cv(p)), numeric(1))
  } else {
    eta_params <- params[params %in% names(iiv)]
    val <- as.numeric(iiv[eta_params])
    if (any(!is.finite(val) | val <= 0)) stop("`iiv` values must be strictly positive.")
    omega <- if (iiv_scale == "cv") .nm_ctl_cv_to_omega(val) else val
  }
  eta_map <- stats::setNames(seq_along(eta_params), eta_params)

  # ---- $SUBROUTINES ------------------------------------------------------------
  advan_line <- NULL; use_ode <- TRUE
  if (dynamics == "analytic" && st$linear) {
    advan_line <- .nm_ctl_advan(st)
    if (!is.null(advan_line)) use_ode <- FALSE
    else advan_line <- paste("$SUBROUTINES", ode_subroutine, "TOL=9 ; ODE: no analytic solution for this model")
  } else {
    # NM-TRAN requires TOL with ADVAN6 and ADVAN13 (error 145)
    advan_line <- trimws(paste("$SUBROUTINES", ode_subroutine, "TOL=9",
                               if (!st$linear) "; ODE required for Michaelis-Menten" else ""))
  }

  # ---- $INPUT / $DATA ------------------------------------------------------------
  headers <- .nm_ctl_headers(input)
  guessed <- if (length(headers)) .nm_ctl_guess_cols(headers) else
    list(id = "ID", time = "TIME", dv = "DV", amt = "AMT", evid = "EVID", mdv = "MDV", cmt = "CMT")
  if (!is.null(cols)) {
    bad <- setdiff(names(cols), c("id", "time", "dv", "amt", "evid", "mdv", "cmt", "rate"))
    if (length(bad)) stop("`cols`: unknown role(s) ", paste(bad, collapse = ", "), ".")
    guessed[names(cols)] <- cols
  }
  cols <- guessed
  if (length(headers)) {
    for (r in names(cols)) if (!is.null(cols[[r]]) && !cols[[r]] %in% headers)
      stop("`cols`: column '", cols[[r]], "' (", r, ") is not in `input`.")
    for (r in c("id", "time", "dv")) if (is.null(cols[[r]]))
      stop("No ", toupper(r), " column found in `input`; give it in `cols`.")
  }
  rate_mapped <- !is.null(cols$rate) && nzchar(cols$rate)
  ignore_char <- .nm_ctl_ignore_char(headers, default = "@")
  auto_drop <- if (route == "EV" && !absorption %in% c("zero-order", "zero-order-d1") && rate_mapped)
    cols$rate else character(0)
  input_line <- .nm_ctl_input_line(headers, cols, ignore_char = ignore_char,
                                   drop_cols = unique(c(drop, auto_drop)))
  data_line <- paste(c(sprintf("$DATA %s IGNORE=%s", data, ignore_char),
                       .nm_ctl_ignore_conditions(ignore)), collapse = " ")

  # ---- $ESTIMATION / $COVARIANCE -------------------------------------------------
  base_cl <- " MAXEVAL=9999 SIGDIGITS=3 PRINT=5"
  est_line <- switch(method,
    # NONMEM knows METHOD=ZERO (or 0), not METHOD=FO (NM-TRAN error 287)
    "FO"       = paste0("$ESTIMATION METHOD=ZERO", if (length(eta_params)) " POSTHOC" else "", base_cl),
    "FOCE"     = paste0("$ESTIMATION METHOD=COND", base_cl),
    "FOCEI"    = paste0("$ESTIMATION METHOD=COND INTERACTION", base_cl),
    "LAPLACE"  = paste0("$ESTIMATION METHOD=COND LAPLACE", base_cl),
    "LAPLACEI" = paste0("$ESTIMATION METHOD=COND LAPLACE INTERACTION", base_cl),
    "IMP"      = "$ESTIMATION METHOD=IMP INTERACTION PRINT=5 ISAMPLE=300 NITER=50",
    "ITS"      = "$ESTIMATION METHOD=ITS INTERACTION PRINT=5 NITER=50",
    "SAEM"     = "$ESTIMATION METHOD=SAEM INTERACTION PRINT=50 NBURN=2000 NITER=1000 ISAMPLE=2",
    "BAYES"    = "$ESTIMATION METHOD=BAYES INTERACTION PRINT=100 NBURN=4000 NITER=10000")
  cov_line <- if (method %in% c("FO", "FOCE", "FOCEI", "LAPLACE", "LAPLACEI"))
    "$COVARIANCE PRINT=E MATRIX=S" else NULL

  # ---- $TABLE (xpose/PsN naming) ---------------------------------------------------
  roles <- c(id = "ID", time = "TIME", amt = "AMT", evid = "EVID", mdv = "MDV", cmt = "CMT", dv = "DV")
  sd_items <- unname(roles[vapply(names(roles), function(r)
    !is.null(cols[[r]]) && nzchar(cols[[r]]), logical(1))])
  sdtab <- paste("$TABLE", paste(c(sd_items, "PRED", "IPRED", "RES", "IRES", "WRES", "IWRES",
                                   "CWRES", "NPDE"), collapse = " "),
                 sprintf("ESAMPLE=%d SEED=%d NOAPPEND NOPRINT ONEHEADER FILE=sdtab",
                         as.integer(npde_esample), as.integer(npde_seed)))
  patab <- if (length(eta_params))
    paste("$TABLE ID", paste0("ETA", seq_along(eta_params), collapse = " "),
          "FIRSTONLY NOPRINT NOAPPEND ONEHEADER FILE=patab") else NULL
  cotab <- if (length(covariates))
    paste("$TABLE ID", paste(unique(covariates), collapse = " "),
          "FIRSTONLY NOPRINT NOAPPEND ONEHEADER FILE=cotab") else NULL

  lines <- c(
    sprintf("$PROBLEM %dcpt %s %s abs=%s elim=%s", cpts, route,
            if (use_ode) "ODE" else "SUB/TRANS", if (route == "IV") "-" else absorption, elimination),
    input_line, data_line, advan_line,
    if (use_ode) .nm_ctl_model_lines(st),
    .nm_ctl_pk_lines(st, use_ode, eta_map, rate_mapped, conc_factor),
    if (use_ode) .nm_ctl_des_lines(st),
    "",
    "$ERROR",
    "  IPRED = F",
    "  W     = SQRT(IPRED**2*SIGMA(1,1) + SIGMA(2,2))",
    "  IF (W.EQ.0) W = 1",
    "  IRES  = DV - IPRED",
    "  IWRES = IRES/W",
    "  Y     = IPRED + IPRED*ERR(1) + ERR(2)",
    "",
    .nm_ctl_theta_lines(params, spec), "",
    .nm_ctl_omega_lines(eta_params, omega, omega_structure, omega_rho),
    if (length(eta_params)) "",
    .nm_ctl_sigma_lines(error, sigma_structure, sigma_init[["prop"]], sigma_init[["add"]],
                        "prop" %in% sigma_fix, "add" %in% sigma_fix), "",
    est_line, cov_line, "",
    sdtab, patab, cotab)

  .nm_ctl_output(lines, file = file, overwrite = overwrite, verbose = verbose)
}
