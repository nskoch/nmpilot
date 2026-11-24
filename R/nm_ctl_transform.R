# =============================================================================
# nm_ctl_transform.R
# Transformations of an existing NONMEM control stream (MU referencing,
# inter-occasion variability, M3 method for BQL data) and static check.
#
# Ported from the Shiny application "Générateur NONMEM .ctl"
# (https://nicolas-simon.shinyapps.io/nonmem_ctl_publi/). Helpers are kept
# self-contained (.nm_ctl_*); merging them with the parsers of
# nm_identify_parameters() / nm_convert_ctl() is left to the interface
# harmonisation pass.
# =============================================================================


# =============================================================================
# Internal helpers: record parsing
# =============================================================================

.nm_ctl_strip_comments <- function(lines) sub(";.*$", "", lines)

# Start and end line of the first record `rec`, or NULL.
.nm_ctl_record_range <- function(lines, rec) {
  i <- grep(paste0("^\\s*\\$", rec, "\\b"), lines, ignore.case = TRUE)
  if (!length(i)) return(NULL)
  nxt <- grep("^\\s*\\$", lines); nxt <- nxt[nxt > i[1]]
  c(i[1], if (length(nxt)) nxt[1] - 1L else length(lines))
}

# Active data items of $INPUT (DROP/SKIP excluded), upper case. A renamed
# item (CONC=DV) is returned under its NONMEM label (DV) and its own name.
.nm_ctl_input_items <- function(lines) {
  r <- .nm_ctl_record_range(lines, "INPUT")
  if (is.null(r)) return(character(0))
  z <- sub("^\\s*\\$INPUT\\w*", "", paste(.nm_ctl_strip_comments(lines[r[1]:r[2]]), collapse = " "),
           ignore.case = TRUE)
  tok <- unlist(strsplit(trimws(z), "[,[:space:]]+"))
  tok <- toupper(tok[nzchar(tok)])
  out <- character(0)
  for (t in tok) {
    p <- strsplit(t, "=", fixed = TRUE)[[1]]
    if (length(p) >= 2L && p[2] %in% c("DROP", "SKIP")) next
    out <- c(out, p)
  }
  unique(out)
}

# Lines of the control stream without the records that hold no abbreviated
# code ($PROBLEM, $INPUT, $THETA, $OMEGA, $TABLE...), so that a title, an
# initial estimate such as "$THETA (18.3 FIXED)" or a $TABLE item is not taken
# for a use in the model code. Record names may be abbreviated (3 letters).
.nm_ctl_code_lines <- function(lines) {
  non_code <- c("PROBLEM", "INPUT", "DATA", "SUBROUTINES", "MODEL", "SIZES",
                "THETA", "THETAI", "THETAR", "THETAP", "THETAPV", "OMEGA", "OMEGAP",
                "OMEGAPD", "SIGMA", "SIGMAP", "SIGMAPD", "PRIOR", "ESTIMATION",
                "COVARIANCE", "TABLE", "SCATTERPLOT", "SIMULATION", "NONPARAMETRIC",
                "ETAS", "PHIS", "LEVEL", "ANNEAL", "MSFI", "SUPER", "WARNINGS",
                "BIND", "FORMAT", "TTDF", "OLKJDF", "OVARF", "SLKJDF", "SVARF",
                "RCOV", "RCOVI", "DEFAULT", "CHAIN", "INCLUDE")
  is_rec <- grepl("^\\s*\\$", lines)
  tok <- toupper(sub("^\\s*\\$([A-Za-z]*).*$", "\\1", lines))
  drop_rec <- vapply(tok, function(t) nchar(t) >= 3L && any(startsWith(non_code, t)),
                     logical(1))
  keep <- logical(length(lines)); cur <- TRUE
  for (k in seq_along(lines)) {
    if (is_rec[k]) cur <- !drop_rec[k]
    keep[k] <- cur
  }
  lines[keep]
}

# Lines of each record whose name starts with `prefix` (record names may be
# abbreviated), comments removed, the record name removed from its first line.
.nm_ctl_records <- function(lines, names) {
  is_rec <- grepl("^\\s*\\$", lines)
  tok <- toupper(sub("^\\s*\\$([A-Za-z]*).*$", "\\1", lines))
  out <- list(); cur <- NULL
  for (k in seq_along(lines)) {
    if (is_rec[k]) {
      if (!is.null(cur)) out[[length(out) + 1L]] <- cur
      cur <- if (tok[k] %in% names) list(name = tok[k], start = k, text = character()) else NULL
      if (is.null(cur)) next
      txt <- sub("^\\s*\\$[A-Za-z]*", "", sub(";.*$", "", lines[k]))
    } else {
      if (is.null(cur)) next
      txt <- sub(";.*$", "", lines[k])
    }
    cur$text <- c(cur$text, txt)
  }
  if (!is.null(cur)) out[[length(out) + 1L]] <- cur
  out
}

# Static checks for errors that NM-TRAN refuses and that were met on control
# streams written by nmpilot: TOL missing with a general ADVAN (error 145),
# unknown METHOD (error 287), THETA initial estimate of 0 not fixed (error
# 260), $TABLE item written after an option of the same line (error 20), Y with EPS/ERR
# defined only in IF blocks testing ICALL.EQ.4 (error 122).
.nm_ctl_nmtran_checks <- function(lines) {
  msg <- list()
  add <- function(cat, txt) msg[[length(msg) + 1L]] <<- c(cat, txt)

  sub_recs <- .nm_ctl_records(lines, c("SUB", "SUBR", "SUBROUTINE", "SUBROUTINES"))
  has_tol_rec <- length(.nm_ctl_records(lines, "TOL")) > 0
  for (r in sub_recs) {
    t <- toupper(paste(r$text, collapse = " "))
    adv <- suppressWarnings(as.integer(sub(".*\\bADVAN\\s*([0-9]+).*", "\\1", t)))
    if (!is.na(adv) && adv %in% c(6L, 8L, 9L, 13L, 14L, 15L, 16L, 17L, 18L) &&
        !grepl("\\bTOL\\s*=", t) && !has_tol_rec)
      add("$SUBROUTINES", paste0("ADVAN", adv, " requires TOL (e.g. TOL=9 in $SUBROUTINES or a $TOL record)."))
  }

  methods <- c("0", "ZERO", "1", "COND", "CONDITIONAL", "HYBRID", "ITS", "IMP", "IMPMAP",
               "SAEM", "BAYES", "DIRECT", "NUTS", "CHAIN")
  for (r in .nm_ctl_records(lines, c("EST", "ESTM", "ESTIMATION", "ESTIMATE"))) {
    t <- toupper(paste(r$text, collapse = " "))
    m <- regmatches(t, gregexpr("\\bMETH(OD)?\\s*=\\s*[A-Z0-9]+", t))[[1]]
    for (v in sub(".*=\\s*", "", m))
      if (!(v %in% methods || (nchar(v) >= 4 && startsWith("CONDITIONAL", v))))
        add("$ESTIMATION", paste0("METHOD=", v, " is not a NONMEM method (FO is written METHOD=ZERO or METHOD=0)."))
  }

  n <- 0L
  for (r in .nm_ctl_records(lines, "THETA")) {
    t <- toupper(paste(r$text, collapse = " "))
    spec <- regmatches(t, gregexpr("\\([^)]*\\)(\\s*FIX(ED)?)?|[-+]?[0-9.]+([ED][-+]?[0-9]+)?(\\s*FIX(ED)?)?", t))[[1]]
    for (sp in spec) {
      n <- n + 1L
      fix <- grepl("FIX", sp)
      p <- trimws(strsplit(gsub("[()]|FIX(ED)?", "", sp), ",")[[1]])
      p <- p[nzchar(p)]
      v <- suppressWarnings(as.numeric(sub("D", "E", if (length(p) >= 2L) p[2] else p[1])))
      if (!fix && is.finite(v) && v == 0)
        add("$THETA", paste0("THETA(", n, "): the initial estimate is 0 and not fixed (refused by NM-TRAN)."))
    }
  }

  # Y with EPS/ERR defined only inside IF blocks that test ICALL.EQ.4: such a
  # block is simulation code for NM-TRAN, which then finds no EPS for the
  # estimation (error 122, "EPS VARIABLES REQUIRED WITH POPULATION DATA").
  # Seen with M3 code written "IF (BQL.EQ.0.OR.ICALL.EQ.4) THEN ... Y = ...
  # ERR ...": accepted once the test of ICALL was moved out of that block.
  for (r in .nm_ctl_records(lines, c("ERROR", "PRED"))) {
    depth <- 0L; icall <- logical(); y_eps_icall <- logical()
    for (ln in toupper(trimws(r$text))) {
      if (!nzchar(ln)) next
      if (grepl("^IF\\s*\\(.*\\)\\s*THEN$", ln)) {
        depth <- depth + 1L
        icall[depth] <- grepl("\\bICALL\\s*\\.EQ\\.\\s*4\\b", ln) || (depth > 1L && icall[depth - 1L])
        next
      }
      if (grepl("^END\\s*IF$", ln)) { depth <- max(0L, depth - 1L); next }
      if (grepl("^ELSE", ln)) {
        if (depth > 0L) icall[depth] <- depth > 1L && icall[depth - 1L]
        next
      }
      one <- grepl("^IF\\s*\\(.*\\)\\s*Y\\s*=", ln)
      if ((one || grepl("^Y\\s*=", ln)) &&
          grepl("\\b(EPS|ERR)\\s*\\(", sub("^.*?\\bY\\s*=", "", ln, perl = TRUE))) {
        in_icall <- (depth > 0L && icall[depth]) ||
          (one && grepl("\\bICALL\\s*\\.EQ\\.\\s*4\\b", sub("\\bY\\s*=.*$", "", ln)))
        y_eps_icall <- c(y_eps_icall, in_icall)
      }
    }
    if (length(y_eps_icall) && all(y_eps_icall))
      add(paste0("$", r$name), paste0("Y with EPS/ERR is defined only in IF blocks testing ICALL.EQ.4, which NM-TRAN ",
                                     "treats as simulation code (error 122): move the test of ICALL out of these blocks."))
  }

  # Checked line by line: NM-TRAN refuses an item written after an option of
  # the same line (A.26: "... FILE=sdtab M3LLOQ"), but accepts items written
  # on a line below the options (bbr run iovmm, NONMEM 7.4.4:
  # "$TABLE NOPRINT ONEHEADER FILE=../10par.tab" then "ID TIME ... ").
  for (r in .nm_ctl_records(lines, c("TAB", "TABLE"))) {
    late <- character()
    for (ln in r$text) {
      t <- gsub("\\s*=\\s*", "=", trimws(ln))
      tok <- strsplit(t, "\\s+")[[1]]; tok <- tok[nzchar(tok)]
      if (!length(tok)) next
      is_opt <- .nm_ctl_table_is_opt(tok)
      after_by <- c(FALSE, toupper(utils::head(tok, -1L)) %in% c("BY", "EXCLUDE_BY"))
      first <- which(is_opt)[1]
      if (!is.na(first)) late <- c(late, tok[seq_along(tok) > first & !is_opt & !after_by])
    }
    if (length(late))
      add("$TABLE", paste0("Item(s) written after an option on the same line of $TABLE (record at line ", r$start, "): ",
                           paste(late, collapse = ", "), ". NM-TRAN reads them as unknown options."))
  }
  msg
}

# Add items to a $TABLE line before its first option: NM-TRAN refuses an
# item written after the options (error 20, "UNKNOWN OPTION").
.nm_ctl_table_opts <- c("PRINT", "NOPRINT", "FILE", "NOHEADER", "ONEHEADER", "ONEHEADERALL",
                        "NOAPPEND", "APPEND", "FIRSTONLY", "FIRSTRECORDONLY", "LASTONLY",
                        "NOFORWARD", "FORWARD", "NOTITLE", "NOLABEL", "FORMAT", "LFORMAT",
                        "RFORMAT", "ESAMPLE", "SEED", "CLOCKSEED", "WRESCHOL", "RANMETHOD",
                        "UNCONDITIONAL", "CONDITIONAL", "OMITTED", "EXCLUDE_BY", "NOSUB",
                        "PARAFILE", "VARCALC", "FIXEDETAS", "NPDTYPE", "IDFORMAT", "BY")
.nm_ctl_table_is_opt <- function(tok) {
  u <- toupper(tok)
  grepl("=", u, fixed = TRUE) | sub("=.*$", "", u) %in% .nm_ctl_table_opts
}
.nm_ctl_table_add_items <- function(line, items) {
  code <- sub(";.*$", "", line); comment <- substring(line, nchar(code) + 1L)
  tok <- strsplit(trimws(code), "\\s+")[[1]]
  opt <- which(.nm_ctl_table_is_opt(tok) & seq_along(tok) > 1L)
  at <- if (length(opt)) opt[1] - 1L else length(tok)
  tok <- append(tok, items, after = at)
  paste0(paste(tok, collapse = " "), if (nzchar(comment)) paste0(" ", comment) else "")
}

# Indices used in the code for THETA(n), ETA(n), EPS(n) or ERR(n), including
# the indices of $ABBR REPLACE: lists ETA(OCC_BOVCL)=ETA(4,6,8) and ranges
# ETA(OCC_CL)=ETA(,3:32) (bbr run acop-iov, NONMEM 7.4).
.nm_ctl_used <- function(lines, token) {
  lines <- .nm_ctl_code_lines(lines)
  ab <- .nm_abbr_replace(lines)
  txt <- toupper(paste(.nm_ctl_strip_comments(lines), collapse = "\n"))
  txt <- gsub("\\$ABBR[^\n]*", "", txt, perl = TRUE)
  pat <- paste0("\\b", token, "\\s*\\(\\s*([0-9]+)")
  hits <- regmatches(txt, gregexpr(pat, txt, perl = TRUE))[[1]]
  used <- if (length(hits)) as.integer(sub(pat, "\\1", hits, perl = TRUE)) else integer()
  for (k in which(ab$token == token)) used <- c(used, ab$index[[k]])
  unique(used)
}

# $ABBR REPLACE records that give a name to THETA, ETA, EPS or ERR elements:
# THETA(CL)=THETA(2), ETA(OCC_BOVCL)=ETA(4,6,8), ETA(OCC_CL)=ETA(,3:32) or
# ETA(3 TO 32 BY 2). Returns a data.frame (token, name, index: list of the
# element numbers, in order).
.nm_abbr_replace <- function(lines) {
  code <- toupper(.nm_ctl_strip_comments(lines))
  rx <- "^\\s*\\$ABBR\\s+REPLACE\\s+(THETA|ETA|EPS|ERR)\\s*\\(\\s*([A-Z_][A-Z0-9_]*)\\s*\\)\\s*=\\s*(THETA|ETA|EPS|ERR)\\s*\\(([^)]*)\\)"
  m <- regmatches(code, regexec(rx, code, perl = TRUE))
  m <- m[lengths(m) == 5L]
  out <- data.frame(token = vapply(m, `[`, "", 2), name = vapply(m, `[`, "", 3),
                    stringsAsFactors = FALSE)
  out$index <- lapply(m, function(x) {
    spec <- gsub("\\s+TO\\s+", ":", gsub("\\s+BY\\s+", "/", trimws(x[5])))
    idx <- integer()
    for (p in trimws(strsplit(spec, ",", fixed = TRUE)[[1]])) {
      if (!nzchar(p)) next
      by <- 1L
      if (grepl("/", p, fixed = TRUE)) {
        by <- as.integer(sub(".*/", "", p)); p <- sub("/.*", "", p)
      }
      ab <- suppressWarnings(as.integer(strsplit(p, ":", fixed = TRUE)[[1]]))
      if (anyNA(ab) || !length(ab)) next
      idx <- c(idx, if (length(ab) == 2L) seq(ab[1], ab[2], by = by) else ab[1])
    }
    idx
  })
  out
}

# Text of every record `rec` (comments removed, record name removed), one
# element per record.
.nm_ctl_records_text <- function(lines, rec) {
  code <- .nm_ctl_strip_comments(lines)
  idx <- grep(paste0("^\\s*\\$", rec, "\\b"), code, ignore.case = TRUE)
  controls <- grep("^\\s*\\$", code)
  vapply(idx, function(i) {
    nxt <- controls[controls > i]
    j <- if (length(nxt)) nxt[1] - 1L else length(code)
    sub(paste0("^\\s*\\$", toupper(rec), "\\w*"), "", toupper(paste(code[i:j], collapse = " ")))
  }, character(1))
}

.nm_ctl_num_pattern <- "(?<![A-Z0-9_.])[-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[EeDd][-+]?[0-9]+)?"

# Number of ETAs (OMEGA) or EPSs (SIGMA) declared over all records: diagonal
# values (several per record allowed, including on the record line),
# BLOCK(n), DIAGONAL(n), BLOCK(n) SAME(k), BLOCK SAME (size of the previous
# block).
.nm_ctl_n_random <- function(lines, rec = c("OMEGA", "SIGMA")) {
  rec <- match.arg(rec)
  total <- 0L; last_dim <- NA_integer_
  for (txt in .nm_ctl_records_text(lines, rec)) {
    g <- regmatches(txt, regexec("\\b(BLOCK|DIAGONAL)\\s*\\(\\s*([0-9]+)\\s*\\)", txt, perl = TRUE))[[1]]
    s <- regmatches(txt, regexec("\\bSAME\\b(\\s*\\(\\s*([0-9]+)\\s*\\))?", txt, perl = TRUE))[[1]]
    mult <- if (length(s) && nzchar(s[3])) as.integer(s[3]) else 1L
    if (length(g)) {
      n <- as.integer(g[3]); total <- total + n * mult; last_dim <- n
    } else if (length(s)) {
      if (!is.na(last_dim)) total <- total + last_dim * mult
    } else {
      z <- gsub("\\b(FIX|FIXED|UNINT|VARIANCE|STANDARD|SD|COVARIANCE|CORRELATION|CORRELATON|CHOLESKY|VALUES)\\b",
                " ", txt, perl = TRUE)
      total <- total + length(regmatches(z, gregexpr(.nm_ctl_num_pattern, z, perl = TRUE))[[1]])
      last_dim <- NA_integer_
    }
  }
  as.integer(total)
}

# Number of THETAs declared over all $THETA records: (low, init, up),
# (low, init), (init), bare values, FIX/FIXED, repetition (...)xN and
# labels NAME=(...).
.nm_ctl_n_theta <- function(lines) {
  total <- 0L
  pat <- "\\(([^()]*)\\)(\\s*(?:FIXED|FIX))?(\\s*X\\s*([0-9]+))?"
  for (txt in .nm_ctl_records_text(lines, "THETA")) {
    txt <- gsub("\\b[A-Z_][A-Z0-9_]*\\s*=", " ", txt, perl = TRUE)
    for (g in regmatches(txt, gregexpr(pat, txt, perl = TRUE))[[1]]) {
      m <- regmatches(g, regexec(pat, g, perl = TRUE))[[1]]
      total <- total + if (nzchar(m[5])) as.integer(m[5]) else 1L
    }
    rest <- gsub("\\b(FIX|FIXED|UNINT|NUMBERPOINTS|ABSOLUTE|NOABORT)\\b", " ",
                 gsub(pat, " ", txt, perl = TRUE), perl = TRUE)
    total <- total + length(regmatches(rest, gregexpr(.nm_ctl_num_pattern, rest, perl = TRUE))[[1]])
  }
  as.integer(total)
}


# =============================================================================
# Internal helpers: line differences
# =============================================================================

# Line differences by longest common subsequence: "  " unchanged, "- " removed,
# "+ " added. An inserted line does not shift the following lines.
.nm_ctl_diff <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  na <- length(a); nb <- length(b)
  L <- matrix(0L, na + 1L, nb + 1L)
  if (na && nb) for (i in na:1) for (j in nb:1)
    L[i, j] <- if (identical(a[i], b[j])) L[i + 1L, j + 1L] + 1L else max(L[i + 1L, j], L[i, j + 1L])
  out <- character(); i <- 1L; j <- 1L
  while (i <= na && j <= nb) {
    if (identical(a[i], b[j])) { out <- c(out, paste0("  ", a[i])); i <- i + 1L; j <- j + 1L }
    else if (L[i + 1L, j] >= L[i, j + 1L]) { out <- c(out, paste0("- ", a[i])); i <- i + 1L }
    else { out <- c(out, paste0("+ ", b[j])); j <- j + 1L }
  }
  while (i <= na) { out <- c(out, paste0("- ", a[i])); i <- i + 1L }
  while (j <= nb) { out <- c(out, paste0("+ ", b[j])); j <- j + 1L }
  out
}


# =============================================================================
# nm_ctl_check
# =============================================================================

#' Static check of a NONMEM control stream
#'
#' Looks for inconsistencies that can be detected without running NM-TRAN:
#' `THETA(n)`, `ETA(n)` and `EPS(n)`/`ERR(n)` used without a matching
#' `$THETA`, `$OMEGA` or `$SIGMA` declaration (and declared parameters that
#' are never used), data items used in the code but absent from `$INPUT`,
#' an active `RATE` with first-order absorption, and the requirements of the
#' M3 method (BQL records not ignored, LLOQ defined, Laplace estimation).
#' It also reports four errors refused by NM-TRAN: `TOL` missing with a
#' general ADVAN (6, 8, 9, 13-18), an unknown `METHOD` (FO is written
#' `METHOD=ZERO`), a `THETA` initial estimate of 0 that is not fixed, a
#' `$TABLE` item written after an option of the same line (items written
#' on a line below the options are accepted), and `Y` with `EPS`/`ERR` defined
#' only in `IF` blocks testing `ICALL.EQ.4` (simulation code for NM-TRAN).
#'
#' Declarations are counted over all `$THETA`, `$OMEGA` and `$SIGMA`
#' records, including values written on the record line, `BLOCK(n)`,
#' `DIAGONAL(n)`, `SAME(k)` and `$ABBR REPLACE` lists and ranges
#' (`ETA(4,6,8)`, `ETA(,3:32)`).
#'
#' This check does not replace NM-TRAN: it catches common mistakes early.
#'
#' @param ctl Path to a control stream, or its lines (character vector or
#'   `nm_ctl` object).
#' @param reference Optional original control stream (same forms as `ctl`).
#'   When given, the line differences from `reference` to `ctl` are added.
#'
#' @return An object of class `nm_ctl_check`: a data frame with columns
#'   `level` (`"ERROR"`, `"WARNING"` or `"OK"`), `category` and `message`.
#'   With `reference`, attribute `"diff"` holds the line differences
#'   (`"  "` unchanged, `"- "` removed, `"+ "` added).
#'
#' @examples
#' ctl <- nm_write_ctl()
#' nm_ctl_check(ctl)
#' nm_ctl_check(nm_ctl_mu(ctl), reference = ctl)
#' @export
nm_ctl_check <- function(ctl, reference = NULL) {
  lines <- .nm_ctl_lines(ctl)
  z <- data.frame(level = character(), category = character(), message = character(),
                  stringsAsFactors = FALSE)
  add <- function(level, cat, msg)
    z <<- rbind(z, data.frame(level = level, category = cat, message = msg, stringsAsFactors = FALSE))

  code <- toupper(paste(.nm_ctl_strip_comments(.nm_ctl_code_lines(lines)), collapse = "\n"))
  full <- toupper(paste(.nm_ctl_strip_comments(lines), collapse = "\n"))
  inp  <- .nm_ctl_input_items(lines)

  nth  <- .nm_ctl_n_theta(lines);           uth  <- .nm_ctl_used(lines, "THETA")
  neta <- .nm_ctl_n_random(lines, "OMEGA"); ueta <- .nm_ctl_used(lines, "ETA")
  nsig <- .nm_ctl_n_random(lines, "SIGMA")
  ueps <- sort(unique(c(.nm_ctl_used(lines, "EPS"), .nm_ctl_used(lines, "ERR"))))

  for (i in sort(uth[uth > nth])) add("ERROR", "Parameters", paste0("THETA(", i, ") is used but not declared."))
  for (i in setdiff(seq_len(nth), uth)) add("WARNING", "Parameters", paste0("THETA(", i, ") is declared but apparently unused."))
  for (i in sort(ueta[ueta > neta])) add("ERROR", "IIV", paste0("ETA(", i, ") is used without a matching OMEGA."))
  for (i in setdiff(seq_len(neta), ueta)) add("WARNING", "IIV", paste0("OMEGA/ETA(", i, ") is declared but ETA(", i, ") is apparently unused."))
  for (i in ueps[ueps > nsig]) add("ERROR", "Residual error", paste0("EPS/ERR(", i, ") is used without a matching SIGMA."))

  # with PREDPP ($SUBROUTINES), NM-TRAN appends the EVID and MDV data items
  # itself when $INPUT has none, so they may then be used in the code
  predpp <- any(grepl("^\\s*\\$SUB", lines, ignore.case = TRUE))
  for (v in c("CMT", if (!predpp) c("EVID", "MDV"), "RATE"))
    if (grepl(paste0("\\b", v, "\\b"), code, perl = TRUE) && !(v %in% inp))
      add("WARNING", "$INPUT", paste0(v, " is used but is not an active $INPUT item."))

  # a zero-order input into any compartment (D1, D2, R1...) uses RATE
  if ("RATE" %in% inp && grepl("\\bKA\\b", code, perl = TRUE) && !grepl("\\b[DR][0-9]+\\b", code, perl = TRUE))
    add("WARNING", "Absorption", "RATE is active in a model with KA; check whether RATE should be RATE=DROP.")

  if (grepl("\\bF_FLAG\\b", code, perl = TRUE) && grepl("\\bPHI\\s*\\(", code, perl = TRUE)) {
    if (grepl("\\$DATA[^\\n]*IGNORE[^\\n]*BQL", full, perl = TRUE))
      add("ERROR", "M3", "IGNORE refers to BQL: censored observations must not be excluded.")
    if (!grepl("\\bM3LLOQ\\b|\\bLLOQ\\b", code, perl = TRUE))
      add("ERROR", "M3", "No explicit censoring threshold (LLOQ) found.")
    if (!grepl("\\bLAPLACE\\b|\\bLAPLACIAN\\b", full, perl = TRUE))
      add("WARNING", "M3", "F_FLAG/PHI found without LAPLACE estimation.")
    if (grepl("\\bBQL\\b", code, perl = TRUE) && !("BQL" %in% inp))
      add("ERROR", "M3", "BQL is used but is not an active $INPUT item.")
  }

  # Errors found by NM-TRAN on control streams written by nmpilot (action 26)
  for (m in .nm_ctl_nmtran_checks(lines)) add("ERROR", m[1], m[2])

  if (!nrow(z)) add("OK", "Check", "No obvious inconsistency found by the static checks.")
  if (!is.null(reference)) attr(z, "diff") <- .nm_ctl_diff(.nm_ctl_lines(reference), lines)
  class(z) <- c("nm_ctl_check", "data.frame")
  z
}

#' @export
print.nm_ctl_check <- function(x, ...) {
  y <- as.data.frame(unclass(x)[c("level", "category", "message")], stringsAsFactors = FALSE)
  cat(sprintf("[%s] %s - %s", y$level, y$category, y$message), sep = "\n")
  d <- attr(x, "diff")
  if (!is.null(d)) {
    ch <- d[!startsWith(d, "  ")]
    cat(sprintf("\nDifferences from reference: %d line(s) removed, %d added\n",
                sum(startsWith(d, "- ")), sum(startsWith(d, "+ "))))
    if (length(ch)) cat(ch, sep = "\n")
  }
  invisible(x)
}


# =============================================================================
# nm_ctl_mu
# =============================================================================

# THETA table of the first $THETA record: one row per value line.
.nm_ctl_theta_table <- function(lines) {
  r <- .nm_ctl_record_range(lines, "THETA")
  if (is.null(r)) return(data.frame())
  z <- list(); n <- 0L
  for (k in seq.int(r[1], r[2])) {
    raw <- lines[k]
    x <- trimws(sub(";.*$", "", raw))
    if (k == r[1]) x <- trimws(sub("^\\s*\\$THETA\\w*", "", x, ignore.case = TRUE))
    if (!nzchar(x)) next
    n <- n + 1L
    fix <- grepl("\\bFIX(ED)?\\b", x, ignore.case = TRUE)
    x <- trimws(gsub("\\bFIX(ED)?\\b", "", x, ignore.case = TRUE))
    lo <- hi <- val <- NA_real_
    f <- function(a) suppressWarnings(as.numeric(a))
    if (grepl("^\\(", x)) {
      p <- trimws(strsplit(sub("^\\((.*)\\).*$", "\\1", x), ",", fixed = TRUE)[[1]])
      if (length(p) == 1L) val <- f(p[1])
      if (length(p) == 2L) { lo <- f(p[1]); val <- f(p[2]) }
      if (length(p) >= 3L) { lo <- f(p[1]); val <- f(p[2]); hi <- f(p[3]) }
    } else val <- f(strsplit(x, "\\s+")[[1]][1])
    com <- if (grepl(";", raw, fixed = TRUE)) trimws(sub("^[^;]*;", "", raw)) else paste0("THETA", n)
    z[[length(z) + 1L]] <- data.frame(theta = n, line = k, value = val, lower = lo, upper = hi,
                                      fixed = fix, comment = com, single = (k != r[1]),
                                      stringsAsFactors = FALSE)
  }
  if (length(z)) do.call(rbind, z) else data.frame()
}

# Parameters written P = THETA(n)*EXP(ETA(m)) or P = TVP*EXP(ETA(m)) with
# TVP = THETA(n).
.nm_ctl_mu_candidates <- function(lines) {
  x <- sub(";.*$", "", lines); a <- list(); tv <- list()
  pt <- "^\\s*([A-Za-z][A-Za-z0-9_]*)\\s*=\\s*THETA\\s*\\(\\s*([0-9]+)\\s*\\)\\s*$"
  for (i in seq_along(x)) {
    g <- regmatches(x[i], regexec(pt, x[i], ignore.case = TRUE, perl = TRUE))[[1]]
    if (length(g)) tv[[toupper(g[2])]] <- list(theta = as.integer(g[3]), line = i)
  }
  pd <- "^\\s*([A-Za-z][A-Za-z0-9_]*)\\s*=\\s*THETA\\s*\\(\\s*([0-9]+)\\s*\\)\\s*\\*\\s*EXP\\s*\\(\\s*ETA\\s*\\(\\s*([0-9]+)\\s*\\)\\s*\\)\\s*$"
  pv <- "^\\s*([A-Za-z][A-Za-z0-9_]*)\\s*=\\s*([A-Za-z][A-Za-z0-9_]*)\\s*\\*\\s*EXP\\s*\\(\\s*ETA\\s*\\(\\s*([0-9]+)\\s*\\)\\s*\\)\\s*$"
  for (i in seq_along(x)) {
    g <- regmatches(x[i], regexec(pd, x[i], ignore.case = TRUE, perl = TRUE))[[1]]
    if (length(g)) a[[length(a) + 1L]] <- data.frame(parameter = g[2], theta = as.integer(g[3]),
                                                     eta = as.integer(g[4]), line = i, tvline = NA_integer_)
    g <- regmatches(x[i], regexec(pv, x[i], ignore.case = TRUE, perl = TRUE))[[1]]
    if (length(g) && !is.null(tv[[toupper(g[3])]])) {
      v <- tv[[toupper(g[3])]]
      a[[length(a) + 1L]] <- data.frame(parameter = g[2], theta = v$theta, eta = as.integer(g[4]),
                                        line = i, tvline = v$line)
    }
  }
  if (length(a)) unique(do.call(rbind, a)) else data.frame()
}

#' MU referencing of a NONMEM control stream
#'
#' Rewrites log-normal individual parameters, `P = THETA(n)*EXP(ETA(m))` or
#' `P = TVP*EXP(ETA(m))` with `TVP = THETA(n)`, as
#' `MU_m = LOG(THETA(n))` and `P = EXP(MU_m + ETA(m))`, as required by the
#' EM and Bayesian methods (SAEM, IMP, BAYES). The former `TVP` line is kept
#' as a comment. Covariate effects on the typical value are not handled:
#' only parameters written exactly in one of the two forms above are
#' transformed.
#'
#' @param ctl Path to a control stream, or its lines (character vector or
#'   `nm_ctl` object).
#' @param parameters Names of the parameters to transform. `NULL` (default):
#'   every parameter written in a supported form.
#' @param mode `"readable"` (`MU_m = LOG(THETA(n))`, THETAs unchanged) or
#'   `"linear"` (`MU_m = THETA(n)`, the THETA initial value and bounds are
#'   log-transformed; requires one value per line in `$THETA`; an initial
#'   value of 1, whose log is 0, is written 0.0001 because NM-TRAN refuses an
#'   initial estimate of 0 unless it is fixed).
#' @param file Output path. If `NULL` (default), nothing is written.
#' @param overwrite Logical. Replace an existing `file`.
#' @param verbose Logical. Print a message when the file is written.
#'
#' @return An object of class `nm_ctl` (control stream lines), invisibly
#'   when `file` is given.
#'
#' @examples
#' ctl <- nm_write_ctl()
#' nm_ctl_mu(ctl)
#' nm_ctl_mu(ctl, parameters = c("CL", "V"), mode = "linear")
#' @export
nm_ctl_mu <- function(ctl, parameters = NULL, mode = c("readable", "linear"),
                      file = NULL, overwrite = FALSE, verbose = TRUE) {
  mode <- match.arg(mode)
  lines <- .nm_ctl_lines(ctl)
  cand <- .nm_ctl_mu_candidates(lines)
  if (!nrow(cand)) stop("MU: no parameter written P = THETA(n)*EXP(ETA(m)) or P = TVP*EXP(ETA(m)).")
  if (!is.null(parameters)) {
    bad <- setdiff(toupper(parameters), toupper(cand$parameter))
    if (length(bad)) stop("MU: parameter(s) not in a supported form: ", paste(bad, collapse = ", "),
                          ". Candidates: ", paste(cand$parameter, collapse = ", "), ".")
    cand <- cand[toupper(cand$parameter) %in% toupper(parameters), , drop = FALSE]
  }
  if (any(duplicated(cand$eta))) stop("MU: the same ETA is used by several parameters.")
  all_txt <- paste(lines, collapse = "\n")
  for (k in cand$eta) if (grepl(paste0("\\bMU_", k, "\\s*="), all_txt, ignore.case = TRUE, perl = TRUE))
    stop("MU: MU_", k, " already exists.")

  if (mode == "linear") {
    th <- .nm_ctl_theta_table(lines)
    for (n in unique(cand$theta)) {
      r <- th[th$theta == n, , drop = FALSE]
      if (!nrow(r) || !r$single)
        stop("MU linear: THETA(", n, ") must be written alone on its own line in the first $THETA record.")
      if (!is.finite(r$value) || r$value <= 0)
        stop("MU linear: the initial value of THETA(", n, ") must be > 0.")
      num <- function(v) format(round(v, 6), trim = TRUE, scientific = FALSE)
      v  <- log(r$value)
      # NM-TRAN refuses an initial estimate of 0 unless fixed (error 260):
      # log(1) = 0 is replaced by a value close to 0
      zero <- !r$fixed && round(v, 6) == 0
      if (zero) v <- 1e-4
      lo <- if (is.finite(r$lower) && r$lower > 0) log(r$lower) else NA_real_
      hi <- if (is.finite(r$upper) && r$upper > 0) log(r$upper) else NA_real_
      sp <- if (is.finite(lo) && is.finite(hi)) paste0("(", num(lo), ", ", num(v), ", ", num(hi), ")")
            else if (is.finite(lo)) paste0("(", num(lo), ", ", num(v), ")")
            else if (is.finite(hi)) paste0("(-INF, ", num(v), ", ", num(hi), ")")
            else num(v)
      if (r$fixed) sp <- paste(sp, "FIX")
      lines[r$line] <- paste0("  ", sp, " ; log(", r$comment, ")",
                              if (zero) " (log(1) = 0 replaced by 0.0001: NONMEM refuses an initial value of 0)" else "")
    }
  }
  for (i in seq_len(nrow(cand))) {
    q <- cand[i, ]
    ind <- sub("^(\\s*).*", "\\1", lines[q$line], perl = TRUE)
    mu <- if (mode == "readable") paste0("MU_", q$eta, " = LOG(THETA(", q$theta, "))")
          else paste0("MU_", q$eta, " = THETA(", q$theta, ")")
    lines[q$line] <- paste0(ind, mu, "\n", ind, q$parameter, " = EXP(MU_", q$eta, " + ETA(", q$eta, "))")
    if (!is.na(q$tvline)) lines[q$tvline] <- paste0("; MU referencing - former definition: ", trimws(lines[q$tvline]))
  }
  .nm_ctl_output(.nm_ctl_lines(lines), file = file, overwrite = overwrite, verbose = verbose)
}


# =============================================================================
# nm_ctl_iov
# =============================================================================

# Parameters that can receive an IOV ETA: P = THETA*EXP(ETA), P = TVP*EXP(ETA),
# P = EXP(MU_x + ETA) (all with an IIV ETA), or P = TVP (no IIV).
.nm_ctl_iov_candidates <- function(lines) {
  code <- sub(";.*$", "", lines); out <- list()
  add <- function(parameter, line, kind)
    out[[length(out) + 1L]] <<- data.frame(parameter = parameter, line = as.integer(line),
                                           kind = kind, stringsAsFactors = FALSE)
  p1 <- "^\\s*([A-Za-z][A-Za-z0-9_]*)\\s*=\\s*THETA\\s*\\(\\s*[0-9]+\\s*\\)\\s*\\*\\s*D?EXP\\s*\\(\\s*ETA\\s*\\(\\s*([0-9]+)\\s*\\)\\s*\\)\\s*$"
  p2 <- "^\\s*([A-Za-z][A-Za-z0-9_]*)\\s*=\\s*([A-Za-z][A-Za-z0-9_]*)\\s*\\*\\s*D?EXP\\s*\\(\\s*ETA\\s*\\(\\s*([0-9]+)\\s*\\)\\s*\\)\\s*$"
  p3 <- "^\\s*([A-Za-z][A-Za-z0-9_]*)\\s*=\\s*D?EXP\\s*\\(\\s*MU_[A-Za-z0-9_]+\\s*\\+\\s*ETA\\s*\\(\\s*([0-9]+)\\s*\\)\\s*\\)\\s*$"
  ptv <- "^\\s*([A-Za-z][A-Za-z0-9_]*)\\s*=\\s*THETA\\s*\\(\\s*[0-9]+\\s*\\)\\s*$"
  ps  <- "^\\s*([A-Za-z][A-Za-z0-9_]*)\\s*=\\s*([A-Za-z][A-Za-z0-9_]*)\\s*$"
  m <- function(p, s) regmatches(s, regexec(p, s, ignore.case = TRUE, perl = TRUE))[[1]]
  tv <- character()
  for (s in code) { g <- m(ptv, s); if (length(g)) tv <- unique(c(tv, toupper(g[2]))) }
  for (i in seq_along(code)) {
    g <- m(p1, code[i]); if (length(g)) add(g[2], i, "iiv")
    g <- m(p2, code[i]); if (length(g) && toupper(g[3]) %in% tv) add(g[2], i, "iiv")
    g <- m(p3, code[i]); if (length(g)) add(g[2], i, "iiv")
  }
  for (i in seq_along(code)) {
    g <- m(ps, code[i])
    if (length(g) && toupper(g[3]) %in% tv && toupper(g[2]) != toupper(g[3])) add(g[2], i, "no_iiv")
  }
  if (!length(out)) return(data.frame())
  z <- do.call(rbind, out)
  z[!duplicated(toupper(z$parameter)), , drop = FALSE]
}

# Check that the occasion column of a dataset holds the integers 1..n_occ.
.nm_ctl_check_occ <- function(data, occ, n_occ) {
  if (is.character(data) && length(data) == 1L) data <- utils::read.csv(data, check.names = FALSE)
  hit <- which(toupper(names(data)) == toupper(occ))
  if (!length(hit)) stop("IOV: column ", occ, " is not in the dataset.")
  x <- trimws(as.character(data[[hit[1]]]))
  x <- x[!is.na(x) & nzchar(x) & x != "."]
  z <- suppressWarnings(as.numeric(x))
  if (any(is.na(z))) stop("IOV: ", occ, " contains non-numeric values.")
  if (any(abs(z - round(z)) > 1e-10)) stop("IOV: ", occ, " must contain integers.")
  z <- unique(as.integer(round(z)))
  bad <- setdiff(z, seq_len(n_occ)); miss <- setdiff(seq_len(n_occ), z)
  if (length(bad)) stop("IOV: ", occ, " values outside 1..", n_occ, ": ", paste(sort(bad), collapse = ", "), ".")
  if (length(miss)) stop("IOV: occasion(s) absent from the dataset: ", paste(miss, collapse = ", "), ".")
  invisible(TRUE)
}

#' Add inter-occasion variability to a NONMEM control stream
#'
#' Adds an inter-occasion (between-occasion) random effect to the selected
#' parameters with `$ABBR REPLACE` (NONMEM 7.3 or later): one ETA per
#' parameter and occasion, selected through the occasion data item, and an
#' `$OMEGA BLOCK(p)` followed by `$OMEGA BLOCK(p) SAME(n_occ - 1)`.
#' A parameter with IIV becomes `P = TVP*EXP(ETA(m)+ETA(OCC_BOVP))`; a
#' parameter without IIV (`P = TVP`) becomes `P = TVP*EXP(ETA(OCC_BOVP))`.
#'
#' @param ctl Path to a control stream, or its lines (character vector or
#'   `nm_ctl` object).
#' @param parameters Names of the parameters receiving IOV.
#' @param n_occ Number of occasions (at least 2).
#' @param occ Name of the occasion data item; it must be active in `$INPUT`
#'   and take the values 1 to `n_occ`. Letters and digits only.
#' @param cv Initial IOV, as CV% (one value, or one per parameter, possibly
#'   named).
#' @param structure `"independent"` (diagonal block) or `"correlated"`
#'   (common initial correlation `rho`).
#' @param rho Initial correlation for `structure = "correlated"`.
#' @param data Optional dataset (data frame or CSV path) used to check that
#'   `occ` holds the integers 1 to `n_occ`.
#' @param file Output path. If `NULL` (default), nothing is written.
#' @param overwrite Logical. Replace an existing `file`.
#' @param verbose Logical. Print a message when the file is written.
#'
#' @return An object of class `nm_ctl` (control stream lines), invisibly
#'   when `file` is given.
#'
#' @examples
#' ctl <- nm_write_ctl(input = c("ID", "TIME", "DV", "AMT", "EVID", "MDV", "CMT", "OCC"))
#' nm_ctl_iov(ctl, parameters = c("CL", "V"), n_occ = 3)
#' @export
nm_ctl_iov <- function(ctl, parameters, n_occ, occ = "OCC", cv = 20,
                       structure = c("independent", "correlated"), rho = 0.1,
                       data = NULL, file = NULL, overwrite = FALSE, verbose = TRUE) {
  structure <- match.arg(structure)
  lines <- .nm_ctl_lines(ctl)
  n_occ <- as.integer(n_occ)
  if (length(n_occ) != 1L || is.na(n_occ) || n_occ < 2L) stop("IOV: `n_occ` must be at least 2.")
  occ <- trimws(occ)
  if (!grepl("^[A-Za-z][A-Za-z0-9]*$", occ))
    stop("IOV: `occ` must contain letters and digits only (e.g. OCC).")
  if (!toupper(occ) %in% .nm_ctl_input_items(lines))
    stop("IOV: ", occ, " must be an active item of $INPUT.")
  if (grepl(paste0("ETA\\s*\\(\\s*", occ, "_"), paste(lines, collapse = "\n"), ignore.case = TRUE, perl = TRUE))
    stop("IOV: an ETA selection by ", occ, " already exists.")
  if (!is.null(data)) .nm_ctl_check_occ(data, occ, n_occ)

  cand <- .nm_ctl_iov_candidates(lines)
  bad <- setdiff(toupper(parameters), toupper(cand$parameter))
  if (length(bad)) stop("IOV: parameter(s) not in a supported form: ", paste(bad, collapse = ", "), ".")
  cand <- cand[match(toupper(parameters), toupper(cand$parameter)), , drop = FALSE]
  p <- nrow(cand)

  if (!is.null(names(cv))) cv <- cv[cand$parameter]
  cv <- rep_len(as.numeric(cv), p)
  if (any(!is.finite(cv)) || any(cv <= 0)) stop("IOV: `cv` must be strictly positive.")
  vars <- log(1 + (cv / 100)^2)
  if (structure == "correlated" && p > 1L && (!is.finite(rho) || rho >= 1 || rho <= -1 / (p - 1)))
    stop(sprintf("IOV: rho must satisfy -1/(p-1) < rho < 1 (here rho > %.3f).", -1 / (p - 1)))

  start <- max(c(0L, .nm_ctl_used(lines, "ETA"), .nm_ctl_n_random(lines, "OMEGA"))) + 1L
  eta_mat <- matrix(seq.int(start, length.out = n_occ * p), nrow = n_occ, ncol = p, byrow = TRUE)
  tags <- paste0(toupper(occ), "_BOV", gsub("[^A-Za-z0-9]", "", toupper(cand$parameter)))
  abbr <- vapply(seq_len(p), function(j)
    paste0("$ABBR REPLACE ETA(", tags[j], ")=ETA(", paste(eta_mat[, j], collapse = ","), ")"), character(1))

  for (j in seq_len(p)) {
    q <- cand[j, ]; ln <- lines[q$line]; code <- sub(";.*$", "", ln)
    com <- if (grepl(";", ln, fixed = TRUE)) paste0(" ;", sub("^[^;]*;", "", ln)) else ""
    indent <- sub("^(\\s*).*", "\\1", code, perl = TRUE)
    rhs <- trimws(sub("^[^=]*=", "", code))
    rhs2 <- if (q$kind == "iiv")
      sub("(ETA\\s*\\(\\s*[0-9]+\\s*\\))", paste0("\\1+ETA(", tags[j], ")"), rhs, ignore.case = TRUE, perl = TRUE)
    else paste0(rhs, "*EXP(ETA(", tags[j], "))")
    lines[q$line] <- paste0(indent, q$parameter, " = ", rhs2, com)
  }

  code_idx <- grep("^\\s*\\$(PK|PRED|ERROR|DES|AES|AESINITIAL)\\b", lines, ignore.case = TRUE, perl = TRUE)
  if (!length(code_idx)) stop("IOV: no $PK, $PRED or $ERROR record found.")
  lines <- append(lines, c("; inter-occasion variability (IOV)", abbr), after = code_idx[1] - 1L)

  om_idx <- grep("^\\s*\\$OMEGA\\b", lines, ignore.case = TRUE)
  if (!length(om_idx)) stop("IOV: no $OMEGA record found.")
  nxt <- grep("^\\s*\\$", lines); nxt <- nxt[nxt > max(om_idx)]
  om_end <- if (length(nxt)) nxt[1] - 1L else length(lines)
  fmt <- function(x) format(round(x, 6), trim = TRUE, scientific = FALSE, nsmall = 0)
  omega <- c(paste0("; IOV: ", n_occ, " occasions, ", structure, " structure"),
             paste0("$OMEGA BLOCK(", p, ")"))
  for (i in seq_len(p)) {
    vals <- vapply(seq_len(i), function(j)
      if (i == j) fmt(vars[i]) else fmt(if (structure == "correlated") rho * sqrt(vars[i] * vars[j]) else 0),
      character(1))
    omega <- c(omega, paste(vals, collapse = " "))
  }
  omega <- c(omega, paste0("$OMEGA BLOCK(", p, ") SAME(", n_occ - 1L, ")"))
  .nm_ctl_output(append(lines, omega, after = om_end), file = file, overwrite = overwrite, verbose = verbose)
}


# =============================================================================
# nm_ctl_m3
# =============================================================================

.nm_ctl_y_exprs <- function(error_lines) {
  z <- .nm_ctl_strip_comments(error_lines)
  m <- grep("(^|\\s)Y\\s*=", z, ignore.case = TRUE)
  if (!length(m)) return(character(0))
  trimws(sub("^.*?Y\\s*=\\s*", "", z[m], ignore.case = TRUE, perl = TRUE))
}

# Residual SD on the DV scale, inferred from usual error models, or NULL.
.nm_ctl_infer_sd <- function(error_lines) {
  ys <- .nm_ctl_y_exprs(error_lines)
  if (!length(ys) || any(grepl("LOG\\(|LOGIT\\(|PHI\\(|F_FLAG", toupper(ys)))) return(NULL)
  y <- gsub("ERR\\(", "EPS(", toupper(gsub("\\s+", "", utils::tail(ys, 1))))
  p <- if (grepl("IPRED", y, fixed = TRUE)) "IPRED" else "F"
  if (y %in% c(paste0(p, "+W*EPS(1)"), paste0("W*EPS(1)+", p))) return("ABS(W)*SQRT(SIGMA(1,1))")
  if (y %in% c(paste0(p, "+EPS(1)"), paste0("EPS(1)+", p))) return("SQRT(SIGMA(1,1))")
  if (y %in% c(paste0(p, "*(1+EPS(1))"), paste0(p, "*(1.0+EPS(1))"), paste0(p, "+", p, "*EPS(1)")))
    return(paste0("ABS(", p, ")*SQRT(SIGMA(1,1))"))
  if (y %in% c(paste0(p, "+", p, "*EPS(1)+EPS(2)"), paste0(p, "+EPS(2)+", p, "*EPS(1)"),
               paste0(p, "*(1+EPS(1))+EPS(2)"), paste0("EPS(2)+", p, "*(1+EPS(1))")))
    return(paste0("SQRT(", p, "**2*SIGMA(1,1)+SIGMA(2,2))"))
  NULL
}

.nm_ctl_pred_expr <- function(error_lines) {
  txt <- toupper(paste(.nm_ctl_strip_comments(error_lines), collapse = "\n"))
  if (grepl("(^|[^A-Z0-9_])IPRED\\s*=", txt, perl = TRUE)) return("IPRED")
  ys <- .nm_ctl_y_exprs(error_lines)
  if (length(ys) && any(grepl("\\bF\\b", toupper(ys)))) return("F")
  NULL
}

#' Implement the M3 method for BQL data in a NONMEM control stream
#'
#' Rewrites `$ERROR` so that quantifiable observations keep their usual
#' likelihood and BQL observations contribute the probability
#' `P(DV < LLOQ) = PHI((LLOQ - IPRED)/SD)` (Beal 2001, method M3), with
#' `F_FLAG`. The dataset is not modified: comments listing its requirements
#' are inserted after `$DATA` (BQL flag coded 0/1, BQL records kept, EVID=0
#' and MDV=0 on BQL records). `LAPLACE` is added to `$ESTIMATION`
#' (EM, Bayesian and FO methods are replaced by FOCE-I Laplace), and
#' `M3LLOQ`, `M3SD`, `M3PBQL` are added to the first `$TABLE`.
#'
#' The LLOQ comes from exactly one of: a constant (`lloq`), a dataset column
#' (`lloq_col`), or a code of the analytical method (`method_col` and
#' `method_lloq`), for studies pooling assays with different LLOQs.
#'
#' @param ctl Path to a control stream, or its lines (character vector or
#'   `nm_ctl` object).
#' @param bql Name of the BQL flag data item (0 = quantifiable, 1 = BQL).
#' @param lloq Constant LLOQ, on the DV scale.
#' @param lloq_col Name of a dataset column holding the LLOQ.
#' @param method_col Name of a dataset column holding a numeric code of the
#'   analytical method.
#' @param method_lloq Named numeric vector of LLOQs by method code, e.g.
#'   `c("1" = 0.1, "2" = 0.05)`.
#' @param residual_sd NONMEM expression of the residual SD on the DV scale.
#'   `NULL` (default): inferred from the `Y =` line for additive,
#'   proportional, combined and `IPRED + W*EPS(1)` models.
#' @param npde_mode Logical. Add `NPDE_MODE` and `DV_LOQ` handling for the
#'   NONMEM NPDE diagnostics.
#' @param file Output path. If `NULL` (default), nothing is written.
#' @param overwrite Logical. Replace an existing `file`.
#' @param verbose Logical. Print a message when the file is written.
#'
#' @return An object of class `nm_ctl` (control stream lines), invisibly
#'   when `file` is given.
#'
#' @references Beal SL. Ways to fit a PK model with some data below the
#'   quantification limit. J Pharmacokinet Pharmacodyn. 2001;28(5):481-504.
#'
#' @examples
#' ctl <- nm_write_ctl(input = c("ID", "TIME", "DV", "AMT", "EVID", "MDV", "CMT", "BQL"),
#'                     error = "combined")
#' nm_ctl_m3(ctl, lloq = 0.1)
#' @export
nm_ctl_m3 <- function(ctl, bql = "BQL", lloq = NULL, lloq_col = NULL,
                      method_col = NULL, method_lloq = NULL, residual_sd = NULL,
                      npde_mode = FALSE, file = NULL, overwrite = FALSE, verbose = TRUE) {
  lines <- .nm_ctl_lines(ctl)
  bql <- toupper(trimws(bql))
  if (!nzchar(bql)) stop("M3: `bql` is empty.")
  src <- c(constant = !is.null(lloq), column = !is.null(lloq_col), method = !is.null(method_col))
  if (sum(src) != 1L) stop("M3: give exactly one of `lloq`, `lloq_col` or `method_col`.")
  source <- names(src)[src]

  lloq_code <- switch(source,
    constant = {
      v <- suppressWarnings(as.numeric(lloq))
      if (length(v) != 1L || !is.finite(v) || v <= 0) stop("M3: `lloq` must be a positive number.")
      sprintf("M3LLOQ = %.8g", v)
    },
    column = paste0("M3LLOQ = ", toupper(lloq_col)),
    method = {
      if (is.null(method_lloq) || !length(method_lloq) || is.null(names(method_lloq)))
        stop("M3: `method_lloq` must be a named vector, e.g. c(\"1\" = 0.1, \"2\" = 0.05).")
      codes <- suppressWarnings(as.numeric(names(method_lloq)))
      vals <- as.numeric(method_lloq)
      if (any(!is.finite(codes))) stop("M3: method codes (names of `method_lloq`) must be numeric.")
      if (any(!is.finite(vals) | vals <= 0)) stop("M3: values of `method_lloq` must be positive.")
      c("M3LLOQ = -1", sprintf("IF (%s.EQ.%.8g) M3LLOQ = %.8g", toupper(method_col), codes, vals))
    })

  er <- .nm_ctl_record_range(lines, "ERROR")
  if (is.null(er)) stop("M3: no $ERROR record found.")
  E <- lines[er[1]:er[2]]
  if (any(grepl("F_FLAG", E, ignore.case = TRUE))) stop("M3: F_FLAG already exists in $ERROR.")
  pred <- .nm_ctl_pred_expr(E)
  if (is.null(pred)) stop("M3: cannot identify IPRED or F on the DV scale.")
  sd <- residual_sd
  if (is.null(sd) || !nzchar(trimws(sd))) sd <- .nm_ctl_infer_sd(E)
  if (is.null(sd)) stop("M3: residual SD cannot be inferred from the Y line; give `residual_sd`.")

  data_ix <- grep("^\\s*\\$DATA\\b", lines, ignore.case = TRUE)
  if (length(data_ix) &&
      grepl(paste0("IGNORE\\s*=\\s*\\([^)]*\\b", bql, "\\b"), toupper(lines[data_ix[1]]), perl = TRUE))
    stop("M3: BQL records must not be excluded by $DATA IGNORE=(...).")

  inp <- .nm_ctl_input_items(lines)
  lloq_item <- switch(source, constant = NULL, column = toupper(lloq_col), method = toupper(method_col))
  db_comments <- c(
    "",
    "; =====================================================================",
    "; M3 BQL - DATASET REQUIREMENTS (the dataset is not modified)",
    paste0("; - ", bql, " must be coded 0 = quantifiable, 1 = BQL."),
    paste0("; - Records with ", bql, " = 1 stay in the analysis: do NOT IGNORE them."),
    "; - BQL records are observations: EVID = 0 and MDV = 0.",
    "; - DV must stay numeric on BQL records; its value is not used as an",
    ";   observed concentration in the censored M3 contribution.",
    switch(source,
      constant = c(paste0("; - Constant LLOQ used in the control stream: ", format(as.numeric(lloq), scientific = FALSE), "."),
                   "; - IMPORTANT: this value enters the likelihood of the BQL records."),
      column = paste0("; - The dataset must contain a positive column ", lloq_item, ", active in $INPUT."),
      method = paste0("; - The dataset must contain column ", lloq_item,
                      "; the LLOQ is assigned in $ERROR from the method code.")),
    if (!bql %in% inp)
      paste0("; *** ACTION REQUIRED: ", bql, " is not in $INPUT. Add it at the position of the column in the dataset. ***")
    else paste0("; - ", bql, " is in $INPUT."),
    if (!is.null(lloq_item) && !lloq_item %in% inp)
      paste0("; *** ACTION REQUIRED: ", lloq_item, " is not in $INPUT. ***"),
    "; =====================================================================",
    "")
  lines <- if (length(data_ix)) append(lines, db_comments, after = data_ix[1]) else c(lines[1], db_comments, lines[-1])

  er <- .nm_ctl_record_range(lines, "ERROR")
  E <- lines[er[1]:er[2]]; header <- E[1]; body <- E[-1]
  y_idx <- grep("(^|\\s)Y\\s*=", .nm_ctl_strip_comments(body), ignore.case = TRUE)
  if (!length(y_idx)) stop("M3: no Y = line found in $ERROR.")
  # One IF/THEN/ELSE around each Y line, the M3 quantities computed just
  # before the first one. The former code put Y with ERR in a block written
  # IF (BQL.EQ.0.OR.ICALL.EQ.4) THEN, which NM-TRAN 7.6 treats as simulation
  # code: it refused the model ("EPS VARIABLES REQUIRED WITH POPULATION
  # DATA", error 122) and gave warning 53. The censoring flag M3CENS keeps the
  # test of ICALL out of the blocks defining Y; simulations (ICALL = 4) use
  # the normal residual error.
  m3 <- c("; --- M3 BQL likelihood ---", lloq_code,
          paste0("M3SD = ", sd), "IF (M3SD.LE.0) M3SD = 1.0E-12",
          paste0("M3Z = (M3LLOQ - ", pred, ")/M3SD"),
          "M3PBQL = PHI(M3Z) + 1e-30",
          "M3CENS = 0",
          paste0("IF (", bql, ".EQ.1) M3CENS = 1"),
          if (isTRUE(npde_mode)) "IF (NPDE_MODE.EQ.1) M3CENS = 0",
          "IF (ICALL.EQ.4) M3CENS = 0",
          if (isTRUE(npde_mode)) c(paste0("IF (", bql, ".EQ.1) DV_LOQ = M3LLOQ"),
                                   "; NPDE_MODE/DV_LOQ: only for the NONMEM NPD/NPDE diagnostics."))
  new_body <- character(0)
  for (k in seq_along(body)) {
    if (k %in% y_idx) {
      indent <- sub("^(\\s*).*", "\\1", body[k])
      if (k == y_idx[1]) new_body <- c(new_body, paste0(indent, m3))
      new_body <- c(new_body,
                    paste0(indent, "IF (M3CENS.EQ.1) THEN"),
                    paste0(indent, "  F_FLAG = 1"),
                    paste0(indent, "  Y = M3PBQL"),
                    paste0(indent, "  MDVRES = 1"),
                    paste0(indent, "ELSE"),
                    paste0(indent, "  F_FLAG = 0"),
                    paste0(indent, "  ", trimws(body[k])),
                    paste0(indent, "ENDIF"))
    } else new_body <- c(new_body, body[k])
  }
  lines <- c(lines[seq_len(er[1] - 1L)], header, new_body,
             if (er[2] < length(lines)) lines[(er[2] + 1L):length(lines)] else character(0))

  # F_FLAG requires a Laplacian estimation
  for (i in grep("^\\s*\\$EST(?:IMATION)?\\b", lines, ignore.case = TRUE, perl = TRUE)) {
    up <- toupper(lines[i])
    if (grepl("METHOD\\s*=\\s*(IMP|ITS|SAEM|BAYES|FO|0)(\\s|$)", up, perl = TRUE)) {
      lines[i] <- "$ESTIMATION METHOD=COND INTERACTION LAPLACE MAXEVAL=9999 SIGDIGITS=3 PRINT=5"
      if (isTRUE(verbose)) message("M3: estimation method replaced by FOCE-I Laplace (line ", i, ").")
    } else if (!grepl("\\bLAPLAC", up)) lines[i] <- paste(lines[i], "LAPLACE")
  }
  ti <- grep("^\\s*\\$TABLE\\b", lines, ignore.case = TRUE)
  if (length(ti)) {
    new <- c("M3LLOQ", "M3SD", "M3PBQL")
    new <- new[!vapply(new, function(x) grepl(paste0("\\b", x, "\\b"), lines[ti[1]], ignore.case = TRUE), logical(1))]
    if (length(new)) lines[ti[1]] <- .nm_ctl_table_add_items(lines[ti[1]], new)
  }

  .nm_ctl_output(lines, file = file, overwrite = overwrite, verbose = verbose)
}
