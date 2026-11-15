# ==============================================================================
# nm_identify_parameters.R
# Semantic parser for NONMEM control streams (.ctl)
#
# Reconstructs the semantic links between PK parameters, THETA and ETA/OMEGA.
# The executable equations are the primary source of truth; comments are kept
# as labels / consistency checks only.
#
# Semantic roles:
#   typical_value     : THETA defining the typical/structural PK value
#   covariate_effect  : THETA controlling a covariate relationship
#   iiv_scale         : THETA multiplying/modulating an ETA effect
#   residual_error    : THETA used in $ERROR
#   iiv               : ETA acting on a PK parameter
#   structural_other  : THETA linked to PK but not classifiable above
#
# Covariate effects are recognised in power form ((WT/70)**THETA(n),
# THETA(n)**SEX), exponential form (EXP(THETA(n)*WT)), linear form
# (THETA(1)*(1 + THETA(n)*(WT-70))) and conditional form
# (IF(SEX.EQ.1) TVCL = TVCL*THETA(n), one-line IF or IF ... THEN blocks).
# For a variable only defined under conditions, the first definition in the
# file is the reference value.
# A covariate is any data item declared in
# $INPUT that is not a NONMEM-defined item, whatever its name; without $INPUT,
# any symbol used in $PK but not defined there (and not a NONMEM item).
#
# $OMEGA BLOCK(n) and DIAGONAL(n) are parsed, including covariance elements.
# Mapping is element-wise (no artificial THETA x ETA Cartesian grid).
#
# Used by nm_extract_runs() to name THETA/ETA after the PK parameters.
# Internal helpers use the prefix `.nm_ip_`.
# ==============================================================================

.nm_ip_trim <- function(x) trimws(x, which = "both")
.nm_ip_strip_inline_comment <- function(x) sub(";.*$", "", x)
.nm_ip_is_full_comment <- function(x) grepl("^\\s*;", x)

.nm_ip_comment <- function(x) {
  if (!grepl(";", x, fixed = TRUE)) return(NA_character_)
  .nm_ip_trim(sub("^[^;]*;", "", x))
}

.nm_ip_label_clean <- function(x) {
  if (length(x) == 0L || is.na(x)) return(NA_character_)
  x <- .nm_ip_trim(x)
  .nm_ip_trim(strsplit(x, ";", fixed = TRUE)[[1]][1])
}

.nm_ip_section_lines <- function(lines, section) {
  sec <- toupper(section)
  starts <- which(grepl(paste0("^\\s*\\", sec, "(?:\\s|$)"),
                         toupper(lines), perl = TRUE))
  if (!length(starts)) return(character())

  out <- character()
  for (s in starts) {
    first <- sub(paste0("^\\s*\\", sec, "\\s*"), "", lines[s],
                 ignore.case = TRUE, perl = TRUE)
    if (nzchar(.nm_ip_trim(first))) out <- c(out, first)

    nxt <- which(seq_along(lines) > s & grepl("^\\s*\\$[A-Za-z]", lines))
    e <- if (length(nxt)) min(nxt) - 1L else length(lines)
    if (e >= s + 1L) out <- c(out, lines[(s + 1L):e])
  }
  out
}

.nm_ip_numbers <- function(x) {
  m <- regmatches(x, gregexpr(
    "[-+]?(?:(?:[0-9]+\\.?[0-9]*)|(?:\\.[0-9]+))(?:[Ee][-+]?[0-9]+)?",
    x, perl = TRUE))[[1]]
  if (!length(m) || identical(m, character(0))) return(numeric())
  suppressWarnings(as.numeric(m))
}

.nm_ip_split_ints <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) return(integer())
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

.nm_ip_split_vars <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) return(character())
  strsplit(x, ",", fixed = TRUE)[[1]]
}

# ------------------------------------------------------------------------------
# THETA
# ------------------------------------------------------------------------------

.nm_ip_parse_theta <- function(lines) {
  z <- .nm_ip_section_lines(lines, "$THETA")
  if (!length(z)) return(data.frame())

  z <- z[!.nm_ip_is_full_comment(z)]
  z <- z[nzchar(.nm_ip_trim(.nm_ip_strip_inline_comment(z)))]
  if (!length(z)) return(data.frame())

  ans <- list(); k <- 0L

  for (ln in z) {
    code <- .nm_ip_trim(.nm_ip_strip_inline_comment(ln))
    if (!nzchar(code)) next

    # Parenthesized THETA records or simple scalar records.
    entries <- regmatches(code, gregexpr(
      "\\([^)]*\\)\\s*(?:FIX(?:ED)?\\b)?|[-+]?(?:(?:[0-9]+\\.?[0-9]*)|(?:\\.[0-9]+))(?:[Ee][-+]?[0-9]+)?\\s*(?:FIX(?:ED)?\\b)?",
      code, perl = TRUE, ignore.case = TRUE))[[1]]

    if (!length(entries) || identical(entries, character(0))) entries <- code
    cm <- .nm_ip_comment(ln)

    for (entry in entries) {
      entry <- .nm_ip_trim(entry)
      fixed <- grepl("\\bFIX(?:ED)?\\b", entry, ignore.case = TRUE, perl = TRUE)
      core <- .nm_ip_trim(gsub("\\bFIX(?:ED)?\\b", "", entry,
                           ignore.case = TRUE, perl = TRUE))
      vals <- .nm_ip_numbers(core)

      if (!length(vals)) next
      if (length(vals) == 1L) {
        lower <- NA_real_; initial <- vals[1]; upper <- NA_real_
      } else if (length(vals) == 2L) {
        lower <- vals[1]; initial <- vals[2]; upper <- NA_real_
      } else {
        lower <- vals[1]; initial <- vals[2]; upper <- vals[3]
      }

      k <- k + 1L
      ans[[k]] <- data.frame(
        theta = k,
        label = .nm_ip_label_clean(cm),
        lower = lower,
        initial = initial,
        upper = upper,
        fixed = fixed,
        raw = ln,
        stringsAsFactors = FALSE
      )
      cm <- NA_character_
    }
  }

  if (!length(ans)) return(data.frame())
  do.call(rbind, ans)
}

# ------------------------------------------------------------------------------
# OMEGA: diagonal, DIAGONAL(n), BLOCK(n)
# ------------------------------------------------------------------------------

.nm_ip_omega_header <- function(code) {
  b <- regexec("\\bBLOCK\\s*\\(\\s*([0-9]+)\\s*\\)", code,
               ignore.case = TRUE, perl = TRUE)
  mb <- regmatches(code, b)[[1]]
  if (length(mb) >= 2L)
    return(list(type = "BLOCK", n = as.integer(mb[2])))

  d <- regexec("\\bDIAGONAL\\s*\\(\\s*([0-9]+)\\s*\\)", code,
               ignore.case = TRUE, perl = TRUE)
  md <- regmatches(code, d)[[1]]
  if (length(md) >= 2L)
    return(list(type = "DIAGONAL", n = as.integer(md[2])))

  list(type = "SCALAR", n = 1L)
}

.nm_ip_remove_omega_header <- function(code) {
  sub("^\\s*(?:BLOCK|DIAGONAL)\\s*\\(\\s*[0-9]+\\s*\\)\\s*", "", code,
      ignore.case = TRUE, perl = TRUE)
}

.nm_ip_parse_omega <- function(lines) {
  z <- .nm_ip_section_lines(lines, "$OMEGA")
  if (!length(z)) return(list(eta = data.frame(), elements = data.frame()))

  # Keep comments because they often label ETA dimensions.
  records <- z[!.nm_ip_is_full_comment(z)]
  records <- records[nzchar(.nm_ip_trim(.nm_ip_strip_inline_comment(records)))]
  if (!length(records)) return(list(eta = data.frame(), elements = data.frame()))

  eta_rows <- list(); element_rows <- list()
  eta_offset <- 0L
  er <- 0L; rr <- 0L
  i <- 1L

  while (i <= length(records)) {
    ln <- records[i]
    code <- .nm_ip_trim(.nm_ip_strip_inline_comment(ln))
    hdr <- .nm_ip_omega_header(code)

    if (hdr$type == "SCALAR") {
      vals <- .nm_ip_numbers(code)
      if (!length(vals)) { i <- i + 1L; next }

      fixed <- grepl("\\bFIX(?:ED)?\\b", code, ignore.case = TRUE, perl = TRUE)
      eta_offset <- eta_offset + 1L
      rr <- rr + 1L
      eta_rows[[rr]] <- data.frame(
        eta = eta_offset,
        omega_label = .nm_ip_label_clean(.nm_ip_comment(ln)),
        omega_initial = vals[1],
        omega_fixed = fixed,
        omega_structure = "DIAGONAL",
        block_id = NA_integer_,
        block_position = 1L,
        omega_raw = ln,
        stringsAsFactors = FALSE
      )
      er <- er + 1L
      element_rows[[er]] <- data.frame(
        eta_i = eta_offset, eta_j = eta_offset,
        element_type = "variance", initial = vals[1], fixed = fixed,
        structure = "DIAGONAL", block_id = NA_integer_, raw = ln,
        stringsAsFactors = FALSE
      )
      i <- i + 1L
      next
    }

    # BLOCK(n) / DIAGONAL(n): values can be on header line and/or following lines.
    n <- hdr$n
    need <- if (hdr$type == "BLOCK") n * (n + 1L) / 2L else n
    vals <- numeric(); raws <- character(); comments <- character(); fixed_flags <- logical()

    first_code <- .nm_ip_remove_omega_header(code)
    if (nzchar(first_code)) {
      vals <- c(vals, .nm_ip_numbers(first_code))
      raws <- c(raws, ln)
      comments <- c(comments, .nm_ip_comment(ln))
      fixed_flags <- c(fixed_flags,
                       grepl("\\bFIX(?:ED)?\\b", first_code,
                             ignore.case = TRUE, perl = TRUE))
    }

    j <- i + 1L
    while (length(vals) < need && j <= length(records)) {
      next_code <- .nm_ip_trim(.nm_ip_strip_inline_comment(records[j]))
      # A new explicit BLOCK/DIAGONAL header starts a new record.
      if (grepl("\\b(?:BLOCK|DIAGONAL)\\s*\\(", next_code,
                ignore.case = TRUE, perl = TRUE)) break

      newvals <- .nm_ip_numbers(next_code)
      if (length(newvals)) {
        vals <- c(vals, newvals)
        raws <- c(raws, records[j])
        comments <- c(comments, .nm_ip_comment(records[j]))
        fixed_flags <- c(fixed_flags,
                         grepl("\\bFIX(?:ED)?\\b", next_code,
                               ignore.case = TRUE, perl = TRUE))
      }
      j <- j + 1L
    }

    # NONMEM SAME is structurally meaningful but has no new numeric initial values.
    if (grepl("\\bSAME\\b", code, ignore.case = TRUE, perl = TRUE)) {
      vals <- rep(NA_real_, need)
    }

    if (length(vals) < need) vals <- c(vals, rep(NA_real_, need - length(vals)))
    if (length(vals) > need) vals <- vals[seq_len(need)]

    block_id <- eta_offset + 1L
    eta_ids <- eta_offset + seq_len(n)

    if (hdr$type == "DIAGONAL") {
      for (a in seq_len(n)) {
        rr <- rr + 1L
        eta_rows[[rr]] <- data.frame(
          eta = eta_ids[a],
          omega_label = if (a <= length(comments)) .nm_ip_label_clean(comments[a]) else NA_character_,
          omega_initial = vals[a],
          omega_fixed = if (length(fixed_flags)) all(fixed_flags) else FALSE,
          omega_structure = "DIAGONAL",
          block_id = block_id,
          block_position = a,
          omega_raw = paste(raws, collapse = " | "),
          stringsAsFactors = FALSE
        )
        er <- er + 1L
        element_rows[[er]] <- data.frame(
          eta_i = eta_ids[a], eta_j = eta_ids[a],
          element_type = "variance", initial = vals[a],
          fixed = if (length(fixed_flags)) all(fixed_flags) else FALSE,
          structure = "DIAGONAL", block_id = block_id,
          raw = paste(raws, collapse = " | "), stringsAsFactors = FALSE
        )
      }
    } else {
      # BLOCK values are lower triangular row-wise:
      # (1,1), (2,1), (2,2), (3,1), (3,2), (3,3), ...
      pos <- 0L
      diag_values <- rep(NA_real_, n)
      diag_fixed <- rep(FALSE, n)

      for (a in seq_len(n)) {
        for (b in seq_len(a)) {
          pos <- pos + 1L
          isdiag <- a == b
          if (isdiag) diag_values[a] <- vals[pos]
          er <- er + 1L
          element_rows[[er]] <- data.frame(
            eta_i = eta_ids[a], eta_j = eta_ids[b],
            element_type = if (isdiag) "variance" else "covariance",
            initial = vals[pos],
            fixed = if (length(fixed_flags)) all(fixed_flags) else FALSE,
            structure = "BLOCK", block_id = block_id,
            raw = paste(raws, collapse = " | "), stringsAsFactors = FALSE
          )
        }
      }

      # Comments in BLOCK records are not standardized. Associate comments from
      # rows preferentially to the corresponding ETA dimension when available.
      block_comments <- unique(stats::na.omit(vapply(comments, .nm_ip_label_clean, character(1))))
      for (a in seq_len(n)) {
        lab <- if (length(block_comments) >= a) block_comments[a] else NA_character_
        rr <- rr + 1L
        eta_rows[[rr]] <- data.frame(
          eta = eta_ids[a], omega_label = lab,
          omega_initial = diag_values[a],
          omega_fixed = if (length(fixed_flags)) all(fixed_flags) else FALSE,
          omega_structure = "BLOCK",
          block_id = block_id,
          block_position = a,
          omega_raw = paste(raws, collapse = " | "),
          stringsAsFactors = FALSE
        )
      }
    }

    eta_offset <- eta_offset + n
    i <- max(j, i + 1L)
  }

  list(
    eta = if (length(eta_rows)) do.call(rbind, eta_rows) else data.frame(),
    elements = if (length(element_rows)) do.call(rbind, element_rows) else data.frame()
  )
}

# ------------------------------------------------------------------------------
# Generic assignment parser for $PK / $ERROR
# ------------------------------------------------------------------------------

# Split a leading IF(condition) from the rest of a line, with balanced
# parentheses (e.g. "IF(WT.GT.80.AND.(SEX.EQ.0)) TVKA = THETA(7)").
.nm_ip_split_if <- function(ln) {
  m <- regexpr("^\\s*IF\\s*\\(", ln, ignore.case = TRUE, perl = TRUE)
  if (m == -1L) return(list(cond = NA_character_, rest = ln))
  ch <- strsplit(ln, "", fixed = TRUE)[[1]]
  op <- m + attr(m, "match.length") - 1L
  depth <- 0L
  for (p in op:length(ch)) {
    if (ch[p] == "(") depth <- depth + 1L
    if (ch[p] == ")") {
      depth <- depth - 1L
      if (depth == 0L)
        return(list(cond = substr(ln, op + 1L, p - 1L),
                    rest = .nm_ip_trim(substr(ln, p + 1L, nchar(ln)))))
    }
  }
  list(cond = NA_character_, rest = ln)
}

# Variable names used in an IF condition (Fortran operators removed).
.nm_ip_cond_vars <- function(cond) {
  if (length(cond) == 0L || is.na(cond) || !nzchar(cond)) return(character())
  v <- toupper(regmatches(cond, gregexpr("\\b[A-Za-z][A-Za-z0-9_]*\\b", cond,
                                         perl = TRUE))[[1]])
  unique(setdiff(v, c("EQ", "NE", "GT", "GE", "LT", "LE", "AND", "OR", "NOT",
                      "TRUE", "FALSE", "EQV", "NEQV")))
}

# Parse the assignments of $PK / $ERROR. Each row records the variables of the
# IF condition(s) under which the assignment is made (one-line IF or
# IF ... THEN / ELSE / ENDIF blocks) in `condition_vars`, and `conditional`.
.nm_ip_parse_assignments <- function(lines, section = "$PK") {
  z <- .nm_ip_section_lines(lines, section)
  if (!length(z)) return(data.frame())

  code <- .nm_ip_trim(.nm_ip_strip_inline_comment(z))
  code <- code[nzchar(code)]
  if (!length(code)) return(data.frame())

  reserved <- c("THETA", "ETA", "ERR", "EPS", "EXP", "LOG", "LOG10", "SQRT", "ABS",
                "MIN", "MAX", "IF", "THEN", "AND", "OR", "NOT", "EQ", "NE",
                "GT", "GE", "LT", "LE")
  rows <- list(); k <- 0L
  block <- list()   # condition variables of the open IF ... THEN blocks

  for (i in seq_along(code)) {
    ln <- code[i]
    up <- toupper(ln)

    # Block IF structure.
    if (grepl("^END\\s*IF\\b", up, perl = TRUE)) {
      if (length(block)) block[[length(block)]] <- NULL
      next
    }
    if (grepl("^ELSE\\s*IF\\b", up, perl = TRUE)) {
      sp <- .nm_ip_split_if(sub("^\\s*ELSE\\s*", "", ln, ignore.case = TRUE, perl = TRUE))
      if (length(block))
        block[[length(block)]] <- unique(c(block[[length(block)]], .nm_ip_cond_vars(sp$cond)))
      next
    }
    if (grepl("^ELSE\\s*$", up, perl = TRUE)) next

    sp <- .nm_ip_split_if(ln)
    if (!is.na(sp$cond) && grepl("^THEN\\s*$", toupper(sp$rest), perl = TRUE)) {
      block[[length(block) + 1L]] <- .nm_ip_cond_vars(sp$cond)
      next
    }

    expr <- sp$rest
    if (!grepl("=", expr, fixed = TRUE)) next
    m <- regexec("^\\s*([A-Za-z][A-Za-z0-9_]*)\\s*=\\s*(.+)$", expr, perl = TRUE)
    r <- regmatches(expr, m)[[1]]
    if (length(r) < 3L) next

    lhs <- toupper(r[2]); rhs <- .nm_ip_trim(r[3])

    th <- regmatches(rhs, gregexpr("THETA\\s*\\(\\s*[0-9]+\\s*\\)", rhs,
                                  ignore.case = TRUE, perl = TRUE))[[1]]
    # (?<![A-Za-z]) so that the "ETA(n)" inside "THETA(n)" is not matched.
    et <- regmatches(rhs, gregexpr("(?<![A-Za-z])ETA\\s*\\(\\s*[0-9]+\\s*\\)", rhs,
                                  ignore.case = TRUE, perl = TRUE))[[1]]
    er <- regmatches(rhs, gregexpr("(?:ERR|EPS)\\s*\\(\\s*[0-9]+\\s*\\)", rhs,
                                  ignore.case = TRUE, perl = TRUE))[[1]]

    thn <- if (length(th)) as.integer(gsub("[^0-9]", "", th)) else integer()
    etn <- if (length(et)) as.integer(gsub("[^0-9]", "", et)) else integer()
    ern <- if (length(er)) as.integer(gsub("[^0-9]", "", er)) else integer()

    cond_vars <- unique(c(.nm_ip_cond_vars(sp$cond), unlist(block)))
    vars <- regmatches(rhs, gregexpr("\\b[A-Za-z][A-Za-z0-9_]*\\b", rhs, perl = TRUE))[[1]]
    vars <- setdiff(unique(toupper(c(vars, cond_vars))), reserved)
    vars <- vars[vars != lhs]

    k <- k + 1L
    rows[[k]] <- data.frame(
      section = toupper(section), line_order = i, lhs = lhs, rhs = rhs,
      theta_direct = paste(sort(unique(thn)), collapse = ","),
      eta_direct = paste(sort(unique(etn)), collapse = ","),
      error_direct = paste(sort(unique(ern)), collapse = ","),
      variables_rhs = paste(vars, collapse = ","),
      condition_vars = paste(cond_vars, collapse = ","),
      conditional = !is.na(sp$cond) || length(block) > 0L,
      stringsAsFactors = FALSE
    )
  }

  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

.nm_ip_parse_pk <- function(lines) .nm_ip_parse_assignments(lines, "$PK")
.nm_ip_parse_error <- function(lines) .nm_ip_parse_assignments(lines, "$ERROR")

# ------------------------------------------------------------------------------
# Dependency resolution
# ------------------------------------------------------------------------------

# $PRED models hold the structural model and the residual error in one
# record. Residual-error statements are those using EPS()/ERR(), plus the
# variables they reference directly (e.g. W in Y = IPRED + W*EPS(1)) that do
# not themselves reference a PK parameter or an ETA. The rest is the
# structural part, handled like $PK.
.nm_ip_split_pred <- function(pred) {
  empty <- data.frame()
  if (!nrow(pred)) return(list(pk = empty, error = empty))
  vars <- strsplit(pred$variables_rhs, ",", fixed = TRUE)
  is_eps <- nzchar(pred$error_direct)
  ref <- unique(toupper(unlist(vars[is_eps])))
  helper <- !is_eps & toupper(pred$lhs) %in% ref & !nzchar(pred$eta_direct) &
    !vapply(vars, function(v) any(vapply(v, .nm_ip_is_pk_parameter, logical(1))),
            logical(1))
  err <- is_eps | helper
  rownames(pred) <- NULL
  list(pk = pred[!err, , drop = FALSE], error = pred[err, , drop = FALSE])
}

.nm_ip_resolve_pk_dependencies <- function(pk) {
  if (!nrow(pk)) return(pk)

  # A variable may be defined several times (e.g. TVCL = THETA(1) then
  # IF(SEX.EQ.1) TVCL = TVCL*THETA(4)): all its definitions are pooled, so that
  # every THETA/ETA that can affect its final value is kept.
  idx_by <- split(seq_len(nrow(pk)), pk$lhs)
  defs <- pk[vapply(idx_by, max, integer(1)), , drop = FALSE]
  pool <- function(x, numeric_sort = FALSE) {
    v <- unique(unlist(strsplit(x[nzchar(x)], ",", fixed = TRUE)))
    if (numeric_sort && length(v)) v <- sort(as.integer(v))
    paste(v, collapse = ",")
  }
  for (q in seq_len(nrow(defs))) {
    rr <- idx_by[[defs$lhs[q]]]
    if (length(rr) > 1L) {
      defs$theta_direct[q]  <- pool(pk$theta_direct[rr], TRUE)
      defs$eta_direct[q]    <- pool(pk$eta_direct[rr], TRUE)
      defs$variables_rhs[q] <- pool(pk$variables_rhs[rr])
    }
  }
  rownames(defs) <- defs$lhs

  memo <- new.env(parent = emptyenv())

  resolve_one <- function(v, stack = character()) {
    if (exists(v, envir = memo, inherits = FALSE)) return(get(v, envir = memo))
    if (!(v %in% rownames(defs)))
      return(list(theta = integer(), eta = integer(), vars = character()))
    if (v %in% stack)
      return(list(theta = integer(), eta = integer(), vars = character()))

    row <- defs[v, , drop = FALSE]
    th <- .nm_ip_split_ints(row$theta_direct)
    et <- .nm_ip_split_ints(row$eta_direct)
    depvars <- .nm_ip_split_vars(row$variables_rhs)
    allvars <- depvars

    for (dv in depvars) {
      if (dv %in% rownames(defs)) {
        rr <- resolve_one(dv, c(stack, v))
        th <- c(th, rr$theta)
        et <- c(et, rr$eta)
        allvars <- c(allvars, rr$vars)
      }
    }

    ans <- list(theta = sort(unique(th)), eta = sort(unique(et)), vars = unique(allvars))
    assign(v, ans, envir = memo)
    ans
  }

  resolved <- lapply(defs$lhs, resolve_one)
  defs$theta_all <- vapply(resolved, function(x) paste(x$theta, collapse = ","), character(1))
  defs$eta_all <- vapply(resolved, function(x) paste(x$eta, collapse = ","), character(1))
  defs$variables_all <- vapply(resolved, function(x) paste(x$vars, collapse = ","), character(1))
  defs
}

# ------------------------------------------------------------------------------
# PK parameter heuristics
# ------------------------------------------------------------------------------

.nm_ip_is_pk_parameter <- function(name) {
  n <- toupper(name)
  patterns <- c(
    "^CL[0-9A-Z_]*$", "^V[0-9A-Z_]*$", "^Q[0-9A-Z_]*$", "^KA[0-9A-Z_]*$",
    "^K[0-9]+$", "^K[0-9]+[0-9]+$", "^KEL$", "^KE$", "^KM$", "^VMAX$",
    "^F[0-9A-Z_]*$", "^ALAG[0-9]+$", "^D[0-9]+$", "^R[0-9]+$",
    "^S[0-9]+$", "^BIO[0-9A-Z_]*$", "^MTT$"
  )
  any(vapply(patterns, function(p) grepl(p, n, perl = TRUE), logical(1))) &&
    !grepl("^TV", n)
}

# Data items declared in $INPUT (aliases kept, DROP/SKIP items removed).
.nm_ip_parse_input <- function(lines) {
  z <- .nm_ip_section_lines(lines, "$INPUT")
  if (!length(z)) return(character())
  z <- .nm_ip_trim(.nm_ip_strip_inline_comment(z[!.nm_ip_is_full_comment(z)]))
  toks <- unlist(strsplit(paste(z, collapse = " "), "[[:space:],]+"))
  toks <- toupper(toks[nzchar(toks)])
  out <- character()
  for (tk in toks) {
    parts <- strsplit(tk, "=", fixed = TRUE)[[1]]
    if (any(parts %in% c("DROP", "SKIP"))) next
    out <- c(out, parts)
  }
  unique(out)
}

# Wings for NONMEM extended syntax: model records may refer to parameters by
# the label written in the $THETA/$OMEGA/$SIGMA comments (e.g. POPCL, eta_cl)
# instead of THETA(n)/ETA(n)/ERR(n). Labels are replaced by their THETA(n),
# ETA(n) or ERR(n) reference before parsing, reusing the nm_convert_ctl()
# mapping. A label is substituted only if it is used in a model record, is
# never assigned there, and is neither a $INPUT item nor a NONMEM item, so
# standard control streams are left unchanged.
.nm_ip_expand_labels <- function(lines) {
  none <- data.frame(label = character(), reference = character(),
                     stringsAsFactors = FALSE)
  lines <- gsub("\r", "", lines, fixed = TRUE)
  map <- c(.nm_build_mapping(lines, "^\\s*\\$THETA", "THETA"),
           .nm_build_eta_mapping(lines),
           .nm_build_err_mapping(lines))
  if (!length(map)) return(list(lines = lines, substitutions = none))
  map <- map[!duplicated(toupper(names(map)))]

  rec <- toupper(sub("^\\s*\\$([A-Za-z]+).*$", "\\1", lines))
  is_rec <- grepl("^\\s*\\$[A-Za-z]", lines)
  cur <- cumsum(is_rec)
  rec_of <- ifelse(cur > 0, rec[is_rec][pmax(cur, 1L)], "")
  model <- rec_of %in% c("PK", "PRED", "ERROR", "DES", "AES", "MIX") &
    !grepl("^\\s*;", lines)
  if (!any(model)) return(list(lines = lines, substitutions = none))

  code <- sub(";.*$", "", lines[model])
  code[is_rec[model]] <- sub("^\\s*\\$[A-Za-z]+", "", code[is_rec[model]])
  lhs <- unlist(regmatches(code, gregexpr("[A-Za-z_][A-Za-z0-9_]*(?=\\s*=(?!=))",
                                          code, perl = TRUE)))
  words <- unlist(regmatches(code, gregexpr("[A-Za-z_][A-Za-z0-9_]*", code)))
  excl <- c(toupper(lhs), .nm_ip_parse_input(lines), .nm_ip_nonmem_items,
            "THETA", "ETA", "ERR", "EPS", "EXP", "LOG", "SQRT")
  keep <- toupper(names(map)) %in% toupper(words) &
    !toupper(names(map)) %in% excl &
    grepl("^[A-Za-z_][A-Za-z0-9_]*$", names(map))
  map <- map[keep]
  if (!length(map)) return(list(lines = lines, substitutions = none))

  idx <- which(model)
  for (k in seq_along(idx)) {
    ln <- lines[idx[k]]
    sc <- regexpr(";", ln, fixed = TRUE)
    cp <- if (sc > 0) substring(ln, 1L, sc - 1L) else ln
    cm <- if (sc > 0) substring(ln, sc) else ""
    for (nm in names(map))
      cp <- gsub(paste0("(?<![A-Za-z0-9_])", nm, "(?![A-Za-z0-9_(])"),
                 map[[nm]], cp, perl = TRUE, ignore.case = TRUE)
    lines[idx[k]] <- paste0(cp, cm)
  }
  list(lines = lines,
       substitutions = data.frame(label = names(map),
                                  reference = unlist(map, use.names = FALSE),
                                  stringsAsFactors = FALSE))
}

# NONMEM-defined data items and system variables: never covariates.
.nm_ip_nonmem_items <- c(
  "ID", "L1", "L2", "TIME", "DATE", "DAT1", "DAT2", "DAT3", "DV", "MDV",
  "EVID", "AMT", "RATE", "SS", "II", "ADDL", "CMT", "PCMT", "CALL", "CONT",
  "XVID1", "XVID2", "XVID3", "XVID4", "XVID5", "NEWIND", "ICALL", "IREP",
  "NIREC", "NDREC", "T", "MTIME", "MNEXT", "MPAST", "MNOW", "NEWL2",
  "F", "A", "IPRED", "PRED", "RES", "WRES")

# Build the covariate predicate used by the role heuristics.
#   - With a $INPUT record: a covariate is any declared data item that is not
#     a NONMEM-defined item (so unusual names such as GENO are recognised).
#   - Without $INPUT: a symbol used in $PK but not defined there must come
#     from the data set, so any such symbol that is not a NONMEM item counts.
.nm_ip_covariate_predicate <- function(input_vars) {
  if (length(input_vars)) {
    covs <- setdiff(input_vars, .nm_ip_nonmem_items)
    function(v) toupper(v) %in% covs
  } else {
    function(v) !(toupper(v) %in% .nm_ip_nonmem_items)
  }
}

.nm_ip_known_data_like <- function(v) {
  # Common NONMEM/data/covariate symbols. Unknown terminal RHS symbols are also
  # treated as data-like by the covariate heuristic below.
  v <- toupper(v)
  v %in% c("AGE", "WT", "WGT", "WEIGHT", "HT", "HEIGHT", "BMI", "IMC", "SEX",
           "RACE", "BSA", "CRCL", "CREAT", "CREA", "SCR", "ALB", "AST", "ALT",
           "BILI", "DOSE", "AMT", "RATE", "TIME", "DAY", "GROUP", "GROUPE",
           "CMT", "EVID", "OCC", "STUDY", "STUD", "COPPM")
}

.nm_ip_theta_direct_context <- function(pk, theta_index) {
  if (!nrow(pk)) return(data.frame())
  hit <- vapply(pk$theta_direct, function(x) theta_index %in% .nm_ip_split_ints(x), logical(1))
  pk[hit, , drop = FALSE]
}

.nm_ip_eta_direct_context <- function(pk, eta_index) {
  if (!nrow(pk)) return(data.frame())
  hit <- vapply(pk$eta_direct, function(x) eta_index %in% .nm_ip_split_ints(x), logical(1))
  pk[hit, , drop = FALSE]
}

# Extract every EXP(...) term of an expression, with balanced parentheses
# (e.g. "EXP(ETA(2)*THETA(5))" is returned whole).
.nm_ip_exp_terms <- function(rhs) {
  st <- gregexpr("(?<![A-Za-z])EXP\\s*\\(", rhs, ignore.case = TRUE, perl = TRUE)[[1]]
  if (st[1] == -1L) return(character())
  len <- attr(st, "match.length")
  ch <- strsplit(rhs, "", fixed = TRUE)[[1]]
  out <- character()
  for (k in seq_along(st)) {
    depth <- 0L
    for (p in (st[k] + len[k] - 1L):length(ch)) {
      if (ch[p] == "(") depth <- depth + 1L
      if (ch[p] == ")") {
        depth <- depth - 1L
        if (depth == 0L) { out <- c(out, substr(rhs, st[k], p)); break }
      }
    }
  }
  out
}

# Split an expression into its top-level additive terms (on + and - outside
# parentheses; the sign of an exponent such as 1E-3 is not a separator).
.nm_ip_additive_terms <- function(x) {
  ch <- strsplit(x, "", fixed = TRUE)[[1]]
  depth <- 0L; cur <- ""; out <- character()
  for (p in seq_along(ch)) {
    c1 <- ch[p]
    if (c1 == "(") depth <- depth + 1L
    if (c1 == ")") depth <- depth - 1L
    is_exp_sign <- p > 2L && toupper(ch[p - 1L]) == "E" && grepl("[0-9.]", ch[p - 2L])
    if (depth == 0L && c1 %in% c("+", "-") && !is_exp_sign) {
      out <- c(out, cur); cur <- ""
    } else {
      cur <- paste0(cur, c1)
    }
  }
  out <- c(out, cur)
  out[nzchar(trimws(out))]
}

.nm_ip_theta_in_iiv_term <- function(rhs, theta_index) {
  # A THETA is an IIV scale when it multiplies an ETA inside an EXP() term,
  # i.e. both appear in the same additive term of the EXP argument.
  # Examples:
  #   CL = THETA(2)*EXP(ETA(1)*THETA(6))           -> THETA(6) = iiv_scale
  #   CL = THETA(1)*EXP(THETA(4)*(WT-70) + ETA(1)) -> THETA(4) is not
  thpat <- paste0("THETA\\s*\\(\\s*", theta_index, "\\s*\\)")
  exps <- .nm_ip_exp_terms(rhs)
  if (!length(exps)) return(FALSE)
  any(vapply(exps, function(e) {
    inner <- sub("^[^(]*\\(", "", sub("\\)\\s*$", "", e))
    any(vapply(.nm_ip_additive_terms(inner), function(t) {
      grepl("(?<![A-Za-z])ETA\\s*\\(", t, ignore.case = TRUE, perl = TRUE) &&
        grepl(thpat, t, ignore.case = TRUE, perl = TRUE)
    }, logical(1)))
  }, logical(1)))
}

# A THETA is an exponential covariate effect when, inside an EXP() term, it
# shares an additive term with a covariate (data-like variable, possibly via
# helper equations) and no ETA. Examples:
#   CL = THETA(1)*EXP(THETA(4)*WT)
#   CL = THETA(1)*EXP(THETA(4)*(WT-70) + ETA(1))      (mu-referencing)
#   WTC = WT-70 ; CL = THETA(1)*EXP(THETA(4)*WTC)       (helper variable)
.nm_ip_theta_exp_covariate <- function(rhs, theta_index, pk_defs, pk_resolved,
                                       is_cov = .nm_ip_known_data_like) {
  thpat <- paste0("THETA\\s*\\(\\s*", theta_index, "\\s*\\)")
  reserved <- c("THETA", "ETA", "ERR", "EPS", "EXP", "LOG", "LOG10", "SQRT",
                "ABS", "MIN", "MAX", "IF", "THEN", "AND", "OR", "NOT", "EQ",
                "NE", "GT", "GE", "LT", "LE")
  for (e in .nm_ip_exp_terms(rhs)) {
    inner <- sub("^[^(]*\\(", "", sub("\\)\\s*$", "", e))
    for (t in .nm_ip_additive_terms(inner)) {
      if (!grepl(thpat, t, ignore.case = TRUE, perl = TRUE)) next
      if (grepl("(?<![A-Za-z])ETA\\s*\\(", t, ignore.case = TRUE, perl = TRUE)) next
      vars <- toupper(regmatches(t, gregexpr("\\b[A-Za-z][A-Za-z0-9_]*\\b", t,
                                             perl = TRUE))[[1]])
      vars <- setdiff(vars, reserved)
      if (!length(vars)) next
      expanded <- vars
      for (v in intersect(vars, pk_resolved$lhs)) {
        rr <- pk_resolved[pk_resolved$lhs == v, , drop = FALSE]
        expanded <- c(expanded, .nm_ip_split_vars(rr$variables_all[1]))
      }
      terminal <- setdiff(unique(expanded), c(pk_defs, "F", "A", "IPRED", "PRED"))
      if (length(terminal) && any(vapply(terminal, is_cov, logical(1))))
        return(TRUE)
    }
  }
  FALSE
}

# Variables reached from helper variables through the right-hand sides of
# their definitions only (IF conditions are not followed). Recursive; the
# pooled definitions of pk_resolved are used.
.nm_ip_rhs_vars <- function(vars, pk_resolved, seen = character()) {
  reserved <- c("THETA", "ETA", "ERR", "EPS", "EXP", "LOG", "LOG10", "SQRT",
                "ABS", "MIN", "MAX")
  out <- character()
  for (v in setdiff(vars, seen)) {
    rr <- pk_resolved[pk_resolved$lhs == v, , drop = FALSE]
    if (!nrow(rr)) next
    rhs <- paste(rr$rhs, collapse = " ")
    w <- toupper(regmatches(rhs, gregexpr("\\b[A-Za-z][A-Za-z0-9_]*\\b", rhs, perl = TRUE))[[1]])
    w <- setdiff(unique(w), reserved)
    seen <- c(seen, v)
    out <- c(out, w, .nm_ip_rhs_vars(intersect(w, pk_resolved$lhs), pk_resolved, seen))
  }
  unique(out)
}

# A THETA is a linear covariate effect when it is multiplied by a covariate
# inside a sum, outside any EXP(). The sum is the innermost parenthesised
# group (or the whole expression) containing the THETA, and it must have at
# least two additive terms. Examples:
#   CL = THETA(1)*(1 + THETA(4)*(WT-70))      -> THETA(4) covariate
#   V  = (THETA(2) + THETA(6)*SEX)            -> THETA(6) covariate,
#                                                THETA(2) typical
#   KA = THETA(3)*(WT/70)**0.75               -> THETA(3) typical (no sum)
.nm_ip_theta_linear_covariate <- function(rhs, theta_index, pk_defs, pk_resolved,
                                          is_cov = .nm_ip_known_data_like) {
  thpat <- paste0("(?<![A-Za-z])THETA\\s*\\(\\s*", theta_index, "\\s*\\)")
  hits <- gregexpr(thpat, rhs, ignore.case = TRUE, perl = TRUE)[[1]]
  if (hits[1] == -1L) return(FALSE)
  ch <- strsplit(rhs, "", fixed = TRUE)[[1]]
  n <- length(ch)
  reserved <- c("THETA", "ETA", "ERR", "EPS", "EXP", "LOG", "LOG10", "SQRT",
                "ABS", "MIN", "MAX", "IF", "THEN", "AND", "OR", "NOT", "EQ",
                "NE", "GT", "GE", "LT", "LE")

  for (h in hits) {
    # Open parentheses enclosing position h, innermost last.
    stack <- integer()
    for (p in seq_len(h - 1L)) {
      if (ch[p] == "(") stack <- c(stack, p)
      if (ch[p] == ")" && length(stack)) stack <- stack[-length(stack)]
    }
    # The sum is the innermost enclosing group with at least two additive
    # terms: a group holding a single product, e.g. (THETA(8)*(CLCR/89)) in
    # ((THETA(1)*FALL)+(THETA(8)*(CLCR/89))), is passed over. Stop at EXP(),
    # which is handled elsewhere.
    found <- FALSE
    for (lev in c(rev(seq_along(stack)), 0L)) {
      if (lev > 0L) {
        op <- stack[lev]
        depth <- 0L; cl <- n
        for (p in op:n) {
          if (ch[p] == "(") depth <- depth + 1L
          if (ch[p] == ")") { depth <- depth - 1L; if (depth == 0L) { cl <- p; break } }
        }
        before <- if (op > 1L) paste(ch[1:(op - 1L)], collapse = "") else ""
        if (grepl("EXP\\s*$", before, ignore.case = TRUE)) break
        g_start <- op + 1L; g_end <- cl - 1L
      } else {
        g_start <- 1L; g_end <- n
      }
      if (g_end < g_start) next

      # Split the group into top-level additive terms, keeping their positions.
      depth <- 0L; starts <- g_start; ends <- integer()
      for (p in g_start:g_end) {
        if (ch[p] == "(") depth <- depth + 1L
        if (ch[p] == ")") depth <- depth - 1L
        is_exp_sign <- p > 2L && toupper(ch[p - 1L]) == "E" && grepl("[0-9.]", ch[p - 2L])
        if (depth == 0L && ch[p] %in% c("+", "-") && !is_exp_sign && p > g_start) {
          ends <- c(ends, p - 1L); starts <- c(starts, p + 1L)
        }
      }
      ends <- c(ends, g_end)
      terms <- mapply(function(a, b) if (b >= a) paste(ch[a:b], collapse = "") else "",
                      starts, ends)
      if (sum(nzchar(trimws(terms))) >= 2L) { found <- TRUE; break }
    }
    if (!found) next

    k <- which(starts <= h & ends >= h)[1]
    if (is.na(k)) next
    t <- terms[k]
    if (grepl("(?<![A-Za-z])ETA\\s*\\(", t, ignore.case = TRUE, perl = TRUE)) next

    vars <- toupper(regmatches(t, gregexpr("\\b[A-Za-z][A-Za-z0-9_]*\\b", t,
                                           perl = TRUE))[[1]])
    vars <- setdiff(vars, reserved)
    if (!length(vars)) next
    # Helper variables are followed through the variables of their
    # right-hand sides only, not through the conditions of IF blocks: in
    # THETA(1)*FALL with IF(CT01.EQ.1) FALL = THETA(7) ELSE FALL = 1, THETA(1)
    # is a typical value (the CT01 effect is THETA(7), found by the
    # conditional rule), whereas in THETA(6)*RATIO with
    # RATIO = C_RES/(THETA(4)+C_RES), THETA(6) is a covariate effect.
    expanded <- c(vars, .nm_ip_rhs_vars(intersect(vars, pk_resolved$lhs), pk_resolved))
    terminal <- setdiff(unique(expanded), pk_defs)
    if (length(terminal) && any(vapply(terminal, is_cov, logical(1)))) return(TRUE)
  }
  FALSE
}

# Saturable (Emax/IC50-type) covariate effect: the THETA is an additive term of
# a denominator whose other terms contain a covariate, e.g.
#   RATIO = C_RES/(THETA(4) + C_RES)      CL = THETA(2)*(1 - RATIO)
#   TVCL  = THETA(1)*(1 + THETA(5)*WT/(THETA(6) + WT))   -> THETA(6)
.nm_ip_theta_saturable_covariate <- function(rhs, theta_index, pk_defs, pk_resolved,
                                             is_cov = .nm_ip_known_data_like) {
  thpat <- paste0("(?<![A-Za-z])THETA\\s*\\(\\s*", theta_index, "\\s*\\)")
  hits <- gregexpr(thpat, rhs, ignore.case = TRUE, perl = TRUE)[[1]]
  if (hits[1] == -1L) return(FALSE)
  ch <- strsplit(rhs, "", fixed = TRUE)[[1]]
  n <- length(ch)
  reserved <- c("THETA", "ETA", "ERR", "EPS", "EXP", "LOG", "LOG10", "SQRT",
                "ABS", "MIN", "MAX")
  for (h in hits) {
    stack <- integer()
    for (p in seq_len(h - 1L)) {
      if (ch[p] == "(") stack <- c(stack, p)
      if (ch[p] == ")" && length(stack)) stack <- stack[-length(stack)]
    }
    if (!length(stack)) next
    op <- stack[length(stack)]
    before <- if (op > 1L) trimws(paste(ch[1:(op - 1L)], collapse = "")) else ""
    if (!grepl("/$", before)) next                     # not a denominator
    depth <- 0L; cl <- n
    for (p in op:n) {
      if (ch[p] == "(") depth <- depth + 1L
      if (ch[p] == ")") { depth <- depth - 1L; if (depth == 0L) { cl <- p; break } }
    }
    g_start <- op + 1L; g_end <- cl - 1L
    if (g_end < g_start) next
    depth <- 0L; starts <- g_start; ends <- integer()
    for (p in g_start:g_end) {
      if (ch[p] == "(") depth <- depth + 1L
      if (ch[p] == ")") depth <- depth - 1L
      is_exp_sign <- p > 2L && toupper(ch[p - 1L]) == "E" && grepl("[0-9.]", ch[p - 2L])
      if (depth == 0L && ch[p] %in% c("+", "-") && !is_exp_sign && p > g_start) {
        ends <- c(ends, p - 1L); starts <- c(starts, p + 1L)
      }
    }
    ends <- c(ends, g_end)
    k <- which(starts <= h & ends >= h)[1]
    if (is.na(k) || length(starts) < 2L) next
    others <- paste(mapply(function(a, b) if (b >= a) paste(ch[a:b], collapse = "") else "",
                           starts[-k], ends[-k]), collapse = " ")
    if (grepl("(?<![A-Za-z])ETA\\s*\\(", others, ignore.case = TRUE, perl = TRUE)) next
    vars <- toupper(regmatches(others, gregexpr("\\b[A-Za-z][A-Za-z0-9_]*\\b", others,
                                                perl = TRUE))[[1]])
    vars <- setdiff(vars, reserved)
    if (!length(vars)) next
    expanded <- vars
    for (v in intersect(vars, pk_resolved$lhs)) {
      rr <- pk_resolved[pk_resolved$lhs == v, , drop = FALSE]
      expanded <- c(expanded, .nm_ip_split_vars(rr$variables_all[1]))
    }
    terminal <- setdiff(unique(expanded), pk_defs)
    if (length(terminal) && any(vapply(terminal, is_cov, logical(1)))) return(TRUE)
  }
  FALSE
}

.nm_ip_theta_role_for_parameter <- function(theta_index, parameter, pk, pk_resolved, error,
                                            is_cov = .nm_ip_known_data_like) {
  err_hit <- if (nrow(error)) any(vapply(error$theta_direct,
                                        function(x) theta_index %in% .nm_ip_split_ints(x), logical(1))) else FALSE

  ctx <- .nm_ip_theta_direct_context(pk, theta_index)
  if (!nrow(ctx)) return(if (err_hit) "residual_error" else "structural_other")

  # Only a THETA located inside an ETA-bearing IIV term is an IIV scale.
  if (any(vapply(ctx$rhs, .nm_ip_theta_in_iiv_term, logical(1), theta_index = theta_index)))
    return("iiv_scale")

  # Exponential covariate effect, e.g. EXP(THETA(4)*WT) or EXP(THETA(4)*(WT-70)).
  if (any(vapply(ctx$rhs, .nm_ip_theta_exp_covariate, logical(1),
                 theta_index = theta_index, pk_defs = unique(pk$lhs),
                 pk_resolved = pk_resolved, is_cov = is_cov)))
    return("covariate_effect")

  # Linear covariate effect, e.g. THETA(1)*(1 + THETA(4)*(WT-70)).
  if (any(vapply(ctx$rhs, .nm_ip_theta_linear_covariate, logical(1),
                 theta_index = theta_index, pk_defs = unique(pk$lhs),
                 pk_resolved = pk_resolved, is_cov = is_cov)))
    return("covariate_effect")

  # Saturable covariate effect, e.g. C_RES/(THETA(4) + C_RES).
  if (any(vapply(ctx$rhs, .nm_ip_theta_saturable_covariate, logical(1),
                 theta_index = theta_index, pk_defs = unique(pk$lhs),
                 pk_resolved = pk_resolved, is_cov = is_cov)))
    return("covariate_effect")

  # Conditional covariate effect: the THETA only appears in assignments made
  # under an IF condition on a covariate, and it is not the reference value of
  # the assigned variable. The reference is the unconditional definition when
  # there is one, otherwise the first definition in the file. Examples:
  #   TVCL = THETA(1)
  #   IF(SEX.EQ.1) TVCL = TVCL*THETA(4)         -> THETA(4) covariate
  #   IF(GENO.EQ.2) THEN / TVV = TVV*THETA(5) / ENDIF
  #   IF(SEX.EQ.1) THEN / CL = THETA(1) / ELSE / CL = THETA(2) / ENDIF
  #                                             -> THETA(1) reference (typical),
  #                                                THETA(2) covariate
  #   IF(CT01.EQ.1) THEN / FALL = THETA(7) / ELSE / FALL = 1 / ENDIF
  #                                             -> THETA(7) covariate (the
  #                                                constant branch is the reference)
  if ("conditional" %in% names(ctx) && all(ctx$conditional)) {
    for (q in seq_len(nrow(ctx))) {
      cv <- .nm_ip_split_vars(ctx$condition_vars[q])
      expanded <- cv
      for (v in intersect(cv, pk_resolved$lhs)) {
        rr <- pk_resolved[pk_resolved$lhs == v, , drop = FALSE]
        expanded <- c(expanded, .nm_ip_split_vars(rr$variables_all[1]))
      }
      terminal <- setdiff(unique(expanded), unique(pk$lhs))
      same_lhs <- pk$lhs == ctx$lhs[q]
      has_unconditional <- any(same_lhs & !pk$conditional)
      is_first_definition <- ctx$line_order[q] == min(pk$line_order[same_lhs])
      # Another branch assigning a plain constant (no THETA, no variable) is
      # the reference: IF(CT01.EQ.1) THEN / FALL = THETA(7) / ELSE / FALL = 1
      constant_branch <- any(same_lhs & pk$line_order != ctx$line_order[q] &
                               grepl("^\\s*[-+]?[0-9.]+([EeDd][-+]?[0-9]+)?\\s*$", pk$rhs))
      not_reference <- has_unconditional || !is_first_definition || constant_branch
      if (not_reference && length(terminal) &&
          any(vapply(terminal, is_cov, logical(1))))
        return("covariate_effect")
    }
  }

  # Recover terminal variables through helper equations (e.g. GRP <- GROUPE),
  # not only variables appearing directly on the THETA line.
  pk_defs <- unique(pk$lhs)
  ctx_vars <- character()
  for (lhs in unique(ctx$lhs)) {
    rr <- pk_resolved[pk_resolved$lhs == lhs, , drop = FALSE]
    if (nrow(rr)) ctx_vars <- c(ctx_vars, .nm_ip_split_vars(rr$variables_all[1]))
  }
  ctx_vars <- unique(c(ctx_vars, unlist(lapply(ctx$variables_rhs, .nm_ip_split_vars))))
  terminal <- setdiff(ctx_vars, pk_defs)
  terminal <- setdiff(terminal, c("F", "A", "IPRED", "PRED"))
  cov_evidence <- length(terminal) > 0L &&
    any(vapply(terminal, is_cov, logical(1)))

  thpat <- paste0("THETA\\s*\\(\\s*", theta_index, "\\s*\\)")

  # Continuous covariate coefficient, e.g. (AGE/62.5)**THETA(7).
  exponent_coefficient <- any(grepl(
    paste0("(?:\\*\\*|\\^)\\s*", thpat),
    ctx$rhs, ignore.case = TRUE, perl = TRUE))
  if (cov_evidence && exponent_coefficient) return("covariate_effect")

  # Categorical/group effect, e.g. THETA(7)**GRP.
  theta_power_indicator <- any(grepl(
    paste0(thpat, "\\s*(?:\\*\\*|\\^)\\s*[A-Za-z]"),
    ctx$rhs, ignore.case = TRUE, perl = TRUE))
  if (cov_evidence && theta_power_indicator) return("covariate_effect")

  "typical_value"
}

# ------------------------------------------------------------------------------
# Element-wise mapping + residual error table
# ------------------------------------------------------------------------------

.nm_ip_label_matches <- function(label, param) {
  if (length(label) == 0L || is.na(label) || !nzchar(label)) return(NA)
  a <- gsub("[^A-Z0-9]", "", toupper(label))
  b <- gsub("[^A-Z0-9]", "", toupper(param))
  a <- sub("^(POP|TV|ETA|OMEGA)", "", a)
  b <- sub("^(POP|TV|ETA|OMEGA)", "", b)
  if (!nzchar(a) || !nzchar(b)) return(NA)
  identical(a, b) || grepl(b, a, fixed = TRUE) || grepl(a, b, fixed = TRUE)
}

# Number of assignment steps between a THETA/ETA and a variable: 0 when it
# appears directly in one of the variable's definitions (CL = THETA(2)*...),
# 1 when it comes through one helper (CL = TVCL*... with TVCL = THETA(2)),
# and so on. Used to name a THETA/ETA after the parameter closest to it rather
# than after a derived one (KE = CL/V).
.nm_ip_dependency_depth <- function(pk, var, index, type = c("THETA", "ETA")) {
  type <- match.arg(type)
  col <- if (type == "THETA") "theta_direct" else "eta_direct"
  lhs_u <- toupper(pk$lhs)
  walk <- function(v, seen) {
    rows <- which(lhs_u == toupper(v))
    if (!length(rows)) return(Inf)
    if (any(vapply(pk[[col]][rows], function(x) index %in% .nm_ip_split_ints(x),
                   logical(1)))) return(0)
    seen <- c(seen, toupper(v))
    nxt <- setdiff(intersect(toupper(unlist(strsplit(pk$variables_rhs[rows], ",",
                                                      fixed = TRUE))), lhs_u), seen)
    if (!length(nxt)) return(Inf)
    1 + min(vapply(nxt, walk, numeric(1), seen = seen))
  }
  walk(var, character())
}

.nm_ip_build_mapping <- function(theta, omega_eta, pk, pk_resolved, error,
                                 is_cov = .nm_ip_known_data_like) {
  if (!nrow(pk_resolved)) return(data.frame())

  rows <- list(); k <- 0L

  # Prefer final PK-facing variables; keep helper variables in detailed mapping
  # but flag likely PK parameters.
  for (i in seq_len(nrow(pk_resolved))) {
    p <- pk_resolved$lhs[i]
    ths <- .nm_ip_split_ints(pk_resolved$theta_all[i])
    ets <- .nm_ip_split_ints(pk_resolved$eta_all[i])

    for (tt in ths) {
      tr <- theta[theta$theta == tt, , drop = FALSE]
      role <- .nm_ip_theta_role_for_parameter(tt, p, pk, pk_resolved, error,
                                              is_cov = is_cov)
      k <- k + 1L
      rows[[k]] <- data.frame(
        pk_parameter = p,
        likely_pk_parameter = .nm_ip_is_pk_parameter(p),
        element_type = "THETA",
        index = tt,
        symbol = paste0("THETA(", tt, ")"),
        role = role,
        label = if (nrow(tr)) tr$label[1] else NA_character_,
        initial = if (nrow(tr)) tr$initial[1] else NA_real_,
        fixed = if (nrow(tr)) tr$fixed[1] else NA,
        comment_matches_pk = if (nrow(tr)) .nm_ip_label_matches(tr$label[1], p) else NA,
        depth = .nm_ip_dependency_depth(pk, p, tt, "THETA"),
        stringsAsFactors = FALSE
      )
    }

    for (ee in ets) {
      or <- omega_eta[omega_eta$eta == ee, , drop = FALSE]
      k <- k + 1L
      rows[[k]] <- data.frame(
        pk_parameter = p,
        likely_pk_parameter = .nm_ip_is_pk_parameter(p),
        element_type = "ETA",
        index = ee,
        symbol = paste0("ETA(", ee, ")"),
        role = "iiv",
        label = if (nrow(or)) or$omega_label[1] else NA_character_,
        initial = if (nrow(or)) or$omega_initial[1] else NA_real_,
        fixed = if (nrow(or)) or$omega_fixed[1] else NA,
        comment_matches_pk = if (nrow(or)) .nm_ip_label_matches(or$omega_label[1], p) else NA,
        depth = .nm_ip_dependency_depth(pk, p, ee, "ETA"),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) return(data.frame())
  unique(do.call(rbind, rows))
}

.nm_ip_residual_error_mapping <- function(theta, error) {
  if (!nrow(error)) return(data.frame())
  th <- sort(unique(unlist(lapply(error$theta_direct, .nm_ip_split_ints))))
  if (!length(th)) return(data.frame())

  out <- lapply(th, function(tt) {
    tr <- theta[theta$theta == tt, , drop = FALSE]
    ctx <- error[vapply(error$theta_direct,
                        function(x) tt %in% .nm_ip_split_ints(x), logical(1)), , drop = FALSE]
    data.frame(
      theta = tt,
      symbol = paste0("THETA(", tt, ")"),
      role = "residual_error",
      label = if (nrow(tr)) tr$label[1] else NA_character_,
      initial = if (nrow(tr)) tr$initial[1] else NA_real_,
      fixed = if (nrow(tr)) tr$fixed[1] else NA,
      error_variables = paste(unique(ctx$lhs), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

.nm_ip_parameter_summary <- function(mapping, pk = NULL) {
  if (!nrow(mapping)) return(data.frame())

  x <- mapping[!grepl("^TV", mapping$pk_parameter), , drop = FALSE]
  if (!nrow(x)) x <- mapping

  # Only present likely PK-facing variables in the concise table.
  xp <- x[x$likely_pk_parameter, , drop = FALSE]
  if (nrow(xp)) x <- xp

  keys <- unique(x$pk_parameter)

  # A parameter whose definition uses other parameters of the table
  # (S2 = V, KE = CL/V) is derived: `derived_from` lists them (NA otherwise).
  derived_from <- function(p) {
    if (is.null(pk) || !nrow(pk) || !"variables_rhs" %in% names(pk)) return(NA_character_)
    v <- unlist(strsplit(as.character(pk$variables_rhs[pk$lhs == p]), ","))
    v <- unique(trimws(v))
    v <- v[v %in% setdiff(keys, p)]
    if (length(v)) paste(v, collapse = ", ") else NA_character_
  }

  out <- lapply(keys, function(p) {
    d <- x[x$pk_parameter == p, , drop = FALSE]

    get_symbols <- function(role = NULL, type = NULL) {
      q <- d
      if (!is.null(role)) q <- q[q$role %in% role, , drop = FALSE]
      if (!is.null(type)) q <- q[q$element_type %in% type, , drop = FALSE]
      s <- unique(q$symbol)
      s <- s[!is.na(s) & nzchar(s)]
      if (length(s)) paste(s, collapse = ", ") else NA_character_
    }

    data.frame(
      PK_parameter = p,
      THETA_typical = get_symbols("typical_value", "THETA"),
      THETA_covariate = get_symbols("covariate_effect", "THETA"),
      THETA_IIV_scale = get_symbols("iiv_scale", "THETA"),
      ETA_IIV = get_symbols("iiv", "ETA"),
      THETA_other = get_symbols("structural_other", "THETA"),
      derived_from = derived_from(p),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

.nm_ip_find_unmapped <- function(theta, omega_eta, mapping, residual_error) {
  mapped_th <- if (nrow(mapping)) unique(mapping$index[mapping$element_type == "THETA"]) else integer()
  err_th <- if (nrow(residual_error)) residual_error$theta else integer()
  mapped_et <- if (nrow(mapping)) unique(mapping$index[mapping$element_type == "ETA"]) else integer()

  mapped_th <- sort(unique(c(mapped_th, err_th)))
  mapped_et <- sort(unique(mapped_et))

  list(
    theta = if (nrow(theta)) theta[!theta$theta %in% mapped_th, , drop = FALSE] else data.frame(),
    omega = if (nrow(omega_eta)) omega_eta[!omega_eta$eta %in% mapped_et, , drop = FALSE] else data.frame()
  )
}

# ------------------------------------------------------------------------------
# Main API
# ------------------------------------------------------------------------------

#' Identify the semantic role of THETA and ETA in a NONMEM control stream
#'
#' Parses `$THETA`, `$OMEGA`, `$PK` and `$ERROR` of a control stream (or
#' `$PRED` when there is no `$PK`) and links
#' each THETA and ETA to the PK parameter(s) it acts on, with a semantic role:
#' `typical_value`, `covariate_effect`, `iiv_scale`, `residual_error`, `iiv`
#' or `structural_other`. Executable equations are the primary source;
#' comments are only used as labels and consistency checks.
#'
#' Covariate effects are recognised in power, exponential, linear,
#' saturable and conditional forms (e.g. `(WT/70)**THETA(4)`,
#' `EXP(THETA(4)*WT)`, `THETA(1)*(1 + THETA(4)*(WT-70))`,
#' `C_RES/(THETA(4) + C_RES)`, `IF(SEX.EQ.1) TVCL = TVCL*THETA(4)`,
#' including `IF ... THEN ... ELSE ... ENDIF` blocks). A variable defined
#' several times keeps the THETA/ETA of all its definitions. When a variable
#' is only defined under conditions (e.g. `IF(SEX.EQ.1) CL = THETA(1)` /
#' `IF(SEX.EQ.0) CL = THETA(2)`), the THETA of its first definition in the
#' file is the reference (`typical_value`) and the others are covariate
#' effects. A covariate is any data item declared
#' in `$INPUT` that is not a NONMEM-defined item (`ID`, `TIME`, `DV`, `AMT`,
#' `MDV`, `EVID`, `CMT`...), whatever its name; items dropped with `=DROP`
#' are ignored. Without `$INPUT`, any symbol used in `$PK` but not defined
#' there is taken as a data item.
#'
#' Wings for NONMEM extended control streams, where model records use the
#' labels of the `$THETA`, `$OMEGA` and `$SIGMA` comments (e.g.
#' `CL = POPCL*EXP(eta_cl)`), are handled: such labels are first replaced by
#' `THETA(n)`, `ETA(n)` or `ERR(n)`, as in [nm_convert_ctl()]. A label is
#' replaced only if it is used but never assigned in a model record and is
#' not a `$INPUT` item.
#'
#' @param file Path to a NONMEM control stream (`.ctl`).
#' @param likely_pk_only If `TRUE` (default), `mapping_detail` keeps only
#'   variables that look like PK parameters (CL, V, Q, KA, F, ALAG...).
#'   `mapping_detail_all` always keeps every variable.
#'
#' @return An object of class `nm_parameter_map`, a list with, among others:
#'   `parameter_table` (one row per PK parameter with its THETA/ETA by role;
#'   `derived_from` names the parameters of the table that a derived
#'   parameter is computed from, e.g. `"V"` for `S2 = V` or `"CL, V"` for
#'   `KE = CL/V`, and is `NA` for the others),
#'   `mapping_detail` and `mapping_detail_all` (one row per PK parameter /
#'   THETA or ETA association; `depth` is the number of assignment steps
#'   between the THETA/ETA and the parameter, 0 when it appears directly in
#'   its definition), `residual_error`, `theta`, `omega`,
#'   `omega_elements` (variances and BLOCK covariances), `unmapped_theta`,
#'   `unmapped_omega`, `missing_theta_declaration`,
#'   `missing_omega_declaration` and `label_substitutions` (the labels
#'   replaced and their `THETA(n)`/`ETA(n)`/`ERR(n)` reference; empty for a
#'   standard control stream).
#'
#' @examples
#' \dontrun{
#' x <- nm_identify_parameters("run.ctl")
#' x$parameter_table
#' x$mapping_detail
#' x$omega_elements
#' }
#' @seealso [nm_identify_directory()], [nm_extract_runs()]
#' @export
nm_identify_parameters <- function(file, likely_pk_only = TRUE) {
  stopifnot(length(file) == 1L, file.exists(file))
  lines <- readLines(file, warn = FALSE)
  lab <- .nm_ip_expand_labels(lines)
  lines <- lab$lines

  theta <- .nm_ip_parse_theta(lines)
  om <- .nm_ip_parse_omega(lines)
  pk <- .nm_ip_parse_pk(lines)
  error <- .nm_ip_parse_error(lines)
  if (!nrow(pk)) {
    sp <- .nm_ip_split_pred(.nm_ip_parse_assignments(lines, "$PRED"))
    pk <- sp$pk
    if (!nrow(error)) error <- sp$error
  }
  pk_resolved <- .nm_ip_resolve_pk_dependencies(pk)

  is_cov <- .nm_ip_covariate_predicate(.nm_ip_parse_input(lines))
  mapping <- .nm_ip_build_mapping(theta, om$eta, pk, pk_resolved, error,
                                  is_cov = is_cov)
  if (likely_pk_only && nrow(mapping))
    mapping_view <- mapping[mapping$likely_pk_parameter, , drop = FALSE]
  else
    mapping_view <- mapping

  summary <- .nm_ip_parameter_summary(mapping, pk)
  residual <- .nm_ip_residual_error_mapping(theta, error)
  unmapped <- .nm_ip_find_unmapped(theta, om$eta, mapping, residual)

  eta_pk <- if (nrow(mapping)) sort(unique(mapping$index[mapping$element_type == "ETA"])) else integer()
  eta_omega <- if (nrow(om$eta)) om$eta$eta else integer()
  th_pk <- if (nrow(mapping)) sort(unique(mapping$index[mapping$element_type == "THETA"])) else integer()
  th_decl <- if (nrow(theta)) theta$theta else integer()

  structure(list(
    file = normalizePath(file, winslash = "/", mustWork = FALSE),
    parameter_table = summary,
    mapping_detail = mapping_view,
    mapping_detail_all = mapping,
    residual_error = residual,
    theta = theta,
    omega = om$eta,
    omega_elements = om$elements,
    pk_assignments = pk,
    error_assignments = error,
    pk_resolved = pk_resolved,
    unmapped_theta = unmapped$theta,
    unmapped_omega = unmapped$omega,
    missing_theta_declaration = setdiff(th_pk, th_decl),
    missing_omega_declaration = setdiff(eta_pk, eta_omega),
    label_substitutions = lab$substitutions
  ), class = "nm_parameter_map")
}

#' @export
print.nm_parameter_map <- function(x, ...) {
  cat("NONMEM semantic parameter identification\n")
  cat("File:", x$file, "\n\n")
  cat("PK parameter convergence table:\n")
  print(x$parameter_table, row.names = FALSE)

  if (nrow(x$residual_error)) {
    cat("\nResidual-error THETA:\n")
    print(x$residual_error, row.names = FALSE)
  }

  if (nrow(x$omega_elements) && any(x$omega_elements$element_type == "covariance")) {
    cat("\nOMEGA covariance elements:\n")
    print(x$omega_elements[x$omega_elements$element_type == "covariance", , drop = FALSE],
          row.names = FALSE)
  }

  if (nrow(x$unmapped_theta)) {
    cat("\nTHETA not classified:\n")
    print(x$unmapped_theta[, c("theta", "label", "initial", "fixed"), drop = FALSE],
          row.names = FALSE)
  }

  if (nrow(x$unmapped_omega)) {
    cat("\nETA/OMEGA not linked to a PK parameter:\n")
    print(x$unmapped_omega[, c("eta", "omega_label", "omega_initial", "omega_fixed",
                               "omega_structure"), drop = FALSE], row.names = FALSE)
  }

  invisible(x)
}

#' Identify THETA/ETA roles for all control streams of a directory
#'
#' Applies [nm_identify_parameters()] to every `.ctl` file found in `path` and
#' stacks the parameter tables. A file that cannot be parsed gives one row with
#' its error message.
#'
#' @param path Directory to search.
#' @param recursive Search subdirectories (default `TRUE`).
#' @param likely_pk_only Passed to [nm_identify_parameters()].
#'
#' @return A data.frame with a `file` column followed by the columns of
#'   `parameter_table` (or an `error` column for files that failed).
#' @export
nm_identify_directory <- function(path, recursive = TRUE, likely_pk_only = TRUE) {
  files <- list.files(path, pattern = "\\.ctl$", recursive = recursive,
                      full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("No .ctl files found in: ", path)

  out <- lapply(files, function(f) {
    z <- tryCatch(nm_identify_parameters(f, likely_pk_only = likely_pk_only),
                  error = function(e) e)
    if (inherits(z, "error"))
      return(data.frame(file = f, error = conditionMessage(z), stringsAsFactors = FALSE))

    tab <- z$parameter_table
    if (!nrow(tab)) return(data.frame(file = f, stringsAsFactors = FALSE))
    cbind(file = f, tab, stringsAsFactors = FALSE)
  })

  cols <- unique(unlist(lapply(out, names)))
  out <- lapply(out, function(d) {
    miss <- setdiff(cols, names(d))
    for (m in miss) d[[m]] <- NA
    d[, cols, drop = FALSE]
  })
  do.call(rbind, out)
}
