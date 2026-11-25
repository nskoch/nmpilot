# Shared readers of the control stream records $INPUT, $THETA, $OMEGA and
# $SIGMA (action 11c). They replace the readers written separately for
# nm_ctl_check and the transformations (0.3.0), nm_identify_parameters,
# nm_convert_ctl and nm_profile, which differed on edge cases (several THETAs
# or ETAs on one line, BLOCK SAME, labels NAME=(...), DROP without label).

# -----------------------------------------------------------------------------
# $INPUT
# -----------------------------------------------------------------------------

# Names of the data items of $INPUT (all records, comments removed),
# upper case, without the dropped ones (ID=DROP, DROP=ID, DROP or SKIP);
# both names of an alias (DV=CONC) are kept.
.nm_ctl_input_items <- function(lines) {
  txt <- paste(.nm_ctl_records_text(lines, "INPUT"), collapse = " ")
  tok <- unlist(strsplit(trimws(txt), "[,[:space:]]+"))
  tok <- toupper(tok[nzchar(tok)])
  out <- character(0)
  for (t in tok) {
    p <- strsplit(t, "=", fixed = TRUE)[[1]]
    if (any(p %in% c("DROP", "SKIP"))) next
    out <- c(out, p[nzchar(p)])
  }
  unique(out)
}

# -----------------------------------------------------------------------------
# $THETA
# -----------------------------------------------------------------------------

.nm_ctl_theta_tok <- paste0(
  "(\\([^)]*\\)|(?<![A-Za-z0-9_.])[-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[EeDd][-+]?[0-9]+)?)",
  "(\\s*FIX(?:ED)?\\b)?(\\s*[xX]\\s*([0-9]+))?")

# One row per THETA, over all $THETA records, in order: index `theta`, line
# number `line`, position of its text in the line (`start`, `end`), the text
# `token`, `lower`, `initial`, `upper` (NA when absent or infinite), `fixed`,
# the comment of its line (`comment`, leading ";" removed, "" if none), the raw
# line `raw`, `record_line` (TRUE on the line of the record name) and
# `n_on_line` (number of THETAs written on the same line). Labels NAME=(...)
# (NONMEM 7.4) and the options UNINT, NUMBERPOINTS=n, ABSOLUTE, NOABORT are
# skipped; (...)xN gives N THETAs.
.nm_ctl_theta_parse <- function(lines) {
  starts <- grep("^\\s*\\$THETA\\b", lines, ignore.case = TRUE, perl = TRUE)
  ctrl <- grep("^\\s*\\$[A-Za-z]", lines)
  num <- function(x) {
    x <- trimws(x)
    if (!nzchar(x) || grepl("INF", x, ignore.case = TRUE)) return(NA_real_)
    suppressWarnings(as.numeric(gsub("[Dd]", "E", x)))
  }
  rows <- list(); k <- 0L
  for (s in starts) {
    nxt <- ctrl[ctrl > s]
    e <- if (length(nxt)) nxt[1] - 1L else length(lines)
    for (li in s:e) {
      line <- lines[li]
      cpos <- regexpr(";", line, fixed = TRUE)
      code <- if (cpos > 0) substr(line, 1L, cpos - 1L) else line
      comment <- if (cpos > 0) trimws(sub("^;+\\s*", "", substr(line, cpos, nchar(line)))) else ""
      # mask (same length, positions kept) the record name, labels and options
      mask <- function(x, pat) {
        m <- gregexpr(pat, x, perl = TRUE, ignore.case = TRUE)[[1]]
        if (m[1] == -1L) return(x)
        for (j in seq_along(m))
          substr(x, m[j], m[j] + attr(m, "match.length")[j] - 1L) <-
            strrep(" ", attr(m, "match.length")[j])
        x
      }
      if (li == s) code <- mask(code, "^\\s*\\$THETA\\w*")
      code <- mask(code, "\\b(NUMBERPOINTS|NUMPOINTS)\\s*=\\s*[0-9]+")
      code <- mask(code, "\\b(UNINT|ABSOLUTE|NOABORT)\\b")
      code <- mask(code, "\\b[A-Za-z_][A-Za-z0-9_]*\\s*=(?!=)")
      g <- gregexpr(.nm_ctl_theta_tok, code, perl = TRUE)[[1]]
      if (g[1] == -1L) next
      len <- attr(g, "match.length")
      cap <- attr(g, "capture.start"); capl <- attr(g, "capture.length")
      line_rows <- list()
      for (t in seq_along(g)) {
        tok <- substr(code, g[t], g[t] + len[t] - 1L)
        body <- substr(code, cap[t, 1], cap[t, 1] + capl[t, 1] - 1L)
        fixed <- capl[t, 2] > 0
        nrep <- if (capl[t, 4] > 0) as.integer(substr(code, cap[t, 4], cap[t, 4] + capl[t, 4] - 1L)) else 1L
        lo <- up <- NA_real_
        if (startsWith(body, "(")) {
          inner <- sub("^\\((.*)\\)$", "\\1", body)
          # FIX may be written inside the parentheses: (0.5 FIX), (0,1,5 FIXED)
          if (grepl("\\bFIX(?:ED)?\\b", inner, ignore.case = TRUE, perl = TRUE)) {
            fixed <- TRUE
            inner <- gsub("\\bFIX(?:ED)?\\b", " ", inner, ignore.case = TRUE, perl = TRUE)
          }
          # values separated by commas and/or blanks: (0,1,5) or (0 1 5)
          p <- strsplit(trimws(inner), "\\s*,\\s*|\\s+", perl = TRUE)[[1]]
          if (length(p) == 1L) ini <- num(p[1])
          else { lo <- num(p[1]); ini <- num(p[2]); if (length(p) >= 3L) up <- num(p[3]) }
        } else ini <- num(body)
        for (r in seq_len(max(1L, nrep)))
          line_rows[[length(line_rows) + 1L]] <- data.frame(
            line = li, start = as.integer(g[t]), end = as.integer(g[t] + len[t] - 1L),
            token = trimws(tok), lower = lo, initial = ini, upper = up, fixed = fixed,
            comment = comment, raw = line, record_line = (li == s), stringsAsFactors = FALSE)
      }
      for (lr in line_rows) {
        k <- k + 1L
        lr$theta <- k
        lr$n_on_line <- length(line_rows)
        rows[[k]] <- lr
      }
    }
  }
  if (!k) return(data.frame(theta = integer(0), line = integer(0), start = integer(0),
                            end = integer(0), token = character(0), lower = numeric(0),
                            initial = numeric(0), upper = numeric(0), fixed = logical(0),
                            comment = character(0), raw = character(0),
                            record_line = logical(0), n_on_line = integer(0),
                            stringsAsFactors = FALSE))
  out <- do.call(rbind, rows)
  out[, c("theta", "line", "start", "end", "token", "lower", "initial", "upper", "fixed",
          "comment", "raw", "record_line", "n_on_line")]
}

# -----------------------------------------------------------------------------
# $OMEGA and $SIGMA
# -----------------------------------------------------------------------------

# Structure of the random effects of $OMEGA (ETAs) or $SIGMA (EPSs), over all
# records: list with `eta` (one row per ETA or EPS: index `eta`, label,
# initial variance, fixed, structure DIAGONAL or BLOCK, block, position in the
# block, raw text, and `omega_comment`, the full comment of the line of its
# diagonal value) and `elements` (one row per variance or covariance).
# Handles several diagonal values on one line, values on the record line,
# DIAGONAL(n), BLOCK(n), BLOCK(n) SAME(k) and BLOCK SAME (size of the
# previous block).
.nm_ctl_random_parse <- function(lines, rec = c("OMEGA", "SIGMA")) {
  rec <- match.arg(rec)
  empty <- list(eta = data.frame(), elements = data.frame())
  starts <- grep(paste0("^\\s*\\$", rec, "\\b"), lines, ignore.case = TRUE, perl = TRUE)
  if (!length(starts)) return(empty)
  ctrl <- grep("^\\s*\\$[A-Za-z]", lines)
  is_fix <- function(x) grepl("\\bFIX(?:ED)?\\b", x, ignore.case = TRUE, perl = TRUE)
  numbers <- function(x) {
    x <- gsub("\\b(BLOCK|DIAGONAL|SAME)\\s*\\(\\s*[0-9]+\\s*\\)", " ", x, ignore.case = TRUE, perl = TRUE)
    m <- regmatches(x, gregexpr("(?<![A-Za-z0-9_.])[-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[EeDd][-+]?[0-9]+)?",
                                x, perl = TRUE))[[1]]
    if (!length(m)) return(numeric())
    suppressWarnings(as.numeric(gsub("[Dd]", "E", m)))
  }
  eta_rows <- list(); element_rows <- list()
  eta_offset <- 0L; last_n <- NA_integer_

  for (s in starts) {
    nxt <- ctrl[ctrl > s]
    e <- if (length(nxt)) nxt[1] - 1L else length(lines)
    rl <- lines[s:e]
    code <- sub(";.*$", "", rl)
    code[1] <- sub(paste0("^\\s*\\$", rec, "\\w*"), "", code[1], ignore.case = TRUE)
    comm <- ifelse(grepl(";", rl, fixed = TRUE), trimws(sub("^[^;]*;", "", rl)), NA_character_)
    head_txt <- paste(code, collapse = " ")
    b <- regmatches(head_txt, regexec("\\b(BLOCK|DIAGONAL)\\s*\\(\\s*([0-9]+)\\s*\\)", head_txt,
                                      ignore.case = TRUE, perl = TRUE))[[1]]
    same <- grepl("\\bSAME\\b", head_txt, ignore.case = TRUE, perl = TRUE)
    sm <- regmatches(head_txt, regexec("\\bSAME\\s*\\(\\s*([0-9]+)\\s*\\)", head_txt,
                                       ignore.case = TRUE, perl = TRUE))[[1]]
    reps <- if (length(sm) == 2L) max(1L, as.integer(sm[2])) else 1L
    type <- if (length(b)) toupper(b[2]) else if (same && grepl("\\bBLOCK\\b", head_txt, ignore.case = TRUE)) "BLOCK" else "SCALAR"
    n <- if (length(b)) as.integer(b[3]) else if (type == "BLOCK") last_n else NA_integer_
    if (type == "BLOCK" && is.na(n)) next
    # values with their line (index in rl)
    vals <- numeric(); vline <- integer()
    for (q in seq_along(code)) {
      v <- numbers(code[q])
      vals <- c(vals, v); vline <- c(vline, rep(q, length(v)))
    }
    # FIX anywhere in a BLOCK/DIAGONAL record fixes the whole record (NONMEM)
    fixed_all <- any(is_fix(code))
    # raw text: record name removed, only lines carrying code; a header line
    # is kept when something follows BLOCK(n)/DIAGONAL(n) (FIX, SAME, values)
    rl_raw <- rl
    rl_raw[1] <- sub(paste0("^\\s*\\$", rec, "\\w*\\s*"), "", rl_raw[1], ignore.case = TRUE, perl = TRUE)
    has_code <- nzchar(trimws(code))
    has_code[1] <- nzchar(trimws(sub("^\\s*(?:BLOCK|DIAGONAL)\\s*\\(\\s*[0-9]+\\s*\\)\\s*", "",
                                     code[1], ignore.case = TRUE, perl = TRUE)))
    raw_join <- paste(rl_raw[has_code], collapse = " | ")
    # comments label ETAs only when they sit on a line carrying values
    # (comment-only lines and bare record headers are ignored)
    lab <- function(q) if (is.na(q) || is.na(comm[q])) NA_character_ else .nm_ip_label_clean(comm[q])

    if (type == "SCALAR") {
      # each value is one ETA (several values may share a line)
      for (a in seq_along(vals)) {
        q <- vline[a]
        first_on_line <- a == 1L || vline[a - 1L] != q
        eta_offset <- eta_offset + 1L
        fx <- is_fix(code[q])
        eta_rows[[length(eta_rows) + 1L]] <- data.frame(
          eta = eta_offset, omega_label = if (first_on_line) lab(q) else NA_character_,
          omega_initial = vals[a], omega_fixed = fx, omega_structure = "DIAGONAL",
          block_id = NA_integer_, block_position = 1L, omega_raw = rl_raw[q],
          omega_comment = if (first_on_line) comm[q] else NA_character_,
          stringsAsFactors = FALSE)
        element_rows[[length(element_rows) + 1L]] <- data.frame(
          eta_i = eta_offset, eta_j = eta_offset, element_type = "variance",
          initial = vals[a], fixed = fx, structure = "DIAGONAL", block_id = NA_integer_,
          raw = rl_raw[q], stringsAsFactors = FALSE)
      }
      last_n <- NA_integer_
      next
    }

    need <- if (type == "BLOCK") n * (n + 1L) / 2L else n
    # SAME: no values; the comment of the record line labels the first ETA
    if (same) { vals <- rep(NA_real_, need); vline <- rep(if (has_code[1]) 1L else NA_integer_, need) }
    if (length(vals) < need) {
      vline <- c(vline, rep(NA_integer_, need - length(vals)))
      vals <- c(vals, rep(NA_real_, need - length(vals)))
    }
    vals <- vals[seq_len(need)]; vline <- vline[seq_len(need)]
    for (rep_k in seq_len(reps)) {
      block_id <- eta_offset + 1L
      ids <- eta_offset + seq_len(n)
      if (type == "DIAGONAL") {
        for (a in seq_len(n)) {
          first_use <- !is.na(vline[a]) && !(vline[a] %in% vline[seq_len(a - 1L)])
          eta_rows[[length(eta_rows) + 1L]] <- data.frame(
            eta = ids[a], omega_label = if (first_use) lab(vline[a]) else NA_character_,
            omega_initial = vals[a],
            omega_fixed = fixed_all, omega_structure = "DIAGONAL", block_id = block_id,
            block_position = a, omega_raw = raw_join,
            omega_comment = if (first_use) comm[vline[a]] else NA_character_,
            stringsAsFactors = FALSE)
          element_rows[[length(element_rows) + 1L]] <- data.frame(
            eta_i = ids[a], eta_j = ids[a], element_type = "variance", initial = vals[a],
            fixed = fixed_all, structure = "DIAGONAL", block_id = block_id,
            raw = raw_join, stringsAsFactors = FALSE)
        }
      } else {
        pos <- 0L; dline <- rep(NA_integer_, n); dval <- rep(NA_real_, n)
        for (a in seq_len(n)) for (bb in seq_len(a)) {
          pos <- pos + 1L
          if (a == bb) { dval[a] <- vals[pos]; dline[a] <- vline[pos] }
          element_rows[[length(element_rows) + 1L]] <- data.frame(
            eta_i = ids[a], eta_j = ids[bb], element_type = if (a == bb) "variance" else "covariance",
            initial = vals[pos], fixed = fixed_all, structure = "BLOCK", block_id = block_id,
            raw = raw_join, stringsAsFactors = FALSE)
        }
        # label of ETA a: comment of the line holding its variance (the
        # diagonal element), used once only (several rows on one line)
        for (a in seq_len(n)) {
          first_use <- !is.na(dline[a]) && !(dline[a] %in% dline[seq_len(a - 1L)])
          eta_rows[[length(eta_rows) + 1L]] <- data.frame(
            eta = ids[a], omega_label = if (first_use) lab(dline[a]) else NA_character_,
            omega_initial = dval[a], omega_fixed = fixed_all, omega_structure = "BLOCK",
            block_id = block_id, block_position = a, omega_raw = raw_join,
            omega_comment = if (first_use) comm[dline[a]] else NA_character_,
            stringsAsFactors = FALSE)
        }
      }
      eta_offset <- eta_offset + n
    }
    last_n <- if (type == "BLOCK") n else NA_integer_
  }
  list(eta = if (length(eta_rows)) do.call(rbind, eta_rows) else data.frame(),
       elements = if (length(element_rows)) do.call(rbind, element_rows) else data.frame())
}
