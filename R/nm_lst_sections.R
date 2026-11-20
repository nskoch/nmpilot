# =============================================================================
# Reading NONMEM outputs with several $PROB, several $EST, $MIX, IOV (SAME),
# empty or inconsistent .ext files, and failed or interrupted runs.
# Internal helpers used by nm_extract_runs().
# =============================================================================

# Line numbers of the "PROBLEM NO.:" headers of a .lst.
.nm_lst_problem_starts <- function(L) {
  grep("^\\s*PROBLEM NO\\.:\\s*[0-9]+", L)
}

# Structure of a .lst: problems, estimation steps ($EST) and the part of the
# file that describes the last estimation step of the last problem that has
# one. A problem without estimation (simulation only) is skipped. Files
# without #METH lines (NONMEM < 7) are read as one step.
.nm_lst_structure <- function(L) {
  ps <- .nm_lst_problem_starts(L)
  if (!length(ps)) ps <- 1L
  pe <- c(ps[-1] - 1L, length(L))
  # text before the first problem (NM-TRAN messages) goes with problem 1
  ps[1] <- 1L

  meth_all <- grep("^\\s*#METH:", L)
  has_est <- vapply(seq_along(ps), function(k)
    any(meth_all >= ps[k] & meth_all <= pe[k]), logical(1))
  k <- if (any(has_est)) max(which(has_est)) else length(ps)

  Lp <- L[ps[k]:pe[k]]
  meth <- grep("^\\s*#METH:", Lp)
  methods <- trimws(sub("^\\s*#METH:\\s*", "", Lp[meth]))
  Ls <- if (length(meth)) Lp[meth[length(meth)]:length(Lp)] else Lp

  list(
    n_problems = length(ps),
    problem = k,
    problems_with_est = sum(has_est),
    n_est = length(meth),
    methods = methods,
    last_method = if (length(methods)) methods[length(methods)] else NA_character_,
    L_problem = Lp,
    L_step = Ls
  )
}

# Short code of an estimation method, from its #METH text.
.nm_est_code <- function(m) {
  if (is.na(m)) return(NA_character_)
  u <- toupper(m)
  if (grepl("ITERATIVE TWO STAGE", u)) return("ITS")
  if (grepl("STOCHASTIC APPROXIMATION", u)) return("SAEM")
  if (grepl("IMPORTANCE SAMPLING", u)) return(if (grepl("EVALUATION|NO ESTIMATION|ONLY", u)) "IMP(eval)" else "IMP")
  if (grepl("MCMC|BAYES", u)) return("BAYES")
  if (grepl("NUTS", u)) return("NUTS")
  if (grepl("DIRECT SAMPLING", u)) return("DIRECT")
  if (grepl("LAPLACIAN", u)) return("LAPLACE")
  if (grepl("FIRST ORDER CONDITIONAL", u)) return(if (grepl("INTERACTION", u)) "FOCEI" else "FOCE")
  if (grepl("FIRST ORDER", u)) return("FO")
  if (grepl("EVALUATION", u)) return("EVAL")
  m
}

# State of the run as a whole (not the convergence of the estimation).
#   "completed"    : the last estimation step reports an objective function;
#   "nmtran_error" : NM-TRAN stopped on the control stream or the data;
#   "terminated"   : NONMEM stopped (PROGRAM TERMINATED) before the end of
#                    the last estimation step;
#   "interrupted"  : the output stops before the end of the last estimation
#                    step (run killed, disk full, still running...).
.nm_lst_run_state <- function(L, st, ofv) {
  if (!length(L)) return("interrupted")
  # NM-TRAN messages come before the first problem header
  ps <- .nm_lst_problem_starts(L)
  head_txt <- if (length(ps)) L[seq_len(ps[1])] else L
  nmtran_err <- any(grepl("AN ERROR WAS FOUND IN THE CONTROL STATEMENTS|\\(DATA ERROR\\)|^\\s*AN ERROR WAS FOUND|NM-TRAN.*ERROR",
                          head_txt, ignore.case = TRUE))
  has_obj <- any(grepl("#OBJV|OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT|MINIMUM VALUE OF OBJECTIVE FUNCTION",
                       st$L_step, ignore.case = TRUE))
  if (!length(grep("^\\s*#METH:", L)) && !length(ps) && nmtran_err) return("nmtran_error")
  if (nmtran_err && !has_obj) return("nmtran_error")
  if (has_obj && is.finite(ofv)) return("completed")
  if (any(grepl("PROGRAM TERMINATED", st$L_step, ignore.case = TRUE))) return("terminated")
  "interrupted"
}


# -----------------------------------------------------------------------------
# Final estimates printed in the .lst (used when the .ext is empty, missing
# or does not match the .lst)
# -----------------------------------------------------------------------------

.nm_lst_num_tokens <- function(x) {
  unlist(regmatches(x, gregexpr("-?[0-9]+\\.[0-9]+E[+-][0-9]+|\\.{5,}", x)))
}

# Values of one printed block ("FINAL PARAMETER ESTIMATE" or "STANDARD ERROR
# OF ESTIMATE") as a named vector THETA1.., OMEGA(i,j), SIGMA(i,j); NA for
# elements printed as dots.
.nm_lst_block_values <- function(B) {
  out <- c()
  tok2num <- function(t) suppressWarnings(ifelse(grepl("^\\.+$", t), NA_real_, as.numeric(t)))

  i_th <- grep("THETA - VECTOR OF FIXED EFFECTS", B)
  i_om <- grep("OMEGA - COV MATRIX FOR RANDOM EFFECTS", B)
  i_si <- grep("SIGMA - COV MATRIX FOR RANDOM EFFECTS", B)
  ends <- sort(c(i_th, i_om, i_si, grep("CORR MATRIX", B), length(B) + 1L))
  seg <- function(i) {
    if (!length(i)) return(character(0))
    e <- ends[ends > i[1]][1]
    B[(i[1] + 1L):(e - 1L)]
  }

  th <- seg(i_th)
  th <- th[!grepl("TH\\s*[0-9]", th)]
  v <- tok2num(.nm_lst_num_tokens(th))
  if (length(v)) out <- c(out, stats::setNames(v, paste0("THETA", seq_along(v))))

  mat <- function(S, lab, pre) {
    rows <- list(); cur <- NA_integer_
    for (ln in S) {
      if (grepl(paste0("^\\s*", lab, "[0-9]+\\s*$"), ln)) {
        cur <- as.integer(sub(paste0("^\\s*", lab), "", trimws(ln)))
        rows[[as.character(cur)]] <- numeric(0)
      } else if (grepl("^\\+", ln) && !is.na(cur)) {
        rows[[as.character(cur)]] <- c(rows[[as.character(cur)]], tok2num(.nm_lst_num_tokens(ln)))
      }
    }
    res <- c()
    for (r in names(rows)) {
      vals <- rows[[r]]
      if (!length(vals)) next
      nm <- sprintf("%s(%s,%d)", pre, r, seq_along(vals))
      res <- c(res, stats::setNames(vals, nm))
    }
    res
  }
  out <- c(out, mat(seg(i_om), "ETA", "OMEGA"), mat(seg(i_si), "EPS", "SIGMA"))
  out
}

# Estimates of the last estimation step read from the .lst, returned in the
# same form as .nm_parse_ext() so that .nm_ext_params() can use them.
.nm_lst_estimates <- function(Ls, ofv = NA_real_) {
  hdr <- function(p) grep(p, Ls)
  i_est <- hdr("FINAL PARAMETER ESTIMATE")
  if (!length(i_est)) return(NULL)
  i_se  <- hdr("STANDARD ERROR OF ESTIMATE")
  i_next <- sort(c(hdr("COVARIANCE MATRIX OF ESTIMATE"), hdr("STANDARD ERROR OF ESTIMATE"),
                   hdr("CORRELATION MATRIX OF ESTIMATE"), hdr("#METH:"), length(Ls) + 1L))
  blk <- function(i) {
    if (!length(i)) return(character(0))
    e <- i_next[i_next > i[1]][1]
    Ls[i[1]:(e - 1L)]
  }
  est <- .nm_lst_block_values(blk(i_est))
  if (!length(est)) return(NULL)
  se <- if (length(i_se)) .nm_lst_block_values(blk(i_se)) else c()
  se_v <- stats::setNames(rep(NA_real_, length(est)), names(est))
  common <- intersect(names(se), names(est))
  se_v[common] <- se[common]

  offdiag <- grepl("^(OMEGA|SIGMA)\\(([0-9]+),([0-9]+)\\)$", names(est)) &
    sub("^.*\\(([0-9]+),([0-9]+)\\)$", "\\1", names(est)) !=
    sub("^.*\\(([0-9]+),([0-9]+)\\)$", "\\2", names(est))
  # structural zeros off the diagonal, zero variances, and (when standard
  # errors were printed) elements without a standard error are fixed
  fixed <- (offdiag & est == 0) | (!offdiag & grepl("^(OMEGA|SIGMA)", names(est)) & est == 0)
  if (length(i_se)) fixed <- fixed | is.na(se_v)

  rse <- ifelse(!fixed & is.finite(se_v) & est != 0, abs(se_v / est) * 100, NA_real_)
  list(
    raw = tibble::tibble(raw_name = names(est), value = as.numeric(est),
                         fixed = as.logical(fixed), se = as.numeric(se_v),
                         rse_pct = as.numeric(rse)),
    ofv = ofv,
    last_iteration = NA_real_,
    final_row_available = TRUE,
    table = NULL,
    source = "lst"
  )
}


# -----------------------------------------------------------------------------
# Control stream: IOV (SAME) and mixture ($MIX)
# -----------------------------------------------------------------------------

# Indices of the ETAs declared by a $OMEGA ... SAME record (inter-occasion
# variability copies), and the number of SAME records.
.nm_ctl_same_etas <- function(lines) {
  same <- integer(0); total <- 0L; last_dim <- NA_integer_; n_rec <- 0L
  for (txt in .nm_ctl_records_text(lines, "OMEGA")) {
    g <- regmatches(txt, regexec("\\b(BLOCK?|DIAG(?:ONAL)?)\\s*\\(\\s*([0-9]+)\\s*\\)", txt, perl = TRUE))[[1]]
    s <- regmatches(txt, regexec("\\bSAME\\b(\\s*\\(\\s*([0-9]+)\\s*\\))?", txt, perl = TRUE))[[1]]
    mult <- if (length(s) && nzchar(s[3])) as.integer(s[3]) else 1L
    if (length(s)) {
      dim <- if (length(g)) as.integer(g[3]) else last_dim
      if (!is.na(dim)) {
        same <- c(same, total + seq_len(dim * mult))
        total <- total + dim * mult
        n_rec <- n_rec + 1L
      }
    } else if (length(g)) {
      last_dim <- as.integer(g[3]); total <- total + last_dim
    } else {
      z <- gsub("\\b(FIX|FIXED|UNINT|VARIANCE|STANDARD|SD|COVARIANCE|CORRELATION|CORRELATON|CHOLESKY|VALUES)\\b",
                " ", txt, perl = TRUE)
      total <- total + length(regmatches(z, gregexpr(.nm_ctl_num_pattern, z, perl = TRUE))[[1]])
      last_dim <- NA_integer_
    }
  }
  list(etas = as.integer(same), n_records = n_rec)
}

# $MIX: number of subpopulations (NSPOP) and THETAs used in P(k).
.nm_ctl_mixture <- function(lines) {
  txt <- .nm_ctl_records_text(lines, "MIX")
  if (!length(txt)) return(list(mix = FALSE, n_subpop = NA_integer_, p_theta = integer(0)))
  t <- paste(txt, collapse = " ")
  ns <- regmatches(t, regexec("NSPOP\\s*=\\s*([0-9]+)", t))[[1]]
  n_subpop <- if (length(ns)) as.integer(ns[2]) else {
    p <- regmatches(t, gregexpr("P\\(\\s*([0-9]+)\\s*\\)", t))[[1]]
    if (length(p)) max(as.integer(gsub("\\D", "", p))) else NA_integer_
  }
  th <- regmatches(t, gregexpr("THETA\\(\\s*([0-9]+)\\s*\\)", t))[[1]]
  list(mix = TRUE, n_subpop = n_subpop,
       p_theta = sort(unique(as.integer(gsub("\\D", "", th)))))
}

# .phm (mixture runs): probability of each subpopulation for each subject;
# the most probable subpopulation is the subject's assignment.
.nm_parse_phm <- function(path) {
  L <- .nm_read_lines_latin1(path)
  if (is.null(L) || !length(L)) return(NULL)
  tb <- grep("^TABLE NO", L)
  if (length(tb)) L <- L[tb[length(tb)]:length(L)]
  h <- grep("^\\s*SUBJECT_NO", L)
  if (!length(h)) return(NULL)
  hdr <- strsplit(trimws(L[h[1]]), "\\s+")[[1]]
  dat <- L[(h[1] + 1L):length(L)]
  dat <- dat[nzchar(trimws(dat)) & !grepl("^TABLE", dat)]
  tok <- strsplit(trimws(dat), "\\s+")
  tok <- tok[lengths(tok) == length(hdr)]
  if (!length(tok)) return(NULL)
  m <- do.call(rbind, tok)
  colnames(m) <- hdr
  if (!all(c("ID", "SUBPOP", "PMIX") %in% hdr)) return(NULL)
  d <- data.frame(ID = suppressWarnings(as.numeric(m[, "ID"])),
                  SUBPOP = suppressWarnings(as.integer(m[, "SUBPOP"])),
                  PMIX = suppressWarnings(as.numeric(m[, "PMIX"])))
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d <- d[order(d$ID, -d$PMIX), , drop = FALSE]
  best <- d[!duplicated(d$ID), , drop = FALSE]
  d$assigned <- d$SUBPOP == best$SUBPOP[match(d$ID, best$ID)]
  tibble::as_tibble(d)
}


# -----------------------------------------------------------------------------
# Status of the last estimation step, shared by nm_run(), nm_profile(),
# nm_bootstrap() and nm_jackknife()
# -----------------------------------------------------------------------------

# code: "successful", "rounding_errors", "terminated", "not_tested" or NA
# (no status line: run interrupted or evaluation only).
.nm_lst_step_status <- function(L) {
  if (!length(L)) return(list(code = NA_character_, line = NA_character_))
  Ls <- .nm_lst_structure(L)$L_step
  pat <- paste0("MINIMIZATION (SUCCESSFUL|TERMINATED|ABORTED)",
                "|OPTIMIZATION (WAS )?(COMPLETED|NOT COMPLETED|NOT TESTED)",
                "|(STOCHASTIC|STATISTICAL|BURN-IN) PORTION (WAS )?(COMPLETED|NOT COMPLETED)",
                "|EXPECTATION ONLY PROCESS (WAS )?COMPLETED")
  idx <- grep(pat, Ls, ignore.case = TRUE)
  if (!length(idx)) return(list(code = NA_character_, line = NA_character_))
  i <- max(idx); ln <- toupper(Ls[i])
  near <- toupper(paste(Ls[i:min(length(Ls), i + 3L)], collapse = " "))
  code <- if (grepl("NOT TESTED", ln)) "not_tested"
    else if (grepl("SUCCESSFUL", ln) || (grepl("COMPLETED", ln) && !grepl("NOT COMPLETED", ln))) "successful"
    else if (grepl("ROUNDING ERRORS", near)) "rounding_errors"
    else "terminated"
  list(code = code, line = trimws(Ls[i]))
}

# Lines of the last estimation step (whole file when there is no #METH).
.nm_lst_step_lines <- function(L) {
  if (!length(L)) return(L)
  .nm_lst_structure(L)$L_step
}
