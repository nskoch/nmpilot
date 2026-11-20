#  ---------------------------------------------------------------------------
#  nm_bootstrap.R  --  Nonparametric case (subject-level) bootstrap for NONMEM
#  Part of the nmtools package (PsN-/WFN-independent).
#
#  Exported :  nm_bootstrap()
#  S3       :  print.nm_bootstrap(), plot.nm_bootstrap()
#  Helpers  :  .nm_*  (dot-prefixed; the $DATA/$INPUT/strip/launch helpers are
#              ported from the *verified* jackknife script -- same logic, kept
#              compatible so they can migrate to R/utils-nm.R unchanged).
#
#  Method (Efron nonparametric bootstrap for popPK):
#    * resampling unit = the subject (ID column), drawn WITH replacement to the
#      same number of subjects N;
#    * EVERY drawn copy receives a NEW, unique ID (1..N) so NONMEM never merges
#      duplicated subjects -- this is the core correctness requirement;
#    * ONE fixed, optimised control stream ($COV/$TABLE stripped, warnings muted,
#      $DATA repointed to a reused replicate file) is run against each replicate;
#    * final THETA/OMEGA/SIGMA are read from each run's .ext (final-estimate row);
#    * confidence intervals are percentile CIs over the SUCCESSFUL runs only.
#
#  Storage strategy (evaluated): draw -> write ONE replicate csv (reused name)
#    -> run -> extract -> clean.  No 500 datasets kept on disk; instead a compact
#    resample_index.rds (seed + sampled IDs per replicate) makes the whole thing
#    reproducible, and bootstrap_estimates.csv is written row-by-row so a crash
#    loses nothing and `resume = TRUE` continues where it stopped.
#  ---------------------------------------------------------------------------

# Named constants ------------------------------------------------------------
.NM_EXT_FINAL_ITER <- -1000000000L   # .ext row holding the final parameter estimates

# ===========================================================================
#  Exported function
# ===========================================================================

#' Nonparametric (case) bootstrap of a NONMEM model
#'
#' @param ctl            Path to the base control stream (.ctl).
#' @param nmfe_path      Path to the nmfe*.bat launcher (e.g. "C:/nm74g64/run/nmfe74.bat").
#' @param n              Number of bootstrap replicates (default 500).
#' @param ci             Confidence level(s) for the percentile CI (default 0.95; e.g. c(0.95,0.90)).
#' @param seed           RNG seed for reproducibility (default: a random seed, stored in the result).
#' @param id_col         Label of the subject ID column, as it appears in $INPUT (default "ID").
#' @param strata         Optional $INPUT column to stratify the resampling on (constant within ID).
#' @param comment_char   Data-file comment/ignore character. NULL = auto-derive from the single-char
#'                       IGNORE in $DATA ("@", "#", a letter...); "" = none.
#' @param sep            Data-file field separator (default ",").
#' @param run_dir        Working directory (default "<ctl_dir>/<base>_bootstrap").
#' @param keep_datasets  Keep each replicate data file? (default FALSE = space-optimal).
#' @param keep_outputs   Keep each replicate's NONMEM outputs? (default FALSE; .ext is read then removed).
#' @param strip_cov      Remove $COVARIANCE (default TRUE; useless for a bootstrap).
#' @param strip_table    Remove $TABLE records (default TRUE; fixed FILE= names would collide).
#' @param add_warnings   Insert a warning-suppression record (default TRUE).
#' @param warnings_line  Exact record inserted (default "$WARN WARN=NONE"; canonical form "$WARNINGS").
#' @param include_rounding Count "minimisation successful, but rounding errors" runs as successful (default FALSE).
#' @param resume         Skip replicates already in bootstrap_estimates.csv (default TRUE).
#' @param run            FALSE = only build the optimised ctl + resample plan (dry run, no NONMEM).
#' @param verbose        Progress messages (default TRUE).
#' @param n_cores Number of replicates run at the same time (default 1). A
#'   whole number,
#'   `"auto"` (physical cores minus one, see [nm_cores()]) or `"ask"`
#'   (interactive sessions: the cores are shown and the number is asked once per R
#'   session, then reused; `options(nmtools.n_cores = NULL)` asks again).
#'   Each run is then started in the background in its own directory; if R
#'   is interrupted, the runs already started go on in the background and a
#'   later call with the same arguments picks up their results. In parallel, each run uses a
#'   working subdirectory `slot<k>/` of `run_dir`; the estimates are appended
#'   to `bootstrap_estimates.csv` as the replicates end (the file is sorted
#'   by replicate when the result is assembled).
#'
#' @return An object of class "nm_bootstrap" (see print()/plot()).
#' @export
nm_bootstrap <- function(ctl,
                         nmfe_path        = NULL,
                         n                = 500L,
                         ci               = 0.95,
                         seed             = NULL,
                         id_col           = "ID",
                         strata           = NULL,
                         comment_char     = NULL,
                         sep              = ",",
                         run_dir          = NULL,
                         keep_datasets    = FALSE,
                         keep_outputs     = FALSE,
                         strip_cov        = TRUE,
                         strip_table      = TRUE,
                         add_warnings     = TRUE,
                         warnings_line    = "$WARN WARN=NONE",
                         include_rounding = FALSE,
                         resume           = TRUE,
                         run              = TRUE,
                         verbose          = TRUE,
                         n_cores          = 1L) {

  stopifnot(file.exists(ctl))
  n  <- as.integer(n)
  ci <- sort(unique(ci))
  if (any(ci <= 0 | ci >= 1)) stop("`ci` must be in (0, 1), e.g. 0.95.")
  if (is.null(seed)) seed <- as.integer(Sys.time()) %% .Machine$integer.max
  set.seed(seed)

  ctl     <- normalizePath(ctl, winslash = "/", mustWork = TRUE)
  ctl_dir <- dirname(ctl)
  base    <- tools::file_path_sans_ext(basename(ctl))
  if (is.null(run_dir)) run_dir <- file.path(ctl_dir, paste0(base, "_bootstrap"))
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

  ctl_lines <- .nm_read_lines(ctl)

  # --- Original data file, ID/strata positions, comment character ----------
  data_path <- .nm_extract_data_path(ctl)                 # absolute, verified helper
  input_nm  <- .nm_extract_input_names(ctl)               # $INPUT labels (positional)
  id_idx    <- .nm_col_index(input_nm, id_col)
  if (is.na(id_idx)) stop(sprintf("ID column '%s' not found in $INPUT: %s",
                                  id_col, paste(input_nm, collapse = " ")))
  strat_idx <- NA_integer_
  if (!is.null(strata)) {
    strat_idx <- .nm_col_index(input_nm, strata)
    if (is.na(strat_idx)) stop(sprintf("strata column '%s' not found in $INPUT.", strata))
  }
  if (is.null(comment_char)) comment_char <- .nm_data_comment_char(ctl_lines)  # "@" / "#" / letter / ""

  # --- Read the data (positional, formatting preserved) --------------------
  raw <- .nm_read_data(data_path, comment_char = comment_char, sep = sep,
                       id_idx = id_idx, strat_idx = strat_idx)
  ids <- unique(raw$id)
  N   <- length(ids)
  if (verbose) message(sprintf("Dataset: %d subjects, %d records (ID = col #%d '%s', comment = '%s').",
                               N, length(raw$lines), id_idx, id_col, comment_char))

  # --- Build the ONE optimised bootstrap control stream --------------------
  boot_csv_name <- paste0(base, "_boot.csv")              # relative; sits in run_dir (NONMEM CWD)
  boot_ctl_name <- paste0(base, "_boot.ctl")
  bl <- ctl_lines
  if (strip_cov)   bl <- .nm_strip_record(bl, "COV")
  if (strip_table) bl <- .nm_strip_record(bl, "TABLE")
  di <- grep("^\\s*\\$DATA", bl, ignore.case = TRUE)[1]
  bl[di] <- .nm_patch_data_line(bl[di], boot_csv_name)    # keep IGNORE/ACCEPT, point to reused csv
  if (add_warnings) bl <- .nm_add_warnings(bl, warnings_line)
  boot_ctl <- file.path(run_dir, boot_ctl_name)
  writeLines(bl, boot_ctl)
  if (verbose) message("Bootstrap control stream: ", boot_ctl)

  # --- Reproducible resample plan ------------------------------------------
  plan <- .nm_resample_plan(raw, ids, n = n, strat_idx = strat_idx)
  saveRDS(list(seed = seed, id_col = id_col, strata = strata, plan = plan),
          file.path(run_dir, "resample_index.rds"))

  est_file <- file.path(run_dir, "bootstrap_estimates.csv")
  boot_csv <- file.path(run_dir, boot_csv_name)

  if (!run) {
    if (verbose) message("run = FALSE: ctl + resample plan generated, no NONMEM launched.")
    return(structure(list(call = match.call(), base_ctl = ctl, boot_ctl = boot_ctl,
                          run_dir = run_dir, n_requested = n, n_run = 0L,
                          n_successful = 0L, success_rate = NA_real_,
                          estimates = NULL, ci = NULL, seed = seed, plan = plan),
                     class = "nm_bootstrap"))
  }
  if (is.null(nmfe_path) || !file.exists(nmfe_path))
    stop("A valid `nmfe_path` is required when run = TRUE.")
  nmfe_win <- normalizePath(nmfe_path, winslash = "\\", mustWork = FALSE)

  # --- Resume --------------------------------------------------------------
  done <- integer(0)
  if (resume && file.exists(est_file)) {
    prev <- utils::read.csv(est_file, check.names = FALSE)
    done <- prev$replicate
    if (verbose) message(sprintf("Resume: %d replicate(s) already done.", length(done)))
  }

  # --- Parallel loop: one working subdirectory per slot ---------------------
  n_cores <- .nm_resolve_cores(n_cores)
  if (n_cores > 1L) {
    todo <- setdiff(seq_len(n), done)
    slot_dir <- file.path(run_dir, paste0("slot", seq_len(n_cores)))
    for (d in slot_dir) { dir.create(d, showWarnings = FALSE); file.copy(boot_ctl, d, overwrite = TRUE) }
    if (verbose) message(sprintf("Running %d replicate(s), %d at a time.", length(todo), n_cores))
    .nm_run_queue(as.list(todo), nmfe_path, n_cores, verbose = FALSE, label = "replicates",
      prepare = function(b, slot) {
        .nm_write_replicate(raw, plan[[b]], file.path(slot_dir[slot], boot_csv_name), sep = sep)
        list(dir = slot_dir[slot], name = sub("\\.ctl$", "", boot_ctl_name))
      },
      finish = function(b, slot, st, secs) {
        sd  <- slot_dir[slot]
        lst <- file.path(sd, sub("\\.ctl$", ".lst", boot_ctl_name))
        status <- .nm_min_status(if (file.exists(lst)) lst else NA_character_)
        ext  <- file.path(sd, paste0(base, "_boot.ext"))
        pars <- if (file.exists(ext)) .nm_ext_final(ext) else NULL
        .nm_append_row(est_file,
                       c(list(replicate = b, status = status, seconds = round(secs, 1)),
                         as.list(pars)))
        if (verbose) message(sprintf("  [%d/%d] %s (%.0fs)", b, n, status, secs))
        if (keep_datasets)
          file.copy(file.path(sd, boot_csv_name),
                    file.path(run_dir, sprintf("%s_boot_%04d.csv", base, b)), overwrite = TRUE)
        unlink(file.path(sd, boot_csv_name))
        if (keep_outputs) {
          for (f in list.files(sd, pattern = paste0("^", base, "_boot\\.(lst|ext)$")))
            file.copy(file.path(sd, f), file.path(run_dir, sub("_boot\\.", sprintf("_boot_%04d.", b), f)),
                      overwrite = TRUE)
        }
        .nm_clean_run_outputs(sd, base, keep_ext = FALSE)
        unlink(file.path(sd, c(paste0(base, "_boot_nmtools.bat"), paste0(base, "_boot_nmtools.sh"),
                               paste0(base, "_boot_nmtools.log"), paste0(base, "_boot_nmtools.err"),
                               paste0(base, "_boot.nmtools_done"))))
      })
    done <- seq_len(n)   # nothing left for the sequential loop
  }

  # --- Main sequential loop ------------------------------------------------
  for (b in seq_len(n)) {
    if (b %in% done) next

    .nm_write_replicate(raw, plan[[b]], boot_csv, sep = sep)   # renumbered IDs 1..N, no header

    t0     <- Sys.time()
    lst    <- .nm_boot_launch(run_dir, boot_ctl_name, nmfe_win)
    secs   <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    status <- .nm_min_status(lst)
    ext    <- file.path(run_dir, paste0(base, "_boot.ext"))
    pars   <- if (file.exists(ext)) .nm_ext_final(ext) else NULL

    .nm_append_row(est_file,
                   c(list(replicate = b, status = status, seconds = round(secs, 1)),
                     as.list(pars)))
    if (verbose) message(sprintf("  [%d/%d] %s (%.0fs)", b, n, status, secs))

    if (!keep_datasets) unlink(boot_csv)
    if (!keep_outputs)  .nm_clean_run_outputs(run_dir, base, keep_ext = FALSE)
  }

  # --- Assemble ------------------------------------------------------------
  est <- utils::read.csv(est_file, check.names = FALSE)
  est <- est[order(est$replicate), , drop = FALSE]
  ok_status <- if (include_rounding) c("successful", "rounding_errors") else "successful"
  ok  <- est$status %in% ok_status
  n_ok <- sum(ok)

  param_names <- setdiff(names(est), c("replicate", "status", "seconds"))
  ci_tab <- .nm_percentile_ci(est[ok, , drop = FALSE], param_names, ci)

  base_ext <- file.path(dirname(ctl), paste0(base, ".ext"))
  original <- if (file.exists(base_ext)) .nm_ext_final(base_ext) else NULL
  if (!is.null(original)) {
    ci_tab$original <- unname(original[match(ci_tab$parameter, names(original))])
    ci_tab <- ci_tab[, c("parameter", "original",
                         setdiff(names(ci_tab), c("parameter", "original")))]
  }

  if (verbose)
    message(sprintf("Bootstrap done: %d/%d successful (%.1f%%).",
                    n_ok, nrow(est), 100 * n_ok / nrow(est)))

  structure(list(call = match.call(), base_ctl = ctl, boot_ctl = boot_ctl, run_dir = run_dir,
                 n_requested = n, n_run = nrow(est), n_successful = n_ok,
                 success_rate = n_ok / nrow(est), estimates = est, ci = ci_tab,
                 original = original, param_names = param_names, ci_levels = ci,
                 seed = seed, include_rounding = include_rounding),
            class = "nm_bootstrap")
}

# ===========================================================================
#  Control-stream helpers  (ported from the verified jackknife script)
#  .nm_read_lines, .nm_strip_record, .nm_extract_data_path and
#  .nm_patch_data_line are shared and live in utils-nm.R.
# ===========================================================================


# $INPUT column labels, positional (continuation lines, LABEL=ALIAS, DROP/SKIP). [verified]
.nm_extract_input_names <- function(ctl) {
  lines <- .nm_read_lines(ctl)
  i <- grep("^\\s*\\$INPUT", lines, ignore.case = TRUE)[1]
  if (is.na(i)) return(character(0))
  block <- lines[i]; k <- i + 1L
  while (k <= length(lines) && !grepl("^\\s*\\$", lines[k])) { block <- c(block, lines[k]); k <- k + 1L }
  txt  <- sub("^\\s*\\$INPUT\\s*", "", paste(block, collapse = " "), ignore.case = TRUE)
  txt  <- strsplit(txt, ";")[[1]][1]
  toks <- strsplit(trimws(txt), "\\s+")[[1]]; toks <- toks[nchar(toks) > 0]
  if (!length(toks)) return(character(0))
  nm <- vapply(toks, function(t) {
    parts <- strsplit(t, "=")[[1]]
    lab   <- parts[!toupper(parts) %in% c("DROP", "SKIP")]
    if (length(lab)) lab[1] else "DROP"
  }, character(1))
  make.unique(nm, sep = "_")
}


# Insert the warning-suppression record just after $PROBLEM (idempotent).
.nm_add_warnings <- function(lines, warnings_line) {
  if (length(grep("^\\s*\\$WARN", lines, ignore.case = TRUE))) return(lines)
  p <- grep("^\\s*\\$PROB", lines, ignore.case = TRUE)
  append(lines, warnings_line, after = if (length(p)) p[1] else 0L)
}

# The single-character IGNORE from $DATA -> the data comment char ("@", "#", ...).
# Conditional IGNORE=(...) are skipped; "" if no single-char ignore is present.
.nm_data_comment_char <- function(lines) {
  i <- grep("^\\s*\\$DATA", lines, ignore.case = TRUE)[1]
  if (is.na(i)) return("")
  dollar <- grep("^\\s*\\$", lines)
  e <- { nx <- dollar[dollar > i]; if (length(nx)) min(nx) - 1L else length(lines) }
  block <- lines[i:e]
  block <- vapply(block, function(l) sub(";.*$", "", l), character(1))   # per-line comment strip
  txt   <- paste(block, collapse = " ")
  m <- regmatches(txt, gregexpr("IGNORE\\s*=\\s*([^( \t])", txt, ignore.case = TRUE))[[1]]
  if (!length(m)) return("")
  ch <- sub(".*IGNORE\\s*=\\s*", "", m[1], ignore.case = TRUE)
  substr(trimws(ch), 1, 1)
}

.nm_col_index <- function(cols, label) {
  idx <- which(toupper(cols) == toupper(label))
  if (length(idx)) idx[1] else NA_integer_
}

# ===========================================================================
#  Data resampling
# ===========================================================================

# Skip predicate matching NONMEM's single-char IGNORE semantics.
.nm_is_comment_line <- function(lines, cc) {
  if (is.null(cc) || !nzchar(cc)) return(rep(FALSE, length(lines)))
  if (cc == "@") !grepl("^\\s*[-+.0-9]", lines)         # first non-blank char not numeric -> ignored
  else substr(lines, 1, 1) == cc                        # column-1 char equals the ignore char
}

# Read data as raw lines; keep conditionally-IGNORE'd rows (NONMEM re-applies them);
# return per-line raw text + ID + stratum, parsed positionally.
.nm_read_data <- function(path, comment_char, sep = ",", id_idx, strat_idx = NA_integer_) {
  all <- .nm_read_lines(path)
  all <- all[nzchar(trimws(all))]
  all <- all[!.nm_is_comment_line(all, comment_char)]
  sp  <- if (sep == ",") "," else if (sep %in% c(" ", "\t")) "[[:space:]]+" else sep
  fld <- strsplit(all, sp)
  id  <- vapply(fld, function(f) trimws(f[id_idx]), character(1))
  strat <- if (!is.na(strat_idx)) vapply(fld, function(f) trimws(f[strat_idx]), character(1)) else NULL
  if (length(id) && suppressWarnings(is.na(as.numeric(id[1]))))
    warning("First data ID '", id[1], "' is non-numeric -- a header row may not be ",
            "matched by comment_char='", comment_char, "'. Check `comment_char`.")
  list(lines = all, id = id, strat = strat, sep = sep, id_idx = id_idx)
}

# For each replicate: a vector of sampled ORIGINAL ids (length N), stratified if asked.
.nm_resample_plan <- function(raw, ids, n, strat_idx = NA_integer_) {
  if (is.na(strat_idx)) {
    N <- length(ids)
    lapply(seq_len(n), function(i) sample(ids, N, replace = TRUE))
  } else {
    strat_of_id <- raw$strat[match(ids, raw$id)]
    groups <- split(ids, strat_of_id)
    lapply(seq_len(n), function(i)
      unlist(lapply(groups, function(g) sample(g, length(g), replace = TRUE)), use.names = FALSE))
  }
}

# Write one replicate: each sampled subject's block, with a fresh sequential ID
# (1..N) so NONMEM never merges duplicated subjects. No header (matches $DATA options).
.nm_write_replicate <- function(raw, sampled_ids, out_path, sep = ",") {
  by_id <- split(seq_along(raw$id), raw$id)
  join  <- if (sep == ",") "," else " "
  sp    <- if (sep == ",") "," else if (sep %in% c(" ", "\t")) "[[:space:]]+" else sep
  con <- file(out_path, open = "wb"); on.exit(close(con))
  new_id <- 0L; buf <- character(0)
  for (oid in sampled_ids) {
    new_id <- new_id + 1L
    block  <- raw$lines[by_id[[as.character(oid)]]]
    block  <- vapply(block, function(ln) {
      f <- strsplit(ln, sp)[[1]]
      if (length(f) >= raw$id_idx) f[raw$id_idx] <- as.character(new_id)
      paste(f, collapse = join)
    }, character(1), USE.NAMES = FALSE)
    buf <- c(buf, block)
    if (length(buf) > 5000L) { writeLines(buf, con, sep = "\n"); buf <- character(0) }
  }
  if (length(buf)) writeLines(buf, con, sep = "\n")
  invisible(out_path)
}

# ===========================================================================
#  NONMEM launch + output parsing
# ===========================================================================

# Verified launcher: setwd into the run dir, call nmfe via cmd.exe /c quoting
# ONLY the executable (nmfe_win uses backslashes; .ctl/.lst have no spaces).
.nm_boot_launch <- function(run_dir, boot_ctl_name, nmfe_win) {
  lst_name <- sub("\\.ctl$", ".lst", boot_ctl_name)
  old <- setwd(run_dir); on.exit(setwd(old))
  suppressWarnings(system2("cmd.exe",
    args = c("/c", shQuote(nmfe_win, type = "cmd"), boot_ctl_name, lst_name)))
  lst <- file.path(run_dir, lst_name)
  if (file.exists(lst)) lst else NA_character_
}

# Final THETA/OMEGA/SIGMA from a .ext (last $EST table, ITERATION == -1e9).
.nm_ext_final <- function(ext_path) {
  lines <- .nm_read_lines(ext_path)
  tb <- grep("^\\s*TABLE NO", lines)
  if (!length(tb)) return(NULL)
  block <- lines[(tb[length(tb)] + 1):length(lines)]
  hdr  <- strsplit(trimws(block[1]), "\\s+")[[1]]
  body <- block[-1]; body <- body[nzchar(trimws(body))]
  mat  <- do.call(rbind, lapply(body, function(l) as.numeric(strsplit(trimws(l), "\\s+")[[1]])))
  colnames(mat) <- hdr
  fin <- mat[mat[, "ITERATION"] == .NM_EXT_FINAL_ITER, , drop = FALSE]
  if (!nrow(fin)) return(NULL)
  v <- fin[1, ]; v <- v[setdiff(names(v), c("ITERATION", "OBJ"))]
  v[v != 0 | !grepl("^(OMEGA|SIGMA)", names(v))]      # drop structural-zero off-diagonals
}

# Minimisation status: "successful" | "rounding_errors" | "terminated" | "nmtran_error".
.nm_min_status <- function(lst_path) {
  if (is.null(lst_path) || is.na(lst_path) || !file.exists(lst_path)) return("nmtran_error")
  txt <- toupper(paste(.nm_read_lines(lst_path), collapse = "\n"))
  if (grepl("AN ERROR WAS FOUND IN THE CONTROL", txt) ||
      grepl("AN ERROR WAS FOUND ON LINE", txt))    return("nmtran_error")
  # status of the last estimation step (several $EST, EM/IMP methods)
  code <- .nm_lst_step_status(.nm_read_lines(lst_path))$code
  if (identical(code, "successful") || identical(code, "rounding_errors")) return(code)
  "terminated"
}

# Remove bulky NONMEM scratch/outputs for one run (optionally keep the .ext).
.nm_clean_run_outputs <- function(run_dir, base, keep_ext = FALSE) {
  keep <- paste0(base, "_boot.ctl")
  if (keep_ext) keep <- c(keep, paste0(base, "_boot.ext"))
  patt <- c(paste0("^", base, "_boot\\.(lst|ext|phi|cov|cor|coi|shk|grd|xml|smr|msf|clt|cnv|fgh|shm)$"),
            "^(FDATA|FCON|FREPORT|FSTREAM|FSUBS.*|PRDERR|OUTPUT|INTER|FMSG|FSIZES|FLIB|GFCOMPILE.*|LINKC.*|nonmem.*|nul|WK.*|trash\\.tmp|fort\\..*|.*\\.(f90|o|lnk|set))$")
  for (f in list.files(run_dir, all.files = TRUE, no.. = TRUE)) {
    if (f %in% keep) next
    if (any(vapply(patt, function(p) grepl(p, f, ignore.case = TRUE), logical(1))))
      unlink(file.path(run_dir, f), recursive = TRUE, force = TRUE)
  }
  invisible(NULL)
}

# ===========================================================================
#  CI + incremental I/O
# ===========================================================================

.nm_append_row <- function(file, row) {
  df <- as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  if (!file.exists(file)) { utils::write.csv(df, file, row.names = FALSE); return(invisible(file)) }
  old <- utils::read.csv(file, check.names = FALSE)
  cols <- union(names(old), names(df))
  for (c in setdiff(cols, names(old))) old[[c]] <- NA
  for (c in setdiff(cols, names(df)))  df[[c]]  <- NA
  utils::write.csv(rbind(old[cols], df[cols]), file, row.names = FALSE)
  invisible(file)
}

.nm_percentile_ci <- function(est_ok, param_names, ci_levels) {
  do.call(rbind, lapply(param_names, function(p) {
    x <- suppressWarnings(as.numeric(est_ok[[p]])); x <- x[is.finite(x)]
    base <- data.frame(parameter = p, n = length(x),
                       median = if (length(x)) stats::median(x) else NA_real_,
                       stringsAsFactors = FALSE)
    for (lv in ci_levels) {
      a <- (1 - lv) / 2
      q <- if (length(x)) stats::quantile(x, c(a, 1 - a), names = FALSE, type = 7) else c(NA, NA)
      base[[sprintf("lower%g", 100 * lv)]] <- q[1]
      base[[sprintf("upper%g", 100 * lv)]] <- q[2]
    }
    base
  }))
}

# ===========================================================================
#  S3 methods
# ===========================================================================

#' @export
print.nm_bootstrap <- function(x, digits = 4, ...) {
  cat("<nm_bootstrap>\n")
  cat(sprintf("  base ctl   : %s\n", basename(x$base_ctl)))
  cat(sprintf("  run dir    : %s\n", x$run_dir))
  if (x$n_run == 0L) { cat("  (dry run: ctl + resample plan only)\n"); return(invisible(x)) }
  cat(sprintf("  replicates : %d requested, %d run\n", x$n_requested, x$n_run))
  cat(sprintf("  successful : %d (%.1f%%)%s\n", x$n_successful, 100 * x$success_rate,
              if (x$include_rounding) "  [rounding-error runs included]" else ""))
  cat(sprintf("  CI level(s): %s\n", paste0(100 * x$ci_levels, "%", collapse = ", ")))
  cat(sprintf("  seed       : %d\n\n", x$seed))
  if (!is.null(x$ci)) {
    tab <- x$ci; num <- vapply(tab, is.numeric, logical(1))
    tab[num] <- lapply(tab[num], round, digits)
    print(tab, row.names = FALSE)
  }
  invisible(x)
}

#' @export
plot.nm_bootstrap <- function(x, type = c("ci", "dist"), level = NULL, ...) {
  type <- match.arg(type)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot.nm_bootstrap() (Suggests).")
  if (is.null(x$ci) || !nrow(x$ci)) stop("No CI table to plot (dry run?).")
  lv <- if (is.null(level)) max(x$ci_levels) else level
  lo <- sprintf("lower%g", 100 * lv); hi <- sprintf("upper%g", 100 * lv)
  if (!all(c(lo, hi) %in% names(x$ci)))
    stop(sprintf("CI level %g%% not available; computed: %s",
                 100 * lv, paste0(100 * x$ci_levels, "%", collapse = ", ")))

  if (type == "ci") {
    d <- x$ci; d$parameter <- factor(d$parameter, levels = rev(d$parameter))
    p <- ggplot2::ggplot(d, ggplot2::aes(y = .data$parameter)) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[lo]], xmax = .data[[hi]]), height = 0.25) +
      ggplot2::geom_point(ggplot2::aes(x = .data$median)) +
      ggplot2::labs(x = sprintf("Bootstrap median + %g%% CI", 100 * lv), y = NULL,
                    title = sprintf("Bootstrap %g%% CI  (%d/%d successful)",
                                    100 * lv, x$n_successful, x$n_run)) +
      ggplot2::theme_bw()
    if (!is.null(x$original))
      p <- p + ggplot2::geom_point(ggplot2::aes(x = .data$original),
                                   shape = 4, size = 2.5, colour = "firebrick")
    return(p)
  }

  ok <- x$estimates$status %in%
    (if (x$include_rounding) c("successful", "rounding_errors") else "successful")
  e  <- x$estimates[ok, , drop = FALSE]
  long <- do.call(rbind, lapply(x$param_names, function(p)
    data.frame(parameter = p, value = suppressWarnings(as.numeric(e[[p]])), stringsAsFactors = FALSE)))
  long <- long[is.finite(long$value), ]
  bounds <- x$ci[, c("parameter", lo, hi)]; names(bounds) <- c("parameter", "lo", "hi")
  ggplot2::ggplot(long, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey80", colour = "grey40") +
    ggplot2::geom_vline(data = bounds, ggplot2::aes(xintercept = .data$lo), linetype = 2) +
    ggplot2::geom_vline(data = bounds, ggplot2::aes(xintercept = .data$hi), linetype = 2) +
    ggplot2::facet_wrap(~ .data$parameter, scales = "free") +
    ggplot2::labs(x = "Bootstrap estimate", y = "count",
                  title = sprintf("Bootstrap distributions (%g%% CI dashed)", 100 * lv)) +
    ggplot2::theme_bw()
}
