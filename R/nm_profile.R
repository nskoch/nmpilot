#-------------------------------------------------------------------------------
# nm_profile.R
# Objective-function profile (profile likelihood) over one model parameter.
# Part of the {nmpilot} package. PsN-independent, Windows-compatible.
#-------------------------------------------------------------------------------

# Signed integer/decimal with optional scientific exponent
.NM_NUM_PAT   <- "-?[0-9]*\\.?[0-9]+(?:[Ee][+-]?[0-9]+)?"
# One $THETA token: a parenthesised group OR a bare number, optional trailing FIX
.NM_THETA_TOK <- paste0("(\\([^)]*\\)|", .NM_NUM_PAT, ")(\\s*FIX(?:ED)?)?")

#-------------------------------------------------------------------------------
# Internal helpers (not exported)
#-------------------------------------------------------------------------------

#' @noRd
.nm_ensure_dir <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
  x
}

#' @noRd
.nm_first_num <- function(x) {
  m <- regexpr(.NM_NUM_PAT, x, perl = TRUE)
  if (m < 0) return(NA_real_)
  as.numeric(regmatches(x, m))
}

#' @noRd
.nm_decimals <- function(x) {
  s <- format(x, scientific = FALSE, trim = TRUE)
  if (!grepl("\\.", s)) return(0L)
  nchar(sub("^[^.]*\\.", "", s))
}

#' @noRd
.nm_fmt_val <- function(x, digits) formatC(x, format = "f", digits = digits)

#' @noRd
.nm_tag_val <- function(x, digits) {
  s <- .nm_fmt_val(x, digits)
  s <- gsub("-", "m", s, fixed = TRUE)   # minus -> m
  gsub(".", "p", s, fixed = TRUE)        # point -> p
}


#' OFV from #OBJV (WITHOUT the additive constant N*LOG(2PI))
#' @noRd
.nm_extract_ofv <- function(lst) {
  lines <- .nm_lst_step_lines(readLines(lst, warn = FALSE))
  idx <- grep("#OBJV", lines)
  if (length(idx) == 0) return(NA_real_)
  .nm_first_num(lines[idx[length(idx)]])
}

#' Minimisation / covariance / NMTRAN outcomes from a .lst
#'
#' Covariance handling copes with modern NONMEM ("COVARIANCE STEP ABORTED") and
#' older formats ("COVARIANCE STEP OMITTED: NO/YES"); the condition-number
#' warning is matched on "EIGENVALUE OF COV" (not "WRESCHOL EIGENVALUE").
#' @noRd
.nm_parse_status <- function(lst) {
  txt <- paste(readLines(lst, warn = FALSE), collapse = "\n")

  nmtran_err <- grepl("AN ERROR WAS FOUND IN THE CONTROL STATEMENTS", txt, ignore.case = TRUE)
  # status of the last estimation step (several $EST, EM methods)
  min_ok     <- identical(.nm_lst_step_status(readLines(lst, warn = FALSE))$code, "successful")

  cov_omitted   <- grepl("COVARIANCE STEP OMITTED:\\s*YES", txt, ignore.case = TRUE)
  cov_requested <- grepl("COVARIANCE STEP OMITTED:\\s*NO", txt, ignore.case = TRUE) |
                   grepl("STANDARD ERROR OF ESTIMATE", txt, ignore.case = TRUE)
  cov_aborted   <- grepl("COVARIANCE STEP ABORTED", txt, ignore.case = TRUE)
  cov_matrix_ok <- grepl("COVARIANCE MATRIX OF ESTIMATE", txt, ignore.case = TRUE)
  cov_singular  <- grepl("MATRIX ALGORITHMICALLY SINGULAR|R MATRIX.*SINGULAR|S MATRIX.*SINGULAR",
                         txt, ignore.case = TRUE)
  cov_warning   <- grepl("EIGENVALUE OF COV", txt, ignore.case = TRUE)  # ill-conditioning

  cov_ok <- cov_requested & cov_matrix_ok & !cov_aborted & !cov_omitted & !cov_singular

  data.frame(min_ok = min_ok, cov_ok = cov_ok, cov_warning = cov_warning,
             nmtran_error = nmtran_err)
}

#' Best-effort initial value of a $THETA token (display only)
#' @noRd
.nm_theta_init_value <- function(tok) {
  body <- sub("\\s*FIX(?:ED)?\\s*$", "", tok, ignore.case = TRUE)
  if (grepl("^\\(", body)) {
    inside <- sub("^\\(\\s*", "", sub("\\s*\\)$", "", body))
    parts  <- trimws(strsplit(inside, ",")[[1]])
    nums   <- suppressWarnings(as.numeric(parts))
    nums_ok <- nums[!is.na(nums)]
    if (length(nums_ok) == 0) return(NA_real_)
    if (length(parts) >= 2)  return(suppressWarnings(as.numeric(parts[2])))  # (low,init[,up])
    return(nums_ok[1])                                                       # (init)
  }
  suppressWarnings(as.numeric(body))
}

#' Enumerate every $THETA initial estimate with its exact character span
#'
#' Handles continuation lines, several $THETA records, several tokens per line,
#' and a token sitting on the "$THETA ..." header line itself.
#' @noRd
.nm_parse_theta_tokens <- function(ctl_lines) {
  in_theta <- FALSE
  rows <- list(); k <- 0L
  for (li in seq_along(ctl_lines)) {
    line      <- ctl_lines[li]
    is_header <- grepl("^\\s*\\$THETA", line, ignore.case = TRUE)
    if (is_header) {
      in_theta <- TRUE
    } else if (grepl("^\\s*\\$", line)) {
      in_theta <- FALSE
    }
    if (!in_theta) next

    cpos    <- regexpr(";", line, fixed = TRUE)
    code    <- if (cpos > 0) substr(line, 1, cpos - 1) else line
    comment <- if (cpos > 0) trimws(sub("^;+\\s*", "", substr(line, cpos, nchar(line)))) else ""

    offset <- 0L
    scan   <- code
    if (is_header) {
      hm     <- regexpr("^\\s*\\$THETA\\b", code, ignore.case = TRUE, perl = TRUE)
      offset <- attr(hm, "match.length")
      scan   <- substring(code, offset + 1)
    }

    g <- gregexpr(.NM_THETA_TOK, scan, perl = TRUE)[[1]]
    if (g[1] == -1) next
    lens <- attr(g, "match.length")
    for (t in seq_along(g)) {
      k <- k + 1L
      abs_start <- offset + g[t]
      abs_end   <- offset + g[t] + lens[t] - 1L
      tok       <- substr(line, abs_start, abs_end)
      rows[[k]] <- data.frame(
        theta_index = k, line = li, start = abs_start, end = abs_end,
        token = tok, value = .nm_theta_init_value(tok),
        fixed = grepl("FIX", tok, ignore.case = TRUE), label = comment,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0)
    return(data.frame(theta_index = integer(0), line = integer(0), start = integer(0),
                      end = integer(0), token = character(0), value = numeric(0),
                      fixed = logical(0), label = character(0), stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

#' Resolve `param` (index, "THETA(n)", "n", or a $THETA label) to one row
#' @noRd
.nm_find_target_theta <- function(tokens, param) {
  if (nrow(tokens) == 0)
    stop("No $THETA record / no initial estimate found in the control file.")

  is_num <- grepl("^\\s*(THETA\\s*\\(?\\s*)?[0-9]+\\s*\\)?\\s*$", param, ignore.case = TRUE)
  if (is_num) {
    idx <- as.integer(sub("^\\s*THETA\\s*\\(?\\s*([0-9]+)\\s*\\)?\\s*$", "\\1",
                          param, ignore.case = TRUE))
    if (idx < 1 || idx > nrow(tokens))
      stop("THETA index ", idx, " out of range (model has ", nrow(tokens), " THETAs).")
    return(tokens[tokens$theta_index == idx, , drop = FALSE])
  }

  esc <- gsub("([.\\^$*+?()\\[\\]{}|])", "\\\\\\1", param, perl = TRUE)
  pat <- paste0("(^|[^A-Za-z0-9_])", esc, "([^A-Za-z0-9_]|$)")
  hit <- grepl(pat, tokens$label, ignore.case = TRUE, perl = TRUE)
  if (sum(hit) == 1) return(tokens[hit, , drop = FALSE])
  if (sum(hit) > 1)
    stop("Parameter name '", param, "' is ambiguous; it matches THETA(",
         paste(tokens$theta_index[hit], collapse = "), THETA("), ").\n",
         "Disambiguate by passing the THETA index instead.")
  stop("Parameter '", param, "' not found among the $THETA labels.\n",
       "Available labels:\n",
       paste0("  THETA(", tokens$theta_index, ")  value=", signif(tokens$value, 4),
              "  ; ", ifelse(nzchar(tokens$label), tokens$label, "(no label)"),
              collapse = "\n"))
}

#' Build one profiling control file (target THETA fixed to `value`)
#' @noRd
.nm_write_profile_ctl <- function(base_ctl, target, value, digits,
                                  data_abs, strip_tables, strip_cov) {
  lines <- base_ctl
  li    <- target$line
  repl  <- paste0(.nm_fmt_val(value, digits), " FIX")
  lines[li] <- paste0(substr(lines[li], 1, target$start - 1), repl,
                      substr(lines[li], target$end + 1, nchar(lines[li])))

  if (strip_tables) lines <- .nm_strip_record(lines, "TABLE")
  if (strip_cov)    lines <- .nm_strip_record(lines, "COV")

  idx <- grep("^\\s*\\$DATA", lines, ignore.case = TRUE)[1]
  if (is.na(idx)) stop("No $DATA record to patch.")
  lines[idx] <- .nm_patch_data_line(lines[idx], data_abs)
  lines
}

#' Launch one NONMEM run in `run_dir` (working directory switched locally)
#'
#' Quotes ONLY the executable so cmd.exe /c keeps exactly two quotes; the
#' generated .ctl/.lst names never contain spaces and are passed unquoted.
#' To reuse the package's own launcher, replace the body with a call to nm_run().
#' @noRd
.nm_launch <- function(nmfe_path, run_dir, run_name) {
  old <- getwd(); on.exit(setwd(old), add = TRUE)
  setwd(run_dir)
  lst <- paste0(run_name, ".lst")
  if (.Platform$OS.type == "windows") {
    nmfe_win <- normalizePath(nmfe_path, winslash = "\\", mustWork = FALSE)
    system2("cmd.exe", args = c("/c", shQuote(nmfe_win, type = "cmd"),
                                paste0(run_name, ".ctl"), lst))
  } else {
    .nm_check_nmfe_exec(nmfe_path)
    system2(nmfe_path, args = c(paste0(run_name, ".ctl"), lst))
  }
}

#' On Linux and macOS, nmfe is run directly: it must be executable
#' @noRd
.nm_check_nmfe_exec <- function(nmfe_path) {
  if (.Platform$OS.type != "windows" && file.exists(nmfe_path) &&
      file.access(nmfe_path, 1L) != 0L)
    stop("The NONMEM launcher is not executable: ", nmfe_path,
         "\n  Make it executable (e.g. chmod +x) or check the path.", call. = FALSE)
  invisible(TRUE)
}

#' Linear-interpolation crossings of a curve through `level`, bracketing `x_min`
#' @noRd
.nm_ci_from_curve <- function(x, y, level, x_min) {
  o <- order(x); x <- x[o]; y <- y[o]
  cross <- numeric(0)
  for (j in seq_len(length(x) - 1)) {
    y1 <- y[j]; y2 <- y[j + 1]
    if (is.na(y1) || is.na(y2)) next
    if ((y1 - level) * (y2 - level) < 0)
      cross <- c(cross, x[j] + (level - y1) / (y2 - y1) * (x[j + 1] - x[j]))
  }
  lo <- cross[cross < x_min]; hi <- cross[cross > x_min]
  c(if (length(lo)) max(lo) else NA_real_,
    if (length(hi)) min(hi) else NA_real_)
}

#-------------------------------------------------------------------------------
# Exported function
#-------------------------------------------------------------------------------

#' Profile the objective function over one model parameter
#'
#' Runs a series of NONMEM estimations in which a single `$THETA` is fixed to
#' each value of a grid while every other parameter is re-estimated, yielding the
#' profile of the objective function value (OFV, *without* the additive
#' constant) for that parameter. The minimum of the curve locates the
#' conditional estimate; the width at `min(OFV) + ci_threshold` gives an
#' approximate likelihood-based confidence interval.
#'
#' Only `$DATA` is remapped (to an absolute path) and `$TABLE`/`$COV` are
#' stripped by default; every other record is preserved byte-for-byte. Each run
#' is written to its own sub-directory under `output_dir/runs/`.
#'
#' @param nmfe_path Path to the NONMEM launcher (e.g. `nmfe76.bat` on Windows,
#'   the `nmfe76` script on Linux or macOS).
#' @param ctl_path Path to the reference control file.
#' @param param The parameter to profile: a `$THETA` label (matched whole-word,
#'   case-insensitive, e.g. `"CL"`) or an index (`3`, `"3"`, or `"THETA(3)"`).
#' @param from,to,by Numeric grid: first value, last value, and increment.
#' @param output_dir Output directory. Defaults to
#'   `<ctl_dir>/<ctl>_<param>_profile`.
#' @param strip_tables,strip_cov Logical; drop `$TABLE` / `$COV` from the
#'   generated controls. Both `TRUE` by default (recommended for a profile).
#' @param overwrite Logical; if `FALSE` (default) a run whose `.lst` already
#'   exists and minimised successfully is skipped (resume after interruption).
#' @param ci_threshold Numeric; delta-OFV cut-off for the interval. Default `3.84`
#'   (chi-square, 1 df, two-sided 95%). Set `NA` to disable. Note: for a
#'   variance component tested at its boundary (\eqn{\omega^2=0}) the appropriate
#'   one-sided cut-off is `2.71`, not `3.84`.
#' @param digits Decimals used to format fixed values and run labels. Auto-
#'   derived from `by` when `NULL`.
#' @param run Logical; if `FALSE`, generate the control files without launching
#'   NONMEM.
#' @param quiet Logical; suppress progress messages.
#' @param n_cores Number of grid values run at the same time (default 1). A
#'   whole number,
#'   `"auto"` (physical cores minus one, see [nm_cores()]) or `"ask"`
#'   (interactive sessions: the cores are shown and the number is asked once per R
#'   session, then reused; `options(nmpilot.n_cores = NULL)` asks again).
#'   Each run is then started in the background in its own directory; if R
#'   is interrupted, the runs already started go on in the background and a
#'   later call with the same arguments picks up their results.
#'
#' @return An object of class `nm_profile`: a list with `results` (a data frame
#'   of `run`, `value`, `ofv`, `dOFV`, `min_ok`, `cov_ok`, `cov_warning`, `dir`,
#'   `lst`), plus `ci`, `min_ofv`, `min_value`, `param`, `param_safe`, `target`,
#'   `ci_threshold`, and `output_dir`. Use [print()] for a summary and [plot()]
#'   (with `type = "ofv"` or `"dofv"`) for figures.
#'
#' @seealso [nm_run()], [nm_extract_runs()], [nm_compare()]
#'
#' @examples
#' \dontrun{
#' prof <- nm_profile(
#'   nmfe_path = "C:/nm76g64/run/nmfe76.bat",
#'   ctl_path  = "C:/nm76g64/run/test/run_test.reg/run_test.ctl",
#'   param     = "CL", from = 0.1, to = 3, by = 0.05
#' )
#' prof                       # summary (print.nm_profile)
#' plot(prof)                 # OFV vs value
#' plot(prof, type = "dofv")  # dOFV vs value, with the 3.84 threshold
#' }
#' @export
nm_profile <- function(nmfe_path, ctl_path, param, from, to, by,
                       output_dir   = NULL,
                       strip_tables = TRUE,
                       strip_cov    = TRUE,
                       overwrite    = FALSE,
                       ci_threshold = 3.84,
                       digits       = NULL,
                       run          = TRUE,
                       quiet        = FALSE,
                       n_cores      = 1L) {

  say <- function(...) if (!quiet) message(...)

  if (run && !file.exists(nmfe_path))
    stop("NONMEM launcher not found: ", nmfe_path)
  if (!file.exists(ctl_path))
    stop("Control file not found: ", ctl_path)
  if (!(is.numeric(from) && is.numeric(to) && is.numeric(by)) || by <= 0)
    stop("`from`, `to` must be numeric and `by` a positive numeric.")
  if (from > to) stop("`from` (", from, ") must be <= `to` (", to, ").")

  ctl_path <- normalizePath(ctl_path, winslash = "/", mustWork = TRUE)
  base     <- tools::file_path_sans_ext(basename(ctl_path))

  if (is.null(digits))
    digits <- max(.nm_decimals(by), .nm_decimals(from), .nm_decimals(to), 1L)
  vals <- unique(round(seq(from, to, by = by), digits))
  if (length(vals) == 0) stop("Empty value grid; check from/to/by.")

  data_abs <- .nm_extract_data_path(ctl_path)
  if (!file.exists(data_abs))
    warning("Data file referenced by $DATA does not exist: ", data_abs)

  base_ctl <- readLines(ctl_path, warn = FALSE)
  tokens   <- .nm_parse_theta_tokens(base_ctl)
  target   <- .nm_find_target_theta(tokens, as.character(param))

  param_safe <- gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", as.character(param)))
  if (!nzchar(param_safe)) param_safe <- paste0("THETA", target$theta_index)

  if (is.null(output_dir))
    output_dir <- file.path(dirname(ctl_path), paste0(base, "_", param_safe, "_profile"))
  .nm_ensure_dir(output_dir)
  runs_dir <- .nm_ensure_dir(file.path(output_dir, "runs"))

  # final estimate from the .ext of the control stream, if it was already run
  est <- NA_real_
  ext <- sub("\\.[^./]*$", ".ext", ctl_path)
  if (file.exists(ext)) {
    e <- tryCatch(.nm_parse_ext(ext), error = function(e) NULL)
    if (!is.null(e) && isTRUE(e$final_row_available)) {
      v <- e$raw$value[e$raw$raw_name == paste0("THETA", target$theta_index)]
      if (length(v)) est <- v[1]
    }
  }
  say("Profiling THETA(", target$theta_index, ")",
      if (nzchar(target$label)) paste0(" [", target$label, "]") else "",
      "  initial value = ", signif(target$value, 4),
      if (!is.na(est)) paste0(", final estimate = ", signif(est, 4)) else "",
      if (isTRUE(target$fixed)) "  (was FIXED)" else "")
  say("Grid: ", length(vals), " values from ", .nm_fmt_val(min(vals), digits),
      " to ", .nm_fmt_val(max(vals), digits), " by ", by)

  runs <- data.frame(
    run   = paste0(base, "_", param_safe, "_", vapply(vals, .nm_tag_val, "", digits = digits)),
    value = vals, stringsAsFactors = FALSE
  )
  runs$dir <- file.path(runs_dir, runs$run)

  for (i in seq_len(nrow(runs))) {
    .nm_ensure_dir(runs$dir[i])
    ctl_i <- .nm_write_profile_ctl(base_ctl, target, runs$value[i], digits,
                                   data_abs, strip_tables, strip_cov)
    writeLines(ctl_i, file.path(runs$dir[i], paste0(runs$run[i], ".ctl")))
  }
  say(nrow(runs), " control files generated in ", runs_dir)

  if (run && (n_cores <- .nm_resolve_cores(n_cores)) > 1L) {
    todo <- lapply(seq_len(nrow(runs)), function(i) i)
    .nm_run_queue(todo, nmfe_path, n_cores, verbose = !quiet, label = "grid values",
      prepare = function(i, slot) {
        lst_i <- file.path(runs$dir[i], paste0(runs$run[i], ".lst"))
        if (!overwrite && file.exists(lst_i) && isTRUE(.nm_parse_status(lst_i)$min_ok)) {
          say("[", i, "/", nrow(runs), "] skip (done): ", runs$run[i]); return(NULL)
        }
        list(dir = runs$dir[i], name = runs$run[i])
      },
      finish = function(i, slot, st, secs) {
        if (!identical(as.integer(st), 0L))
          warning("Non-zero exit status (", st, ") for ", runs$run[i])
      })
  } else if (run) {
    for (i in seq_len(nrow(runs))) {
      lst_i <- file.path(runs$dir[i], paste0(runs$run[i], ".lst"))
      if (!overwrite && file.exists(lst_i) && isTRUE(.nm_parse_status(lst_i)$min_ok)) {
        say("[", i, "/", nrow(runs), "] skip (done): ", runs$run[i]); next
      }
      say("[", i, "/", nrow(runs), "] running: ", runs$run[i],
          "  (", param_safe, " = ", .nm_fmt_val(runs$value[i], digits), ")")
      st <- .nm_launch(nmfe_path, runs$dir[i], runs$run[i])
      if (!identical(as.integer(st), 0L))
        warning("Non-zero exit status (", st, ") for ", runs$run[i])
    }
  } else {
    say("run = FALSE: control files generated but NONMEM not launched.")
  }

  # Collect results
  runs$lst         <- file.path(runs$dir, paste0(runs$run, ".lst"))
  runs$ofv         <- NA_real_
  runs$min_ok      <- NA
  runs$cov_ok      <- NA
  runs$cov_warning <- NA
  for (i in seq_len(nrow(runs))) {
    if (!file.exists(runs$lst[i])) next   # missing .lst => terminal / NMTRAN error
    runs$ofv[i] <- .nm_extract_ofv(runs$lst[i])
    st <- .nm_parse_status(runs$lst[i])
    runs$min_ok[i]      <- st$min_ok
    runs$cov_ok[i]      <- st$cov_ok
    runs$cov_warning[i] <- st$cov_warning
  }

  ok        <- which(runs$min_ok %in% TRUE & is.finite(runs$ofv))
  min_ofv   <- if (length(ok)) min(runs$ofv[ok]) else NA_real_
  min_value <- if (length(ok)) runs$value[ok][which.min(runs$ofv[ok])] else NA_real_
  runs$dOFV <- runs$ofv - min_ofv

  runs <- runs[order(runs$value),
               c("run", "value", "ofv", "dOFV", "min_ok", "cov_ok", "cov_warning", "dir", "lst")]
  rownames(runs) <- NULL

  ci <- c(NA_real_, NA_real_)
  if (!is.na(ci_threshold) && length(ok) >= 2 && is.finite(min_ofv)) {
    s <- runs[runs$min_ok %in% TRUE & is.finite(runs$dOFV), ]
    ci <- .nm_ci_from_curve(s$value, s$dOFV, ci_threshold, min_value)
  }

  structure(
    list(results = runs, ci = ci, min_ofv = min_ofv, min_value = min_value,
         param = param, param_safe = param_safe, target = target,
         ci_threshold = ci_threshold, digits = digits, output_dir = output_dir),
    class = "nm_profile"
  )
}

#-------------------------------------------------------------------------------
# S3 methods
#-------------------------------------------------------------------------------

#' @param x An `nm_profile` object.
#' @param ... Ignored.
#' @rdname nm_profile
#' @export
print.nm_profile <- function(x, ...) {
  n_ok <- sum(x$results$min_ok %in% TRUE)
  cat("<nm_profile>\n")
  cat("  Parameter : THETA(", x$target$theta_index, ")",
      if (nzchar(x$target$label)) paste0(" [", x$target$label, "]") else "", "\n", sep = "")
  cat("  Grid      : ", nrow(x$results), " values [",
      .nm_fmt_val(min(x$results$value), x$digits), ", ",
      .nm_fmt_val(max(x$results$value), x$digits), "]  (",
      n_ok, " minimised)\n", sep = "")
  if (is.finite(x$min_ofv))
    cat("  Minimum   : OFV = ", formatC(x$min_ofv, format = "f", digits = 3),
        " at ", x$param_safe, " = ", .nm_fmt_val(x$min_value, x$digits), "\n", sep = "")
  if (!is.na(x$ci_threshold) && all(is.finite(x$ci)))
    cat("  ~95% CI   : [", .nm_fmt_val(x$ci[1], x$digits + 1), ", ",
        .nm_fmt_val(x$ci[2], x$digits + 1), "]  (dOFV = ",
        x$ci_threshold, ")\n", sep = "")
  else if (!is.na(x$ci_threshold))
    cat("  ~95% CI   : not bracketed within the grid (widen the range)\n")
  cat("  Output    : ", x$output_dir, "\n", sep = "")
  invisible(x)
}

#' @param type Which figure to draw: `"ofv"` (OFV vs value) or `"dofv"`
#'   (delta-OFV vs value, with the `ci_threshold` line).
#' @rdname nm_profile
#' @export
plot.nm_profile <- function(x, type = c("ofv", "dofv"), ...) {
  type <- match.arg(type)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot.nm_profile(); install it to draw figures.")

  d      <- x$results
  ok_df  <- d[d$min_ok %in% TRUE & is.finite(d$ofv), ]
  if (nrow(ok_df) == 0)
    stop("No successfully minimised run to plot; inspect the .lst files.")
  xlab   <- paste0(x$param_safe, " (fixed value)")

  if (type == "ofv") {
    d2 <- d[is.finite(d$ofv), ]
    d2$status <- ifelse(d2$min_ok %in% TRUE, "minimised", "not minimised")
    p <- ggplot2::ggplot(d2, ggplot2::aes(x = .data$value, y = .data$ofv)) +
      ggplot2::geom_line(data = ok_df,
                         ggplot2::aes(x = .data$value, y = .data$ofv),
                         colour = "steelblue") +
      ggplot2::geom_point(ggplot2::aes(shape = .data$status, colour = .data$status),
                          size = 2.4) +
      ggplot2::scale_colour_manual(values = c("minimised" = "steelblue",
                                              "not minimised" = "firebrick")) +
      ggplot2::scale_shape_manual(values = c("minimised" = 16, "not minimised" = 4)) +
      ggplot2::geom_vline(xintercept = x$min_value, linetype = "dotted", colour = "grey40") +
      ggplot2::theme_bw() +
      ggplot2::labs(x = xlab, y = "OFV (without constant)",
                    title = paste0("OFV profile - ", x$param_safe),
                    shape = NULL, colour = NULL)
    return(p)
  }

  # type == "dofv"
  p <- ggplot2::ggplot(ok_df, ggplot2::aes(x = .data$value, y = .data$dOFV)) +
    ggplot2::geom_line(colour = "steelblue") +
    ggplot2::geom_point(size = 2.2, colour = "steelblue") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = xlab, y = expression(Delta * "OFV"),
                  title = paste0("Profile likelihood - ", x$param_safe))
  if (!is.na(x$ci_threshold)) {
    p <- p +
      ggplot2::geom_hline(yintercept = x$ci_threshold, linetype = "dashed", colour = "firebrick") +
      ggplot2::annotate("text", x = min(ok_df$value), y = x$ci_threshold,
                        label = sprintf("dOFV = %.2f", x$ci_threshold),
                        hjust = 0, vjust = -0.4, size = 3, colour = "firebrick")
    if (all(is.finite(x$ci)))
      p <- p + ggplot2::geom_vline(xintercept = x$ci, linetype = "dotted", colour = "grey40")
  }
  p
}
