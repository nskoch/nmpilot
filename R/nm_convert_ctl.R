# =============================================================================
# nm_convert_ctl.R
# Convert a Wings for NONMEM (WFN) extended control stream to standard NONMEM.
#
# WFN extended control streams use parameter names (e.g. POP_CL, PPV_CL)
# instead of THETA(n)/ETA(n)/ERR(n) references in model records ($PK, $PRED,
# $ERROR, $DES, $AES). This function parses the parameter names from $THETA,
# $OMEGA and $SIGMA comment fields and substitutes them back into standard
# NONMEM syntax.
#
# Reference: https://wfn.sourceforge.net/wfncs.htm#control_streams
# =============================================================================


# =============================================================================
# Internal helpers
# =============================================================================

# Extract parameter name from a $THETA/$OMEGA/$SIGMA comment line.
# Rules (from WFN docs):
#   - Name = first non-numeric word after ";"
#   - Name ends before "(", "{", "[", ";" or whitespace
.nm_extract_param_name <- function(line) {
  # Find comment
  m <- regexpr(";\\s*", line)
  if (m == -1L) return(NA_character_)
  after <- substring(line, m + attr(m, "match.length"))
  # Skip leading numeric tokens (e.g. "1 POP_E0" -> skip "1")
  after <- trimws(after)
  tokens <- strsplit(after, "\\s+")[[1]]
  tokens <- tokens[tokens != ""]
  if (length(tokens) == 0L) return(NA_character_)
  # First token that is not purely numeric
  for (tok in tokens) {
    tok_clean <- sub("[({\\[;].*$", "", tok)  # strip trailing bracket
    if (!grepl("^[0-9\\.\\-\\+eE]+$", tok_clean) && nchar(tok_clean) > 0L)
      return(tok_clean)
  }
  NA_character_
}



# Mappings name -> "THETA(n)", "ETA(n)", "ERR(n)", built on the shared readers
# of $THETA, $OMEGA and $SIGMA (nm_ctl_read.R, action 11c). The name is read in
# the comment of the line of each parameter (for a BLOCK, the line of its
# diagonal value); a line holding several THETAs names the first one.
.nm_theta_name_mapping <- function(ctl_lines) {
  th <- .nm_ctl_theta_parse(ctl_lines)
  mapping <- list()
  first <- !duplicated(th$line)
  for (i in which(first)) {
    name <- .nm_extract_param_name(th$raw[i])
    if (!is.na(name) && nchar(name) > 0L) mapping[[name]] <- sprintf("THETA(%d)", th$theta[i])
  }
  mapping
}

.nm_random_name_mapping <- function(ctl_lines, rec = c("OMEGA", "SIGMA")) {
  rec <- match.arg(rec)
  prefix <- if (rec == "OMEGA") "ETA" else "ERR"
  et <- .nm_ctl_random_parse(gsub("\r", "", ctl_lines, fixed = TRUE), rec)$eta
  mapping <- list()
  for (i in seq_len(nrow(et))) {
    cm <- et$omega_comment[i]
    if (is.na(cm)) next
    name <- .nm_extract_param_name(paste0(";", cm))
    if (!is.na(name) && nchar(name) > 0L) mapping[[name]] <- sprintf("%s(%d)", prefix, et$eta[i])
  }
  mapping
}

# Apply mapping to model records ($PK, $PRED, $ERROR, $DES, $AES)
# WFN rule: once a name appears on the LHS of an assignment, further RHS
# occurrences in subsequent lines are NOT substituted.
.nm_substitute_model <- function(ctl_lines, mapping) {
  model_records <- c("PK", "PRED", "ERROR", "DES", "AES", "MIX")
  pattern_model <- paste0("^\\s*\\$(", paste(model_records, collapse="|"), ")\\b")

  in_model  <- FALSE
  lhs_seen  <- character(0)  # names that have appeared on LHS
  out_lines <- character(length(ctl_lines))

  for (i in seq_along(ctl_lines)) {
    ln <- ctl_lines[i]

    # Detect model record start
    if (grepl(pattern_model, ln, ignore.case = TRUE)) {
      in_model <- TRUE
      # lhs_seen is kept across model records: a variable assigned in $PK
      # is still a variable in $DES and $ERROR
      out_lines[i] <- ln
      next
    }
    # Detect end of model record
    if (in_model && grepl("^\\s*\\$[A-Za-z]", ln)) {
      in_model <- FALSE
    }

    if (!in_model) {
      out_lines[i] <- ln
      next
    }

    # Substitute: split line into code and comment parts
    comment_part <- ""
    code_part    <- ln
    sc_pos <- regexpr(";", ln, fixed = TRUE)
    if (sc_pos > 0L) {
      code_part    <- substring(ln, 1L, sc_pos - 1L)
      comment_part <- substring(ln, sc_pos)
    }

    # Detect LHS: word before "=" (ignoring spaces), also after a leading
    # IF(condition): "IF(PATIENT.EQ.1) SCALAR = ..."
    lhs_txt <- .nm_ip_split_if(code_part)$rest
    lhs_m <- regmatches(lhs_txt,
                        regexpr("^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*(\\([^)]*\\))?\\s*=(?!=)",
                                lhs_txt, perl = TRUE))
    if (length(lhs_m) > 0L) {
      # name only, without an index: MTIME(1) = ... -> MTIME
      lhs_name <- sub("^\\s*([A-Za-z_][A-Za-z0-9_]*).*$", "\\1", lhs_m)
      lhs_seen <- unique(c(lhs_seen, lhs_name))
    }

    # Substitute names on RHS (names not in lhs_seen)
    for (nm in names(mapping)) {
      if (nm %in% lhs_seen) next
      # Word-boundary replacement
      code_part <- gsub(paste0("\\b", nm, "\\b"), mapping[[nm]], code_part)
    }

    out_lines[i] <- paste0(code_part, comment_part)
  }

  out_lines
}


# =============================================================================
# Main exported function
# =============================================================================

#' Convert a Wings for NONMEM extended control stream to standard NONMEM
#'
#' Wings for NONMEM (WFN) supports an extended control stream format where
#' parameter names (e.g. \code{POP_CL}, \code{PPV_CL}, \code{RUV_PROP}) are
#' used directly in model records (\code{$PK}, \code{$PRED}, \code{$ERROR},
#' \code{$DES}) instead of \code{THETA(n)}, \code{ETA(n)}, \code{ERR(n)}.
#' This function parses the parameter names from \code{$THETA}, \code{$OMEGA}
#' and \code{$SIGMA} comment fields and substitutes them back into standard
#' NONMEM syntax, producing a control stream compatible with direct
#' \code{nmfe} execution.
#'
#' @param ctl Path to the Wings extended control stream (\code{.ctl}).
#'   (Formerly \code{ctl_file}, still accepted with a warning.)
#' @param output_file Output path for the converted control stream. If
#'   \code{NULL} (default), the output file is named
#'   \code{<run>_converted.ctl} in the same directory as \code{ctl}. (Formerly `out_file`, still accepted with a warning.)
#' @param overwrite Logical. If \code{TRUE}, overwrite \code{output_file} if it
#'   already exists (default \code{FALSE}).
#' @param verbose Logical. If \code{TRUE} (default), print a summary of
#'   substitutions made.
#'
#' @param ... Former argument name, accepted with a warning until version
#'   1.0: \code{ctl_file} (now \code{ctl}).
#' @return Invisibly returns the path to the converted control stream.
#'
#' @details
#' The conversion follows the WFN specification:
#' \itemize{
#'   \item Parameter names are extracted from the first non-numeric word
#'     after the \code{;} comment delimiter in each \code{$THETA},
#'     \code{$OMEGA} and \code{$SIGMA} line.
#'   \item Names are substituted by \code{THETA(n)}, \code{ETA(n)} and
#'     \code{ERR(n)} in model records (\code{$PK}, \code{$PRED},
#'     \code{$ERROR}, \code{$DES}, \code{$AES}).
#'   \item Once a name appears on the left-hand side of an assignment, it
#'     is no longer substituted on subsequent right-hand sides (WFN rule).
#'   \item \code{$OMEGA BLOCK} structures are handled correctly.
#' }
#'
#' @examples
#' \dontrun{
#' # Convert a WFN extended control stream
#' out <- nm_convert_ctl("C:/analyses/run001_wfn.ctl")
#'
#' # Then run with nmpilot
#' log <- nm_run(
#'   ctl          = out,
#'   nmfe_path    = "C:/nm76g64/run/nmfe76.bat",
#'   analysis_dir = "C:/analyses"
#' )
#' }
#'
#' @seealso \code{\link{nm_run}}
#'
#' @export
nm_convert_ctl <- function(ctl,
                           output_file = NULL,
                           overwrite = FALSE,
                           verbose   = TRUE,
                           ...) {
  dep <- .nm_deprecated_args(list(...), c(ctl_file = "ctl", out_file = "output_file"),
                             "nm_convert_ctl")
  if (!is.null(dep$ctl)) ctl <- dep$ctl
  if (!is.null(dep$output_file)) output_file <- dep$output_file
  out_file <- output_file
  ctl_file <- ctl

  if (!file.exists(ctl_file))
    stop("Control stream not found: ", ctl_file)

  # Default output path
  if (is.null(out_file)) {
    base    <- tools::file_path_sans_ext(ctl_file)
    ext     <- tools::file_ext(ctl_file)
    out_file <- paste0(base, "_converted.", ext)
  }

  if (file.exists(out_file) && !overwrite)
    stop("Output file already exists (use overwrite = TRUE): ", out_file)

  ctl_lines <- readLines(ctl_file, warn = FALSE)
  # Normalize line endings (remove Windows \r)
  ctl_lines <- gsub("\r", "", ctl_lines, fixed = TRUE)

  # --- Build mappings ---------------------------------------------------------
  theta_map <- .nm_theta_name_mapping(ctl_lines)
  eta_map   <- .nm_random_name_mapping(ctl_lines, "OMEGA")
  err_map   <- .nm_random_name_mapping(ctl_lines, "SIGMA")
  mapping   <- c(theta_map, eta_map, err_map)

  if (verbose) {
    message("=== WFN parameter mapping ===")
    if (length(mapping) == 0L) {
      message("  No parameter names found - file may already be standard NONMEM.")
    } else {
      for (nm in names(mapping))
        message(sprintf("  %-20s -> %s", nm, mapping[[nm]]))
    }
  }

  # --- Substitute in model records --------------------------------------------
  out_lines <- .nm_substitute_model(ctl_lines, mapping)

  # --- Write output -----------------------------------------------------------
  writeLines(out_lines, out_file)
  if (verbose) message("Written: ", out_file)

  invisible(out_file)
}
