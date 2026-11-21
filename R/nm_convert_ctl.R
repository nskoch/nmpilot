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

# Parse $OMEGA BLOCK lines: each ROW gets the name from its diagonal comment.
# Returns a vector of names in ETA order (one per diagonal element).
.nm_parse_omega_block <- function(block_lines) {
  # Find BLOCK dimension
  dim_m <- regmatches(block_lines[1],
                      regexpr("BLOCK\\s*\\(\\s*([0-9]+)\\s*\\)",
                              block_lines[1], ignore.case = TRUE))
  if (length(dim_m) == 0L) return(character(0))
  dim_n <- as.integer(gsub("[^0-9]", "", dim_m))

  names_out <- character(dim_n)
  row_idx   <- 0L
  for (ln in block_lines[-1]) {
    if (grepl("^\\s*;;", ln)) next   # SAME convention comment
    # Count numeric tokens on this line (before comment)
    code_part <- sub(";.*$", "", ln)
    nums <- strsplit(trimws(code_part), "\\s+")[[1]]
    nums <- nums[nums != "" & grepl("^[0-9\\.\\-\\+eE]+$", nums)]
    if (length(nums) == 0L) next
    row_idx <- row_idx + 1L
    # Diagonal element = last numeric on row = name from comment
    name <- .nm_extract_param_name(ln)
    if (!is.na(name)) names_out[row_idx] <- name
    if (row_idx >= dim_n) break
  }
  names_out
}

# Build mapping list from a record type
# Returns list: name -> substitution string e.g. "POP_CL" -> "THETA(1)"
.nm_build_mapping <- function(ctl_lines, record_pattern,
                           sub_prefix, start_idx = 1L) {
  mapping <- list()
  in_rec  <- FALSE
  counter <- start_idx - 1L

  for (ln in ctl_lines) {
    # Detect record start ($THETA, $SIGMA)
    if (grepl(record_pattern, ln, ignore.case = TRUE)) {
      if (grepl("BLOCK", ln, ignore.case = TRUE) && sub_prefix == "ETA") {
        next  # handled by .nm_build_eta_mapping
      }
      in_rec <- TRUE

      # Format B (inline): "$THETA (0.1,10) ;POPV2 ;L" - value and name on same line
      # Extract: everything after the $THETA keyword
      m <- regexpr("\\$[A-Za-z]+\\s*", ln, ignore.case = TRUE)
      if (m > 0) rest <- substring(ln, m + attr(m, "match.length"))
      rest <- trimws(rest)

      # Check if this line has a numeric value (not just $THETA alone)
      code_part <- sub(";.*$", "", rest)
      if (grepl("[0-9]", code_part) && grepl(";", rest)) {
        counter <- counter + 1L
        name <- .nm_extract_param_name(rest)
        if (!is.na(name) && nchar(name) > 0L)
          mapping[[name]] <- sprintf("%s(%d)", sub_prefix, counter)
      }
      next
    }

    # Detect next record (stop)
    if (in_rec && grepl("^\\s*\\$[A-Za-z]", ln)) {
      in_rec <- FALSE
    }
    if (!in_rec) next

    # Skip pure comment or empty lines
    if (grepl("^\\s*;", ln) || grepl("^\\s*$", ln)) next
    # Skip SAME lines
    if (grepl("SAME", ln, ignore.case = TRUE)) next

    # Format A: value on its own line (multi-line $THETA block)
    code_part <- sub(";.*$", "", ln)
    if (!grepl("[0-9]", code_part)) next

    counter <- counter + 1L
    name    <- .nm_extract_param_name(ln)
    if (!is.na(name) && nchar(name) > 0L) {
      mapping[[name]] <- sprintf("%s(%d)", sub_prefix, counter)
    }
  }
  mapping
}
# Build ETA mapping handling both diagonal and BLOCK $OMEGA
.nm_build_eta_mapping <- function(ctl_lines) {
  mapping     <- list()
  eta_counter <- 0L
  i           <- 1L

  while (i <= length(ctl_lines)) {
    ln <- gsub("\r", "", ctl_lines[i], fixed = TRUE)

    if (grepl("^\\s*\\$OMEGA", ln, ignore.case = TRUE)) {

      if (grepl("SAME", ln, ignore.case = TRUE)) {
        i <- i + 1L
        next
      }

      if (grepl("BLOCK", ln, ignore.case = TRUE)) {
        # Collect block lines until next $record
        block <- ln
        j <- i + 1L
        while (j <= length(ctl_lines) &&
               !grepl("^\\s*\\$[A-Za-z]", ctl_lines[j])) {
          block <- c(block, gsub("\r", "", ctl_lines[j], fixed = TRUE))
          j <- j + 1L
        }
        names_block <- .nm_parse_omega_block(block)
        for (nm in names_block) {
          eta_counter <- eta_counter + 1L
          if (!is.na(nm) && nchar(nm) > 0L)
            mapping[[nm]] <- sprintf("ETA(%d)", eta_counter)
        }
        i <- j
        next
      }

      # Format B (inline): "$OMEGA 0.5 ; ETA_CL" - value and name on same line
      m <- regexpr("\\$OMEGA\\s*", ln, ignore.case = TRUE)
      rest <- if (m > 0) trimws(substring(ln, m + attr(m, "match.length"))) else ""
      code_part <- sub(";.*$", "", rest)

      if (grepl("[0-9]", code_part) && grepl(";", rest)) {
        # Inline: one parameter on this line; further values of the same
        # record on the following lines are read below (Format A).
        eta_counter <- eta_counter + 1L
        name <- .nm_extract_param_name(rest)
        if (!is.na(name) && nchar(name) > 0L)
          mapping[[name]] <- sprintf("ETA(%d)", eta_counter)
      }

      # Format A: multi-line $OMEGA block (diagonal, non-BLOCK)
      j <- i + 1L
      while (j <= length(ctl_lines) &&
             !grepl("^\\s*\\$[A-Za-z]", ctl_lines[j])) {
        cln <- gsub("\r", "", ctl_lines[j], fixed = TRUE)
        if (!grepl("^\\s*;", cln) && !grepl("^\\s*$", cln) &&
            !grepl("SAME", cln, ignore.case = TRUE)) {
          code_j <- sub(";.*$", "", cln)
          if (grepl("[0-9]", code_j)) {
            eta_counter <- eta_counter + 1L
            name <- .nm_extract_param_name(cln)
            if (!is.na(name) && nchar(name) > 0L)
              mapping[[name]] <- sprintf("ETA(%d)", eta_counter)
          }
        }
        j <- j + 1L
      }
      i <- j
      next
    }
    i <- i + 1L
  }
  mapping
}
# Build ERR mapping from $SIGMA
.nm_build_err_mapping <- function(ctl_lines) {
  .nm_build_mapping(ctl_lines, "^\\s*\\$SIGMA", "ERR")
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
      lhs_seen <- character(0)
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

    # Detect LHS: word before "=" (ignoring spaces)
    lhs_m <- regmatches(code_part,
                        regexpr("^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*=",
                                code_part))
    if (length(lhs_m) > 0L) {
      lhs_name <- trimws(sub("\\s*=\\s*$", "", trimws(lhs_m)))
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
#' @param ctl_file Path to the Wings extended control stream (\code{.ctl}).
#' @param out_file Output path for the converted control stream. If
#'   \code{NULL} (default), the output file is named
#'   \code{<run>_converted.ctl} in the same directory as \code{ctl_file}.
#' @param overwrite Logical. If \code{TRUE}, overwrite \code{out_file} if it
#'   already exists (default \code{FALSE}).
#' @param verbose Logical. If \code{TRUE} (default), print a summary of
#'   substitutions made.
#'
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
#'   nmfe_path    = "C:/nm76g64/run/nmfe76.bat",
#'   ctl_files    = out,
#'   analysis_dir = "C:/analyses"
#' )
#' }
#'
#' @seealso \code{\link{nm_run}}
#'
#' @export
nm_convert_ctl <- function(ctl_file,
                           out_file  = NULL,
                           overwrite = FALSE,
                           verbose   = TRUE) {

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
  theta_map <- .nm_build_mapping(ctl_lines, "^\\s*\\$THETA", "THETA")
  eta_map   <- .nm_build_eta_mapping(ctl_lines)
  err_map   <- .nm_build_err_mapping(ctl_lines)
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
