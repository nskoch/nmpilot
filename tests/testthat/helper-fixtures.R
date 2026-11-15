# Helpers that write small synthetic NONMEM files in a temporary directory.
# Fixtures are created at test time (no binary or NONMEM output stored in the
# repository).

write_ctl <- function(...) {
  path <- tempfile(fileext = ".ctl")
  writeLines(c(...), path)
  path
}

# Parameter table row for one PK parameter.
pt_row <- function(x, param) {
  pt <- x$parameter_table
  pt[pt$PK_parameter == param, , drop = FALSE]
}

# Role of THETA(index) for a given parameter in the full mapping.
theta_role <- function(x, index, param) {
  m <- x$mapping_detail_all
  unique(m$role[m$element_type == "THETA" & m$index == index &
                  m$pk_parameter == param])
}

write_ext <- function(path, header, rows) {
  writeLines(c(
    paste("TABLE NO.     1: First Order Conditional Estimation with Interaction:",
          "Goal Function=MINIMUM VALUE OF OBJECTIVE FUNCTION: Problem=1",
          "Subproblem=0 Superproblem1=0 Iteration1=0 Superproblem2=0 Iteration2=0"),
    header, rows), path)
}

# Small NONMEM-style data set (2 subjects, 1 dose and 3 observations each).
write_data_csv <- function() {
  f <- tempfile(fileext = ".csv")
  writeLines(c(
    "ID,TIME,AMT,DV,EVID,MDV,CMT,WT",
    "1,0,100,.,1,1,1,70", "1,1,.,5.2,0,0,1,70", "1,4,.,3.1,0,0,1,70",
    "1,12,.,0.8,0,0,1,70",
    "2,0,100,.,1,1,1,85", "2,1,.,4.1,0,0,1,85", "2,4,.,2.6,0,0,1,85",
    "2,12,.,0.05,0,0,1,85"), f)
  f
}

# NONMEM sdtab with the usual diagnostic columns; `shift` is added to CWRES.
write_sdtab <- function(path = tempfile("sdtab"), shift = 0) {
  rows <- rbind(
    c(1, 1, 5.2, 5.0, 5.1, 0.2, 0.1, 0.15, 0, 5.0, 0),
    c(1, 4, 3.1, 3.3, 3.2, -0.3, -0.1, -0.2, 0, 3.3, 0),
    c(1, 12, 0.8, 0.9, 0.85, -0.1, -0.05, 0.1, 0, 0.9, 0),
    c(2, 1, 4.1, 4.5, 4.2, -0.4, -0.2, -0.3, 0, 4.5, 0),
    c(2, 4, 2.6, 2.5, 2.6, 0.1, 0, 0.05, 0, 2.5, 0),
    c(2, 12, 0.05, 0.1, 0.06, -0.5, -0.1, -0.4, 0, 0.1, 0))
  rows[, 6] <- rows[, 6] + shift
  hdr <- c("ID", "TIME", "DV", "PRED", "IPRED", "CWRES", "IWRES", "NPDE",
           "MDV", "CPRED", "EVID")
  writeLines(c("TABLE NO.  1", paste(sprintf("%-12s", hdr), collapse = ""),
               apply(rows, 1, function(r) paste(sprintf("%12.4E", r), collapse = ""))),
             path)
  path
}

# Fixed-width numeric row as written in a NONMEM .ext file.
ext_row <- function(iteration, values) {
  paste0(sprintf("%13s", format(iteration, scientific = FALSE)),
         paste(sprintf("%13.5E", values), collapse = ""))
}
