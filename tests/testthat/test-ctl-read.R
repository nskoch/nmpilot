# Shared readers of $INPUT, $THETA, $OMEGA, $SIGMA (action 11c)

test_that("$INPUT: DROP alone, ID=DROP, DROP=ID and aliases", {
  l <- c("$PROBLEM x",
         "$INPUT C ID TIME DROP AMT=DOSE DV=CONC WT=DROP SKIP=SEX ; comment",
         "       CMT",
         "$DATA d.csv")
  it <- nmpilot:::.nm_ctl_input_items(l)
  expect_equal(it, c("C", "ID", "TIME", "AMT", "DOSE", "DV", "CONC", "CMT"))
  expect_false(any(c("DROP", "WT", "SEX", "SKIP") %in% it))
})

test_that("$THETA: several per line, FIX inside or outside, xN, NAME=, INF", {
  l <- c("$THETA (0,1,10) 2 FIX ; CL",
         "  (0.5 FIX) ; KA",
         "  (0 1 5)   ; blanks as separators",
         "  (-INF,3,INF) ; V",
         "  (2.0)x3",
         "$THETA CL2=(0, 4) NUMBERPOINTS=3 ; named",
         "$THETAP (1 FIX) ; prior, not a THETA",
         "$OMEGA 0.1")
  th <- nmpilot:::.nm_ctl_theta_parse(l)
  expect_equal(nrow(th), 9L)
  expect_equal(th$initial, c(1, 2, 0.5, 1, 3, 2, 2, 2, 4))
  expect_equal(th$lower[c(1, 4, 5, 9)], c(0, 0, NA, 0))
  expect_equal(th$upper[c(1, 4, 5)], c(10, 5, NA))
  expect_equal(th$fixed, c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(th$n_on_line[1:2], c(2L, 2L))
  expect_true(th$record_line[1]); expect_false(th$record_line[3])
  expect_equal(th$comment[1], "CL")
})

test_that("identify: THETA label only on the first THETA of a line, raw without record name", {
  l <- c("$THETA (0,1,10) 2 ; CL", "$OMEGA 0.1")
  th <- nmpilot:::.nm_ip_parse_theta(l)
  expect_equal(th$label, c("CL", NA))
  expect_equal(th$raw[1], "(0,1,10) 2 ; CL")
})

test_that("$OMEGA: scalars on one line, BLOCK FIX, labels from the diagonal line", {
  l <- c("$OMEGA 0.1 0.2 ; CL",
         "$OMEGA BLOCK(3) FIX",
         "  0.38          ; V2",
         " -0.2 0.45      ; V4",
         "  0  0.5 1.8    ; V3",
         "; comment-only line, not a label",
         "$SIGMA 1 FIX")
  om <- nmpilot:::.nm_ctl_random_parse(l, "OMEGA")$eta
  expect_equal(nrow(om), 5L)
  expect_equal(om$omega_initial, c(0.1, 0.2, 0.38, 0.45, 1.8))
  expect_equal(om$omega_label, c("CL", NA, "V2", "V4", "V3"))
  expect_equal(om$omega_fixed, c(FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(om$omega_structure, c("DIAGONAL", "DIAGONAL", "BLOCK", "BLOCK", "BLOCK"))
  el <- nmpilot:::.nm_ctl_random_parse(l, "OMEGA")$elements
  expect_equal(nrow(el), 2L + 6L)
})

test_that("$OMEGA: one value per line in a BLOCK gives labels of the variances", {
  l <- c("$OMEGA BLOCK(2)", "0.15 ; [P] CL", "0.01 ; [F] cov", "0.15 ; [P] V")
  om <- nmpilot:::.nm_ctl_random_parse(l, "OMEGA")$eta
  expect_equal(om$omega_label, c("[P] CL", "[P] V"))
})

test_that("$OMEGA: BLOCK(1) SAME, BLOCK SAME without size, SAME(k), DIAGONAL(n)", {
  l <- c("$OMEGA .25 ; ETACL",
         "$OMEGA BLOCK(1) .25 ; IOV1",
         "$OMEGA BLOCK(1) SAME ; IOV2",
         "$OMEGA .01 ; ETAD1",
         "$OMEGA BLOCK(2) 0.1 0.01 0.1",
         "$OMEGA BLOCK SAME",
         "$OMEGA BLOCK(1) 0.2 ; OCC",
         "$OMEGA BLOCK(1) SAME(3)",
         "$OMEGA DIAGONAL(2) 0.3 0.4 FIX")
  om <- nmpilot:::.nm_ctl_random_parse(l, "OMEGA")$eta
  expect_equal(nrow(om), 1L + 1L + 1L + 1L + 2L + 2L + 1L + 3L + 2L)
  expect_equal(om$omega_label[1:4], c("ETACL", "IOV1", "IOV2", "ETAD1"))
  expect_equal(om$omega_initial[4], 0.01)
  expect_true(all(is.na(om$omega_initial[7:8])))
  expect_equal(tail(om$omega_fixed, 2), c(TRUE, TRUE))
  expect_equal(nmpilot:::.nm_ctl_n_random(l, "OMEGA"), nrow(om))
})

test_that("$SIGMA is read the same way", {
  l <- c("$OMEGA 0.1", "$SIGMA 0.04 ; prop", "  0.1 FIX ; add")
  sg <- nmpilot:::.nm_ctl_random_parse(l, "SIGMA")$eta
  expect_equal(sg$omega_initial, c(0.04, 0.1))
  expect_equal(sg$omega_fixed, c(FALSE, TRUE))
  expect_equal(sg$omega_label, c("prop", "add"))
})

test_that("nm_convert_ctl: no substitution of a variable assigned after IF(...) or in an earlier record", {
  ctl <- c("$PROBLEM wfn",
           "$INPUT ID TIME AMT DV PAT",
           "$DATA d.csv",
           "$SUBROUTINES ADVAN13 TOL=6",
           "$MODEL COMP=(C1)",
           "$PK",
           "IF(PAT.EQ.1) SCAL = POPS",
           "IF(PAT.EQ.0) SCAL = 0",
           "MTIM(1) = POPM",
           "$DES",
           "DADT(1) = -SCAL*A(1)",
           "$ERROR",
           "Y = F + EPS(1)",
           "$THETA (0,1) ; POPS",
           "$THETA 2 ; POPM",
           "$THETA 3 ; SCAL",
           "$THETA 4 ; MTIM",
           "$OMEGA 0.1",
           "$SIGMA 0.1")
  f <- tempfile(fileext = ".ctl"); writeLines(ctl, f)
  out <- suppressMessages(nm_convert_ctl(f, output_file = tempfile(fileext = ".ctl"), verbose = FALSE))
  txt <- readLines(out)
  expect_true(any(grepl("IF\\(PAT.EQ.1\\) SCAL = THETA\\(1\\)", txt)))
  expect_true(any(grepl("MTIM\\(1\\) = THETA\\(2\\)", txt)))
  expect_true(any(grepl("DADT\\(1\\) = -SCAL\\*A\\(1\\)", txt)))
})
