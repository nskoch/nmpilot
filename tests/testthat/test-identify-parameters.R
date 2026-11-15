test_that("a standard control stream is mapped", {
  f <- write_ctl(
    "$INPUT ID TIME AMT DV WT",
    "$THETA (0,5) ; CL", " (0,50) ; V", " (0,0.2) ; PROP",
    "$OMEGA 0.1 ; CL", " 0.1 ; V",
    "$SIGMA 1 FIX",
    "$PK",
    "  TVCL = THETA(1)",
    "  CL = TVCL*EXP(ETA(1))",
    "  V  = THETA(2)*EXP(ETA(2))",
    "$ERROR",
    "  Y = F + F*THETA(3)*EPS(1)")
  x <- nm_identify_parameters(f)
  expect_s3_class(x, "nm_parameter_map")
  expect_equal(pt_row(x, "CL")$THETA_typical, "THETA(1)")
  expect_equal(pt_row(x, "CL")$ETA_IIV, "ETA(1)")
  expect_equal(pt_row(x, "V")$THETA_typical, "THETA(2)")
  expect_equal(x$residual_error$theta, 3L)
  expect_equal(nrow(x$unmapped_theta), 0L)
  expect_equal(nrow(x$unmapped_omega), 0L)
  expect_equal(nrow(x$label_substitutions), 0L)
})

test_that("ETA( inside THETA( is not read as an ETA", {
  f <- write_ctl(
    "$THETA (0,5)", "$OMEGA 0.1",
    "$PK", "  CL = THETA(1)*EXP(ETA(1))")
  x <- nm_identify_parameters(f)
  m <- x$mapping_detail_all
  expect_setequal(m$index[m$element_type == "ETA"], 1L)
})

test_that("power, exponential, linear and saturable covariate effects", {
  f <- write_ctl(
    "$INPUT ID TIME AMT DV WT AGE GENO CRES",
    "$THETA (0,5)", " 0.75", " (0,50)", " 0.01", " (0,1)", " 0.01", " (0,2)",
    "$OMEGA 0.1 0.1",
    "$PK",
    "  TVCL = THETA(1)*(WT/70)**THETA(2)*EXP(THETA(4)*(AGE-40))",
    "  RATIO = CRES/(THETA(7)+CRES)",
    "  CL = TVCL*(1-RATIO)*EXP(ETA(1))",
    "  V  = THETA(3)*(1 + THETA(6)*(GENO-1))*EXP(ETA(2))",
    "  KA = THETA(5)")
  x <- nm_identify_parameters(f)
  expect_equal(theta_role(x, 1, "CL"), "typical_value")
  expect_equal(theta_role(x, 2, "CL"), "covariate_effect")   # power
  expect_equal(theta_role(x, 4, "CL"), "covariate_effect")   # exponential
  expect_equal(theta_role(x, 7, "CL"), "covariate_effect")   # saturable
  expect_equal(theta_role(x, 3, "V"), "typical_value")
  expect_equal(theta_role(x, 6, "V"), "covariate_effect")    # linear, unusual name
})

test_that("conditional covariate effects", {
  f <- write_ctl(
    "$INPUT ID TIME AMT DV SEX",
    "$THETA (0,5)", " (0,1)", " (0,50)", " (0,60)",
    "$OMEGA 0.1",
    "$PK",
    "  TVCL = THETA(1)",
    "  IF (SEX.EQ.1) TVCL = TVCL*THETA(2)",
    "  CL = TVCL*EXP(ETA(1))",
    "  IF (SEX.EQ.1) THEN",
    "    V = THETA(3)",
    "  ELSE",
    "    V = THETA(4)",
    "  ENDIF")
  x <- nm_identify_parameters(f)
  expect_equal(theta_role(x, 1, "CL"), "typical_value")
  expect_equal(theta_role(x, 2, "CL"), "covariate_effect")
  # Only conditional definitions: the first one is the reference.
  expect_equal(theta_role(x, 3, "V"), "typical_value")
  expect_equal(theta_role(x, 4, "V"), "covariate_effect")
})

test_that("a THETA scaling an ETA is an IIV scale", {
  f <- write_ctl(
    "$THETA (0,5)", " (0,1)", "$OMEGA 0.1",
    "$PK", "  CL = THETA(1)*EXP(ETA(1)*THETA(2))")
  x <- nm_identify_parameters(f)
  expect_equal(theta_role(x, 2, "CL"), "iiv_scale")
})

test_that("Wings label syntax is converted before parsing", {
  f <- write_ctl(
    "$INPUT ID TIME AMT DV C_RES",
    "$THETA", " (0.31) FIXED ;POPKA", " (0.001,.1) ;POPCL", " (0.001,10) ;POPV",
    " (0.0001,0.01) ;COV1", " (.01) ;RUV_P",
    "$OMEGA 0 FIXED ; eta_ka",
    ";$OMEGA 0.1 ; eta_old",
    "$OMEGA 0.1 ; eta_v",
    "$SIGMA 1 FIXED ; EPS1",
    "$PK",
    "  KA = POPKA*EXP(eta_ka)",
    "  RATIO = C_RES/(COV1+C_RES)",
    "  CL = POPCL*(1-RATIO)",
    "  V  = POPV*EXP(eta_v)",
    "$ERROR",
    "  Y = F + F*RUV_P*ERR(1)")
  x <- nm_identify_parameters(f)
  subs <- stats::setNames(x$label_substitutions$reference, x$label_substitutions$label)
  expect_equal(unname(subs[c("POPKA", "POPCL", "POPV", "COV1", "eta_ka", "eta_v")]),
               c("THETA(1)", "THETA(2)", "THETA(3)", "THETA(4)", "ETA(1)", "ETA(2)"))
  expect_equal(pt_row(x, "KA")$THETA_typical, "THETA(1)")
  expect_equal(pt_row(x, "V")$ETA_IIV, "ETA(2)")
  expect_equal(theta_role(x, 4, "CL"), "covariate_effect")
  expect_equal(x$residual_error$theta, 5L)
  expect_equal(nrow(x$unmapped_theta), 0L)
})

test_that("$INPUT items are never treated as labels", {
  f <- write_ctl(
    "$INPUT ID TIME AMT DV WT",
    "$THETA (0,5) ; WT", "$OMEGA 0.1",
    "$PK", "  CL = THETA(1)*(WT/70)*EXP(ETA(1))")
  x <- nm_identify_parameters(f)
  expect_equal(nrow(x$label_substitutions), 0L)
})

test_that("$PRED models are split into structural and residual parts", {
  f <- write_ctl(
    "$INPUT ID TIME AMT DV DOSE",
    "$THETA (0,1)", " (0,0.2)", " (0,10)", " (0,0.1)", " (0,0.5)",
    "$OMEGA 0.1 0.1 0.1", "$SIGMA 1 FIX",
    "$PRED",
    "  KA = THETA(1)*EXP(ETA(1))",
    "  KE = THETA(2)*EXP(ETA(2))",
    "  V  = THETA(3)*EXP(ETA(3))",
    "  IPRED = DOSE/V*KA/(KA-KE)*(EXP(-KE*TIME)-EXP(-KA*TIME))",
    "  W = SQRT(THETA(4)**2*IPRED**2 + THETA(5)**2)",
    "  Y = IPRED + W*EPS(1)")
  x <- nm_identify_parameters(f)
  expect_setequal(pt_row(x, "KE")$THETA_typical, "THETA(2)")
  expect_setequal(x$residual_error$theta, c(4L, 5L))
  expect_equal(nrow(x$unmapped_theta), 0L)
})

test_that("derived parameters are kept and depth points to the direct one", {
  f <- write_ctl(
    "$THETA (0,2) ; MTT", " (0,20) ; CL", " (0,60) ; V",
    "$OMEGA 0.1 0.1 0.1",
    "$PK",
    "  MTT = THETA(1)*EXP(ETA(1))",
    "  CL = THETA(2)*EXP(ETA(2))",
    "  V  = THETA(3)*EXP(ETA(3))",
    "  KE = CL/V")
  x <- nm_identify_parameters(f)
  expect_true(all(c("MTT", "CL", "V", "KE") %in% x$parameter_table$PK_parameter))
  m <- x$mapping_detail_all
  d <- m$depth[m$element_type == "THETA" & m$index == 2]
  names(d) <- m$pk_parameter[m$element_type == "THETA" & m$index == 2]
  expect_equal(d[["CL"]], 0)
  expect_equal(d[["KE"]], 1)
})

test_that("nm_identify_directory returns one row per control stream", {
  d <- file.path(tempdir(), "nmtools-identify-dir")
  dir.create(d, showWarnings = FALSE)
  writeLines(c("$THETA (0,5)", "$OMEGA 0.1", "$PK", "  CL = THETA(1)*EXP(ETA(1))"),
             file.path(d, "a.ctl"))
  writeLines(c("$THETA (0,50)", "$OMEGA 0.1", "$PK", "  V = THETA(1)*EXP(ETA(1))"),
             file.path(d, "b.ctl"))
  res <- nm_identify_directory(d, recursive = FALSE)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2L)
  expect_setequal(res$PK_parameter, c("CL", "V"))
})

test_that("a THETA under a condition whose other branch is a constant is a covariate effect", {
  f <- write_ctl(
    "$INPUT ID TIME AMT DV CLCR CT01",
    "$THETA (0,4)", " (0,20)", " (0,0.7)", " (0,0.1)",
    "$OMEGA 0.4 0.4",
    "$PK",
    "IF (CT01.EQ.1) THEN",
    "  FALL = THETA(3)",
    "ELSE",
    "  FALL = 1",
    "ENDIF",
    "  CL = ((THETA(1)*FALL)+(THETA(4)*(CLCR/89)))*(EXP(ETA(1)))",
    "  V1 = THETA(2)*FALL*(EXP(ETA(2)))")
  x <- nm_identify_parameters(f, likely_pk_only = FALSE)
  expect_equal(theta_role(x, 1, "CL"), "typical_value")
  expect_equal(theta_role(x, 3, "CL"), "covariate_effect")   # FALL = THETA(3) vs 1
  expect_equal(theta_role(x, 4, "CL"), "covariate_effect")   # additive THETA(4)*CLCR/89
  expect_equal(theta_role(x, 2, "V1"), "typical_value")
})

test_that("a THETA multiplying a covariate-dependent helper stays a covariate effect", {
  f <- write_ctl(
    "$INPUT ID TIME AMT DV C_RES",
    "$THETA (0,0.1)", " (0.057) FIX", " (0.9) FIX", "$OMEGA 0.1",
    "$PK",
    "  RATIO = C_RES/(THETA(2)+C_RES)",
    "  CL = THETA(1)*(1-THETA(3)*RATIO)*EXP(ETA(1))")
  x <- nm_identify_parameters(f, likely_pk_only = FALSE)
  expect_equal(theta_role(x, 1, "CL"), "typical_value")
  expect_equal(theta_role(x, 3, "CL"), "covariate_effect")
})
