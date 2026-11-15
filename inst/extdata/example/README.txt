nmtools example project - SIMULATED DATA, FICTITIOUS SUBJECTS
=============================================================

example_data.csv was simulated in R (script sim_example_data.R, seed 20261115).
It does not come from any patient or study.

Design: 60 subjects, single oral dose of 100 mg at time 0, 9 samples
(0.5, 1, 2, 3, 4, 6, 8, 12, 24 h). Covariates: body weight WT (kg) and a
binary genotype GENO (0/1).

True model used for the simulation (one-compartment, first-order absorption):
  CL = 5 * (WT/70)^0.75 * 1.6^GENO * exp(eta_CL)   omega2 = 0.09
  V  = 50 * (WT/70) * exp(eta_V)                   omega2 = 0.04
  KA = 1.2 * exp(eta_KA)                           omega2 = 0.16
  residual error: proportional 15 % + additive 0.02 mg/L

Models:
  run1.ctl  BASE  - no covariate
  run2.ctl  FINAL - WT on CL (estimated exponent) and V (fixed 1), GENO on CL

To reproduce the NONMEM outputs:
  nm_run(nmfe_path = "<path>/nmfe76.bat", ctl_files = c("run1.ctl", "run2.ctl"),
         analysis_dir = "<this folder>")
