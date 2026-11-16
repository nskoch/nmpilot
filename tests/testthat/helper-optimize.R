# Emulates the derivatives table of run1 (NONMEM unavailable here): analytic
# 1-cpt oral model, derivatives by finite differences at the EBEs of run1.
.emul_deriv <- function(run_dir, data_file, out) {
  ext <- nmtools:::.nm_vpc_read_ext_final(file.path(run_dir, "run1.ext"))
  th <- as.numeric(ext[1, grep("^THETA", names(ext))])
  phi <- nmtools:::.nm_w_read_phi(file.path(run_dir, "run1.phi"))
  d <- read.csv(data_file, na.strings = "."); names(d)[1] <- "ID"
  d <- d[d$EVID == 0, ]
  conc <- function(eta, t, dose = 100) {
    CL <- th[1]*exp(eta[1]); V <- th[2]*exp(eta[2]); KA <- th[3]*exp(eta[3]); k <- CL/V
    dose*KA/(V*(KA-k))*(exp(-k*t)-exp(-KA*t))
  }
  rows <- lapply(split(d, d$ID), function(s) {
    e0 <- as.numeric(phi[phi$ID == s$ID[1], c("ETA(1)","ETA(2)","ETA(3)")])
    f <- conc(e0, s$TIME)
    G <- sapply(1:3, function(k) { h <- 1e-6; e <- e0; e[k] <- e[k]+h; (conc(e, s$TIME)-f)/h })
    W <- sqrt((th[4]*f)^2 + th[5]^2)
    D <- th[4]^2 * f * G / W
    data.frame(ID = s$ID, TIME = s$TIME, MDV = 0, DV = s$DV, CIPREDI = f,
               ETA1 = e0[1], ETA2 = e0[2], ETA3 = e0[3],
               G_ETA1 = G[,1], G_ETA2 = G[,2], G_ETA3 = G[,3], H_EPS1 = W,
               D_EPS1_ETA1 = D[,1], D_EPS1_ETA2 = D[,2], D_EPS1_ETA3 = D[,3])
  })
  x <- do.call(rbind, rows)
  con <- file(out, "w")
  writeLines("TABLE NO.  1", con)
  writeLines(paste(" ", paste(names(x), collapse = " ")), con)
  write.table(format(x, scientific = TRUE, digits = 9), con, row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(con); invisible(out)
}
