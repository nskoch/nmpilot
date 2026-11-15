# Simulated example data for nmtools (fictitious subjects).
set.seed(20261115)
n <- 60
id <- seq_len(n)
wt <- round(pmin(pmax(rnorm(n, 75, 14), 45), 120), 1)
geno <- rbinom(n, 1, 0.3)
eta <- cbind(rnorm(n, 0, sqrt(0.09)), rnorm(n, 0, sqrt(0.04)), rnorm(n, 0, sqrt(0.16)))
cl <- 5 * (wt / 70)^0.75 * 1.6^geno * exp(eta[, 1])
v  <- 50 * (wt / 70) * exp(eta[, 2])
ka <- 1.2 * exp(eta[, 3])
dose <- 100
times <- c(0.5, 1, 2, 3, 4, 6, 8, 12, 24)
rows <- list()
for (i in id) {
  k <- cl[i] / v[i]
  ipred <- dose * ka[i] / (v[i] * (ka[i] - k)) * (exp(-k * times) - exp(-ka[i] * times))
  dv <- ipred * (1 + rnorm(length(times), 0, 0.15)) + rnorm(length(times), 0, 0.02)
  dv <- pmax(dv, 0.001)
  rows[[i]] <- rbind(
    data.frame(ID = i, TIME = 0, AMT = dose, DV = ".", EVID = 1, MDV = 1, CMT = 1, WT = wt[i], GENO = geno[i]),
    data.frame(ID = i, TIME = times, AMT = ".", DV = sprintf("%.4f", dv), EVID = 0, MDV = 0, CMT = 2, WT = wt[i], GENO = geno[i]))
}
d <- do.call(rbind, rows)
names(d)[1] <- "#ID"
write.csv(d, "example_data.csv", row.names = FALSE, quote = FALSE)
cat("records:", nrow(d), " subjects:", n, " GENO=1:", sum(geno), " WT range:", range(wt), "\n")
