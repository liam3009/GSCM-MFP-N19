source("prep.R"); set.seed(SEED)
williams <- function(rjk, rjh, rkh, n) { detR <- 1 - rjk^2 - rjh^2 - rkh^2 + 2 * rjk * rjh * rkh; rbar <- (rjk + rjh) / 2
  t <- (rjk - rjh) * sqrt(((n - 1) * (1 + rkh)) / (2 * ((n - 1) / (n - 3)) * detR + rbar^2 * (1 - rkh)^3)); c(t = t, p = 2 * pt(-abs(t), n - 3)) }
o <- subset(full, mah_out); n_o <- subset(full, !mah_out)
cat("outliers:", nrow(o), " also careless/dup:", sum(o$careless | o$dup), "\n")
cat("in clean sample:", sum(clean$mah_out), "\n")
ms <- function(x) round(c(sapply(x[c("GSCM","LEAN","MP","FP")], mean), sd_items = mean(x$sd_items), diffMF = mean(x$MP - x$FP)), 2)
print(rbind(outliers = ms(o), others = ms(n_o)))
print(round(cor(o[c("GSCM","LEAN","MP","FP")]), 2))
# within the clean sample without outliers
cn <- subset(clean, !mah_out); R <- cor(cn[c("GSCM","LEAN","MP","FP")])
cat("clean & no outliers n =", nrow(cn), "\n"); print(williams(R["GSCM","MP"], R["GSCM","FP"], R["MP","FP"], nrow(cn)))
# rank-based (Spearman) Williams test on the full sample, all cases kept
S <- cor(full[c("GSCM","LEAN","MP","FP")], method = "spearman"); print(round(S, 3))
cat("Spearman Williams H1:\n"); print(williams(S["GSCM","MP"], S["GSCM","FP"], S["MP","FP"], nrow(full)))
cat("Spearman Williams H2:\n"); print(williams(S["LEAN","FP"], S["GSCM","FP"], S["GSCM","LEAN"], nrow(full)))
# robust (MM) regression composite H1 with all cases
suppressPackageStartupMessages(library(MASS))
z <- z_df(full); f <- "~ zGSCM + zLEAN + zCSR + size + industry + position"
a <- coef(rlm(as.formula(paste("zMP", f)), z, method = "MM")); b <- coef(rlm(as.formula(paste("zFP", f)), z, method = "MM"))
cat("MM-robust composite H1:", round(a[["zGSCM"]] - b[["zGSCM"]], 3), " H2:", round(b[["zLEAN"]] - b[["zGSCM"]], 3), "\n")
bt <- replicate(1000, { i <- sample(nrow(z), replace = TRUE); x <- z[i, ]
  a <- coef(rlm(as.formula(paste("zMP", f)), x, method = "MM", maxit = 100)); b <- coef(rlm(as.formula(paste("zFP", f)), x, method = "MM", maxit = 100)); a[["zGSCM"]] - b[["zGSCM"]] })
print(round(quantile(bt, c(.025, .975)), 3))
# how many outliers have MP > FP pattern etc: H1 within outliers
Ro <- cor(o[c("GSCM","MP","FP")]); cat("outliers r(G,MP) r(G,FP):", round(Ro[1,2:3], 3), "\n")
