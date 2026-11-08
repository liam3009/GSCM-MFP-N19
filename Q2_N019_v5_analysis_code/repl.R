source("prep.R")
# original model (Alkandi et al., 2025): second-order GSCM and lean, CSR; business performance = MP + FP items
m <- paste("GD =~ GD1 + GD2 + GD3 + GD4", "GP =~ GP1 + GP2 + GP3 + GP4", "IEM =~ IEM1 + IEM2 + IEM3 + IEM4",
  "LQ =~ LQ1 + LQ2 + LQ3", "JIT =~ JIT1 + JIT2 + JIT3 + JIT4", "TPM =~ TPM1 + TPM2 + TPM3 + TPM4",
  "CSR =~ CSR1 + CSR2 + CSR3 + CSR4", "BP =~ MP1 + MP2 + MP3 + MP4 + FP1 + FP2 + FP3",
  "GSCM =~ GD + GP + IEM", "LEAN =~ LQ + JIT + TPM", "LEAN ~ GSCM", "CSR ~ GSCM", "BP ~ GSCM + LEAN + CSR", sep = "\n")
for (nm in c("all345", "consent342")) {
  d <- if (nm == "all345") d0 else full
  f <- sem(m, data = d[codes], estimator = "ML"); ss <- standardizedSolution(f)
  r <- ss[ss$op == "~", c("lhs", "rhs", "est.std", "pvalue")]; cat("\n", nm, "n =", nrow(d), "\n"); print(r, row.names = FALSE, digits = 3) }
cat("\nconsent table:\n"); print(table(d0$consent))
