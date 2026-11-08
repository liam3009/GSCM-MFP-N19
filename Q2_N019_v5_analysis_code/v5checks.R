# ---- v5checks.R: checks added for manuscript v5 (pre-submission review) ----
# M1 discriminant validity of market vs financial performance: CI of the CFA correlation and chi2(cut) tests (Ronkko & Cho, 2022)
# M7 bootstrap: delta-method (Wald) CI for H1/H2 without discarding draws; reasons for inadmissible draws
# M3 Holm-adjusted p for the H1/H2 Williams tests; M4 H1/H2 without the 48 multivariate outliers; M6 RQ2 without CSR
source("prep.R"); set.seed(SEED)
out <- list()
m_rq1 <- paste(
  "GD =~ GD1 + GD2 + GD3 + GD4", "GP =~ GP1 + GP2 + GP3 + GP4", "IEM =~ IEM1 + IEM2 + IEM3 + IEM4",
  "LQ =~ LQ1 + LQ2 + LQ3", "JIT =~ JIT1 + JIT2 + JIT3 + JIT4", "TPM =~ TPM1 + TPM2 + TPM3 + TPM4",
  "MP =~ MP1 + MP2 + MP3 + MP4", "FP =~ FP1 + FP2 + FP3", "GSCM =~ GD + GP + IEM", "LEAN =~ LQ + JIT + TPM",
  paste("MP ~ a*GSCM + c*LEAN +", paste(ctrl_names, collapse = " + ")),
  paste("FP ~ b*GSCM + d*LEAN +", paste(ctrl_names, collapse = " + ")), "MP ~~ FP", "H1 := a - b", "H2 := d - b", sep = "\n")

## M1
m1 <- list()
for (s in c("full", "clean")) {
  d <- get(s)
  f0 <- cfa(m_first, data = d[codes], estimator = "MLR", std.lv = TRUE)
  pe <- parameterEstimates(f0); r <- pe[pe$lhs == "MP" & pe$op == "~~" & pe$rhs == "FP", ]
  t1 <- cfa(paste(m_first, "MP ~~ 1*FP", sep = "\n"), data = d[codes], estimator = "MLR", std.lv = TRUE)
  t9 <- cfa(paste(m_first, "MP ~~ 0.9*FP", sep = "\n"), data = d[codes], estimator = "MLR", std.lv = TRUE)
  x1 <- lavTestLRT(t1, f0, method = "satorra.2000"); x9 <- lavTestLRT(t9, f0, method = "satorra.2000")
  # all other factor correlations: largest upper limits
  cc <- pe[pe$op == "~~" & pe$lhs != pe$rhs, ]; cc <- cc[order(-cc$ci.upper), c("lhs", "rhs", "est", "ci.lower", "ci.upper")]
  m1[[s]] <- data.frame(Sample = s, r_MP_FP = r$est, ci_low = r$ci.lower, ci_high = r$ci.upper,
                        chi_cut1 = x1[2, "Chisq diff"], p_cut1 = x1[2, "Pr(>Chisq)"], chi_cut90 = x9[2, "Chisq diff"], p_cut90 = x9[2, "Pr(>Chisq)"])
  cat("\n", s, "\n"); print(round(m1[[s]][-1], 4)); print(head(cc, 6), row.names = FALSE, digits = 3)
}
out$M1 <- do.call(rbind, m1)

## M7a delta-method CI for H1, H2 (standardized)
m7 <- list()
for (s in c("full", "clean")) {
  dd <- add_dummies(get(s)); f <- sem(m_rq1, data = dd[c(codes, ctrl_names)], estimator = "ML")
  ss <- standardizedSolution(f); q <- ss[ss$op == ":=", c("lhs", "est.std", "se", "ci.lower", "ci.upper", "pvalue")]
  m7[[s]] <- cbind(Sample = s, q); print(m7[[s]], digits = 3)
}
out$M7a <- do.call(rbind, m7)

## M3 Holm for the zero-order tests (p values from Table 7)
lev <- read.csv("X13_levels_of_evidence.csv")
holm <- do.call(rbind, lapply(split(lev, lev$Sample), function(x) data.frame(Sample = x$Sample, p_H1 = x$W_H1_p, p_H2 = x$W_H2_p,
  holm_H1 = p.adjust(c(x$W_H1_p, x$W_H2_p), "holm")[1], holm_H2 = p.adjust(c(x$W_H1_p, x$W_H2_p), "holm")[2])))
out$M3 <- holm; print(holm)

## M4 without multivariate outliers
williams <- function(rjk, rjh, rkh, n) { detR <- 1 - rjk^2 - rjh^2 - rkh^2 + 2 * rjk * rjh * rkh; rbar <- (rjk + rjh) / 2
  t <- (rjk - rjh) * sqrt(((n - 1) * (1 + rkh)) / (2 * ((n - 1) / (n - 3)) * detR + rbar^2 * (1 - rkh)^3)); c(t = t, p = 2 * pt(-abs(t), n - 3)) }
comp <- function(x) { x <- z_df(droplevels(x)); f <- "~ zGSCM + zLEAN + zCSR + size + industry + position"
  a <- coef(lm(as.formula(paste("zMP", f)), x)); b <- coef(lm(as.formula(paste("zFP", f)), x))
  c(H1 = a[["zGSCM"]] - b[["zGSCM"]], H2 = b[["zLEAN"]] - b[["zGSCM"]]) }
no <- subset(full, !mah_out); n <- nrow(no); R <- cor(no[c("GSCM", "LEAN", "MP", "FP")])
w1 <- williams(R["GSCM", "MP"], R["GSCM", "FP"], R["MP", "FP"], n); w2 <- williams(R["LEAN", "FP"], R["GSCM", "FP"], R["GSCM", "LEAN"], n)
ce <- comp(no); cb <- t(replicate(B_COMP, comp(no[sample(n, replace = TRUE), ])))
dd <- add_dummies(no); fo <- sem(m_rq1, data = dd[c(codes, ctrl_names)], estimator = "ML")
so <- standardizedSolution(fo); so <- so[so$op == ":=", ]
fun <- function(x) { s2 <- standardizedSolution(x, se = FALSE); s2$est.std[s2$op == ":="] }
bt <- bootstrapLavaan(fo, R = 1000, FUN = fun, parallel = "multicore", ncpus = 2)
bad <- attr(bt, "nonadmissible"); if (is.null(bad)) bad <- integer(0); ok <- setdiff(seq_len(nrow(bt)), bad)
lci <- apply(bt[ok, ], 2, quantile, c(.025, .975))
out$M4 <- data.frame(n = n, r_GM = R["GSCM", "MP"], r_GF = R["GSCM", "FP"], r_LF = R["LEAN", "FP"], W_H1_t = w1[["t"]], W_H1_p = w1[["p"]], W_H2_t = w2[["t"]], W_H2_p = w2[["p"]],
  comp_H1 = ce[["H1"]], comp_H1_lo = quantile(cb[, "H1"], .025), comp_H1_hi = quantile(cb[, "H1"], .975),
  comp_H2 = ce[["H2"]], comp_H2_lo = quantile(cb[, "H2"], .025), comp_H2_hi = quantile(cb[, "H2"], .975),
  lat_H1 = so$est.std[1], lat_H1_lo = lci[1, 1], lat_H1_hi = lci[2, 1], lat_H2 = so$est.std[2], lat_H2_lo = lci[1, 2], lat_H2_hi = lci[2, 2],
  boot_ok = length(ok), boot_total = nrow(bt))
print(t(round(out$M4, 3)))

## M6 RQ2 without CSR
tost <- function(b, se, dfr, bound) max(pt((b + bound) / se, dfr, lower.tail = FALSE), pt((b - bound) / se, dfr))
out$M6 <- do.call(rbind, lapply(c("full", "clean"), function(s) {
  d <- z_df(get(s)); m <- lm(zPERF ~ zGSCM * zLEAN + size + industry + position, d)
  ct <- coeftest(m, vcov = vcovHC(m, type = "HC3")); b <- ct["zGSCM:zLEAN", 1]; se <- ct["zGSCM:zLEAN", 2]; dfr <- df.residual(m)
  data.frame(Sample = s, beta = b, se = se, p = ct["zGSCM:zLEAN", 4], lo = b - qt(.975, dfr) * se, hi = b + qt(.975, dfr) * se,
             tost10 = tost(b, se, dfr, .10), tost05 = tost(b, se, dfr, .05)) }))
print(out$M6, digits = 3)
saveRDS(out, "v5checks.rds"); cat("\nV5 CHECKS DONE\n")
