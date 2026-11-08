source("prep.R"); set.seed(SEED); suppressPackageStartupMessages(library(MASS))
z <- z_df(full); f <- "~ zGSCM + zLEAN + zCSR + size + industry + position"
g <- function(x) { a <- coef(rlm(as.formula(paste("zMP", f)), x, method = "MM", maxit = 100)); b <- coef(rlm(as.formula(paste("zFP", f)), x, method = "MM", maxit = 100))
  c(H1 = a[["zGSCM"]] - b[["zGSCM"]], H2 = b[["zLEAN"]] - b[["zGSCM"]]) }
e <- g(z); bt <- t(replicate(2000, g(z[sample(nrow(z), replace = TRUE), ])))
print(round(rbind(est = e, apply(bt, 2, quantile, c(.025, .975))), 3))
