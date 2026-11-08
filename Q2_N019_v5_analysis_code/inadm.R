source("prep.R"); set.seed(SEED + 7)
library(parallel)
m_rq1 <- paste("GD =~ GD1 + GD2 + GD3 + GD4", "GP =~ GP1 + GP2 + GP3 + GP4", "IEM =~ IEM1 + IEM2 + IEM3 + IEM4",
  "LQ =~ LQ1 + LQ2 + LQ3", "JIT =~ JIT1 + JIT2 + JIT3 + JIT4", "TPM =~ TPM1 + TPM2 + TPM3 + TPM4",
  "MP =~ MP1 + MP2 + MP3 + MP4", "FP =~ FP1 + FP2 + FP3", "GSCM =~ GD + GP + IEM", "LEAN =~ LQ + JIT + TPM",
  paste("MP ~ GSCM + LEAN +", paste(ctrl_names, collapse = " + ")), paste("FP ~ GSCM + LEAN +", paste(ctrl_names, collapse = " + ")), "MP ~~ FP", sep = "\n")
dd <- add_dummies(full)[c(codes, ctrl_names)]; n <- nrow(dd); B <- as.integer(Sys.getenv("BB", "1000"))
idx <- replicate(B, sample(n, replace = TRUE), simplify = FALSE)
one <- function(i) { f <- tryCatch(suppressWarnings(sem(m_rq1, data = dd[i, ], estimator = "ML")), error = function(e) NULL)
  if (is.null(f) || !lavInspect(f, "converged")) return(c(conv = 0, negvar = NA, npd = NA, std_gt1 = NA))
  pe <- parameterEstimates(f); v <- pe$est[pe$op == "~~" & pe$lhs == pe$rhs]
  ev <- tryCatch(min(eigen(lavInspect(f, "cov.lv"), only.values = TRUE)$values), error = function(e) -1)
  ss <- standardizedSolution(f, se = FALSE); s1 <- any(abs(ss$est.std[ss$op %in% c("~", "=~")]) > 1, na.rm = TRUE)
  c(conv = 1, negvar = as.numeric(any(v < 0)), npd = as.numeric(ev <= 0), std_gt1 = as.numeric(s1)) }
r <- do.call(rbind, mclapply(idx, one, mc.cores = 2)); r <- as.data.frame(r)
cat("draws:", B, "\nnon-converged:", sum(r$conv == 0), "\n")
cv <- r[r$conv == 1, ]; bad <- cv$negvar == 1 | cv$npd == 1
cat("converged:", nrow(cv), " improper:", sum(bad), sprintf("(%.1f%%)", 100 * mean(c(r$conv == 0) | c(r$conv == 1 & (r$negvar == 1 | r$npd == 1)), na.rm = TRUE)), "\n")
cat("negative variance:", sum(cv$negvar), " non-PD latent cov:", sum(cv$npd), " both:", sum(cv$negvar == 1 & cv$npd == 1), " |std|>1:", sum(cv$std_gt1), "\n")
saveRDS(r, "inadm.rds")
whichneg <- function(i) { f <- suppressWarnings(sem(m_rq1, data = dd[i, ], estimator = "ML")); pe <- parameterEstimates(f)
  v <- pe[pe$op == "~~" & pe$lhs == pe$rhs & pe$est < 0, "lhs"]; if (length(v)) paste(v, collapse = "+") else "none" }
print(sort(table(unlist(mclapply(idx[1:300], whichneg, mc.cores = 2))), decreasing = TRUE))
