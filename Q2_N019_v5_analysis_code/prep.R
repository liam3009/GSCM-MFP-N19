# ---- prep.R: settings, packages, data, coding, screening (sourced by every analysis script) ----
source("settings.R")
suppressPackageStartupMessages({ library(lavaan); library(semTools); library(sandwich); library(lmtest) })
dir.create("output", showWarnings = FALSE)

raw <- read.csv(DATA_FILE, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
stopifnot(ncol(raw) == 44)
codes <- c(paste0("GD", 1:4), paste0("GP", 1:4), paste0("IEM", 1:4), paste0("MP", 1:4), paste0("FP", 1:3),
           paste0("LQ", 1:3), paste0("JIT", 1:4), paste0("TPM", 1:4), paste0("CSR", 1:4))
d0 <- raw
names(d0)[1:10]  <- c("timestamp", "consent", "age", "edu", "gender", "industry", "tot_exp", "org_exp", "position", "size")
names(d0)[11:44] <- codes
d0[codes] <- lapply(d0[codes], as.numeric)

dims <- list(GD = codes[1:4], GP = codes[5:8], IEM = codes[9:12], MP = codes[13:16], FP = codes[17:19],
             LQ = codes[20:22], JIT = codes[23:26], TPM = codes[27:30], CSR = codes[31:34])
dims$GSCM <- c(dims$GD, dims$GP, dims$IEM); dims$LEAN <- c(dims$LQ, dims$JIT, dims$TPM); dims$PERF <- c(dims$MP, dims$FP)
for (k in names(dims)) d0[[k]] <- rowMeans(d0[dims[[k]]])

# screening (uses the original industry labels)
longstring  <- function(x) max(rle(as.numeric(x))$lengths)
d0$consent_ok <- d0$consent == "Agree"
d0$sd_items   <- apply(d0[codes], 1, sd)
d0$longstring <- apply(d0[codes], 1, longstring)
d0$careless   <- d0$sd_items == 0 | d0$longstring >= LONGSTR
d0$dup        <- duplicated(d0[c("age", "edu", "gender", "industry", "tot_exp", "org_exp", "position", "size", codes)])
mah           <- mahalanobis(d0[codes], colMeans(d0[codes]), cov(d0[codes]))
d0$mah_out    <- pchisq(mah, df = length(codes), lower.tail = FALSE) < .001

# controls: industries with < 20 consenting respondents merged into "Other"
d0$industry_raw <- d0$industry
small <- names(which(table(d0$industry[d0$consent_ok]) < 20))
d0$industry[d0$industry %in% small] <- "Other"
d0$size     <- relevel(factor(d0$size), ref = "More than 249")
d0$industry <- relevel(factor(d0$industry), ref = "Other")
d0$position <- relevel(factor(d0$position), ref = "Lower Manger")

full  <- subset(d0, consent_ok)
clean <- subset(d0, consent_ok & !careless & !dup)

# shared helpers and model strings
z_df <- function(d) { v <- c("GSCM", "LEAN", "CSR", "MP", "FP", "PERF"); d[paste0("z", v)] <- scale(d[v]); d }
add_dummies <- function(d) { mm <- model.matrix(~ size + industry + position, d)[, -1, drop = FALSE]
  colnames(mm) <- make.names(colnames(mm)); cbind(d, mm) }
ctrl_names <- setdiff(colnames(add_dummies(full)), colnames(full))
m_first  <- paste(sapply(names(dims)[1:9], function(k) paste(k, "=~", paste(dims[[k]], collapse = " + "))), collapse = "\n")
fitcfa   <- function(model, d) cfa(model, data = d[codes], estimator = "MLR")
