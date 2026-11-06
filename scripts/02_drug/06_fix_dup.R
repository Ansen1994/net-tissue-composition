## 06_fix_dup.R —— 重复筛选药物改用原稿批次做主分析，另一批次留作可重复性检验
suppressPackageStartupMessages(library(glmnet))
setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))
set.seed(42)

G <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")
tcga <- readRDS("tcga_log2tpm.rds"); expr_t <- tcga$expr; clin <- tcga$clin
gdsc_e <- as.matrix(readRDS(file.path(G, "GDSC2_Expr (RMA Normalized and Log Transformed).rds")))
gdsc_r <- readRDS(file.path(G, "GDSC2_Res.rds"))
genes <- intersect(rownames(gdsc_e), rownames(expr_t))
gx <- t(gdsc_e[genes, ]); tx <- t(expr_t[genes, ])
txz <- scale(tx); txz[is.na(txz)] <- 0

res  <- read.csv("out/drug_results.csv", stringsAsFactors = FALSE)
pred <- read.csv("out/pred_matrix.csv", row.names = 1, check.names = FALSE)

cohen_d <- function(a, b) {
  sp <- sqrt(((length(a)-1)*var(a) + (length(b)-1)*var(b))/(length(a)+length(b)-2))
  (mean(a)-mean(b))/sp
}
fitp <- function(col, z = FALSE) {
  y <- gdsc_r[, col]; ok <- !is.na(y)
  Xtr <- gx[rownames(gdsc_r)[ok], , drop = FALSE]
  if (z) { Xtr <- scale(Xtr); Xtr[is.na(Xtr)] <- 0 }
  cv <- cv.glmnet(Xtr, as.numeric(y[ok]), alpha = 0, nfolds = 5)
  as.numeric(predict(cv, s = "lambda.min", newx = if (z) txz else tx))
}

prim <- c(Fulvestrant = "Fulvestrant_1200", Docetaxel = "Docetaxel_1819")
for (d in names(prim)) {
  cn <- prim[[d]]
  p  <- fitp(cn); pz <- fitp(cn, TRUE)
  hi <- p[clin$RiskGroup == "High"]; lo <- p[clin$RiskGroup == "Low"]
  set.seed(42)
  bd <- replicate(1000, mean(sample(hi, length(hi), TRUE)) - mean(sample(lo, length(lo), TRUE)))
  i <- which(res$drug == d)
  res$gdsc_col[i]  <- cn
  res$n_cell_lines[i] <- sum(!is.na(gdsc_r[, cn]))
  tr <- gdsc_r[, cn]; res$train_lo[i] <- round(min(tr, na.rm = TRUE), 2)
  res$train_hi[i]  <- round(max(tr, na.rm = TRUE), 2)
  res$pred_min[i]  <- round(min(p), 2); res$pred_max[i] <- round(max(p), 2)
  res$med_Low[i]   <- round(median(lo), 3); res$med_High[i] <- round(median(hi), 3)
  res$diff[i]      <- round(mean(hi) - mean(lo), 3)
  res$ci_lo[i]     <- round(quantile(bd, .025), 3); res$ci_hi[i] <- round(quantile(bd, .975), 3)
  res$p_wilcox[i]  <- wilcox.test(hi, lo)$p.value
  res$cohen_d[i]   <- round(cohen_d(hi, lo), 3)
  r <- suppressWarnings(cor.test(p, clin$RiskScore, method = "spearman"))
  res$rho_risk[i]  <- round(r$estimate, 3); res$p_rho[i] <- r$p.value
  res$d_zscore[i]  <- round(cohen_d(pz[clin$RiskGroup=="High"], pz[clin$RiskGroup=="Low"]), 3)
  pred[, d] <- as.numeric(scale(p))
  cat(sprintf("%-14s %-16s diff=%.3f p=%.3g d=%.3f\n", d, cn, res$diff[i], res$p_wilcox[i], res$cohen_d[i]))
}

res$fdr <- p.adjust(res$p_wilcox, "BH")
res$direction <- ifelse(res$diff > 0, "High-risk resistant", "High-risk sensitive")
res$robust <- ifelse(sign(res$cohen_d) == sign(res$d_zscore), "consistent", "DISCORDANT")
res$fdr_within <- ave(res$p_wilcox, res$category, FUN = function(p) p.adjust(p, "BH"))
res <- res[order(res$category, res$p_wilcox), ]

write.csv(res, "out/drug_results.csv", row.names = FALSE)
write.csv(pred, "out/pred_matrix.csv")
cat("\nFDR<0.05 的药物:\n")
print(res[res$fdr < 0.05, c("drug","category","diff","p_wilcox","fdr","cohen_d")], row.names = FALSE)
