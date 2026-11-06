## 05_dupcheck.R —— 同一药物的两个 GDSC2 筛选批次分别建模，检验可重复性
suppressPackageStartupMessages(library(glmnet))
setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))
set.seed(42)

G <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")
tcga <- readRDS("tcga_log2tpm.rds"); expr_t <- tcga$expr; clin <- tcga$clin
gdsc_e <- as.matrix(readRDS(file.path(G, "GDSC2_Expr (RMA Normalized and Log Transformed).rds")))
gdsc_r <- readRDS(file.path(G, "GDSC2_Res.rds"))
genes <- intersect(rownames(gdsc_e), rownames(expr_t))
gx <- t(gdsc_e[genes, ]); tx <- t(expr_t[genes, ])

run <- function(col) {
  y <- gdsc_r[, col]
  ok <- !is.na(y) & rownames(gdsc_r) %in% rownames(gx)
  Xtr <- gx[rownames(gdsc_r)[ok], ]
  yv <- as.numeric(y[ok])
  set.seed(42)
  cv <- cv.glmnet(Xtr, yv, alpha = 0, nfolds = 5)
  p <- as.numeric(predict(cv, s = "lambda.min", newx = tx))
  hi <- p[clin$RiskGroup == "High"]; lo <- p[clin$RiskGroup == "Low"]
  sp <- sqrt(((length(hi)-1)*var(hi) + (length(lo)-1)*var(lo))/(length(hi)+length(lo)-2))
  data.frame(col = col, n = length(yv), med_lo = round(median(lo),3), med_hi = round(median(hi),3),
             diff = round(mean(hi)-mean(lo),3),
             p = signif(wilcox.test(hi, lo)$p.value,3),
             d = round((mean(hi)-mean(lo))/sp,3),
             cv_rmse = round(sqrt(min(cv$cvm)),3))
}

pairs <- list(c("Fulvestrant_1200","Fulvestrant_1816"),
              c("Docetaxel_1819","Docetaxel_1007"))
out <- do.call(rbind, lapply(pairs, function(pp) do.call(rbind, lapply(pp, run))))
print(out, row.names = FALSE)
write.csv(out, "out/dupcheck.csv", row.names = FALSE)
