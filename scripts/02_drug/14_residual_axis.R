## 14_residual_axis.R
## 决定性检验：扣除共同轴后，各药自身信号是否仍然存在
## 做法：pred ~ 群体均值(共同轴) 取残差，对残差再检验组间差异
## 若原始显著但残差不显著 → 该结论由共同轴驱动，非药物特异
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(org.Hs.eg.db); library(glmnet)
})
G <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")
setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))
set.seed(42)

load(file.path(G, "brca.gdc_2024.rda"))
cnt <- assay(expquery2, "unstranded"); tpm <- assay(expquery2, "tpm_unstrand")
sym_of <- function(m) {
  s <- mapIds(org.Hs.eg.db, keys = gsub("\\..*", "", rownames(m)), column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
  k <- !is.na(s); m <- m[k, ]; s <- s[k]
  o <- order(rowMeans(m), decreasing = TRUE); m <- m[o, ]; s <- s[o]; d <- duplicated(s)
  m <- m[!d, ]; rownames(m) <- s[!d]; m
}
cnt <- sym_of(cnt); tpm <- sym_of(tpm)
colnames(cnt) <- substr(gsub("\\.", "-", colnames(cnt)), 1, 12); colnames(tpm) <- colnames(cnt)
rm(expquery2); gc(verbose = FALSE)
clin <- read.csv(file.path(G, "final_data.csv"), stringsAsFactors = FALSE)
clin$RiskGroup <- factor(clin$RiskGroup, levels = c("Low","High"))
cm <- intersect(colnames(cnt), clin$Patient_ID)
cnt <- cnt[, cm]; tpm <- tpm[, cm]; clin <- clin[match(cm, clin$Patient_ID), ]

gdsc_e <- as.matrix(readRDS(file.path(G, "GDSC2_Expr (RMA Normalized and Log Transformed).rds")))
gdsc_r <- readRDS(file.path(G, "GDSC2_Res.rds"))
genes <- intersect(rownames(gdsc_e), rownames(cnt))
gx <- t(gdsc_e[genes, , drop = FALSE]); rm(gdsc_e); gc(verbose = FALSE)
gd_q <- apply(gx, 2, function(v) quantile(v, probs = seq(0,1,0.01), names = FALSE))
qn_map <- function(m) {
  q <- t(apply(m, 1, function(v) rank(v, ties.method="average")/(length(v)+1)))
  out <- sapply(seq_len(ncol(q)), function(j)
    approx(seq(0,1,length.out=nrow(gd_q)), gd_q[,j], xout=q[,j], rule=2)$y)
  dimnames(out) <- dimnames(m); out
}
rk <- function(m) t(apply(m, 1, function(v) rank(v, ties.method="average")/(length(v)+1)))
zs <- function(m) { z <- scale(m); z[is.na(z)] <- 0; as.matrix(z) }

cnt_t <- t(cnt[genes, , drop=FALSE]); tpm_t <- t(tpm[genes, , drop=FALSE])
rm(cnt, tpm); gc(verbose = FALSE)
TESTS <- list(counts = cnt_t, log2tpm = log2(tpm_t+1), zscore = zs(log2(tpm_t+1)),
              qn = qn_map(tpm_t), rank = rk(tpm_t))
TRAINS <- list(counts = gx, log2tpm = gx, zscore = zs(gx), qn = gx, rank = rk(gx))
rm(cnt_t, tpm_t); gc(verbose = FALSE)

## 6 药：原稿 4 药 + 紫杉醇、他莫昔芬作共同轴参照
panel <- data.frame(
  drug = c("Tamoxifen","Fulvestrant","Epirubicin","5-Fluorouracil","Docetaxel","Paclitaxel"),
  col  = c("Tamoxifen_1199","Fulvestrant_1200","Epirubicin_1511","5-Fluorouracil_1073",
           "Docetaxel_1819","Paclitaxel_1080"), stringsAsFactors = FALSE)
core <- c("Epirubicin","5-Fluorouracil","Docetaxel","Fulvestrant")

cohen_d <- function(a,b){n1<-length(a);n2<-length(b)
  sp<-sqrt(((n1-1)*var(a)+(n2-1)*var(b))/(n1+n2-2)); (mean(a)-mean(b))/sp}

out <- list()
for (md in names(TESTS)) {
  cat("\n##########", md, "##########\n")
  Xte <- TESTS[[md]][, genes, drop=FALSE]; Xtr <- TRAINS[[md]][, genes, drop=FALSE]
  P <- matrix(NA, nrow(clin), 0)
  for (i in seq_len(nrow(panel))) {
    y <- gdsc_r[, panel$col[i]]; ok <- !is.na(y)
    set.seed(42)
    fit <- cv.glmnet(Xtr[rownames(gdsc_r)[ok], , drop=FALSE], as.numeric(y[ok]),
                     alpha=0, nfolds=5, standardize=TRUE)
    pr <- as.numeric(predict(fit, s="lambda.min", newx=Xte))
    P <- cbind(P, pr); colnames(P)[ncol(P)] <- panel$drug[i]
  }
  avg_all <- rowMeans(scale(P))                   # 全部 6 药均值（有自身循环）
  for (d in core) {
    pr <- P[, d]
    ## 留一法：共同轴只用其余 5 药构建，避免扣掉被检验药自身的信号
    others <- setdiff(panel$drug, d)
    avg_loo <- rowMeans(scale(P[, others, drop = FALSE]))
    r  <- resid(lm(pr ~ avg_loo))
    hi <- pr[clin$RiskGroup=="High"]; lo <- pr[clin$RiskGroup=="Low"]
    rh <- r[clin$RiskGroup=="High"];  rl <- r[clin$RiskGroup=="Low"]
    out[[length(out)+1]] <- data.frame(
      mode = md, drug = d,
      d_raw = round(cohen_d(hi,lo),3), p_raw = signif(wilcox.test(hi,lo)$p.value,3),
      d_resid = round(cohen_d(rh,rl),3), p_resid = signif(wilcox.test(rh,rl)$p.value,3),
      rho_with_axis = round(cor(pr, avg_loo, method="spearman"),3),
      stringsAsFactors = FALSE)
    cat(sprintf("  %-15s 原始 d=%6.3f p=%.3g | 留一残差 d=%6.3f p=%.3g | 与该轴相关 %.2f\n",
                d, cohen_d(hi,lo), wilcox.test(hi,lo)$p.value,
                cohen_d(rh,rl), wilcox.test(rh,rl)$p.value,
                cor(pr, avg_loo, method="spearman")))
  }
  saveRDS(P, paste0("out/P_", md, ".rds"))
}
res <- do.call(rbind, out)
res$survives <- ifelse(res$p_resid < 0.05, "YES", "NO")
write.csv(res, "out/residual_axis.csv", row.names=FALSE)
cat("\n\n===== 扣除共同轴后的存活情况 =====\n")
options(width=200); print(res, row.names=FALSE)
cat("\n各药原稿显著方案中，扣除共同轴后仍显著的比例:\n")
for (d in core) {
  s <- res[res$drug==d, ]
  n <- sum(s$p_raw < 0.05); k <- sum(s$p_raw < 0.05 & s$p_resid < 0.05)
  cat(sprintf("  %-15s %d/%d\n", d, k, n))
}
sink("out/sessionInfo_resid.txt"); print(sessionInfo()); sink()
