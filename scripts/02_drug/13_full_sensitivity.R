## 13_full_sensitivity.R
## 完整敏感性矩阵：12 药（跨 6 类）× 5 种输入预处理
## 每个组合报告：Cohen's d、p、域外比例、CV R2
## 每种处理再算共同轴判据：12药均值 vs 风险评分、首PC方差占比、药间相关中位数
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(org.Hs.eg.db); library(glmnet)
})
G <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")
setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))
set.seed(42)

## ---------- 表达矩阵 ----------
cat("== 载入 TCGA ==\n")
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
clin$RiskGroup <- factor(clin$RiskGroup, levels = c("Low", "High"))
cm <- intersect(colnames(cnt), clin$Patient_ID)
cnt <- cnt[, cm]; tpm <- tpm[, cm]; clin <- clin[match(cm, clin$Patient_ID), ]
cat("TCGA:", dim(cnt), " High/Low:", paste(table(clin$RiskGroup), collapse = "/"), "\n")

gdsc_e <- as.matrix(readRDS(file.path(G, "GDSC2_Expr (RMA Normalized and Log Transformed).rds")))
gdsc_r <- readRDS(file.path(G, "GDSC2_Res.rds"))
genes <- intersect(rownames(gdsc_e), rownames(cnt))
gx <- t(gdsc_e[genes, , drop = FALSE])       # 细胞系 x 基因
rm(gdsc_e); gc(verbose = FALSE)
cat("共同基因:", length(genes), " GDSC:", dim(gx), "\n")

## ---------- 5 种预处理 ----------
gd_q <- apply(gx, 2, function(v) quantile(v, probs = seq(0, 1, 0.01), names = FALSE))

qn_map <- function(m) {
  q <- t(apply(m, 1, function(v) rank(v, ties.method = "average") / (length(v) + 1)))
  out <- sapply(seq_len(ncol(q)), function(j)
    approx(seq(0, 1, length.out = nrow(gd_q)), gd_q[, j], xout = q[, j], rule = 2)$y)
  dimnames(out) <- dimnames(m); out
}
rk <- function(m) t(apply(m, 1, function(v) rank(v, ties.method = "average") / (length(v) + 1)))
zs <- function(m) { z <- scale(m); z[is.na(z)] <- 0; as.matrix(z) }

cat("== 构建 5 种预处理矩阵 ==\n")
cnt_t <- t(cnt[genes, , drop = FALSE])
tpm_t <- t(tpm[genes, , drop = FALSE])
rm(cnt, tpm); gc(verbose = FALSE)

TESTS <- list(
  counts  = cnt_t,
  log2tpm = log2(tpm_t + 1),
  zscore  = zs(log2(tpm_t + 1)),
  qn      = qn_map(tpm_t),
  rank    = rk(tpm_t)
)
rm(cnt_t, tpm_t); gc(verbose = FALSE)

TRAINS <- list(
  counts  = gx,
  log2tpm = gx,
  zscore  = zs(gx),
  qn      = gx,
  rank    = rk(gx)
)
for (nm in names(TESTS)) cat(sprintf("  %-8s test %.2f~%.2f  train %.2f~%.2f\n", nm,
        min(TESTS[[nm]]), max(TESTS[[nm]]), min(TRAINS[[nm]]), max(TRAINS[[nm]])))

## ---------- 药物 ----------
panel <- data.frame(
  drug = c("Tamoxifen","Fulvestrant","Epirubicin","5-Fluorouracil","Docetaxel","Paclitaxel",
           "Lapatinib","Sapitinib","Alpelisib","Pictilisib","Olaparib","Palbociclib"),
  col = c("Tamoxifen_1199","Fulvestrant_1200","Epirubicin_1511","5-Fluorouracil_1073",
          "Docetaxel_1819","Paclitaxel_1080","Lapatinib_1558","Sapitinib_1549",
          "Alpelisib_1560","Pictilisib_1058","Olaparib_1017","Palbociclib_1054"),
  cat = c("Endocrine","Endocrine","Chemo","Chemo","Chemo","Chemo",
          "HER2","HER2","PI3K","PI3K","PARP","CDK4-6"),
  stringsAsFactors = FALSE)

cohen_d <- function(a, b) { n1 <- length(a); n2 <- length(b)
  sp <- sqrt(((n1-1)*var(a) + (n2-1)*var(b))/(n1+n2-2)); (mean(a)-mean(b))/sp }

rows <- list(); mode_summary <- list()

for (md in names(TESTS)) {
  cat("\n########## 方案:", md, "##########\n")
  Xte <- TESTS[[md]][, genes, drop = FALSE]
  Xtr_all <- TRAINS[[md]][, genes, drop = FALSE]
  P <- matrix(NA, nrow(clin), 0)
  for (i in seq_len(nrow(panel))) {
    cn <- panel$col[i]
    y <- gdsc_r[, cn]; ok <- !is.na(y)
    yv <- as.numeric(y[ok]); cells <- rownames(gdsc_r)[ok]
    set.seed(42)
    fit <- cv.glmnet(Xtr_all[cells, , drop = FALSE], yv, alpha = 0, nfolds = 5, standardize = TRUE)
    pr <- as.numeric(predict(fit, s = "lambda.min", newx = Xte))
    hi <- pr[clin$RiskGroup == "High"]; lo <- pr[clin$RiskGroup == "Low"]
    p <- wilcox.test(hi, lo)$p.value
    rows[[length(rows)+1]] <- data.frame(
      drug = panel$drug[i], category = panel$cat[i], mode = md, n_cell = length(yv),
      outside_frac = round(mean(pr < min(yv) | pr > max(yv)), 3),
      cv_r2 = round(1 - min(fit$cvm)/var(yv), 3),
      diff = round(mean(hi) - mean(lo), 4),
      cohen_d = round(cohen_d(hi, lo), 3),
      p = signif(p, 3), sig = ifelse(p < 0.05, "YES", "-"),
      stringsAsFactors = FALSE)
    P <- cbind(P, pr); colnames(P)[ncol(P)] <- panel$drug[i]
    cat(sprintf("  %-15s d=%6.3f p=%.3g 域外%3.0f%%\n", panel$drug[i], cohen_d(hi,lo), p,
                mean(pr < min(yv) | pr > max(yv))*100))
  }
  ## 共同轴判据
  cmx <- cor(P, method = "spearman"); off <- cmx[upper.tri(cmx)]
  ev <- svd(scale(P))$d^2; ev <- ev / sum(ev)
  avg <- rowMeans(scale(P))
  rho <- suppressWarnings(cor(avg, clin$RiskScore, method = "spearman"))
  p_rho <- suppressWarnings(cor.test(avg, clin$RiskScore, method = "spearman")$p.value)
  mode_summary[[md]] <- data.frame(
    mode = md, n_sig = sum(sapply(rows, function(r) r$mode == md && r$sig == "YES")),
    rho_avg_risk = round(rho, 3), p_rho = signif(p_rho, 3),
    pc1_pct = round(ev[1]*100, 1),
    inter_drug_med = round(median(off), 3),
    mean_outside = round(mean(sapply(rows, function(r) if (r$mode == md) r$outside_frac else NA), na.rm = TRUE), 3),
    stringsAsFactors = FALSE)
  cat(sprintf("  >> 共同轴: rho=%+.3f (p=%.2g)  首PC=%.1f%%  药间相关=%.3f\n", rho, p_rho, ev[1]*100, median(off)))
}

res <- do.call(rbind, rows)
ms <- do.call(rbind, mode_summary)
write.csv(res, "out/sensitivity_matrix.csv", row.names = FALSE)
write.csv(ms, "out/sensitivity_by_mode.csv", row.names = FALSE)

## ---------- 每药稳健性汇总 ----------
wide <- reshape(res[, c("drug","category","mode","cohen_d","p","sig")],
                idvar = c("drug","category"), timevar = "mode", direction = "wide")
agg <- do.call(rbind, lapply(split(res, res$drug), function(x) {
  data.frame(drug = x$drug[1], category = x$category[1],
             n_sig = sum(x$sig == "YES"),
             dir_pos = sum(x$cohen_d > 0), dir_neg = sum(x$cohen_d < 0),
             d_min = min(x$cohen_d), d_max = max(x$cohen_d),
             d_med = median(x$cohen_d), stringsAsFactors = FALSE)
}))
agg$direction_consistent <- ifelse(agg$dir_pos == 5 | agg$dir_neg == 5, "YES", "NO")
agg$verdict <- ifelse(agg$n_sig >= 4 & agg$direction_consistent == "YES", "robust",
               ifelse(agg$n_sig >= 3 & agg$direction_consistent == "YES", "moderate", "unstable"))
agg <- agg[order(-agg$n_sig, -abs(agg$d_med)), ]
write.csv(agg, "out/sensitivity_verdict.csv", row.names = FALSE)

cat("\n\n===== 方案层面判据 =====\n"); print(ms, row.names = FALSE)
cat("\n===== 药物稳健性排序 =====\n"); print(agg, row.names = FALSE)

## ---------- 热图 ----------
suppressPackageStartupMessages(library(ggplot2))
d <- res
d$mode <- factor(d$mode, levels = c("counts","log2tpm","zscore","qn","rank"))
d$drug <- factor(d$drug, levels = rev(agg$drug))
d$star <- ifelse(d$p < 0.05, "*", "")
p <- ggplot(d, aes(mode, drug, fill = cohen_d)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = star), size = 6, vjust = 0.78, color = "black") +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                       midpoint = 0, name = "Cohen's d") +
  labs(x = "输入预处理", y = NULL,
       title = "药敏效应量对输入预处理的敏感性",
       subtitle = "* p < 0.05 (Wilcoxon)。蓝=高危敏感，红=高危耐药") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 0),
        plot.title = element_text(face = "bold"))
ggsave("fig/sensitivity_heatmap.pdf", p, width = 7.5, height = 5.5)
ggsave("fig/sensitivity_heatmap.png", p, width = 7.5, height = 5.5, dpi = 200)
cat("\n图已保存 fig/sensitivity_heatmap.pdf / .png\n")
sink("out/sessionInfo_sens.txt"); print(sessionInfo()); sink()
