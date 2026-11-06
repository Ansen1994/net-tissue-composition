## 16_targeted_residual.R
## 精简补齐：29 药面板中「被写进结论」的 11 个药 × 5 种预处理
## 双共同轴留一法对照：
##   (A) 11 药内部轴：其余 10 药标准化均值
##   (B) 29 药全量轴：已有 pred_matrix*.csv 中其余 28 药标准化均值（仅 log2tpm / qn 可用）
## 目的：判断 29 药版的分类结论（含各阴性类）在扣除共同轴后是否成立
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(org.Hs.eg.db); library(glmnet)
})
G <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")
setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))
set.seed(42)
t0 <- Sys.time()

## ---------- 表达矩阵（与 13_full_sensitivity.R 完全一致） ----------
cat("== Load TCGA ==\n")
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
gx <- t(gdsc_e[genes, , drop = FALSE])
rm(gdsc_e); gc(verbose = FALSE)
cat("Common genes:", length(genes), " GDSC:", dim(gx), "\n")

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

cnt_t <- t(cnt[genes, , drop = FALSE]); tpm_t <- t(tpm[genes, , drop = FALSE])
rm(cnt, tpm); gc(verbose = FALSE)
TESTS <- list(counts = cnt_t, log2tpm = log2(tpm_t + 1),
              zscore = zs(log2(tpm_t + 1)), qn = qn_map(tpm_t), rank = rk(tpm_t))
rm(cnt_t, tpm_t); gc(verbose = FALSE)
TRAINS <- list(counts = gx, log2tpm = gx, zscore = zs(gx), qn = gx, rank = rk(gx))

## ---------- 目标药（覆盖 6 类，含全部被写进结论的药） ----------
panel <- data.frame(
  drug = c("Docetaxel", "5-Fluorouracil", "Fulvestrant", "Epirubicin",
           "Sapitinib", "Lapatinib", "Tamoxifen", "Paclitaxel",
           "Palbociclib", "Olaparib", "Alpelisib"),
  col  = c("Docetaxel_1819", "5-Fluorouracil_1073", "Fulvestrant_1200", "Epirubicin_1511",
           "Sapitinib_1549", "Lapatinib_1558", "Tamoxifen_1199", "Paclitaxel_1080",
           "Palbociclib_1054", "Olaparib_1017", "Alpelisib_1560"),
  cat  = c("Chemo", "Chemo", "Endocrine", "Chemo", "HER2", "HER2", "Endocrine",
           "Chemo", "CDK4-6", "PARP", "PI3K-AKT-mTOR"),
  prior = c("resistant", "resistant", "sensitive", "REJECTED(1/4)",
            "resistant(strongest)", "no-diff", "trend-sensitive", "no-diff",
            "no-diff", "no-diff", "no-diff"),
  stringsAsFactors = FALSE)

## 29 药预测矩阵（外部轴），按患者 ID 对齐
EXT <- list(
  log2tpm = read.csv("out/pred_matrix.csv", row.names = 1, check.names = FALSE),
  qn      = read.csv("out/pred_matrix_qn.csv", row.names = 1, check.names = FALSE)
)
for (nm in names(EXT)) {
  m <- as.matrix(EXT[[nm]])
  EXT[[nm]] <- m[match(clin$Patient_ID, rownames(m)), , drop = FALSE]
  cat("External 29-drug axis [", nm, "]:", dim(EXT[[nm]]),
      " NA:", sum(is.na(EXT[[nm]])), "\n")
}

cohen_d <- function(a, b) { n1 <- length(a); n2 <- length(b)
  sp <- sqrt(((n1-1)*var(a) + (n2-1)*var(b))/(n1+n2-2)); (mean(a)-mean(b))/sp }
axis_from <- function(M, exclude) {
  keep <- setdiff(colnames(M), exclude)
  rowMeans(scale(M[, keep, drop = FALSE]))
}
test_one <- function(pr, ax) {
  r <- resid(lm(pr ~ ax))
  hi <- clin$RiskGroup == "High"; lo <- !hi
  c(d_raw = cohen_d(pr[hi], pr[lo]), p_raw = wilcox.test(pr[hi], pr[lo])$p.value,
    d_res = cohen_d(r[hi], r[lo]), p_res = wilcox.test(r[hi], r[lo])$p.value)
}

out <- list()
for (md in names(TESTS)) {
  cat("\n##########", md, "##########\n")
  Xte <- TESTS[[md]][, genes, drop = FALSE]
  Xtr <- TRAINS[[md]][, genes, drop = FALSE]
  P <- matrix(NA, nrow(clin), 0)
  for (i in seq_len(nrow(panel))) {
    cn <- panel$col[i]
    y <- gdsc_r[, cn]; ok <- !is.na(y)
    yv <- as.numeric(y[ok]); cells <- rownames(gdsc_r)[ok]
    set.seed(42)
    fit <- cv.glmnet(Xtr[cells, , drop = FALSE], yv, alpha = 0, nfolds = 5, standardize = TRUE)
    pr <- as.numeric(predict(fit, s = "lambda.min", newx = Xte))
    P <- cbind(P, pr); colnames(P)[ncol(P)] <- panel$drug[i]
  }
  ## 方案层面的 11 药共同轴质量
  avg11 <- rowMeans(scale(P))
  rho11 <- suppressWarnings(cor(avg11, clin$RiskScore, method = "spearman"))
  cat(sprintf("  11-drug axis vs RiskScore: rho=%.3f\n", rho11))
  if (md %in% names(EXT)) {
    E <- EXT[[md]]
    rho29 <- suppressWarnings(cor(rowMeans(scale(E)), clin$RiskScore, method = "spearman"))
    cat(sprintf("  29-drug axis vs RiskScore: rho=%.3f\n", rho29))
    ## 两轴一致性
    cat(sprintf("  axis-axis Spearman: %.3f\n",
                suppressWarnings(cor(avg11, rowMeans(scale(E)), method = "spearman"))))
  }

  for (d in panel$drug) {
    pr <- P[, d]
    a <- test_one(pr, axis_from(P, d))                       # (A) 11 药内部轴
    b <- if (md %in% names(EXT)) test_one(pr, axis_from(EXT[[md]], d)) else rep(NA_real_, 4)
    out[[length(out)+1]] <- data.frame(
      mode = md, drug = d, category = panel$cat[panel$drug == d],
      prior = panel$prior[panel$drug == d],
      d_raw = round(a[1], 3), p_raw = signif(a[2], 3),
      d_res11 = round(a[3], 3), p_res11 = signif(a[4], 3),
      d_res29 = round(b[3], 3), p_res29 = signif(b[4], 3),
      n_axis11 = ncol(P) - 1, n_axis29 = if (md %in% names(EXT)) ncol(EXT[[md]]) - 1 else NA,
      stringsAsFactors = FALSE)
    cat(sprintf("  %-15s raw d=%6.3f p=%.3g | LOO11 d=%6.3f p=%.3g | LOO29 d=%s p=%s\n",
                d, a[1], a[2], a[3], a[4],
                ifelse(is.na(b[3]), "  NA ", sprintf("%6.3f", b[3])),
                ifelse(is.na(b[4]), " NA ", sprintf("%.3g", b[4]))))
  }
  saveRDS(P, paste0("out/P11_", md, ".rds"))
}
res <- do.call(rbind, out)

## ---------- 逐药最终判定 ----------
verdict <- do.call(rbind, lapply(panel$drug, function(d) {
  s <- res[res$drug == d, ]
  sig <- sum(s$p_res11 < 0.05, na.rm = TRUE)
  pos <- sum(s$d_res11 > 0); neg <- sum(s$d_res11 < 0)
  dir <- if (pos == 5) "all-positive" else if (neg == 5) "all-negative" else "MIXED"
  ## 29 药轴（仅 2 方案）方向对照
  s29 <- s[!is.na(s$p_res29), ]
  agree29 <- if (nrow(s29) == 0) NA else
    all(sign(s29$d_res11) == sign(s29$d_res29))
  cls <- if (sig >= 4 && dir != "MIXED") "ROBUST-SIGNAL" else
         if (sig <= 1 && mean(abs(s$d_res11)) < 0.15) "ROBUST-NULL" else "UNCERTAIN"
  data.frame(drug = d, category = panel$cat[panel$drug == d],
             prior_claim = panel$prior[panel$drug == d],
             n_sig_loo11 = sig, direction = dir,
             median_d = round(median(s$d_res11), 3),
             min_d = round(min(s$d_res11), 3), max_d = round(max(s$d_res11), 3),
             axis29_direction_agree = agree29,
             verdict = cls, stringsAsFactors = FALSE)
}))

write.csv(res, "out/targeted_residual.csv", row.names = FALSE)
write.csv(verdict, "out/targeted_verdict.csv", row.names = FALSE)

cat("\n\n===== FINAL VERDICT (11 target drugs x 5 modes, LOO-11 axis) =====\n")
options(width = 220); print(verdict, row.names = FALSE)
cat("\n===== 29-drug axis cross-check (log2tpm / qn only) =====\n")
print(res[!is.na(res$p_res29),
          c("mode","drug","d_res11","p_res11","d_res29","p_res29")], row.names = FALSE)
cat("\nElapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
cat("sessionInfo saved\n")
writeLines(capture.output(sessionInfo()), "out/sessionInfo_targeted.txt")
