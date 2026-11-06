## 02_drug_panel.R
## 目的：扩展药物面板，按乳腺癌术后辅助治疗类别分组，统一用 log2(TPM+1) 输入重跑药敏预测
## 关键修正：输入尺度 = log2(TPM+1)（原代码误用 unstranded 原始 counts）

suppressPackageStartupMessages({
  library(glmnet)
  library(parallel)
})

G   <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")
OUT <- Sys.getenv("WB_OUTPUT_DIR", unset = "out")
setwd(OUT)
set.seed(42)

## ============ 1. 药物面板定义 ============
panel <- data.frame(
  drug = c("Tamoxifen", "Fulvestrant",
           "Epirubicin", "Mitoxantrone", "Paclitaxel", "Docetaxel",
           "Cyclophosphamide", "5-Fluorouracil", "Gemcitabine", "Cisplatin",
           "Vinorelbine", "Vinblastine",
           "Lapatinib", "Afatinib", "Sapitinib",
           "Alpelisib", "Buparlisib", "Pictilisib", "Taselisib",
           "Ipatasertib", "MK-2206", "Rapamycin", "AZD2014",
           "Olaparib", "Talazoparib", "Niraparib",
           "Palbociclib", "Ribociclib", "Dinaciclib"),
  col = c("Tamoxifen_1199", "Fulvestrant_1200",
          "Epirubicin_1511", "Mitoxantrone_1810", "Paclitaxel_1080", "Docetaxel_1819",
          "Cyclophosphamide_1512", "5-Fluorouracil_1073", "Gemcitabine_1190", "Cisplatin_1005",
          "Vinorelbine_2048", "Vinblastine_1004",
          "Lapatinib_1558", "Afatinib_1032", "Sapitinib_1549",
          "Alpelisib_1560", "Buparlisib_1873", "Pictilisib_1058", "Taselisib_1561",
          "Ipatasertib_1924", "MK-2206_1053", "Rapamycin_1084", "AZD2014_1441",
          "Olaparib_1017", "Talazoparib_1259", "Niraparib_1177",
          "Palbociclib_1054", "Ribociclib_1632", "Dinaciclib_1180"),
  cat = c("Endocrine", "Endocrine",
          "Chemo", "Chemo", "Chemo", "Chemo",
          "Chemo", "Chemo", "Chemo", "Chemo",
          "Chemo", "Chemo",
          "HER2", "HER2", "HER2",
          "PI3K-AKT-mTOR", "PI3K-AKT-mTOR", "PI3K-AKT-mTOR", "PI3K-AKT-mTOR",
          "PI3K-AKT-mTOR", "PI3K-AKT-mTOR", "PI3K-AKT-mTOR", "PI3K-AKT-mTOR",
          "PARP", "PARP", "PARP",
          "CDK4-6", "CDK4-6", "CDK4-6-control"),
  stringsAsFactors = FALSE
)
## 重复筛选的同一药物（不同化合物编号）作为技术重复，取两筛平均
dup_extra <- list(Docetaxel = "Docetaxel_1007", Fulvestrant = "Fulvestrant_1816")

## ============ 2. 数据 ============
cat("== 载入数据 ==\n")
tcga <- readRDS("tcga_log2tpm.rds")
expr_t <- tcga$expr; clin <- tcga$clin
cat("TCGA:", dim(expr_t), "\n")

gdsc_e <- readRDS(file.path(G, "GDSC2_Expr (RMA Normalized and Log Transformed).rds"))
gdsc_r <- readRDS(file.path(G, "GDSC2_Res.rds"))
gdsc_e <- as.matrix(gdsc_e)
cat("GDSC expr:", dim(gdsc_e), " resid:", dim(gdsc_r), "\n")

genes <- intersect(rownames(gdsc_e), rownames(expr_t))
cat("共同基因:", length(genes), "\n")

gx <- t(gdsc_e[genes, , drop = FALSE])
## 判断 GDSC 列名是否为细胞系（与 Res 列名对应）
cat("GDSC expr 列名示例:", paste(head(colnames(gdsc_e), 3), collapse = ", "), "\n")
cat("GDSC res  列名示例:", paste(head(colnames(gdsc_r), 3), collapse = ", "), "\n")
cat("细胞系交集:", length(intersect(colnames(gdsc_e), rownames(gdsc_r))), "\n")

## ============ 3. 标准化方式 ============
## mode = "raw"    : 各自 log2 原值直接对接（主分析，与原文一致的建模方式）
## mode = "zscore" : 各数据集内逐基因 z 标准化（消除平台间整体偏移，稳健性检验）
prep_matrix <- function(m, mode) {
  if (mode == "raw") return(m)
  z <- scale(m)
  z[is.na(z)] <- 0
  as.matrix(z)
}

tx <- prep_matrix(t(expr_t[genes, , drop = FALSE]), "raw")

## ============ 4. 建模与预测 ============
fit_one <- function(drug, colname, mode = "raw") {
  y <- gdsc_r[, colname]
  ok <- !is.na(y) & rownames(gdsc_r) %in% rownames(gx)
  yv <- as.numeric(y[ok])
  cells <- rownames(gdsc_r)[ok]

  Xtr <- if (mode == "raw") gx[cells, genes, drop = FALSE] else prep_matrix(gx[cells, genes, drop = FALSE], "zscore")
  Xte <- if (mode == "raw") tx else prep_matrix(tx, "zscore")
  ## 保证测试集列顺序与训练集一致
  Xte <- Xte[, genes, drop = FALSE]

  set.seed(42)
  cvfit <- cv.glmnet(Xtr, yv, alpha = 0, nfolds = 5, standardize = TRUE)
  pred  <- as.numeric(predict(cvfit, s = "lambda.min", newx = Xte))
  list(drug = drug, pred = pred, n_cell = length(yv),
       train_lo = min(yv), train_hi = max(yv),
       train_med = median(yv), lambda = cvfit$lambda.min,
       cv_rmse = sqrt(min(cvfit$cvm)))
}

## ============ 5. 统计 ============
cohen_d <- function(a, b) {           # a = High, b = Low
  n1 <- length(a); n2 <- length(b)
  sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(b)) / (n1 + n2 - 2))
  (mean(a) - mean(b)) / sp
}
boot_diff <- function(a, b, B = 1000) {
  set.seed(42)
  d <- replicate(B, mean(sample(a, replace = TRUE)) - mean(sample(b, replace = TRUE)))
  quantile(d, c(0.025, 0.975))
}

res_list <- list(); pred_mat <- matrix(NA, nrow = nrow(clin), ncol = 0)

for (i in seq_len(nrow(panel))) {
  d <- panel$drug[i]; cn <- panel$col[i]
  cat(sprintf("[%2d/%d] %-16s %s\n", i, nrow(panel), d, cn))

  f1 <- fit_one(d, cn, "raw")
  ## 技术重复：若存在第二个化合物编号，取两筛预测值平均（同一 ln(IC50) 量纲）
  preds <- list(f1$pred); dup_r <- NA_real_
  if (!is.null(dup_extra[[d]])) {
    f2 <- fit_one(d, dup_extra[[d]], "raw")
    preds[[2]] <- f2$pred
    dup_r <- cor(f1$pred, f2$pred)
    cat("      重复筛选相关性 r =", round(dup_r, 3), "\n")
  }
  pred <- if (length(preds) == 1) preds[[1]] else rowMeans(do.call(cbind, preds))

  grp <- clin$RiskGroup
  hi <- pred[grp == "High"]; lo <- pred[grp == "Low"]
  wt <- wilcox.test(hi, lo)
  ci <- boot_diff(hi, lo)
  rho <- suppressWarnings(cor.test(pred, clin$RiskScore, method = "spearman"))

  ## 稳健性：z-score 标准化方式下方向是否一致
  fz <- fit_one(d, cn, "zscore")
  dz <- cohen_d(fz$pred[grp == "High"], fz$pred[grp == "Low"])

  res_list[[i]] <- data.frame(
    drug = d, category = panel$cat[i], gdsc_col = cn,
    n_cell_lines = f1$n_cell,
    dup_screen_r = round(dup_r, 3),
    train_lo = round(f1$train_lo, 2), train_hi = round(f1$train_hi, 2),
    pred_min = round(min(pred), 2), pred_max = round(max(pred), 2),
    med_Low = round(median(lo), 3), med_High = round(median(hi), 3),
    diff = round(mean(hi) - mean(lo), 3),
    ci_lo = round(ci[1], 3), ci_hi = round(ci[2], 3),
    p_wilcox = wt$p.value,
    cohen_d = round(cohen_d(hi, lo), 3),
    rho_risk = round(rho$estimate, 3), p_rho = rho$p.value,
    d_zscore = round(dz, 3),
    cv_rmse = round(f1$cv_rmse, 3),
    stringsAsFactors = FALSE
  )
  pred_mat <- cbind(pred_mat, scale(pred)[, 1])
  colnames(pred_mat)[ncol(pred_mat)] <- d
}

res <- do.call(rbind, res_list)
res$fdr      <- p.adjust(res$p_wilcox, method = "BH")
res$direction <- ifelse(res$diff > 0, "High-risk resistant", "High-risk sensitive")
res$robust   <- ifelse(sign(res$cohen_d) == sign(res$d_zscore), "consistent", "DISCORDANT")

## 可靠性标记：训练细胞系过少 / 体外无活性的前药
res$reliability <- ifelse(res$n_cell_lines < 100, "low (n<100 cell lines)", "ok")
res$reliability[res$drug == "Cyclophosphamide"] <- "not evaluable in vitro (prodrug needs CYP activation)"
res$reliability[res$drug == "Ribociclib"] <- "low (n=46 cell lines only)"

## 类别内 FDR
res$fdr_within <- ave(res$p_wilcox, res$category, FUN = function(p) p.adjust(p, "BH"))

res <- res[order(res$category, res$p_wilcox), ]
write.csv(res, "out/drug_results.csv", row.names = FALSE)
rownames(pred_mat) <- clin$Patient_ID
write.csv(pred_mat, "out/pred_matrix.csv")

cat("\n== 结果 ==\n")
print(res[, c("drug", "category", "diff", "p_wilcox", "fdr", "cohen_d", "rho_risk", "direction", "robust")],
      row.names = FALSE)

cat("\n== 输入尺度检查 ==\n")
cat("TCGA log2TPM 范围:", round(range(tx), 2), " GDSC RMA 范围:", round(range(gx[, genes]), 2), "\n")
