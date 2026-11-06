## 01_prep.R
## 目的：构建 log2(TPM+1) 基因符号矩阵，与 risk 表对齐，保存为 RDS 供后续复用
## 说明：原 figure6.R 使用 assay(expquery2) 即 unstranded 原始 counts，属尺度错误；
##       此处统一改为 tpm_unstrand + log2(x+1)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(org.Hs.eg.db)
})

G   <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")
OUT <- Sys.getenv("WB_OUTPUT_DIR", unset = "out")

setwd(OUT)
set.seed(42)

cat("== 1. 载入 TCGA 数据 ==\n")
load(file.path(G, "brca.gdc_2024.rda"))
cat("assayNames:", paste(assayNames(expquery2), collapse = " | "), "\n")

use_assay <- if ("tpm_unstrand" %in% assayNames(expquery2)) "tpm_unstrand" else "unstranded"
cat("使用 assay:", use_assay, "\n")
mat <- assay(expquery2, use_assay)
cat("原始维度:", dim(mat), "\n")
cat("数值范围:", round(range(mat[, 1:5], na.rm = TRUE), 1), " 中位数:",
    round(median(mat[, 1], na.rm = TRUE), 2), "\n")

colnames(mat) <- substr(gsub("\\.", "-", colnames(mat)), 1, 12)
mat <- mat[, !duplicated(colnames(mat))]

## log2 变换（TPM 已是标准化值，取 log2(x+1)）
mat <- log2(mat + 1)
cat("log2 后范围:", round(range(mat[, 1:5], na.rm = TRUE), 2), "\n")
rm(expquery2); gc(verbose = FALSE)

cat("== 2. Ensembl -> Symbol ==\n")
ens <- gsub("\\..*", "", rownames(mat))
sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
keep <- !is.na(sym)
mat  <- mat[keep, ]
sym  <- sym[keep]

## 同一符号多转录本：保留平均表达最高的那个
ord <- order(rowMeans(mat), decreasing = TRUE)
mat <- mat[ord, ]
sym <- sym[ord]
dup <- duplicated(sym)
mat <- mat[!dup, ]
rownames(mat) <- sym[!dup]
cat("基因数:", nrow(mat), " 样本数:", ncol(mat), "\n")

cat("== 3. 对齐 risk 表 ==\n")
clin <- read.csv(file.path(G, "final_data.csv"), stringsAsFactors = FALSE)
clin$RiskGroup <- factor(clin$RiskGroup, levels = c("Low", "High"))
common <- intersect(colnames(mat), clin$Patient_ID)
cat("交集样本:", length(common), "\n")

mat  <- mat[, common]
clin <- clin[match(common, clin$Patient_ID), ]
stopifnot(all(clin$Patient_ID == colnames(mat)))

cat("RiskGroup 分布:\n"); print(table(clin$RiskGroup))
cat("Stage 分布:\n"); print(table(clin$Stage, useNA = "ifany"))

saveRDS(list(expr = mat, clin = clin), file.path(OUT, "tcga_log2tpm.rds"))
cat("已保存 tcga_log2tpm.rds\n")

sink("out/sessionInfo.txt"); print(sessionInfo()); sink()
cat("== 完成 ==\n")
