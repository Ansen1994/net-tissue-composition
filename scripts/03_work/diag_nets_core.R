# 诊断 v2：用 z-score 均值替代 GSVA，避开 ID 匹配问题
suppressMessages({
  library(SummarizedExperiment); library(org.Hs.eg.db); library(survival)
})

tab <- read.csv("Table_S2_Target_NETs_Matrix.csv", row.names = 1, check.names = FALSE)
cand <- rownames(tab); cand <- cand[cand != "" & !is.na(cand)]

load("../raw/TCGA/brca.gdc_2024.rda")
expr <- assay(expquery2)
colnames(expr) <- substr(gsub("\\.", "-", colnames(expr)), 1, 12)
ens <- gsub("\\..*", "", rownames(expr))
sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(expr) <- sym
expr <- expr[!is.na(rownames(expr)), ]
expr <- expr[!duplicated(rownames(expr)), ]
cat("表达矩阵基因数:", nrow(expr), "| 例子:", paste(head(rownames(expr), 5), collapse=", "), "\n")

fd <- read.csv("Table_S7_Final_Clinical_Data.csv", stringsAsFactors = FALSE)
s <- intersect(colnames(expr), fd$Patient_ID)
E <- expr[, s, drop = FALSE]
rs <- fd$RiskScore[match(s, fd$Patient_ID)]
cat("样本:", length(s), "\n\n")

zscore_mean <- function(genes, mat) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 2) return(NULL)
  m <- mat[g, , drop = FALSE]
  z <- t(scale(t(m)))
  colMeans(z, na.rm = TRUE)
}

sets <- list(
  `NETosis machinery (8)` = c("PADI4","ELANE","MPO","CTSG","PRTN3","AZU1","DNASE1L3","GSDMD"),
  `Neutrophil granule (10)` = c("ELANE","MPO","PRTN3","CTSG","AZU1","LTF","CAMP","BPI","DEFA3","DEFA4"),
  `Neutrophil markers (8)` = c("FCGR3B","CSF3R","CXCR2","FPR1","S100A8","S100A9","S100A12","ITGAM"),
  `4 hub genes` = c("PLG","CD79A","MMP7","ANGPT1"),
  `All candidates` = cand
)

cat("=== 各基因集 z-score 均值 vs risk score (Spearman) ===\n")
sc <- list()
for (nm in names(sets)) {
  v <- zscore_mean(sets[[nm]], E)
  if (is.null(v)) { cat(sprintf("%-24s 匹配基因不足\n", nm)); next }
  sc[[nm]] <- v
  ct <- suppressWarnings(cor.test(rs, v, method = "spearman"))
  cat(sprintf("%-24s 匹配 %2d 个  rho=%+.3f  p=%.2e\n", nm, length(intersect(sets[[nm]], rownames(E))), ct$estimate, ct$p.value))
}

cat("\n=== Cox 单变量（每 SD 的 HR）===\n")
d <- data.frame(OS = fd$OS[match(s, fd$Patient_ID)], time = fd$OS.time[match(s, fd$Patient_ID)])
for (nm in names(sc)) {
  d$v <- as.numeric(sc[[nm]])
  f <- try(coxph(Surv(time, OS) ~ v, data = d), silent = TRUE)
  if (inherits(f, "try-error")) next
  sm <- summary(f)$coef
  cat(sprintf("%-24s HR=%.3f  p=%.2e\n", nm, exp(sm[1,1]), sm[1,5]))
}

cat("\n=== 联合模型：risk score + NETosis machinery ===\n")
d$rs <- rs; d$net <- as.numeric(sc[[`NETosis machinery (8)`]])
print(round(summary(coxph(Surv(time, OS) ~ rs + net, data = d))$coef, 4))
