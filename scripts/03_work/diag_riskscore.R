# 诊断：risk score 与增殖 / 受体 / 免疫 / 纯度的相关性
suppressMessages({
  library(SummarizedExperiment)
  library(org.Hs.eg.db)
})

load("../raw/TCGA/brca.gdc_2024.rda")
expr <- assay(expquery2)
colnames(expr) <- substr(gsub("\\.", "-", colnames(expr)), 1, 12)
ens <- gsub("\\..*", "", rownames(expr))
sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(expr) <- sym
expr <- expr[!is.na(rownames(expr)), ]
expr <- expr[!duplicated(rownames(expr)), ]

fd <- read.csv("Table_S7_Final_Clinical_Data.csv", stringsAsFactors = FALSE)
s <- intersect(colnames(expr), fd$Patient_ID)
cat("匹配样本:", length(s), " 风险表:", nrow(fd), "\n")

genes <- c("ESR1","PGR","ERBB2","MKI67","AURKA","AURKB","BIRC5","CCNB1",
           "CD8A","CD68","FOXP3","PTPRC",
           "PADI4","ELANE","MPO","HMGB1","S100A8","S100A9","CXCL8","TGFB1",
           "PLG","CD79A","MMP7","ANGPT1")
genes <- intersect(genes, rownames(expr))
m <- t(expr[genes, s])
df <- data.frame(Patient_ID = s, m, check.names = FALSE)
df <- merge(df, fd[, c("Patient_ID","RiskScore","Age","Stage","OS")], by = "Patient_ID")

rs <- df$RiskScore
cat("\n=== risk score vs 各基因 (Spearman rho) ===\n")
res <- data.frame()
for (g in genes) {
  ct <- suppressWarnings(cor.test(rs, df[[g]], method = "spearman"))
  res <- rbind(res, data.frame(gene = g, rho = round(ct$estimate, 3), p = signif(ct$p.value, 3)))
}
res <- res[order(-abs(res$rho)), ]
print(res, row.names = FALSE)

cat("\n=== risk score vs 临床变量 ===\n")
cat("Age  rho =", round(cor(rs, df$Age, method = "spearman", use = "complete.obs"), 3), "\n")
st <- df$Stage[!is.na(df$Stage)]
rs2 <- rs[!is.na(df$Stage)]
cat("Stage rho =", round(cor(rs2, as.numeric(factor(st)), method = "spearman"), 3), "\n")
cat("Stage 分组均值:\n"); print(tapply(rs2, st, mean))

cat("\n=== 4 个 hub 基因之间的相关 ===\n")
print(round(cor(df[, c("PLG","CD79A","MMP7","ANGPT1")], method = "spearman"), 3))

cat("\n=== 技术混杂：risk score vs 文库大小 ===\n")
lib <- colSums(expr[, s])
cat("library size rho =", round(cor(rs, lib, method = "spearman"), 3), "\n")

cat("\n=== risk score 四分位 vs ESR1/MKI67 ===\n")
q <- cut(rs, quantile(rs, 0:4/4), include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4"))
print(round(sapply(df[, c("ESR1","MKI67","ERBB2","PADI4","ELANE","CD8A")], function(x) tapply(x, q, median)), 2))
