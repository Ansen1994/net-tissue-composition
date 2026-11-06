# 可行性验证 · 队列二：GSE96058 (n = 2,840)
suppressMessages({ library(data.table); library(survival) })
cat("读取表达矩阵...\n")
e <- fread("../raw/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",
           data.table = FALSE)
rownames(e) <- e[[1]]; e <- as.matrix(e[, -1])
cat("矩阵:", dim(e), "\n")

fd <- read.csv("Table_S6_GSE96058_Validation_Result.csv", stringsAsFactors = FALSE)
s <- intersect(colnames(e), fd$title)
cat("匹配样本:", length(s), "\n")

prolif <- c("ANLN","BIRC5","BUB1","CCNB1","CDC20","CDC25B","CEP55","CKS2","EXO1",
            "KPNA2","MELK","MKI67","MYBL2","NUF2","ORC6L","PTTG1","RRM2","TYMS","UBE2C")
immune <- c("CD8A","CD8B","GZMA","GZMB","PRF1","CXCL9","CXCL10","PDCD1","CTLA4",
            "FOXP3","MS4A1","CD79A","CD79B","TNFRSF17","IGHM")

sc <- function(genes, samples) {
  g <- intersect(genes, rownames(e))
  m <- e[g, samples, drop = FALSE]
  z <- t(scale(t(m))); z[is.na(z)] <- 0
  colMeans(z)
}
df <- data.frame(title = s, Prolif = sc(prolif, s), Immune = sc(immune, s))
df <- merge(df, fd[, c("title","RiskScore","OS","OS.time")], by = "title")
df <- df[!is.na(df$OS) & !is.na(df$OS.time) & df$OS.time > 0, ]
cat("纳入:", nrow(df), " 事件:", sum(df$OS == 1), "\n\n")

ct <- suppressWarnings(cor.test(df$RiskScore, df$Prolif, method = "spearman"))
cat(sprintf("risk vs 增殖模块: rho = %.3f, p = %.3g\n", ct$estimate, ct$p.value))
ct2 <- suppressWarnings(cor.test(df$RiskScore, df$Immune, method = "spearman"))
cat(sprintf("risk vs 免疫模块: rho = %.3f, p = %.3g\n", ct2$estimate, ct2$p.value))

# 随机对照
set.seed(42); allg <- rownames(e); null <- numeric(300)
for (i in 1:300) {
  g <- sample(allg, length(prolif)); v <- sc(g, df$title)
  null[i] <- suppressWarnings(cor(v, df$RiskScore, method = "spearman"))
}
cat(sprintf("随机集 rho: 中位 %.3f, 95%% [%.3f, %.3f]; 实际 %.3f, 百分位 %.1f%%\n",
            median(null), quantile(null,.025), quantile(null,.975),
            ct$estimate, mean(null < ct$estimate)*100))

df$ProlifGrp <- ifelse(df$Prolif > median(df$Prolif), "High", "Low")
cat("\n按增殖分层的 risk 组 HR:\n")
for (lv in c("Low","High")) {
  sub <- df[df$ProlifGrp == lv, ]
  sub$RG <- as.numeric(sub$RiskScore > median(sub$RiskScore))
  f <- summary(coxph(Surv(OS.time, OS) ~ RG, data = sub))
  cat(sprintf("  增殖%s (n=%d, 事件=%d): HR = %.2f (%.2f-%.2f), p = %.4f\n",
              lv, nrow(sub), sum(sub$OS), f$conf.int[1], f$conf.int[3], f$conf.int[4], f$coefficients[1,5]))
}
f2 <- summary(coxph(Surv(OS.time, OS) ~ RiskScore + Prolif, data = df))
cat(sprintf("\n多因素: RiskScore HR %.2f (p=%.4g) | Prolif HR %.2f (p=%.4g)\n",
            f2$conf.int[1,1], f2$coefficients[1,5], f2$conf.int[2,1], f2$coefficients[2,5]))
