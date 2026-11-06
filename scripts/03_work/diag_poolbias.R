s1 <- read.csv("Table_S1_All_DEGs.csv", stringsAsFactors = FALSE, check.names = FALSE)
s2 <- read.csv("Table_S2_Target_NETs_Matrix.csv", row.names = 1, check.names = FALSE,
               nrows = 200, header = TRUE)
genes132 <- rownames(s2)
cat("132 集基因数:", length(genes132), "\n\n")

idx <- match(genes132, s1$Real_Symbol)
cat("匹配到 DEG 表的:", sum(!is.na(idx)), "\n\n")

d <- data.frame(gene = genes132,
                logFC = s1$log2FoldChange[idx],
                padj  = s1$padj[idx],
                stringsAsFactors = FALSE)
d <- d[!is.na(d$logFC), ]

cat("=== 1. 132 基因的 logFC 方向 ===\n")
cat("在 ER+ 中下调 (logFC<0):", sum(d$logFC < 0), " (", round(100*mean(d$logFC<0),1), "%)\n", sep="")
cat("在 ER+ 中上调 (logFC>0):", sum(d$logFC > 0), " (", round(100*mean(d$logFC>0),1), "%)\n", sep="")
cat("logFC 中位数:", round(median(d$logFC), 3), "  范围:",
    round(min(d$logFC),2), "~", round(max(d$logFC),2), "\n")
cat("全部 12852 DEG 的 logFC 中位数:", round(median(s1$log2FoldChange, na.rm=TRUE),3), "\n")
cat("132 集 |logFC| 中位数:", round(median(abs(d$logFC)),3),
    " vs 全部 DEG:", round(median(abs(s1$log2FoldChange), na.rm=TRUE),3), "\n\n")

cat("=== 2. 四个 hub 基因在 ER+ vs ER- 中的表现 ===\n")
hub <- c("PLG","CD79A","MMP7","ANGPT1")
for (g in hub) {
  r <- d[d$gene == g, ]
  if (nrow(r) == 0) { cat(sprintf("  %-7s 不在 132 集内\n", g)); next }
  cat(sprintf("  %-7s logFC=%+7.3f  padj=%.2e  %s\n", g, r$logFC, r$padj,
              ifelse(r$logFC < 0, "ER+ 中低表达", "ER+ 中高表达")))
}

cat("\n=== 3. 增殖相关基因是否构成候选池 ===\n")
prolif <- c("MKI67","CCNB1","CCNB2","TOP2A","AURKA","BIRC5","CDK1","RRM2",
            "TYMS","UBE2C","CCNE1","CDC20","PTTG1","ANLN","KIF2C","MYBL2")
inp <- prolif[prolif %in% d$gene]
cat("132 集中包含的增殖基因:", ifelse(length(inp)==0, "无", paste(inp, collapse=", ")), "\n")
for (g in inp) {
  r <- d[d$gene == g, ]
  cat(sprintf("  %-7s logFC=%+7.3f\n", g, r$logFC))
}

cat("\n=== 4. 其他关键基因是否在候选池 ===\n")
for (g in c("ESR1","PGR","ERBB2","EPCAM","PTPRC","COL1A1","ACTA2","PECAM1","CDH1")) {
  cat(sprintf("  %-7s %s\n", g, ifelse(g %in% d$gene, "在 132 集内", "不在")))
}

cat("\n=== 5. 132 基因在 ER 分组判别上的强度 ===\n")
cat("|logFC| > 2 的占比:", round(100*mean(abs(d$logFC) > 2),1), "%\n")
cat("|logFC| > 3 的占比:", round(100*mean(abs(d$logFC) > 3),1), "%\n")
cat("padj < 1e-10 的占比:", round(100*mean(d$padj < 1e-10),1), "%\n")

write.csv(d[order(-abs(d$logFC)), ], "Table_CandidatePool_ERcontrast.csv", row.names = FALSE)
cat("\n已写出 Table_CandidatePool_ERcontrast.csv\n")
