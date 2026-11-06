# 复核：686 例（无 merge 膨胀）下的单基因相关性
suppressMessages({library(SummarizedExperiment); library(org.Hs.eg.db)})
load("../raw/TCGA/brca.gdc_2024.rda")
expr <- assay(expquery2)
colnames(expr) <- substr(gsub("\\.", "-", colnames(expr)), 1, 12)
ens <- gsub("\\..*", "", rownames(expr))
sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(expr) <- sym
expr <- expr[!is.na(rownames(expr)), ]; expr <- expr[!duplicated(rownames(expr)), ]

fd <- read.csv("Table_S7_Final_Clinical_Data.csv", stringsAsFactors = FALSE)
stopifnot(!any(duplicated(fd$Patient_ID)))
s <- intersect(colnames(expr), fd$Patient_ID)
cat("样本数:", length(s), "| 风险表行数:", nrow(fd), "\n")
rs <- fd$RiskScore[match(s, fd$Patient_ID)]

genes <- c("PADI4","ELANE","MPO","CTSG","PRTN3","AZU1","DNASE1L3","GSDMD",
           "FCGR3B","CSF3R","CXCR2","FPR1","S100A8","S100A9","S100A12","ITGAM",
           "PLG","CD79A","MMP7","ANGPT1",
           "MKI67","AURKA","BIRC5","CCNB1","AURKB",
           "ESR1","PGR","ERBB2","CD8A","FOXP3","CD68","PTPRC")
genes <- intersect(genes, rownames(expr))
out <- data.frame()
for (g in genes) {
  ct <- suppressWarnings(cor.test(rs, as.numeric(expr[g, s]), method = "spearman"))
  out <- rbind(out, data.frame(gene = g, rho = round(unname(ct$estimate), 3), p = signif(ct$p.value, 3)))
}
out <- out[order(out$rho), ]
print(out, row.names = FALSE)
cat("\n中性粒相关基因的平均 rho:", round(mean(out$rho[out$gene %in% c("PADI4","ELANE","MPO","CTSG","PRTN3","AZU1","DNASE1L3","GSDMD","FCGR3B","CSF3R","CXCR2","FPR1","S100A8","S100A9","S100A12","ITGAM")]), 3), "\n")
cat("增殖相关基因的平均 rho:", round(mean(out$rho[out$gene %in% c("MKI67","AURKA","BIRC5","CCNB1","AURKB")]), 3), "\n")
