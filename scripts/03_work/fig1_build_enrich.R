# GO / KEGG 富集（132 个候选基因），结果存 CSV 供绘图复用
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db)})
d <- read.csv("fig1d_cell_origin_final.csv", stringsAsFactors = FALSE)
g <- unique(d$gene)
cat("输入基因:", length(g), "\n")

ids <- bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go <- enrichGO(gene = g, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
               ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
               readable = TRUE)
go_df <- as.data.frame(go)
go_df <- go_df[go_df$pvalue < 0.05, ]
go_df$ONTOLOGY <- factor(go_df$ONTOLOGY, levels = c("BP","CC","MF"))
go_df <- go_df[order(go_df$ONTOLOGY, go_df$pvalue), ]
write.csv(go_df, "fig1_go_enrich.csv", row.names = FALSE)
cat("GO 显著条目:", nrow(go_df), "\n"); print(table(go_df$ONTOLOGY))

kg <- enrichKEGG(ids$ENTREZID, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)
kg_df <- as.data.frame(setReadable(kg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
kg_df <- kg_df[kg_df$pvalue < 0.05, ]
kg_df <- kg_df[order(kg_df$pvalue), ]
write.csv(kg_df, "fig1_kegg_enrich.csv", row.names = FALSE)
cat("KEGG 显著条目:", nrow(kg_df), "\n")
cat("\nTop10 KEGG:/n"); print(head(kg_df[, c("Description","pvalue","Count","GeneRatio")], 10), row.names = FALSE)
cat("\n每本体 top3 GO:/n")
for (o in c("BP","CC","MF")) {
  s <- head(go_df[go_df$ONTOLOGY == o, c("Description","pvalue","GeneRatio")], 3)
  cat("--", o, "--\n"); print(s, row.names = FALSE)
}
