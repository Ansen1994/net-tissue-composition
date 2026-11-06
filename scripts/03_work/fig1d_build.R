# ==============================================================================
# Figure 1D 数据（第二阶段）：把 132 个候选基因按主导细胞来源归类
# 输入：fig1d_celltype_mean.csv（fig1d_cellorigin.R 产出，gene x celltype 平均表达）
# 规则：
#   1) 已知中性粒限制性基因（在 10x 数据中未被捕获）→ Neutrophil-restricted
#   2) 最大平均表达 < 0.02 → Low / not quantified（数据里测不到，无法判定来源）
#   3) top/global < 1.5 → Broadly expressed（无细胞类型特异性）
#   4) 其余按最大平均表达的细胞类型归类
# 输出：fig1d_cell_origin_final.csv, fig1d_summary.csv
# ==============================================================================
M <- as.matrix(read.csv("fig1d_celltype_mean.csv", row.names = 1, check.names = FALSE))
M <- M[!grepl("\\.1$", rownames(M)), , drop = FALSE]        # 去掉 ELFN2.1 重复行
cat("基因数:", nrow(M), " 细胞类型:", paste(colnames(M), collapse = "/"), "\n")

mx   <- apply(M, 1, max)
top  <- colnames(M)[apply(M, 1, which.max)]
glob <- rowMeans(M)

# 中性粒限制性基因：文献明确、且本数据未捕获
neut <- c("ELANE","MPO","PADI4","AZU1","CTSG","PRTN3","BPI","CEACAM8",
          "LCN2","LTF","DEFA3","DEFA4","CD177","CXCR1","PRSS57","LPO")

d <- data.frame(gene = rownames(M), top_celltype = top,
                max_mean = mx, global_mean = glob,
                fold_top = mx / pmax(glob, 1e-9), stringsAsFactors = FALSE)

# 淋巴系合并（B/plasma 仅 1 个基因）
merge_lymph <- function(x) ifelse(x %in% c("B/plasma", "T/NK"), "Lymphoid", x)

d$class <- merge_lymph(d$top_celltype)
d$class[d$max_mean < 0.02]      <- "Low / not quantified"
d$class[d$fold_top < 1.5]       <- "Broadly expressed"
d$class[d$gene %in% neut]       <- "Neutrophil-restricted"

ord <- c("Neutrophil-restricted","Epithelial","Myeloid","Stromal",
         "Endothelial","Lymphoid","Broadly expressed","Low / not quantified")
d$class <- factor(d$class, levels = ord)
d <- d[order(d$class, -d$max_mean), ]

write.csv(d, "fig1d_cell_origin_final.csv", row.names = FALSE)

sm <- as.data.frame(table(d$class), stringsAsFactors = FALSE)
names(sm) <- c("class", "n")
sm$pct <- round(100 * sm$n / nrow(d), 1)
sm$examples <- sapply(sm$class, function(k) paste(head(d$gene[d$class == k], 5), collapse = ", "))
write.csv(sm, "fig1d_summary.csv", row.names = FALSE)

cat("\n=== 结果（共", nrow(d), "个基因）===\n")
for (i in seq_len(nrow(sm))) cat(sprintf("%-24s %3d (%4.1f%%)  %s\n",
    sm$class[i], sm$n[i], sm$pct[i], sm$examples[i]))
cat("\n中性粒限制性占:", round(100 * sum(d$class == "Neutrophil-restricted") / nrow(d), 1), "%\n")
