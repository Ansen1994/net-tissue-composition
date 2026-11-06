# ==============================================================================
# Figure 1D 数据：132 个候选基因由哪些细胞表达？
# 数据：GSE176078 单细胞（26 例，metadata 自带 celltype_major）
# 逻辑：每个样本内按细胞类型算平均表达（CP10K + log1p），跨样本取均值；
#       取表达最高的细胞类型作为该基因的「主导细胞来源」。
#       检出的基因里若没有主导类型 → 标注为未检出；再看是否属于中性粒颗粒基因。
# 输出：fig1d_cell_origin.csv, fig1d_celltype_mean.csv
# 注意：CP10K 的分母用「全转录组」counts，不是取子集后的 counts。
# ==============================================================================
suppressMessages({ library(Matrix) })

base <- "../raw/scRNA/sc"
dirs <- list.dirs(base, recursive = FALSE)
dirs <- dirs[file.exists(file.path(dirs, "metadata.csv"))]
cat("样本数:", length(dirs), "\n")

# 候选基因：与 Table_S2 的唯一 symbol 一致（132）
s2 <- read.csv("Table_S2_Target_NETs_Matrix.csv", row.names = 1, check.names = FALSE)
genes <- unique(rownames(s2))
cat("候选基因:", length(genes), "\n")

# 中性粒颗粒/限制性基因（用于给「未检出」的基因定性，不参与表达计算）
neut_restricted <- c("ELANE","MPO","PADI4","AZU1","CTSG","PRTN3","BPI","CAMP",
                     "LCN2","DEFA1","DEFA1B","DEFA3","DEFA4","CEACAM8","LTF",
                     "MS4A3","FCGR3B","CXCR2","S100A8","S100A9","S100A12",
                     "MMP8","MMP9","ARG1","CSF3R","CXCL8","GSDMD")

# 大类归并
grp_map <- c("Cancer Epithelial" = "Epithelial", "Normal Epithelial" = "Epithelial",
             "CAFs" = "Stromal", "PVL" = "Stromal",
             "Endothelial" = "Endothelial",
             "B-cells" = "B/plasma", "Plasmablasts" = "B/plasma",
             "T-cells" = "T/NK", "Myeloid" = "Myeloid")

acc <- list()   # 每个样本一张 gene x celltype 均值表

for (d in dirs) {
  sid <- basename(d)
  meta <- read.csv(file.path(d, "metadata.csv"), row.names = 1, stringsAsFactors = FALSE)
  gn   <- read.table(file.path(d, "count_matrix_genes.tsv"), stringsAsFactors = FALSE)[, 1]
  mat  <- readMM(file.path(d, "count_matrix_sparse.mtx"))
  if (nrow(mat) != length(gn)) mat <- t(mat)

  # barcode 对齐
  bc <- read.table(file.path(d, "count_matrix_barcodes.tsv"), stringsAsFactors = FALSE)[, 1]
  key_meta <- sub("^[^_]+_", "", rownames(meta)); key_meta <- sub("-1$", "", key_meta)
  key_bc   <- sub("-1$", "", bc)
  idx <- match(key_bc, key_meta)
  if (sum(!is.na(idx)) < 0.9 * length(bc)) idx <- seq_len(min(nrow(meta), ncol(mat)))
  ok <- !is.na(idx)
  mat <- mat[, ok, drop = FALSE]
  ct  <- meta$celltype_major[idx[ok]]

  # 分母：全转录组 counts
  cp <- Matrix::colSums(mat); cp[cp == 0] <- 1

  keep <- intersect(genes, gn)
  sub  <- mat[match(keep, gn), , drop = FALSE]
  Mn   <- log1p(sweep(as.matrix(sub), 2, cp, "/") * 1e4)
  rownames(Mn) <- keep

  grp <- grp_map[ct]; grp[is.na(grp)] <- "Other"
  tb <- sapply(split(seq_along(grp), grp), function(ix) rowMeans(Mn[, ix, drop = FALSE]))
  tb[is.na(tb)] <- 0
  # 补全缺失基因/类型
  full <- matrix(0, length(genes), length(unique(grp_map)),
                 dimnames = list(genes, sort(unique(grp_map))))
  full[rownames(tb), colnames(tb)] <- tb
  acc[[sid]] <- full
  cat(sprintf("  %s  cells=%d  types=%d\n", sid, ncol(Mn), length(unique(grp))))
}

celltypes <- colnames(acc[[1]])
M <- Reduce(`+`, acc) / length(acc)     # 跨样本均值
write.csv(round(M, 5), "fig1d_celltype_mean.csv")

# 主导细胞类型
mx <- apply(M, 1, max)
top <- colnames(M)[apply(M, 1, which.max)]
glob <- rowMeans(M)

out <- data.frame(
  gene = genes,
  top_celltype = top,
  max_mean = round(mx, 4),
  global_mean = round(glob, 4),
  fold_top = round(mx / pmax(glob, 1e-9), 2),
  detected = mx >= 0.10,
  stringsAsFactors = FALSE
)
# 检出但极弱 → 也算未检出
out$detected <- out$detected & out$global_mean >= 0.10

# 分类
out$class <- ifelse(!out$detected, "Not detected", out$top_celltype)
out$class[out$class == "Not detected" & out$gene %in% neut_restricted] <- "Neutrophil (not captured)"
out$class[out$class == "Not detected" & !(out$gene %in% neut_restricted)] <- "Not detected"

# 主导类型但表达很弱的，也降级（避免把噪声当来源）
weak <- out$detected & out$max_mean < 0.5
out$class[weak] <- paste0(out$top_celltype[weak], " (weak)")

out <- out[order(out$class, -out$max_mean), ]
write.csv(out, "fig1d_cell_origin.csv", row.names = FALSE)

cat("\n=== 分类汇总 ===\n")
print(table(out$class))
cat("\n=== 各类代表基因（按表达强度）===\n")
for (k in unique(out$class)) {
  g <- head(out$gene[out$class == k], 8)
  cat(sprintf("%-26s %s\n", k, paste(g, collapse = ", ")))
}
cat("\n未检出基因数:", sum(out$class %in% c("Not detected","Neutrophil (not captured)")), "\n")
