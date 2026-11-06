# ==============================================================================
# 判定 3：四个 hub 基因在单细胞层面到底由谁表达？
# 数据：GSE176078（Wu et al. 2021 Nat Genet），26 例乳腺肿瘤，metadata 自带细胞类型
# 关键输出：
#   (a) 4 个 hub 基因的主导表达细胞类型
#   (b) NETs 执行基因（MPO/ELANE/PADI4/GSDMD）在肿瘤组织中的检出率
#       —— 若检出率极低，则「NETs 相关基因集」在实体瘤转录组中无对应细胞来源
# 输出：j3_celltype_expression.csv, j3_nets_detection.csv,
#       j3_annotation_audit.csv, j3_marker_verification.csv
# ==============================================================================
suppressMessages({ library(Matrix) })

base <- "../raw/scRNA/sc"
dirs <- list.dirs(base, recursive = FALSE)
dirs <- dirs[file.exists(file.path(dirs, "metadata.csv"))]
cat("样本目录数:", length(dirs), "\n")

genes_need <- c("PLG","CD79A","MMP7","ANGPT1",
                "EPCAM","KRT8","CD3D","MS4A1","LYZ","PECAM1","COL1A1","MKI67",
                "PADI4","ELANE","MPO","GSDMD","FCGR3B","CXCR2","S100A8","S100A9")

all_cells <- list()
for (d in dirs) {
  sid <- basename(d)
  mf <- file.path(d, "metadata.csv")
  gf <- file.path(d, "count_matrix_genes.tsv")
  xf <- file.path(d, "count_matrix_sparse.mtx")
  bf <- file.path(d, "count_matrix_barcodes.tsv")
  if (!all(file.exists(c(mf, gf, xf)))) { cat("跳过", sid, "\n"); next }

  meta <- read.csv(mf, row.names = 1, stringsAsFactors = FALSE)
  gn <- read.table(gf, stringsAsFactors = FALSE)[, 1]
  mat <- readMM(xf)
  keep <- intersect(genes_need, gn)
  mat <- mat[match(keep, gn), , drop = FALSE]
  rownames(mat) <- keep

  # barcode 对齐：metadata 行名形如 CID3941_AAAGTAGTCCTTGCCA
  if (file.exists(bf)) {
    bc <- read.table(bf, stringsAsFactors = FALSE)[, 1]
    if (length(bc) == ncol(mat)) {
      colnames(mat) <- bc
      key_meta <- sub("^[^_]+_", "", rownames(meta))
      key_meta <- sub("-1$", "", key_meta)
      key_bc   <- sub("-1$", "", bc)
      idx <- match(key_bc, key_meta)
      if (sum(!is.na(idx)) < 0.9 * nrow(meta)) {
        # 顺序一致时直接对齐
        idx <- seq_len(min(nrow(meta), ncol(mat)))
      }
    } else idx <- seq_len(ncol(mat))
  } else idx <- seq_len(ncol(mat))

  M <- as.matrix(mat)
  # 文库大小必须用全转录组 UMI 总数（metadata$nCount_RNA），不能用上面裁出的
  # 20 基因子集的 colSums()：后者随细胞类型变化（髓系 S100A8/A9 高、分母大），
  # 会让 log1p(count/cp*1e4) 不再是 CP10K，并放大/扭曲各区室的相对均值。
  # 与 12_可复用资源/sc_umap 的 Seurat::NormalizeData(LogNormalize, 1e4) 同口径。
  cp <- as.numeric(meta$nCount_RNA[idx])
  cp[is.na(cp) | cp <= 0] <- 1
  Mn <- log1p(sweep(M, 2, cp, "/") * 1e4)   # CP10K + log1p（全转录组文库）

  df <- data.frame(
    cell   = colnames(M),
    sample = sid,
    subtype = meta$subtype[idx],
    major   = meta$celltype_major[idx],
    minor   = meta$celltype_minor[idx],
    nCount  = meta$nCount_RNA[idx],
    stringsAsFactors = FALSE)
  for (g in keep) { df[[g]] <- as.numeric(M[g, ]); df[[paste0(g, "_n")]] <- as.numeric(Mn[g, ]) }
  all_cells[[sid]] <- df
  cat(sprintf("  %-12s cells=%-6d genes=%d  subtype=%s\n",
              sid, ncol(M), length(keep), paste(unique(meta$subtype[idx]), collapse="/")))
  rm(mat, M, Mn); gc(verbose = FALSE)
}

A <- do.call(rbind, all_cells)
cat("\n总细胞数:", nrow(A), "  ER+ 细胞:", sum(A$subtype == "ER+", na.rm = TRUE), "\n")
cat("细胞类型大类:\n"); print(table(A$major, useNA = "ifany"))

# ------------------------------------------------------------------ (a) hub 基因表达定位
hub <- c("PLG","CD79A","MMP7","ANGPT1")
cat("\n\n================ 一、4 个 hub 基因的主要表达细胞类型 ================\n")
cat("（基于 ER+ 样本；数值 = CP10K-log 均值，括号内为检出率）\n\n")
E <- A[!is.na(A$subtype) & A$subtype == "ER+", ]
cat("ER+ 细胞数:", nrow(E), "\n\n")

out <- list()
for (g in hub) {
  agg <- aggregate(E[[paste0(g, "_n")]], by = list(major = E$major), FUN = mean)
  names(agg)[2] <- "mean_expr"
  det <- aggregate(E[[g]] > 0, by = list(major = E$major), FUN = mean)
  agg$detect <- det$x
  agg$n <- as.numeric(table(E$major)[agg$major])
  agg <- agg[order(-agg$mean_expr), ]
  cat(sprintf("---- %s ----\n", g))
  for (i in seq_len(nrow(agg)))
    cat(sprintf("   %-18s %.3f  (%.0f%% 检出)  n=%d\n",
                agg$major[i], agg$mean_expr[i], 100*agg$detect[i], agg$n[i]))
  cat(sprintf("   最高/次高 比值 = %.1f\n\n", agg$mean_expr[1]/max(agg$mean_expr[2],1e-6)))
  agg$gene <- g
  out[[g]] <- agg
}
O <- do.call(rbind, out)
write.csv(O, "j3_celltype_expression.csv", row.names = FALSE)

# ------------------------------------------------------------------ (b) NETs 基因检出率
cat("================ 二、NETs 执行基因在肿瘤组织中的检出情况 ================\n")
cat("（若这些基因在全部细胞中检出率极低，说明肿瘤组织缺乏中性粒细胞来源）\n\n")
nets <- c("PADI4","ELANE","MPO","GSDMD","FCGR3B","CXCR2","S100A8","S100A9")
nd <- data.frame(gene = nets,
                 detect_all = sapply(nets, function(g) mean(A[[g]] > 0, na.rm = TRUE)),
                 mean_all   = sapply(nets, function(g) mean(A[[paste0(g,"_n")]], na.rm = TRUE)),
                 detect_ERp = sapply(nets, function(g) mean(E[[g]] > 0, na.rm = TRUE)),
                 mean_ERp   = sapply(nets, function(g) mean(E[[paste0(g,"_n")]], na.rm = TRUE)),
                 stringsAsFactors = FALSE)
nd <- nd[order(-nd$detect_all), ]
for (i in seq_len(nrow(nd)))
  cat(sprintf("  %-8s 全样本检出 %5.1f%%   ER+ 检出 %5.1f%%   均值 %.3f\n",
              nd$gene[i], 100*nd$detect_all[i], 100*nd$detect_ERp[i], nd$mean_all[i]))
write.csv(nd, "j3_nets_detection.csv", row.names = FALSE)

# ------------------------------------------------------------------ (c) 对照：谱系 marker
cat("\n================ 三、谱系 marker 对照（注释可信度检查）================\n")
mk <- c("EPCAM","CD3D","MS4A1","LYZ","PECAM1","COL1A1","MKI67")
mg <- data.frame(gene = mk,
                 top_celltype = sapply(mk, function(g) {
                   t <- aggregate(E[[paste0(g,"_n")]], by=list(m=E$major), FUN=mean)
                   t$m[which.max(t$x)] }),
                 top_expr = sapply(mk, function(g) {
                   t <- aggregate(E[[paste0(g,"_n")]], by=list(m=E$major), FUN=mean)
                   round(max(t$x), 2) }),
                 stringsAsFactors = FALSE)
for (i in seq_len(nrow(mg)))
  cat(sprintf("  %-8s 最高表达于 %-18s (%.2f)\n", mg$gene[i], mg$top_celltype[i], as.numeric(mg$top_expr[i])))

saveRDS(A, "j3_allcells.rds")

# ------------------------------------------ (d) 注释体系审计（供补充材料 / 后续出图）
# 本段复用作者发布的 curated 标签，未做 QC 过滤 / 重新聚类 / 重新注释；
# 展示用 UMAP 由 12_可复用资源/sc_umap/ 单独重算，不参与本段统计。
# 这里只记录标签体系本身，特别是「有无中性粒细胞类」。
cat("================ 四、注释体系审计 ================\n")
mk_tab <- function(x, lev) data.frame(level = lev, label = names(x),
                                      n = as.integer(x), stringsAsFactors = FALSE)
audit <- rbind(mk_tab(table(A$major), "major"),
               mk_tab(table(A$minor), "minor"),
               mk_tab(table(A$subtype), "subtype"))
write.csv(audit, "j3_annotation_audit.csv", row.names = FALSE)
cat(sprintf("  major=%d minor=%d subtype=%d; neutrophil/granulocyte labels=%d\n",
            sum(audit$level == "major"), sum(audit$level == "minor"),
            sum(audit$level == "subtype"),
            sum(grepl("neutro|granulo|pmn", audit$label[audit$level != "subtype"],
                      ignore.case = TRUE))))

mk_ver <- do.call(rbind, lapply(
  c("EPCAM", "CD3D", "MS4A1", "LYZ", "PECAM1", "COL1A1", "MKI67"),
  function(g) {
    t <- aggregate(E[[paste0(g, "_n")]], by = list(major = E$major), FUN = mean)
    data.frame(marker = g, top_major = as.character(t$major[which.max(t$x)]),
               mean_expr = round(max(t$x), 3), stringsAsFactors = FALSE)
  }))
write.csv(mk_ver, "j3_marker_verification.csv", row.names = FALSE)
print(mk_ver)

cat("\n完成。\n")
