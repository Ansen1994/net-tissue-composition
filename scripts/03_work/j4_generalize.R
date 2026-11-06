# ==============================================================================
# 泛化验证：7 个已发表 NETs 签名 + 本研究签名，在 bulk 层面到底测的是什么？
#
# 问题：如果这些签名真的反映中性粒活性，评分应与中性粒/NETosis 模块强相关；
#       如果反映组织构成，则应与基质/上皮比例强相关。
#
# 数据：TCGA-BRCA（brca.gdc_2024.rda，tpm_unstrand，Primary Tumor）
# 方法：基因 z-score → 模块内均值；签名与 9 个参照模块的 Spearman 相关；
#       并以 lm 的 R² 量化各签名被中性粒与基质分别解释的方差份额。
# 输出：j4_signature_purity.csv, j4_module_correlations.csv
# ==============================================================================
suppressMessages({ library(SummarizedExperiment) })

# ---- 1. 已发表签名（基因列表均取自可核实的文献正文/摘要） ----
sigs <- list(
  BRCA_A = c("LTF","LCP1","AZU1","ENO1"),
  BRCA_B = c("F2RL2","AZU1","IL33","ELANE"),
  BRCA_C = c("ZMYND10","IL12B","CXCL13","TFF1","LTB","SPIB"),
  GC_A   = c("PDE4B","CD93","CTSG","IL6","KCNJ15","CRISPLD2","ELANE"),
  GC_B   = c("MPO","BST1","IL6","KCNJ15","SELP","ELANE"),
  CRC_A  = c("AKT2","ATG7","CLEC6A","CRISPLD2","ENTPD4","IL17A","ITGAM","MAPK1","SELPLG","SIGLEC14"),
  Ours   = c("PLG","CD79A","MMP7","ANGPT1")
)
sig_src <- c(
  BRCA_A = "PMC10250592 (Front Genet, TCGA BRCA)",
  BRCA_B = "PMC11985209 (Transl Cancer Res, TCGA BRCA)",
  BRCA_C = "TCGA-BRCA NETs+inflammation, 6 genes",
  GC_A   = "PMC10438473 (gastric cancer)",
  GC_B   = "AJCO 2025 (gastric cancer)",
  CRC_A  = "Springer 2025 (colorectal cancer)",
  Ours   = "This study (TCGA BRCA ER+)"
)

# 公认的中性粒限制/颗粒基因（HPA + NETs 文献），用于判定签名内中性粒基因占比
neutrophil_ref <- c("AZU1","ELANE","MPO","CTSG","LTF","PRTN3","DEFA3","DEFA4",
                    "CAMP","LCN2","CEACAM8","FCGR3B","CXCR2","S100A8","S100A9",
                    "PADI4","CSF3R","FPR1")

# ---- 2. 参照模块 ----
modules <- list(
  Neutrophil    = c("FCGR3B","CXCR2","S100A8","S100A9","CSF3R","FPR1","CAMP","LCN2","MMP9"),
  NETosis       = c("PADI4","ELANE","MPO","AZU1","CTSG","PRTN3"),
  Stroma_ECM    = c("COL1A1","COL1A2","COL3A1","FAP","PDGFRB","ACTA2","POSTN","FN1","THY1"),
  Endothelium   = c("PECAM1","VWF","CDH5","CLDN5","KDR"),
  Proliferation = c("MKI67","CCNB1","TOP2A","BIRC5","AURKA","CDK1","RRM2","TYMS","UBE2C","CCNE1"),
  B_cells       = c("CD19","MS4A1","CD79A","CD79B","PAX5"),
  T_cells       = c("CD3D","CD3E","CD8A","CD4"),
  Myeloid       = c("LYZ","CD68","ITGAM","CSF1R"),
  Epithelial    = c("EPCAM","KRT8","KRT18","CDH1")
)

all_need <- unique(c(unlist(sigs), unlist(modules), neutrophil_ref))

# ---- 3. 取 TCGA 表达 ----
load("../raw/TCGA/brca.gdc_2024.rda")
x <- expquery2
keep_col <- which(as.character(colData(x)$definition) == "Primary solid Tumor")
x <- x[, keep_col]
cat("原发肿瘤样本数:", ncol(x), "\n")

gn <- as.character(rowData(x)$gene_name)
ok <- !is.na(gn) & gn != ""
x  <- x[ok, ]; gn <- gn[ok]

need_rows <- which(gn %in% all_need)
cat("命中目标基因行数:", length(need_rows), " 目标基因数:", length(all_need), "\n")

m <- assay(x, "tpm_unstrand")[need_rows, , drop = FALSE]
m <- log2(as.matrix(m) + 1)
rg <- gn[need_rows]

# 多行映射到同一符号：保留表达均值最高的一行
o <- order(rowMeans(m), decreasing = TRUE)
m <- m[o, ]; rg <- rg[o]
dup <- duplicated(rg)
m <- m[!dup, ]; rg <- rg[!dup]
rownames(m) <- rg
cat("去重后基因数:", nrow(m), "\n")

missing <- setdiff(all_need, rownames(m))
if (length(missing)) cat("缺失基因:", paste(missing, collapse = ", "), "\n")

# ---- 4. 模块评分 ----
z <- t(scale(t(m)))
module_score <- sapply(modules, function(g) {
  g <- intersect(g, rownames(z))
  if (!length(g)) return(rep(NA, ncol(z)))
  colMeans(z[g, , drop = FALSE], na.rm = TRUE)
})
sig_score <- sapply(sigs, function(g) {
  g <- intersect(g, rownames(z))
  colMeans(z[g, , drop = FALSE], na.rm = TRUE)
})
cat("模块评分 dim:", paste(dim(module_score), collapse = "x"), "\n")

# ---- 5. 相关矩阵 ----
res <- NULL
for (s in colnames(sig_score)) {
  row <- data.frame(signature = s, stringsAsFactors = FALSE)
  for (md in colnames(module_score)) {
    ct <- suppressWarnings(cor.test(sig_score[, s], module_score[, md],
                                    method = "spearman", exact = FALSE))
    row[[paste0("rho_", md)]] <- round(unname(ct$estimate), 3)
    row[[paste0("p_", md)]]   <- signif(ct$p.value, 3)
  }
  res <- rbind(res, row)
}

# ---- 6. 方差分解：签名 ~ Neutrophil + Stroma + Prolif ----
dec <- t(sapply(colnames(sig_score), function(s) {
  y  <- sig_score[, s]
  r2 <- function(v) summary(lm(y ~ module_score[, v]))$r.squared
  full <- summary(lm(y ~ module_score[, "Neutrophil"] +
                       module_score[, "Stroma_ECM"] +
                       module_score[, "Proliferation"]))$r.squared
  c(R2_Neutrophil = r2("Neutrophil"),
    R2_NETosis    = r2("NETosis"),
    R2_Stroma     = r2("Stroma_ECM"),
    R2_Prolif     = r2("Proliferation"),
    R2_joint      = full)
}))

out <- data.frame(
  signature   = res$signature,
  source      = sig_src[res$signature],
  n_genes     = sapply(sigs[res$signature], length),
  n_neutrophil= sapply(sigs[res$signature], function(g) sum(g %in% neutrophil_ref)),
  res[, grep("^rho_", names(res)), drop = FALSE],
  round(dec[res$signature, ], 3),
  check.names = FALSE, row.names = NULL
)
out$pct_neutrophil <- round(100 * out$n_neutrophil / out$n_genes, 1)

# 主导成分 = 除 NETosis 外相关最强（按 |rho|）的模块
rho_cols <- grep("^rho_", names(out), value = TRUE)
out$dominant <- apply(out[, rho_cols, drop = FALSE], 1, function(v) {
  sub("^rho_", "", rho_cols[which.max(abs(v))])
})

write.csv(out, "j4_signature_purity.csv", row.names = FALSE)
write.csv(as.data.frame(module_score), "j4_module_scores.csv", row.names = TRUE)

cat("\n==== 汇总 ====\n")
print(out[, c("signature","n_genes","n_neutrophil","pct_neutrophil",
              "rho_Neutrophil","rho_NETosis","rho_Stroma_ECM",
              "rho_Proliferation","rho_Epithelial","dominant")], row.names = FALSE)

cat("\n==== 方差分解（R²） ====\n")
print(out[, c("signature","R2_Neutrophil","R2_NETosis","R2_Stroma","R2_Prolif","R2_joint")],
      row.names = FALSE)

cat("\n中性粒基因总占比:", sum(out$n_neutrophil), "/", sum(out$n_genes), "=",
    round(100 * sum(out$n_neutrophil) / sum(out$n_genes), 1), "%\n")

cat("\n==== 参照模块内部一致性（与 Neutrophil 模块的相关） ====\n")
for (md in c("NETosis","Myeloid","Stroma_ECM")) {
  ct <- suppressWarnings(cor.test(module_score[, "Neutrophil"], module_score[, md],
                                  method = "spearman", exact = FALSE))
  cat(sprintf("  Neutrophil vs %-11s rho=%+.3f p=%.2e\n", md, ct$estimate, ct$p.value))
}
