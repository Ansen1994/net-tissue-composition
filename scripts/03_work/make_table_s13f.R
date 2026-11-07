# ==============================================================================
# Table S13f —— 签名纯度画像（与 Table S13e / Figure 7D 同源重算）
#
# 背景：旧 j4_signature_purity.csv 由 j4_generalize.R 产生，签名一律等权；
#       而 Table S13e 与 Figure 7D 由 j4c_generalize.R 产生，BRCA_A 与本研究
#       签名按原 Cox 系数加权。两种评分方式下同一签名的相关符号相反
#       （BRCA_A × Neutrophil +0.523 vs -0.310；Ours × B_cells +0.695 vs -0.298），
#       并列在一份补充材料里无法解释。本脚本沿用 j4c 的评分与模块定义重算，
#       使 S13f 的 rho 列与 S13e 逐值相等。
#
# 输入：../raw/TCGA/brca.gdc_2024.rda
# 输出：j4c_signature_purity.csv
# 用法：在 work/ 目录下 Rscript make_table_s13f.R
# ==============================================================================
suppressMessages({ library(SummarizedExperiment) })

# ---- 1. 签名（与 j4c_generalize.R 完全一致） ----
sigs <- list(
  BRCA_A = c("LTF","LCP1","AZU1","ENO1"),
  BRCA_B = c("F2RL2","AZU1","IL33","ELANE"),
  BRCA_C = c("ZMYND10","IL12B","CXCL13","TFF1","LTB","SPIB"),
  GC_A   = c("PDE4B","CD93","CTSG","IL6","KCNJ15","CRISPLD2","ELANE"),
  GC_B   = c("MPO","BST1","IL6","KCNJ15","SELP","ELANE"),
  CRC_A  = c("AKT2","ATG7","CLEC6A","CRISPLD2","ENTPD4","IL17A","ITGAM","MAPK1","SELPLG","SIGLEC14"),
  Ours   = c("PLG","CD79A","MMP7","ANGPT1")
)
# 已知已发表系数的签名按原公式加权，其余等权（同 j4c）
w_sigs <- list(
  BRCA_A = c(ENO1 = 0.012, AZU1 = -0.059, LCP1 = -0.073, LTF = -0.033),
  Ours   = c(ANGPT1 = log(1.258), PLG = log(1.117), CD79A = log(0.918), MMP7 = log(0.918))
)
sig_src <- c(
  BRCA_A = "PMC10250592 (Front Mol Biosci, TCGA BRCA)",
  BRCA_B = "PMC11985209 (Transl Cancer Res, TCGA BRCA)",
  BRCA_C = "TCGA-BRCA NETs+inflammation, 6 genes",
  GC_A   = "PMC10438473 (gastric cancer)",
  GC_B   = "Asia Pac J Clin Oncol 2026 (gastric cancer)",
  CRC_A  = "Springer 2025 (colorectal cancer)",
  Ours   = "This study (TCGA BRCA ER+)"
)

# 16 基因共识（Methods 2.9 / Figure 1D）
neutrophil_ref <- c("ELANE","MPO","PADI4","AZU1","CTSG","PRTN3","BPI","CEACAM8",
                    "LCN2","LTF","DEFA3","DEFA4","CD177","CXCR1","PRSS57","LPO")

# ---- 2. 九个泛化模块（同 j4c） ----
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
mod_order <- names(modules)

all_need <- unique(c(unlist(sigs), unlist(modules), neutrophil_ref))

# ---- 3. 取 TCGA 原发瘤表达 ----
load("../raw/TCGA/brca.gdc_2024.rda")
x <- expquery2
x <- x[, which(as.character(colData(x)$definition) == "Primary solid Tumor")]
gn <- as.character(rowData(x)$gene_name); ok <- !is.na(gn) & gn != ""
x <- x[ok, ]; gn <- gn[ok]
need_rows <- which(gn %in% all_need)
m <- log2(as.matrix(assay(x, "tpm_unstrand")[need_rows, , drop = FALSE]) + 1)
rg <- gn[need_rows]
o <- order(rowMeans(m), decreasing = TRUE); m <- m[o, ]; rg <- rg[o]
m <- m[!duplicated(rg), ]; rg <- rg[!duplicated(rg)]; rownames(m) <- rg
z <- t(scale(t(m)))
cat("样本:", ncol(z), " 基因:", nrow(z), "\n\n")

# ---- 4. 签名评分：与 j4c 完全相同的加权规则 ----
sig_score <- sapply(names(sigs), function(s) {
  g <- intersect(sigs[[s]], rownames(z))
  if (!is.null(w_sigs[[s]])) {
    w <- w_sigs[[s]][g]; w[is.na(w)] <- 0
    if (all(w == 0)) return(colMeans(z[g, , drop = FALSE]))
    as.numeric(crossprod(w, z[g, , drop = FALSE]) / sum(abs(w)))
  } else colMeans(z[g, , drop = FALSE])
})

# 完整模块评分（方差分解用；模块代表整个区室，不剔除共享基因）
full_score <- sapply(modules, function(g) {
  colMeans(z[intersect(g, rownames(z)), , drop = FALSE])
})

# ---- 5. 相关：模块内与签名共享的基因先剔除（与 S13e 同一规则） ----
rho <- matrix(NA_real_, nrow = length(sigs), ncol = length(mod_order),
              dimnames = list(names(sigs), mod_order))
n_module <- rho
for (s in names(sigs)) for (md in mod_order) {
  g_use <- intersect(setdiff(modules[[md]], sigs[[s]]), rownames(z))
  n_module[s, md] <- length(g_use)
  if (length(g_use) < 3) next          # 少于 3 个基因不计算，与 j4c 一致
  sc <- colMeans(z[g_use, , drop = FALSE])
  ct <- suppressWarnings(cor.test(sig_score[, s], sc, method = "spearman", exact = FALSE))
  rho[s, md] <- round(unname(ct$estimate), 3)
}

# 与 S13e 的一致性自检
if (file.exists("j4c_signature_module_correlations.csv")) {
  e <- read.csv("j4c_signature_module_correlations.csv")
  mism <- 0L
  for (i in seq_len(nrow(e))) {
    v <- rho[e$signature[i], e$module[i]]
    if (is.na(v) || abs(v - e$rho[i]) > 1e-9) mism <- mism + 1L
  }
  cat("与 Table S13e 逐值比对:", nrow(e), "格, 不一致", mism, "格\n\n")
}

# ---- 6. 方差分解：签名 ~ 各区室（单变量 R²）与 Neutrophil+Stroma+Prolif 联合 ----
dec <- t(sapply(names(sigs), function(s) {
  y   <- sig_score[, s]
  r2  <- function(v) summary(lm(y ~ full_score[, v]))$r.squared
  joint <- summary(lm(y ~ full_score[, "Neutrophil"] +
                         full_score[, "Stroma_ECM"] +
                         full_score[, "Proliferation"]))$r.squared
  c(R2_Neutrophil = r2("Neutrophil"), R2_NETosis = r2("NETosis"),
    R2_Stroma = r2("Stroma_ECM"), R2_Prolif = r2("Proliferation"),
    R2_joint = joint)
}))

# ---- 7. 汇总 ----
out <- data.frame(
  signature    = names(sigs),
  source       = sig_src[names(sigs)],
  n_genes      = sapply(sigs, length),
  n_neutrophil = sapply(sigs, function(g) sum(g %in% neutrophil_ref)),
  stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE
)
for (md in mod_order) out[[paste0("rho_", md)]] <- as.numeric(rho[, md])
out <- cbind(out, round(dec[names(sigs), ], 3))
out$pct_neutrophil <- round(100 * out$n_neutrophil / out$n_genes, 1)
out$dominant <- apply(rho, 1, function(v) mod_order[which.max(abs(v))])

write.csv(out, "j4c_signature_purity.csv", row.names = FALSE)

cat("==== Table S13f（重算） ====\n")
print(out[, c("signature","n_genes","n_neutrophil","pct_neutrophil","dominant")], row.names = FALSE)
cat("\n==== rho ====\n"); print(rho, quote = FALSE)
cat("\n==== 主导成分 ====\n")
for (i in seq_len(nrow(out)))
  cat(sprintf("  %-7s -> %-12s rho=%+.3f\n", out$signature[i], out$dominant[i],
              rho[i, out$dominant[i]]))
cat("\n==== 方差分解 ====\n")
print(out[, c("signature","R2_Neutrophil","R2_NETosis","R2_Stroma","R2_Prolif","R2_joint")],
      row.names = FALSE)
