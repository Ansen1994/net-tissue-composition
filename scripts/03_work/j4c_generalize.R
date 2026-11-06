# ==============================================================================
# 泛化验证（最终版）
# 修正 (1) 主导成分判定漏列 bug；(2) 本研究签名使用原 Cox 系数（保留正负号）
#        (3) 输出模块间相关矩阵以说明共变程度；(4) 报告 NETs 执行分子的可测性
# ==============================================================================
suppressMessages({ library(SummarizedExperiment) })

sigs <- list(
  BRCA_A = c("LTF","LCP1","AZU1","ENO1"),
  BRCA_B = c("F2RL2","AZU1","IL33","ELANE"),
  BRCA_C = c("ZMYND10","IL12B","CXCL13","TFF1","LTB","SPIB"),
  GC_A   = c("PDE4B","CD93","CTSG","IL6","KCNJ15","CRISPLD2","ELANE"),
  GC_B   = c("MPO","BST1","IL6","KCNJ15","SELP","ELANE"),
  CRC_A  = c("AKT2","ATG7","CLEC6A","CRISPLD2","ENTPD4","IL17A","ITGAM","MAPK1","SELPLG","SIGLEC14"),
  Ours   = c("PLG","CD79A","MMP7","ANGPT1")
)
# 已知系数的签名按原公式加权（其余等权）
w_sigs <- list(
  BRCA_A = c(ENO1 = 0.012, AZU1 = -0.059, LCP1 = -0.073, LTF = -0.033),
  Ours   = c(ANGPT1 = log(1.258), PLG = log(1.117), CD79A = log(0.918), MMP7 = log(0.918))
)

# 16-gene neutrophil-restricted / granule consensus used in Methods 2.9 and in
# fig1d_build.R for Figure 1D. Replaced the earlier 18-gene list on 2026-11-05
# after an audit: both lists give identical counts for all seven signatures
# (2/2/0/2/2/0/0; 8 of 41 loci, 19.5%). The seven genes removed are all already
# present in modules$Neutrophil, so all_need and every published result are
# unchanged. See Table S13b of the consolidated supplementary workbook.
neutrophil_ref <- c("ELANE","MPO","PADI4","AZU1","CTSG","PRTN3","BPI","CEACAM8",
                    "LCN2","LTF","DEFA3","DEFA4","CD177","CXCR1","PRSS57","LPO")
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

sig_score <- sapply(names(sigs), function(s) {
  g <- intersect(sigs[[s]], rownames(z))
  if (!is.null(w_sigs[[s]])) {
    w <- w_sigs[[s]][g]
    w[is.na(w)] <- 0
    if (all(w == 0)) return(colMeans(z[g, , drop = FALSE]))
    as.numeric(crossprod(w, z[g, , drop = FALSE]) / sum(abs(w)))
  } else colMeans(z[g, , drop = FALSE])
})
epi <- colMeans(z[intersect(modules$Epithelial, rownames(z)), , drop = FALSE])

pscor <- function(a, b, ctrl) {
  ra <- residuals(lm(rank(a) ~ rank(ctrl)))
  rb <- residuals(lm(rank(b) ~ rank(ctrl)))
  ct <- suppressWarnings(cor.test(ra, rb, method = "pearson"))
  c(rho = unname(ct$estimate), p = ct$p.value)
}

rows <- list()
for (s in names(sigs)) for (md in names(modules)) {
  g_use <- intersect(setdiff(modules[[md]], sigs[[s]]), rownames(z))
  if (length(g_use) < 3) next
  sc <- colMeans(z[g_use, , drop = FALSE])
  ct <- suppressWarnings(cor.test(sig_score[, s], sc, method = "spearman", exact = FALSE))
  pc <- pscor(sig_score[, s], sc, epi)
  rows[[length(rows)+1]] <- data.frame(signature = s, module = md, n_module = length(g_use),
    rho = unname(ct$estimate), p = ct$p.value, prho = pc["rho"], p_partial = pc["p"])
}
R_ <- do.call(rbind, rows)
R_$rho <- round(R_$rho, 3); R_$prho <- round(R_$prho, 3)
write.csv(R_, "j4c_signature_module_correlations.csv", row.names = FALSE)

key <- c("Neutrophil","NETosis","Stroma_ECM","Proliferation","Epithelial",
         "B_cells","T_cells","Myeloid","Endothelium")
W <- reshape(R_[, c("signature","module","rho")], idvar = "signature", timevar = "module", direction = "wide")
names(W) <- sub("rho\\.", "", names(W)); W <- W[, c("signature", key)]
cat("==== Spearman rho（已剔除共享基因） ====\n"); print(W, row.names = FALSE)

W2 <- reshape(R_[, c("signature","module","prho")], idvar = "signature", timevar = "module", direction = "wide")
names(W2) <- sub("prho\\.", "", names(W2))
cat("\n==== 偏相关（控制肿瘤上皮含量） ====\n")
print(W2[, c("signature","Neutrophil","NETosis","Stroma_ECM","Proliferation")], row.names = FALSE)

cat("\n==== 主导成分（全部 9 个模块中 |rho| 最大） ====\n")
for (i in seq_len(nrow(W))) {
  v <- as.numeric(W[i, key]); names(v) <- key
  cat(sprintf("  %-7s -> %-12s rho=%+.3f\n", W$signature[i],
              names(v)[which.max(abs(v))], v[which.max(abs(v))]))
}

cat("\n==== 模块间共变（Spearman rho） ====\n")
ms <- sapply(modules, function(g) colMeans(z[intersect(g, rownames(z)), , drop = FALSE]))
CM <- round(cor(ms, method = "spearman"), 2)
print(CM[c("Neutrophil","NETosis","Stroma_ECM","Epithelial","Proliferation"),
         c("Neutrophil","NETosis","Stroma_ECM","Epithelial","Proliferation","B_cells","T_cells","Myeloid")])

cat("\n==== NETs 执行分子在 TCGA 肿瘤中的可测性 ====\n")
exec <- c("PADI4","MPO","ELANE","AZU1","CTSG","PRTN3","FCGR3B","CXCR2")
det <- data.frame(gene = exec,
  median_TPM = sapply(exec, function(g) round(median(2^m[g, ] - 1), 2)),
  pct_detected = sapply(exec, function(g) round(100 * mean((2^m[g, ] - 1) > 1), 1)))
print(det, row.names = FALSE)
write.csv(det, "j4c_nets_detectability.csv", row.names = FALSE)

np <- sapply(sigs, function(g) sum(g %in% neutrophil_ref))
cat("\n==== 签名内中性粒特异基因含量 ====\n")
print(data.frame(signature = names(sigs), n = sapply(sigs, length), n_neut = np,
                 pct = round(100*np/sapply(sigs, length), 1)), row.names = FALSE)
cat("合计:", sum(np), "/", sum(sapply(sigs, length)), "=",
    round(100*sum(np)/sum(sapply(sigs, length)), 1), "%\n")
write.csv(data.frame(signature = names(sigs), n_genes = sapply(sigs, length),
                     n_neutrophil = np), "j4c_nets_content.csv", row.names = FALSE)
