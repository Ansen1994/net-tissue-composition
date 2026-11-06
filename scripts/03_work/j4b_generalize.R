# ==============================================================================
# 泛化验证（修正版）：剔除签名与参照模块的共享基因，排除自我相关；
# 以偏相关控制肿瘤上皮含量（纯度），分离"组织构成"与"中性粒活性"。
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
neutrophil_ref <- c("AZU1","ELANE","MPO","CTSG","LTF","PRTN3","DEFA3","DEFA4",
                    "CAMP","LCN2","CEACAM8","FCGR3B","CXCR2","S100A8","S100A9",
                    "PADI4","CSF3R","FPR1")

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
gn <- as.character(rowData(x)$gene_name)
ok <- !is.na(gn) & gn != ""
x <- x[ok, ]; gn <- gn[ok]
need_rows <- which(gn %in% all_need)
m <- log2(as.matrix(assay(x, "tpm_unstrand")[need_rows, , drop = FALSE]) + 1)
rg <- gn[need_rows]
o <- order(rowMeans(m), decreasing = TRUE); m <- m[o, ]; rg <- rg[o]
m <- m[!duplicated(rg), ]; rg <- rg[!duplicated(rg)]; rownames(m) <- rg
z <- t(scale(t(m)))
cat("样本:", ncol(z), " 基因:", nrow(z), "\n\n")

# ---- 中性粒相关基因在 bulk 中的可测性 ----
cat("==== 中性粒基因在 TCGA 肿瘤中的表达水平 ====\n")
chk <- intersect(c("FCGR3B","CXCR2","CSF3R","FPR1","S100A8","S100A9","MPO","ELANE","AZU1","PADI4"), rownames(m))
for (g in chk) {
  v <- 2^m[g, ] - 1
  cat(sprintf("  %-7s 中位TPM=%7.2f  零值比例=%5.1f%%  检出(TPM>1)=%5.1f%%\n",
              g, median(v), 100*mean(v < 0.01), 100*mean(v > 1)))
}

# ---- 偏 Spearman：以秩残差实现 ----
pscor <- function(a, b, ctrl) {
  ra <- residuals(lm(rank(a) ~ rank(ctrl)))
  rb <- residuals(lm(rank(b) ~ rank(ctrl)))
  ct <- suppressWarnings(cor.test(ra, rb, method = "pearson"))
  c(rho = unname(ct$estimate), p = ct$p.value)
}

sig_score <- sapply(sigs, function(g) colMeans(z[intersect(g, rownames(z)), , drop = FALSE]))
epi <- colMeans(z[intersect(modules$Epithelial, rownames(z)), , drop = FALSE])

rows <- list()
for (s in names(sigs)) {
  for (md in names(modules)) {
    shared <- intersect(sigs[[s]], modules[[md]])
    g_use  <- setdiff(modules[[md]], sigs[[s]])
    g_use  <- intersect(g_use, rownames(z))
    if (length(g_use) < 3) {
      rows[[length(rows)+1]] <- data.frame(signature = s, module = md,
        n_module = length(g_use), n_shared = length(shared),
        rho = NA, p = NA, prho_ctrlEpi = NA, p_partial = NA)
      next
    }
    sc <- colMeans(z[g_use, , drop = FALSE])
    ct <- suppressWarnings(cor.test(sig_score[, s], sc, method = "spearman", exact = FALSE))
    pc <- pscor(sig_score[, s], sc, epi)
    rows[[length(rows)+1]] <- data.frame(signature = s, module = md,
      n_module = length(g_use), n_shared = length(shared),
      rho = round(unname(ct$estimate), 3), p = signif(ct$p.value, 3),
      prho_ctrlEpi = round(pc["rho"], 3), p_partial = signif(pc["p"], 3))
  }
}
R <- do.call(rbind, rows)
write.csv(R, "j4b_signature_module_correlations.csv", row.names = FALSE)

key <- c("Neutrophil","NETosis","Stroma_ECM","Proliferation","Epithelial","B_cells","Endothelium","T_cells","Myeloid")
W <- reshape(R[, c("signature","module","rho")], idvar = "signature",
             timevar = "module", direction = "wide")
names(W) <- sub("rho\\.", "", names(W))
W <- W[, c("signature", key)]
cat("\n==== Spearman rho（已剔除共享基因） ====\n")
print(W, row.names = FALSE)

W2 <- reshape(R[, c("signature","module","prho_ctrlEpi")], idvar = "signature",
              timevar = "module", direction = "wide")
names(W2) <- sub("prho_ctrlEpi\\.", "", names(W2))
cat("\n==== 偏相关：控制肿瘤上皮含量后 ====\n")
print(W2[, c("signature", "Neutrophil", "NETosis", "Stroma_ECM", "Proliferation")], row.names = FALSE)

cat("\n==== 各签名的主导成分（|rho| 最大） ====\n")
for (s in W$signature) {
  v <- unlist(W[W$signature == s, key[-1], drop = TRUE])
  cat(sprintf("  %-7s -> %-12s (rho=%+.3f)\n", s, names(v)[which.max(abs(v))], v[which.max(abs(v))]))
}

cat("\n==== 摘要统计 ====\n")
np <- sapply(sigs, function(g) sum(g %in% neutrophil_ref))
cat("中性粒基因位点:", sum(np), "/", sum(sapply(sigs, length)),
    "=", round(100*sum(np)/sum(sapply(sigs,length)), 1), "%\n")
cat("含 >=1 中性粒基因的签名:", sum(np > 0), "/", length(sigs), "\n")
write.csv(data.frame(signature = names(sigs), n_genes = sapply(sigs, length),
                     n_neutrophil = np,
                     pct = round(100*np/sapply(sigs, length), 1)),
          "j4b_nets_content.csv", row.names = FALSE)
