# ==============================================================================
# Figure 7  A–D  六个已发表 NET 签名 + 本研究签名的泛化检验
#   合计 7 个签名 / 41 位点（已发表 37 + 本文 4）；8/41 = 19.5%，8/37 = 21.6%
#   签名与参考文献编号的对应见 Table S13a 的 Reference 列：
#     BRCA_A [14], BRCA_B [15], BRCA_C [16], GC_A [17], GC_B [18], CRC_A [19]
#   A 各签名中中性粒限制基因的比例
#   B 签名基因构成（中性粒基因标红）
#   C NETs 执行分子在 1,111 例 TCGA-BRCA 中的可测性
#   D 各签名与 9 个泛化模块的相关（已剔除共享基因；不同于 Methods 2.5 的 11 模块）
# ==============================================================================
source("figs_v2/style.R")
suppressMessages(library(ggrepel))

cnt <- read.csv("j4c_nets_content.csv")
det <- read.csv("j4c_nets_detectability.csv")
cor <- read.csv("j4c_signature_module_correlations.csv")

sigs <- list(
  BRCA_A = c("LTF","LCP1","AZU1","ENO1"),
  BRCA_B = c("F2RL2","AZU1","IL33","ELANE"),
  BRCA_C = c("ZMYND10","IL12B","CXCL13","TFF1","LTB","SPIB"),
  GC_A   = c("PDE4B","CD93","CTSG","IL6","KCNJ15","CRISPLD2","ELANE"),
  GC_B   = c("MPO","BST1","IL6","KCNJ15","SELP","ELANE"),
  CRC_A  = c("AKT2","ATG7","CLEC6A","CRISPLD2","ENTPD4","IL17A","ITGAM",
             "MAPK1","SELPLG","SIGLEC14"),
  Ours   = c("PLG","CD79A","MMP7","ANGPT1"))
neut <- c("AZU1","ELANE","MPO","CTSG","LTF","PRTN3","DEFA3","DEFA4","CAMP","LCN2",
          "CEACAM8","FCGR3B","CXCR2","S100A8","S100A9","PADI4","CSF3R","FPR1")
sig_lab <- c(BRCA_A = "BRCA_A", BRCA_B = "BRCA_B", BRCA_C = "BRCA_C",
             GC_A = "GC_A", GC_B = "GC_B", CRC_A = "CRC_A", Ours = "This study")

cnt$frac <- cnt$n_neutrophil / cnt$n_genes * 100
cnt$signature <- factor(cnt$signature, levels = rev(names(sigs)))

# ---------------- A ----------------
pA <- ggplot(cnt, aes(frac, signature)) +
  geom_col(aes(fill = signature == "Ours"), width = 0.66,
           colour = "black", linewidth = 0.2, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d/%d", n_neutrophil, n_genes)),
            hjust = -0.25, size = 2.1) +
  scale_fill_manual(values = c("FALSE" = "#BDBDBD", "TRUE" = "#B2182B")) +
  scale_y_discrete(labels = rev(sig_lab)) +
  scale_x_continuous(limits = c(0, 68), breaks = c(0, 25, 50),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Neutrophil-restricted genes (%)", y = NULL,
       title = "8 of 41 loci are neutrophil-restricted",
       subtitle = "6 published + this study; 8/37 = 21.6% published-only") +
  theme_pub(8) +
  theme(axis.text.y = element_text(size = 6.8))

# ---------------- B ----------------
gm <- do.call(rbind, lapply(names(sigs), function(s) {
  g <- sigs[[s]]
  data.frame(signature = s, idx = seq_along(g), gene = g,
             cls = ifelse(g %in% neut, "Neutrophil-restricted", "Other"))
}))
gm$signature <- factor(gm$signature, levels = rev(names(sigs)))

pB <- ggplot(gm, aes(idx, signature)) +
  geom_tile(aes(fill = cls), colour = "white", linewidth = 0.7, width = 0.94, height = 0.86) +
  geom_text(aes(label = gene, colour = cls), size = 2.1) +
  scale_fill_manual(values = c("Neutrophil-restricted" = "#F4C7C3",
                               "Other" = "#EDEDED"), name = NULL) +
  scale_colour_manual(values = c("Neutrophil-restricted" = "#8C1D18",
                                 "Other" = "#333333"), guide = "none") +
  scale_y_discrete(labels = rev(sig_lab)) +
  scale_x_continuous(breaks = 1:10, expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = NULL, y = NULL, title = "Signature composition") +
  theme_pub(8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.text.y = element_text(size = 6.8),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(2.6, "mm"),
        legend.margin = margin(t = -4, b = 0))

# ---------------- C ----------------
det$highlight <- det$gene %in% c("PADI4","MPO","ELANE")
pC <- ggplot(det, aes(median_TPM, pct_detected)) +
  geom_hline(yintercept = 10, linetype = "dotted", linewidth = 0.3, colour = "grey60") +
  geom_point(aes(colour = highlight, size = highlight), show.legend = FALSE) +
  geom_text_repel(aes(label = gene, colour = highlight), size = 2.0,
                  show.legend = FALSE, fontface = "italic",
                  box.padding = 0.28, point.padding = 0.18,
                  min.segment.length = 0.2, segment.size = 0.2,
                  segment.colour = "grey55", max.overlaps = 20, seed = 1) +
  scale_colour_manual(values = c("FALSE" = "grey45", "TRUE" = "#B2182B")) +
  scale_size_manual(values = c("FALSE" = 1.6, "TRUE" = 2.4)) +
  scale_x_log10(limits = c(0.04, 3.2), breaks = c(0.05, 0.1, 0.5, 1, 3),
                labels = c("0.05","0.1","0.5","1","3")) +
  scale_y_continuous(limits = c(0, 78), breaks = seq(0, 60, 20)) +
  labs(x = "Median TPM (log scale)", y = "Tumours above 1 TPM (%)",
       title = "Seven of eight effectors fall below 1 TPM",
       subtitle = "TCGA-BRCA, 1,111 primary tumours; CTSG is the exception") +
  theme_pub(8)

# ---------------- D ----------------
mod_lab <- c(Neutrophil = "Neutrophil", NETosis = "NETosis", Stroma_ECM = "Stroma",
             Endothelium = "Endothelium", Epithelial = "Epithelial",
             Proliferation = "Proliferation", T_cells = "T cells",
             B_cells = "B cells", Myeloid = "Myeloid")
cor$module <- factor(mod_lab[cor$module], levels = mod_lab)
cor$signature <- factor(cor$signature, levels = rev(names(sigs)))
cor$lab <- sprintf("%.2f", cor$rho)

pD <- ggplot(cor, aes(module, signature, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = lab, colour = abs(rho) > 0.55), size = 2.0, show.legend = FALSE) +
  scale_fill_gradientn(colours = div_cols, limits = c(-0.8, 0.8),
                       breaks = seq(-0.8, 0.8, 0.4), name = "Spearman rho") +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white")) +
  scale_y_discrete(labels = rev(sig_lab)) +
  labs(x = NULL, y = NULL,
       title = "No shared biological reading across signatures",
       subtitle = "correlation with compartment modules; genes shared between a signature and a module removed") +
  theme_pub(8) +
  theme(axis.text.y = element_text(size = 7.0),
        axis.text.x = element_text(size = 6.8, angle = 30, hjust = 1),
        legend.key.height = unit(7, "mm"),
        legend.key.width = unit(2.6, "mm"))

# ---------------- 拼合 ----------------
pA <- tag_panel(pA, "A"); pB <- tag_panel(pB, "B")
pC <- tag_panel(pC, "C"); pD <- tag_panel(pD, "D")

row3 <- plot_spacer() + pD + plot_spacer() +
  plot_layout(widths = c(0.14, 1, 0.14))

fig <- (pA + pC + plot_layout(widths = c(0.78, 1.02))) /
  pB / row3 + plot_layout(heights = c(1, 0.92, 1.2))

save_fig(fig, "Figure7_ABCD", W2, 7.2)
