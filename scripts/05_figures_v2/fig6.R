# ==============================================================================
# Figure 6  A–F
#   A 四个 hub 基因在 9 类细胞中的表达与检出率（点图）
#   B 肿瘤上皮 vs 主导细胞区室的表达对比
#   C 各基因表达量在细胞区室间的来源占比（哪一类细胞贡献了 bulk 信号）
#   D NETs 执行分子在 ER+ 细胞中的检出率
#   E 评分与基质 / 增殖 / 上皮纯度模块的相关（四队列）
#   F 扣除 CD79A 前后的免疫模块相关
# ==============================================================================
source("figs_v2/style.R")

ce <- read.csv("j3_celltype_expression.csv", stringsAsFactors = FALSE)
nd <- read.csv("j3_nets_detection.csv")
ica <- read.csv("Table_ImmuneAxis_AllCorrelations.csv")
scd <- read.csv("Table_ImmuneAxis_ScoreDecomposition.csv")

ct_order <- c("Cancer Epithelial","Normal Epithelial","CAFs","PVL","Endothelial",
              "Myeloid","T-cells","B-cells","Plasmablasts")
ct_lab <- c("Cancer epithelial","Normal epithelial","CAFs","Perivascular",
            "Endothelial","Myeloid","T cells","B cells","Plasmablasts")
genes <- c("PLG","CD79A","MMP7","ANGPT1")

ce$major <- factor(ce$major, levels = rev(ct_order))
ce$gene  <- factor(ce$gene, levels = genes)
ce$det_pct <- ce$detect * 100
n_tot <- sum(ce %>% distinct(major, n) %>% pull(n))

# ---------------- A ----------------
pA <- ggplot(ce, aes(gene, major)) +
  geom_point(aes(size = det_pct, colour = mean_expr)) +
  geom_text(data = ce %>% group_by(gene) %>% slice_max(det_pct, n = 1) %>%
              ungroup() %>% filter(det_pct >= 20),
            aes(label = sprintf("%.0f%%", det_pct)), size = 1.9,
            colour = "black", hjust = -0.55) +
  scale_size_area(max_size = 4.0, limits = c(0, 70), breaks = c(1, 10, 40),
                  name = "Detected (%)") +
  scale_colour_gradientn(colours = c("#EFEFEF", "#9ECAE1", "#2171B5", "#08306B"),
                         limits = c(0, 1.7), name = "Mean expr") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(labels = rev(ct_lab)) +
  labs(x = NULL, y = NULL,
       title = "Cellular origin of the four hub genes",
       subtitle = sprintf("GSE176078, %s ER+ cells", format(n_tot, big.mark = ","))) +
  theme_pub(8) +
  theme(axis.text.x.top = element_text(face = "italic", size = 7.0),
        axis.text.y = element_text(size = 6.8),
        legend.position = "right", legend.box = "vertical",
        legend.spacing.y = unit(0.5, "mm"),
        legend.key.height = unit(3.6, "mm"),
        legend.key.width = unit(2.2, "mm"))

# ---------------- B ----------------
dom <- ce %>% filter(major != "Cancer Epithelial") %>%
  group_by(gene) %>% slice_max(mean_expr, n = 1) %>% ungroup() %>%
  transmute(gene, compartment = as.character(major), expr = mean_expr)
tum <- ce %>% filter(major == "Cancer Epithelial") %>%
  transmute(gene, compartment = "Cancer epithelium", expr = mean_expr)
bb <- rbind(tum, dom)
bb$gene <- factor(bb$gene, levels = rev(genes))
bb$grp  <- factor(ifelse(bb$compartment == "Cancer epithelium", "Tumour", "Dominant"),
                  levels = c("Tumour","Dominant"))
fold <- bb %>% group_by(gene) %>% summarise(f = max(expr) / min(expr), .groups = "drop")

pB <- ggplot(bb, aes(expr, gene)) +
  geom_line(aes(group = gene), colour = "grey70", linewidth = 0.35) +
  geom_point(aes(colour = grp, shape = grp), size = 2.0) +

  scale_x_log10(limits = c(6e-5, 3), breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1),
                labels = c("0.0001","0.001","0.01","0.1","1")) +
  scale_y_discrete(labels = sprintf("%s  (%.0fx)", rev(genes), fold$f[match(rev(genes), fold$gene)])) +
  scale_colour_manual(values = c(Tumour = "#B2182B", Dominant = "#2166AC"), name = NULL) +
  scale_shape_manual(values = c(Tumour = 17, Dominant = 16), name = NULL) +
  labs(x = "Mean expression (log scale)", y = NULL,
       title = "Tumour vs dominant",
       subtitle = "fold difference in brackets") +
  theme_pub(8) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(2.6, "mm"),
        legend.margin = margin(t = -4, b = 0),
        axis.text.y = element_text(size = 7.0, face = "italic"))

# ---------------- C ----------------
mass <- ce %>%
  mutate(group = dplyr::case_when(
    major == "Cancer Epithelial" ~ "Tumour",
    major == "Normal Epithelial" ~ "Normal epithelium",
    major %in% c("CAFs","PVL","Endothelial") ~ "Stroma / vasculature",
    TRUE ~ "Immune")) %>%
  group_by(gene, group) %>% summarise(m = sum(mean_expr * n), .groups = "drop") %>%
  group_by(gene) %>% mutate(frac = m / sum(m) * 100) %>% ungroup()
mass$group <- factor(mass$group, levels = c("Tumour","Normal epithelium",
                                            "Stroma / vasculature","Immune"))
# PLG is near zero in every compartment: percentages of a near-zero total are
# not interpretable, so grey out its bar and label it N/A instead of 12%.
tum_frac <- mass %>% filter(group == "Tumour") %>%
  mutate(interpretable = gene != "PLG",
         lab = ifelse(interpretable, sprintf("%.0f%%", frac), "N/A"),
         x = ifelse(interpretable, 101.5, 50),
         hj = ifelse(interpretable, 0, 0.5),
         tc = ifelse(interpretable, "black", "grey25"))

pC <- ggplot(mass, aes(frac, gene, fill = group)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.4) +
  geom_rect(data = mass %>% filter(gene == "PLG"),
            aes(xmin = 0, xmax = 100,
                ymin = as.numeric(gene) - 0.34,
                ymax = as.numeric(gene) + 0.34),
            inherit.aes = FALSE, fill = "grey90", colour = "grey55",
            linewidth = 0.3) +
  geom_text(data = tum_frac %>% filter(interpretable),
            aes(x = x, y = gene, label = lab, hjust = hj,
                colour = tc), inherit.aes = FALSE,
            size = 2.2, fontface = "bold", show.legend = FALSE) +
  geom_text(data = tum_frac %>% filter(!interpretable),
            aes(x = x, y = gene, label = lab, hjust = hj), inherit.aes = FALSE,
            size = 2.2, colour = "grey25", show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_manual(values = c(Tumour = "#B2182B", "Normal epithelium" = "#4393C3",
                               "Stroma / vasculature" = "#B8A55D",
                               Immune = "#7B7B7B"), name = NULL) +
  scale_x_continuous(limits = c(0, 112), breaks = c(0, 50, 100),
                     expand = expansion(mult = c(0, 0.005))) +
  labs(x = "Contribution to total expression (%)", y = NULL,
       title = "Source of the bulk signal",
       subtitle = "weighted by cell number") +
  theme_pub(8) +
  theme(axis.text.y = element_text(size = 7.0, face = "italic"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 5.8),
        legend.margin = margin(t = -4, b = 0))

# ---------------- D ----------------
nd$gene <- factor(nd$gene, levels = nd$gene[order(nd$detect_ERp)])
nd$pct  <- nd$detect_ERp * 100
pD <- ggplot(nd, aes(pct, gene)) +
  geom_segment(aes(x = 0.04, xend = pct, yend = gene), colour = "grey70", linewidth = 0.4) +
  geom_point(size = 2.1, colour = "#6A51A3") +
  geom_text(aes(x = 45, label = sprintf("%.2f", pct)), hjust = 0, size = 2.0) +
  scale_x_log10(limits = c(0.04, 90), breaks = c(0.1, 1, 10),
                labels = c("0.1","1","10")) +
  labs(x = "Detected cells in ER+ tumours (%, log scale)", y = NULL,
       title = "Most effectors sit at the detection floor",
       subtitle = "PADI4, MPO, ELANE and CXCR2 below 0.25%") +
  theme_pub(8) +
  theme(axis.text.y = element_text(size = 7.0, face = "italic"))

# ---------------- E ----------------
e5 <- subset(ica, analysis == "primary" & variable %in% c("Stroma","Prolif","Purity_epi"))
e5$variable <- factor(e5$variable, levels = c("Stroma","Prolif","Purity_epi"),
                      labels = c("Stroma","Proliferation","Epithelial purity"))
e5$cohort <- factor(e5$cohort, levels = rev(c("TCGA","GSE20685","METABRIC","GSE96058")))
e5$sig <- ifelse(e5$p_BH < 0.05, "p < 0.05", "n.s.")

pE <- ggplot(e5, aes(rho, cohort, colour = variable, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_point(size = 2.0, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(Stroma = "#8C6BB1", Proliferation = "#D7301F",
                                 "Epithelial purity" = "#2166AC"), name = NULL) +
  scale_shape_manual(values = c("p < 0.05" = 16, "n.s." = 1), name = NULL) +
  ## 显著性用点形区分，图例已由图注承担（open symbols = n.s. after BH），
  ## 面板宽度放不下第二个图例，去掉以免与变量名重叠
  guides(shape = "none") +
  scale_x_continuous(limits = c(-0.45, 0.38), breaks = seq(-0.4, 0.2, 0.2)) +
  labs(x = "Spearman rho with risk score", y = NULL,
       title = "What the score tracks",
       subtitle = "same direction in 12 of 12 pairs") +
  theme_pub(8) +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.direction = "horizontal",
        legend.key.size = unit(2.6, "mm"),
        legend.margin = margin(t = -4, b = 0),
        legend.spacing.y = unit(0.4, "mm"),
        axis.text.y = element_text(size = 6.8))

# ---------------- F ----------------
mm <- c("Bcell_noCD79A","CD8_T","Treg","Cytotoxic")
f6 <- subset(scd, module %in% mm)
f6$module <- factor(f6$module, levels = mm,
                    labels = c("B cells","CD8 T cells","Treg","Cytotoxic"))
f6 <- f6 %>% select(cohort, module, rho_full, rho_noCD79A) %>%
  pivot_longer(c(rho_full, rho_noCD79A), names_to = "version", values_to = "rho") %>%
  mutate(version = factor(version, levels = c("rho_full","rho_noCD79A"),
                          labels = c("Full","no CD79A")),
         cohort = factor(cohort, levels = c("TCGA","GSE20685","METABRIC","GSE96058")))

pF <- ggplot(f6, aes(version, rho, colour = cohort, group = cohort)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_line(linewidth = 0.45, alpha = 0.85) +
  geom_point(size = 1.7) +
  facet_wrap(~ module, nrow = 1) +
  scale_colour_manual(values = COH, name = NULL) +
  scale_y_continuous(limits = c(-0.45, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  labs(x = NULL, y = "Spearman rho with module",
       title = "The immune-cold signal is arithmetic",
       subtitle = "removing CD79A removes the negative association in every cohort") +
  theme_pub(8) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(2.6, "mm"),
        legend.margin = margin(t = -4, b = 0),
        axis.text.x = element_text(size = 5.9, angle = 45, hjust = 1),
        strip.text = element_text(size = 7.2),
        panel.spacing = unit(3, "mm"))

# ---------------- 拼合 ----------------
pA <- tag_panel(pA, "A"); pB <- tag_panel(pB, "B"); pC <- tag_panel(pC, "C")
pD <- tag_panel(pD, "D"); pE <- tag_panel(pE, "E"); pF <- tag_panel(pF, "F")

fig <- ((pA + pB + plot_layout(widths = c(1.25, 1))) /
  (pC + pD + plot_layout(widths = c(1, 1))) /
  (pE + pF + plot_layout(widths = c(0.82, 1.35)))) +
  plot_layout(heights = c(1.18, 1, 1.05))

save_fig(fig, "Figure6_ABCDEF", W2, 7.0)
