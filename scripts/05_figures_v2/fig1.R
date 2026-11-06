# ==============================================================================
# Figure 1 — 候选池的构建（v2 重排版）
# A 火山图 | B 韦恩 | C 热图 + logFC | D 细胞来源构成 | E GO | F KEGG
# 配色统一 Okabe-Ito，尺寸双栏 170 mm
#
# 2026-11-06 改版要点：
#   1) A 面板 26 个超出纵轴上限的候选基因此前未被任何图层绘制（圆圈层与三角
#      层筛选条件互斥），现统一画在 y = YMAX 处，132 个候选全部可见；
#   2) A 面板补 Up / Down / n.s. / NET candidate 图例；
#   3) C 面板 ER 色条内标注分组与例数；
#   4) D 面板收紧 x 轴、放大基因名、修正浅色条内数字的颜色；
#   5) E 面板分面标签加底色框；F 面板改为按基因数排序（颜色仍表显著性）；
#   6) 面板字母与 y 轴标题改为锚定面板边缘的标注，避免被拼版列宽拉飞。
# ==============================================================================
source("figs_v2/style.R")
suppressMessages({ library(ggrepel) })

BLU <- "#0072B2"; GRN <- "#009E73"; ORG <- "#D55E00"
MAG <- "#CC79A7"; SKY <- "#56B4E9"; YEL <- "#E69F00"; GRY <- "#999999"

# ---------------------------------------------------------------- 数据
deg <- read.csv("Table_S1_All_DEGs.csv", stringsAsFactors = FALSE)
deg <- deg[!is.na(deg$padj) & deg$Real_Symbol != "", ]
deg$nlp <- -log10(pmax(deg$padj, 1e-300))

pool <- read.csv("fig1d_cell_origin_final.csv", stringsAsFactors = FALSE)
pool_genes <- pool$gene

# 面板字母锚定在「面板」左上角，而不是绘图区左上角。
# 用 x = -Inf / y = Inf 的标注 + clip = "off"：位置随面板走，不受左侧标签宽度影响。
tag_at_panel <- function(p, lab, use_coord = TRUE, size = 3.7) {
  if (use_coord) p <- p + coord_cartesian(clip = "off")
  p +
    annotate("text", x = -Inf, y = Inf, label = lab, hjust = 1, vjust = -0.5,
             size = size, fontface = "bold") +
    theme(plot.title = element_text(margin = ggplot2::margin(l = 11, b = 3)))
}

# 分面图专用：不带分面变量的 annotate 会在每个分面各画一次，这里把字母限定到 BP 分面
tag_at_panel_facet <- function(p, lab, ont = "BP", size = 3.7) {
  p + coord_cartesian(clip = "off") +
    geom_text(data = data.frame(ONTOLOGY = factor(ont, levels = c("BP", "CC", "MF")),
                                x = -Inf, y = Inf, lab = lab),
              aes(x = x, y = y, label = lab), inherit.aes = FALSE,
              hjust = 1, vjust = -0.5, size = size, fontface = "bold") +
    theme(plot.title = element_text(margin = ggplot2::margin(l = 11, b = 3)))
}

# ---------------------------------------------------------------- A 火山图
YMAX <- 110
deg$y_plot <- pmin(deg$nlp, YMAX)
deg$above  <- deg$nlp > YMAX
deg$grp <- ifelse(deg$change == "Not", "NS", deg$change)
deg$is_pool <- deg$Real_Symbol %in% pool_genes

pA <- ggplot(deg, aes(log2FoldChange, y_plot)) +
  geom_vline(xintercept = c(-1, 1), linetype = 3, linewidth = 0.3, colour = "grey55") +
  geom_hline(yintercept = -log10(0.05), linetype = 3, linewidth = 0.3, colour = "grey55") +
  # 非候选基因：上调 / 下调 / 不显著
  geom_point(data = deg[!deg$is_pool, ], aes(colour = grp),
             size = 0.55, stroke = 0, alpha = 0.55) +
  # 超出纵轴上限的非候选基因：灰色三角
  geom_point(data = deg[deg$above & !deg$is_pool, ], shape = 17, size = 0.7,
             colour = "grey55") +
  # 132 个候选基因：全部空心黑圈；超出上限者落在 y = YMAX 的同一条带上
  geom_point(data = deg[deg$is_pool, ], aes(shape = "Candidate"),
             size = 1.5, fill = NA, colour = "black", stroke = 0.35) +
  annotate("segment", x = c(5.95, 1.45, -4.15, -2.28), y = c(108, 50.5, 108, 93.5),
           xend = c(6.32, 2.15, -5.32, -1.60), yend = c(101.5, 57, 101.5, 87.5),
           linewidth = 0.2, colour = "grey45") +
  annotate("label", x = c(6.40, 2.22, -5.42, -1.52), y = c(100, 58.5, 100, 86),
           label = c("ESR1", "ERBB2", "MMP7", "ANGPT1"),
           size = 2.25, fontface = "italic", hjust = 0,
           fill = "white", alpha = 0.82, linewidth = 0, label.padding = unit(0.5, "mm")) +
  scale_colour_manual(values = c(Up = ORG, Down = BLU, NS = "#D8D8D8"),
                      breaks = c("Up", "Down", "NS"),
                      labels = c("Up in ER+", "Down in ER\u2212", "n.s."),
                      name = NULL) +
  scale_shape_manual(values = c(Candidate = 21), labels = c(Candidate = "NET candidate"),
                     name = NULL) +
  guides(colour = guide_legend(order = 1, nrow = 1,
                               override.aes = list(size = 1.2, stroke = 0, alpha = 1)),
         shape  = guide_legend(order = 2, nrow = 1,
                               override.aes = list(size = 1.5, stroke = 0.35, fill = NA))) +
  scale_x_continuous(limits = c(-7.2, 8.2), breaks = seq(-6, 8, 2)) +
  scale_y_continuous(limits = c(0, YMAX + 4), breaks = seq(0, 100, 25)) +
  # y 轴标题同样锚定面板左侧，避免被拼版列宽拉离坐标轴；
  # vjust = -1.85 与左内边距 6 pt 配合，把标题推到刻度标签左侧约 4 pt，
  # 原 -1.1 时标题右缘与「25」「75」只隔 1.4 pt，印刷后会贴在一起
  annotate("text", x = -Inf, y = 55,
           label = expression(-log[10]~"adjusted"~italic(P)),
           angle = 90, vjust = -1.85, size = 2.65) +
  labs(x = "log2 fold change (ER+ vs ER\u2212)",
       y = NULL,
       title = "Differential expression, ER+ vs ER\u2212") +
  theme_pub(7.6) +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.direction = "horizontal", legend.spacing.x = unit(1.6, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.key.height = unit(3, "mm"),
        plot.margin = margin(2, 2, 2, 6))

# ---------------------------------------------------------------- B 韦恩
circ <- function(cx, cy, r, n = 240) {
  th <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(th), y = cy + r * sin(th))
}
c1 <- circ(-0.50, 0, 0.92); c1$set <- "TCGA DEGs"
c2 <- circ( 0.50, 0, 0.92); c2$set <- "GeneCards NET"
nA <- length(unique(deg$Real_Symbol[deg$change != "Not"]))   # 12984
nB <- 509
nAB <- length(pool_genes)                                     # 132

pB <- ggplot() +
  geom_polygon(data = c1, aes(x, y, fill = set), alpha = 0.55, colour = BLU, linewidth = 0.4) +
  geom_polygon(data = c2, aes(x, y, fill = set), alpha = 0.55, colour = ORG, linewidth = 0.4) +
  scale_fill_manual(values = c("TCGA DEGs" = SKY, "GeneCards NET" = YEL), guide = "none") +
  annotate("text", x = -1.02, y = 0, label = format(nA - nAB, big.mark = ","), size = 2.9) +
  annotate("text", x =  0.00, y = 0, label = nAB, size = 3.2, fontface = "bold") +
  annotate("text", x =  1.02, y = 0, label = nB - nAB, size = 2.9) +
  annotate("text", x = -0.62, y = 1.16, label = sprintf("TCGA DEGs\nn = %s", format(nA, big.mark = ",")),
           size = 2.1, colour = BLU, lineheight = 0.92, vjust = 0) +
  annotate("text", x =  0.62, y = 1.16, label = sprintf("GeneCards NET\nn = %d", nB),
           size = 2.1, colour = "#B87A00", lineheight = 0.92, vjust = 0) +
  coord_equal(xlim = c(-1.55, 1.55), ylim = c(-1.02, 1.62), clip = "off") +
  labs(title = "Intersection with NET-related entries") +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 8.5, hjust = 0,
                                  margin = margin(l = 13, b = 2)),
        plot.title.position = "plot",
        plot.tag.position = c(0, 1),
        plot.margin = margin(2, 2, 2, 2))

# ---------------------------------------------------------------- C 热图
hm <- as.matrix(read.csv("Table_S2_Target_NETs_Matrix.csv", row.names = 1, check.names = FALSE))
hm <- hm[!grepl("\\.1$", rownames(hm)), , drop = FALSE]
er <- read.csv("fig1c_ER_status.csv", stringsAsFactors = FALSE)
er <- er[match(colnames(hm), er$sample), ]
ord_col <- order(factor(er$ER, levels = c("Negative", "Positive")))
hm <- hm[, ord_col, drop = FALSE]
er <- er[ord_col, ]

n_neg <- sum(er$ER == "Negative")
n_pos <- sum(er$ER == "Positive")

z <- t(scale(t(hm)))
z[z >  2] <-  2; z[z < -2] <- -2
hc <- hclust(dist(z), method = "complete")
row_ord <- hc$order
z <- z[row_ord, , drop = FALSE]

hm_long <- data.frame(
  gene = factor(rownames(z), levels = rev(rownames(z))),
  col  = rep(seq_len(ncol(z)), each = nrow(z)),
  val  = as.vector(z)
)
# 纵轴改用连续数值（行号），不再用离散因子：色条要在热图上方占掉一段真实的
# 纵轴空间，右侧柱状图必须用同一套 y 映射，132 行才会逐行对齐。
hm_long$y <- as.numeric(hm_long$gene)

# ER 分组色条画在热图自己的坐标系里，占纵轴最上方的 11 个单位（约 6.4 pt）。
# 早先把它做成独立子图纵向堆叠，热图顶部为对齐柱状图必须留白，色条就浮在
# 热图上方 5.2 pt 处，像断掉的表头。
ER_TOP  <- 143.5     # 纵轴上限：132 行 + 色条 11 个单位
ER_CUT  <- n_neg + 0.5
ER_BAR_LO <- 133.0
ER_LAB_Y  <- 138.2

pC_heat <- ggplot(hm_long, aes(col, y, fill = val)) +
  geom_raster() +
  # 常量填充，不干扰 z-score 的连续色标
  annotate("rect", xmin = -Inf, xmax = ER_CUT, ymin = ER_BAR_LO, ymax = ER_TOP,
           fill = SKY, colour = NA) +
  annotate("rect", xmin = ER_CUT, xmax = Inf, ymin = ER_BAR_LO, ymax = ER_TOP,
           fill = MAG, colour = NA) +
  # 色条内直接标注分组与例数，取代原先缺失的图例
  annotate("text", x = 6, y = ER_LAB_Y, label = sprintf("ER\u2212 (n = %d)", n_neg),
           hjust = 0, size = 1.9, colour = "white", fontface = "bold") +
  annotate("text", x = n_neg + 8, y = ER_LAB_Y, label = sprintf("ER+ (n = %d)", n_pos),
           hjust = 0, size = 1.9, colour = "white", fontface = "bold") +
  scale_fill_gradient2(low = BLU, mid = "#F7F7F7", high = ORG, midpoint = 0,
                       limits = c(-2, 2), name = "Row z-score",
                       breaks = c(-2, 0, 2), labels = c("\u2264-2", "0", "\u22652")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, ER_TOP), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_pub(7) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "right",
        legend.key.width = unit(2.2, "mm"), legend.key.height = unit(5, "mm"),
        # 底部 4.2 pt 抵掉右侧柱状图 x 轴刻度占的高度，两侧面板才等高
        plot.margin = margin(0, 2, 4.2, 2))

# 右侧 logFC 条
lfc <- deg[!duplicated(deg$Real_Symbol), c("Real_Symbol", "log2FoldChange")]
names(lfc)[1] <- "gene"
lfc <- lfc[lfc$gene %in% rownames(z), ]
lfc$gene <- factor(lfc$gene, levels = levels(hm_long$gene))
lfc$y <- as.numeric(lfc$gene)
pC_lfc <- ggplot(lfc, aes(log2FoldChange, y)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
  # 纵轴为连续值后 ggplot 无法自动判定方向，必须显式 orientation = "y"
  geom_col(aes(fill = log2FoldChange > 0), width = 0.9, orientation = "y") +
  scale_fill_manual(values = c("TRUE" = ORG, "FALSE" = BLU), guide = "none") +
  # 上限放到 ±6：最极端的 ESR1 为 +5.72，原来的 ±4.2 会让 5 根条（ESR1、CXCL5、
  # BPI、COL11A2、S100A8）整根消失，读者只会看到空白行
  scale_x_continuous(limits = c(-6, 6), breaks = c(-4, 0, 4)) +
  # 与热图同一套 y 映射（上限含色条占位），132 根条才与热图 132 行严格对齐
  scale_y_continuous(limits = c(0.5, ER_TOP), expand = c(0, 0)) +
  labs(x = expression(log[2]~FC), y = NULL) +
  # 标题改用不占版面高度的标注（clip = "off"），否则柱状图面板会比热图矮 12 pt
  annotate("text", x = 0, y = Inf, label = "ER+ vs ER\u2212",
           vjust = -0.2, size = 2.3, colour = "black") +
  coord_cartesian(clip = "off") +
  theme_pub(6.6) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = 6.3, margin = margin(t = 1.2)),
        plot.margin = margin(0, 2, 0, 2))

pC_inner <- pC_heat | pC_lfc
pC_inner <- pC_inner + plot_layout(widths = c(1, 0.13))
pC <- wrap_elements(full = pC_inner) + labs(tag = "C") +
  # tag 的 npc 是相对面板的：1.02 只抬高约 3 pt，字母下半截仍落在面板内，会压住
  # ER− 色条。抬到 1.065，字母底边距色条约 4.5 pt，与 A 面板的留白接近。
  theme(plot.tag = element_text(face = "bold", size = 11),
        plot.tag.position = c(0, 1.065))

# ---------------------------------------------------------------- D 细胞来源
cls_col <- c("Neutrophil-restricted" = ORG, "Epithelial" = BLU, "Myeloid" = GRN,
             "Stromal" = MAG, "Endothelial" = SKY, "Lymphoid" = YEL,
             "Broadly expressed" = GRY, "Low / not quantified" = "#CFCFCF")
cls_lab <- c("Neutrophil-restricted", "Epithelial", "Myeloid", "Stromal",
             "Endothelial", "Lymphoid", "Broadly expressed", "Not quantified")
d_cnt <- as.data.frame(table(pool$class), stringsAsFactors = FALSE)
names(d_cnt) <- c("cls", "n")
d_cnt$cls <- factor(d_cnt$cls, levels = rev(names(cls_col)))
d_cnt$lab <- cls_lab[match(as.character(d_cnt$cls), names(cls_col))]
d_cnt$lab <- factor(d_cnt$lab, levels = rev(cls_lab))
ex <- sapply(as.character(d_cnt$cls), function(k) {
  g <- pool$gene[as.character(pool$class) == k]
  if (k == "Neutrophil-restricted") g <- c("ELANE","MPO","PADI4")
  paste(head(g, 3), collapse = ", ")
})

pD <- ggplot(d_cnt, aes(n, lab, fill = cls)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = n, x = n - 0.55), hjust = 1, size = 2.4,
            # 浅色条（天蓝、黄、灰）上用深色字：白字压这些底色对比度只有 2.0–2.8:1
            colour = ifelse(d_cnt$cls %in% c("Low / not quantified", "Broadly expressed",
                                             "Endothelial", "Lymphoid"),
                            "grey15", "white")) +
  geom_text(aes(label = ex, x = n + 0.8), hjust = 0, size = 2.3,
            fontface = "italic", colour = "grey30") +
  scale_fill_manual(values = cls_col, guide = "none") +
  # 原上限 46 让最长的 29 只占约 2/3 画幅，收紧到 41 同时保证基因名不出界
  scale_x_continuous(limits = c(0, 41), expand = c(0, 0)) +
  labs(x = "Number of candidate genes", y = NULL,
       title = "Cell of origin of the 132 candidates (single-cell mapping)") +
  theme_pub(7.4) +
  # 行独立排版后 D 的标签会贴到纸张左缘，留出一点左边距
  theme(panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey92"),
        plot.margin = margin(4, 6, 4, 10))

# ---------------------------------------------------------------- E GO
go <- read.csv("fig1_go_enrich.csv", stringsAsFactors = FALSE)
go$ONTOLOGY <- factor(go$ONTOLOGY, levels = c("BP", "CC", "MF"))
go <- do.call(rbind, lapply(c("BP", "CC", "MF"), function(o) head(go[go$ONTOLOGY == o, ], 5)))
go$nlp <- -log10(go$pvalue)
go$ratio <- sapply(strsplit(go$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
wrap <- function(x, n = 44) sapply(strwrap(x, n, simplify = FALSE), paste, collapse = "\n")
go$desc <- wrap(go$Description)
go$desc <- factor(go$desc, levels = rev(unique(go$desc[order(go$ONTOLOGY, go$ratio)])))

pE <- ggplot(go, aes(ratio, desc)) +
  geom_point(aes(size = Count, colour = nlp)) +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradient(low = SKY, high = ORG, name = expression(-log[10]~italic(P))) +
  scale_size_continuous(range = c(0.6, 3.0), name = "Count") +
  scale_x_continuous(limits = c(0.10, 0.30), breaks = c(0.15, 0.20, 0.25),
                     labels = c("0.15", "0.2", "0.25")) +
  labs(x = "Gene ratio", y = NULL, title = "GO enrichment") +
  theme_pub(7.2) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 6.8),
        # 分组标签加浅底框，避免读成某条通路的数据标签
        strip.background = element_rect(fill = "grey94", colour = NA),
        panel.spacing.y = unit(2.2, "mm"),
        legend.position = "right",
        legend.key.height = unit(3.2, "mm"))

# ---------------------------------------------------------------- F KEGG
kg <- read.csv("fig1_kegg_enrich.csv", stringsAsFactors = FALSE)
kg <- head(kg, 10)
kg$nlp <- -log10(kg$pvalue)
kg$desc <- wrap(kg$Description, 42)
# 条长表示基因数，就按基因数排序；显著性仍由颜色表示，避免锯齿状长条
kg <- kg[order(-kg$Count, kg$pvalue), ]
kg$desc <- factor(kg$desc, levels = rev(kg$desc))

pF <- ggplot(kg, aes(Count, desc)) +
  geom_col(aes(fill = nlp), width = 0.72) +
  geom_text(aes(label = Count), hjust = -0.3, size = 2.3) +
  scale_fill_gradient(low = SKY, high = ORG, name = expression(-log[10]~italic(P))) +
  scale_x_continuous(limits = c(0, 33), expand = c(0, 0)) +
  labs(x = "Gene count", y = NULL, title = "KEGG enrichment") +
  theme_pub(7.2) +
  theme(axis.text.y = element_text(size = 6.2, lineheight = 0.9),
        legend.key.height = unit(3.2, "mm"))

# ---------------------------------------------------------------- 拼版
# 每行用 wrap_elements 独立排版：否则 patchwork 会把四行面板左边缘强行对齐，
# 由 E 面板最长的 GO 名称撑出约 150 pt 的空白，A / D 的轴标题与字母被拉飞。
row1 <- wrap_elements(full = (tag_at_panel(pA, "A") |
                                tag_panel(pB, "B")) +
                        plot_layout(widths = c(1.06, 1)))
row2 <- pC
row3 <- wrap_elements(full = tag_at_panel(pD, "D"))
row4 <- wrap_elements(full = (tag_at_panel_facet(pE, "E") |
                                tag_at_panel(pF, "F")) +
                        plot_layout(widths = c(1.18, 1)))

final <- row1 / row2 / row3 / row4 +
  plot_layout(heights = c(1.42, 1.10, 0.95, 1.62))

save_fig(final, "Figure1_ABCDEF", W2, 8.95)
