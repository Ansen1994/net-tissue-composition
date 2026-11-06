setwd('G:\\Project\\nets_breast canner_1')
# =====================================================================
# Figure 1 完整分析流水线：差异分析、GeneCards交集与功能聚类
# 逻辑流：Fig 1A(火山图) -> Fig 1B(韦恩图) -> Fig 1C(热图) -> Fig 1D(GO) -> Fig 1E(KEGG)
# =====================================================================

# 1. 核心包加载
setwd('G:\\Project\\nets_breast canner_1')
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ggVennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)

# 2. 数据载入与临床分组 (ER+ vs ER-)
load("brca.gdc_2024.rda") 
exp_mat <- assay(expquery2)
clin_info <- as.data.frame(colData(expquery2))

clin_er <- clin_info %>%
  filter(paper_BRCA_Subtype_PAM50 %in% c("LumA", "LumB", "Basal")) %>%
  mutate(ER_status = case_when(
    paper_BRCA_Subtype_PAM50 %in% c("LumA", "LumB") ~ "Positive",
    paper_BRCA_Subtype_PAM50 == "Basal" ~ "Negative"
  )) %>%
  mutate(ER_status = factor(ER_status, levels = c("Negative", "Positive"))) %>%
  filter(!is.na(sample_submitter_id))

shared_samples <- intersect(colnames(exp_mat), rownames(clin_er))
exp_mat_er <- exp_mat[, shared_samples]
clin_er <- clin_er[shared_samples, ]

# 3. 运行 DESeq2 差异分析
dds <- DESeqDataSetFromMatrix(countData = round(exp_mat_er),
                              colData = clin_er,
                              design = ~ ER_status)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = c("ER_status", "Positive", "Negative"))

res_df <- as.data.frame(res) %>%
  rownames_to_column("GeneSymbol") %>% 
  filter(!is.na(padj)) %>%
  mutate(change = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "Not"
  ))

# ID 转换 (Ensembl -> Symbol)
gene_info <- as.data.frame(rowData(expquery2))
res_df_mapped <- res_df %>%
  mutate(Real_Symbol = gene_info$gene_name[match(GeneSymbol, rownames(gene_info))]) %>%
  filter(!is.na(Real_Symbol))
DEGs_real <- res_df_mapped %>% filter(change != "Not") %>% pull(Real_Symbol)

# 4. 导入 GeneCards 数据并取交集
genecards_data <- read.csv("GeneCards_NETs.csv", header = FALSE)
nets_genes_new <- as.character(genecards_data$V1)
target_NETs_genes <- intersect(DEGs_real, nets_genes_new)

# 输出辅助表格
write.csv(res_df_mapped, file = "Table_S1_All_DEGs.csv", row.names = FALSE)
target_data <- res_df_mapped[res_df_mapped$Real_Symbol %in% target_NETs_genes, ]
heat_mat <- as.matrix(log2(exp_mat_er[target_data$GeneSymbol, ] + 1))
rownames(heat_mat) <- make.names(target_data$Real_Symbol, unique = TRUE)
write.csv(heat_mat, file = "Table_S2_Target_NETs_Matrix.csv", row.names = TRUE)

# ================= 绘图模块 =================

# ==============================================================================
# Figure 1 顶刊马卡龙色系全套绘图代码 (按 1A-1E 顺序严密排版)
# ==============================================================================

# ==============================================================================
# Figure 1 顶刊排版级绘图代码 (完整保留原图逻辑，全局统一 Sans 字体)
# ==============================================================================

library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ggVennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
library(circlize)
library(ComplexHeatmap)
library(grid)
library(RColorBrewer)

# 统一提取 6 色标准色值
macaque_colors <- c(
  yellow = "#f9eba7", orange = "#fbd6be", pink   = "#f6c0ce", 
  purple = "#cac0e3", blue   = "#abdaec", green  = "#97d1a0"
)

# ------------------------------------------------------------------------------
# 1A: 火山图 Volcano Plot 
# ------------------------------------------------------------------------------
print("正在绘制 1A 火山图...")
p_volcano_pure <- ggplot(res_df_mapped, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey60", alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60", alpha = 0.8) +
  geom_point(aes(color = change), alpha = 0.75, size = 1.3, stroke = 0) +
  scale_color_manual(values = c("Up" = "#f6c0ce", "Down" = "#abdaec", "Not" = "#E0E0E0")) +
  theme_classic() +
  labs(title = "Volcano Plot (ER+ vs ER-)", 
       x = "log2(Fold Change)", 
       y = "-log10(Adjusted P-value)") +
  # 🌟 字体升级：强制使用 sans
  theme(text = element_text(family = "sans", color = "black"),
        axis.text = element_text(family = "sans", size = 12, color = "black"),
        axis.title = element_text(family = "sans", size = 14, face = "bold"),
        plot.title = element_text(family = "sans", hjust = 0.5, face = "bold", size = 16), 
        legend.position = "none") 

ggsave("Figure1A_Volcano_Pure.pdf", plot = p_volcano_pure, width = 6, height = 6)

# ------------------------------------------------------------------------------
# 1B: 韦恩图 Venn Diagram 
# ------------------------------------------------------------------------------
print("正在绘制 1B 韦恩图...")
venn_list <- list("TCGA BRCA DEGs" = DEGs_real, 
                  "GeneCards NETs" = nets_genes_new)

p_venn_perfect <- ggvenn(
  venn_list, 
  fill_color = c("#abdaec", "#fbd6be"),  
  fill_alpha = 0.7,                      
  stroke_size = 0.3,                     
  stroke_color = "white",                
  set_name_color = c("#5c92b8", "#d69c73"), 
  set_name_size = 5,                     
  text_size = 6,                         
  text_color = "black",
  show_percentage = FALSE                
) +
  labs(title = "NETs-related DEGs") +
  theme(
    # 🌟 字体升级：强制使用 sans
    text = element_text(family = "sans", color = "black"),
    plot.title = element_text(family = "sans", hjust = 0.5, face = "bold", size = 15, margin = margin(b = 20))
  )

ggsave("Figure1B_Venn_Perfect.pdf", plot = p_venn_perfect, width = 7, height = 5)


# ------------------------------------------------------------------------------
# 1C: 表达热图 Heatmap 
# ------------------------------------------------------------------------------
print("正在绘制 1C 表达热图...")
target_data <- res_df_mapped[res_df_mapped$Real_Symbol %in% target_NETs_genes, ]
heat_mat <- as.matrix(log2(exp_mat_er[target_data$GeneSymbol, ] + 1))
rownames(heat_mat) <- make.names(target_data$Real_Symbol, unique = TRUE)

clin_er_sorted <- clin_er[order(clin_er$ER_status), , drop = FALSE]
heat_mat_sorted <- heat_mat[, rownames(clin_er_sorted)]

annotation_col_sorted <- data.frame(ER_Status = clin_er_sorted$ER_status)
rownames(annotation_col_sorted) <- rownames(clin_er_sorted)
ann_colors <- list(ER_Status = c(Negative = "#abdaec", Positive = "#f6c0ce"))

pdf(file = "Figure1C_Heatmap_Optimized.pdf", width = 8.5, height = 6)
pheatmap(heat_mat_sorted,
         scale = "row",
         cluster_cols = FALSE,  
         cluster_rows = TRUE,
         treeheight_row = 35,          
         show_colnames = FALSE,
         show_rownames = FALSE,     
         border_color = NA,         
         annotation_col = annotation_col_sorted,
         annotation_colors = ann_colors,
         annotation_names_col = FALSE, 
         color = colorRampPalette(c("#abdaec", "#FFFFFF", "#fbd6be"))(100), 
         breaks = seq(-2, 2, length.out = 101),
         fontsize = 11,  
         fontfamily = "sans", # 🌟 字体升级：pheatmap 专用参数
         main = "Expression Profiling of Candidate NETs Genes")
dev.off()


# ------------------------------------------------------------------------------
# 1D: GO 气泡图 (恢复 BP/CC/MF 分面，统一字体)
# ------------------------------------------------------------------------------
print("正在绘制 1D GO气泡图...")
p_go <- dotplot(go_res, showCategory = 5, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scale = "free") +
  scale_color_gradient(low = "#f6c0ce", high = "#abdaec") + 
  scale_fill_gradient(low = "#f6c0ce", high = "#abdaec") + 
  theme_bw() + 
  labs(title = "GO Enrichment Analysis") +
  # 🌟 字体升级：统一坐标轴和右侧 BP/CC/MF 标签的字体
  theme(
    text = element_text(family = "sans", color = "black"),
    plot.title = element_text(family = "sans", hjust = 0.5, face = "bold", size = 15),
    axis.text = element_text(family = "sans", size = 10, color = "black"),
    strip.text = element_text(family = "sans", size = 11, face = "bold") # 分割面的字体
  )

ggsave("Figure1D_GO_Bubble_Fixed.pdf", plot = p_go, width = 7, height = 8)


# ------------------------------------------------------------------------------
# 1D: GO 圈图 (完整保留代码，强行注入字体)
# ------------------------------------------------------------------------------
print("正在绘制 1D GO圈图...")
GO <- as.data.frame(go_res)
ontology.col = c("#abdaec", "#f6c0ce", "#97d1a0") 
data = GO[order(GO$pvalue),]
datasig = data[data$pvalue < 0.05, , drop=F]

BP = head(datasig[datasig$ONTOLOGY == "BP", , drop=F], 6)
CC = head(datasig[datasig$ONTOLOGY == "CC", , drop=F], 6)
MF = head(datasig[datasig$ONTOLOGY == "MF", , drop=F], 6)
data = rbind(BP, CC, MF)

data$ONTOLOGY <- factor(data$ONTOLOGY, levels = c("BP", "CC", "MF"))
main.col = ontology.col[as.numeric(data$ONTOLOGY)]

BgGene = as.numeric(sapply(strsplit(data$BgRatio, "/"), '[', 1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio, '/'), '[', 1))
ratio = Gene / BgGene
logpvalue = -log(data$pvalue, 10)

logpvalue.col = colorRampPalette(c("#fbd6be", "#f6c0ce", "#cac0e3"))(6) 
f = colorRamp2(breaks = c(0, 2, 4, 6, 8, 10), colors = logpvalue.col)
BgGene.col = f(logpvalue)

df = data.frame(GO = data$ID, start = 1, end = max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO = data$ID, start = 1, end = BgGene, BgGene = BgGene, BgGene.col = BgGene.col)
bed3 = data.frame(GO = data$ID, start = 1, end = Gene, BgGene = Gene)
bed4 = data.frame(GO = data$ID, start = 1, end = max(BgGene), ratio = ratio, col = main.col)
bed4$ratio = bed4$ratio / max(bed4$ratio) * 9.5

pdf(file = "Figure1D_GO_Circlize_Final.pdf", width = 10, height = 10)
# 🌟 字体升级：底层 par 参数设置
par(family = "sans", omi = c(0.1, 0.1, 0.1, 1.5))
circos.par(track.margin = c(0.01, 0.01))
circos.genomicInitialize(df, plotType = "none")

circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA, bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si, track.index = 1,
              major.at = seq(0, max(BgGene), by = 100), labels.facing = "clockwise")
}

circos.genomicTrack(bed2, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj = 0, cex = 0.8, ...)
                    })

circos.genomicTrack(bed3, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#cac0e3', border = NA, ...) 
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex = 0.9, adj = 0, ...)
                    })

circos.genomicTrack(bed4, ylim = c(0, 10), track.height = 0.35, bg.border = "white", bg.col = "grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2] - cell.ylim[1]) / 10 * j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], border = NA, ...)
                    })
circos.clear()

# 🌟 字体升级：ComplexHeatmap 图例字体设置 (gpar)
middle.legend = Legend(
  labels = c('Number of Genes', 'Number of Select', 'Rich Factor(0-1)'),
  type = "points", pch = c(15, 15, 17), 
  legend_gp = gpar(col = c('#f6c0ce', '#cac0e3', ontology.col[1])),
  labels_gp = gpar(fontfamily = "sans"), # 字体
  title = "", nrow = 3, size = unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend, x = circle_size * 0.42)

main.legend = Legend(
  labels = c("Biological Process", "Cellular Component", "Molecular Function"), type = "points", pch = 15,
  legend_gp = gpar(col = ontology.col), title_position = "topcenter",
  labels_gp = gpar(fontfamily = "sans"), title_gp = gpar(fontfamily = "sans", fontface = "bold"), # 字体
  title = "ONTOLOGY", nrow = 3, size = unit(3, "mm"), grid_height = unit(5, "mm"), grid_width = unit(5, "mm")
)

logp.legend = Legend(
  labels = c('(0,2]', '(2,4]', '(4,6]', '(6,8]', '(8,10]', '>=10'),
  type = "points", pch = 16, legend_gp = gpar(col = logpvalue.col), title = "-log10(Pvalue)",
  title_position = "topcenter", 
  labels_gp = gpar(fontfamily = "sans"), title_gp = gpar(fontfamily = "sans", fontface = "bold"), # 字体
  grid_height = unit(5, "mm"), grid_width = unit(5, "mm"), size = unit(3, "mm")
)

lgd = packLegend(main.legend, logp.legend)
draw(lgd, x = circle_size * 0.85, y = circle_size * 0.55, just = "left")
dev.off()


# ------------------------------------------------------------------------------
# 1E: KEGG 气泡图与 KEGG 柱状图 
# ------------------------------------------------------------------------------
if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
  
  print("正在绘制 1E KEGG气泡图与柱状图...")
  p_kegg_dot <- dotplot(kegg_res, showCategory = 10, color = "pvalue") +
    scale_color_gradient(low = "#f6c0ce", high = "#abdaec") + 
    scale_fill_gradient(low = "#f6c0ce", high = "#abdaec") + 
    theme_bw() + 
    labs(title = "KEGG Pathway Enrichment") +
    # 🌟 字体升级
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(family = "sans", hjust = 0.5, face = "bold", size = 15),
      axis.text = element_text(family = "sans", size = 10, color = "black")
    )
  
  ggsave("Figure1E_KEGG_Dotplot_Fixed.pdf", plot = p_kegg_dot, width = 7, height = 5)
  
  kegg_df <- as.data.frame(kegg_res)
  kegg_df <- kegg_df %>%
    mutate(Ratio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
    arrange(desc(Ratio_num)) 
  
  top10_kegg <- head(kegg_df, 10)
  top10_kegg$Description <- factor(top10_kegg$Description, levels = rev(top10_kegg$Description))
  
  p_bar_matched <- ggplot(top10_kegg, aes(x = Count, y = Description, fill = pvalue)) +
    geom_col(width = 0.7) +
    scale_fill_gradient(low = "#f6c0ce", high = "#abdaec") + 
    theme_bw() +
    labs(title = "KEGG Pathway Enrichment", x = "Count", y = "") +
    # 🌟 字体升级
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(family = "sans", hjust = 0.5, face = "bold", size = 15),
      axis.text.y = element_text(family = "sans", size = 10, color = "black"),
      axis.text.x = element_text(family = "sans", size = 10, color = "black")
    )
  
  ggsave("Figure1E_KEGG_Bar_Matched.pdf", plot = p_bar_matched, width = 8, height = 5)  
  
  print("✅ Figure 1 (1A-1E) 全部图表生成完毕！完美保留原有细节，全剧终入顶刊 Sans 字体！")
}