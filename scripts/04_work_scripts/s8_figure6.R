# [patched] setwd removed; cwd is work/
# ================= 0. 加载所有依赖包 =================
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(GSVA)
library(glmnet)
library(estimate)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db", ask = FALSE)
}
library(org.Hs.eg.db)

# 全局高级配色设定（马卡龙蓝/莫兰迪粉）
my_pal <- c("Low" = "#abdaec", "High" = "#f6c0ce")

# 🌟 统一顶刊字体主题 (强制 Sans 无衬线体，AI 打开即为 Arial/Helvetica)
journal_font <- theme(
  text = element_text(family = "sans", color = "black"),
  plot.title = element_text(family = "sans", hjust = 0.5, face = "bold", size = 14),
  axis.title = element_text(family = "sans", face = "bold"),
  axis.text = element_text(family = "sans", color = "black"),
  legend.text = element_text(family = "sans"),
  legend.title = element_text(family = "sans", face = "bold")
)

# ================= 1. 数据准备与基因映射 =================
print("正在提取并匹配表达矩阵...")
load("../raw/TCGA/brca.gdc_2024.rda")
expr_matrix <- assay(expquery2)
colnames(expr_matrix) <- substr(gsub("\\.", "-", colnames(expr_matrix)), 1, 12)

# [patched] 改读 figure5 新生成的 Table S7（与 final_data.csv 同结构，避免用到旧副本）
final_data <- read.csv("Table_S7_Final_Clinical_Data.csv", stringsAsFactors = FALSE)
final_data$RiskGroup <- factor(final_data$RiskGroup, levels = c("Low", "High"))
common_samples <- intersect(colnames(expr_matrix), final_data$Patient_ID)
expr_matched <- expr_matrix[, common_samples]

risk_info <- final_data %>%
  dplyr::filter(Patient_ID %in% common_samples) %>%
  dplyr::select(Patient_ID, RiskScore, RiskGroup) %>%
  dplyr::arrange(match(Patient_ID, common_samples))

ensembl_ids <- gsub("\\..*", "", rownames(expr_matched))
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

expr_symbol <- expr_matched
rownames(expr_symbol) <- gene_symbols
expr_symbol <- expr_symbol[!is.na(rownames(expr_symbol)), ]
expr_symbol <- expr_symbol[!duplicated(rownames(expr_symbol)), ]


# ================= Figure 6A: 28 种免疫细胞全景图 =================
charoentong_28 <- list(
  `Activated B cell` = c("CD79A", "CD79B", "BLNK", "CD19", "MS4A1"),
  `Activated CD4 T cell` = c("CD4", "IL2RA", "CD69", "IL2", "CD40LG"),
  `Activated CD8 T cell` = c("CD8A", "CD8B", "GZMB", "PRF1", "IFNG"),
  `Central memory CD4 T cell` = c("SELL", "CCR7", "CD27", "CD28", "IL7R"),
  `Central memory CD8 T cell` = c("SELL", "CCR7", "CD27", "CD28", "EOMES"),
  `Effector memory CD4 T cell` = c("IL7R", "FAS", "CCR7", "CD27"),
  `Effector memory CD8 T cell` = c("GZMB", "PRF1", "FAS", "CCR7"),
  `Gamma delta T cell` = c("TRDC", "TRGC1", "TRGC2", "TRDV1"),
  `Immature B cell` = c("CD19", "MS4A1", "CD38", "IL4R"),
  `Memory B cell` = c("CD27", "CD38", "IGHM", "IGHD"),
  `Regulatory T cell` = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18"),
  `T follicular helper cell` = c("CXCR5", "BCL6", "IL21", "PDCD1"),
  `Type 1 T helper cell` = c("TBX21", "IFNG", "CXCR3", "IL12RB2"),
  `Type 17 T helper cell` = c("RORC", "IL17A", "IL17F", "CCR6"),
  `Type 2 T helper cell` = c("GATA3", "IL4", "IL5", "IL13"),
  `Activated dendritic cell` = c("CD83", "CD86", "HLA-DRA", "CD40"),
  `CD56bright natural killer cell` = c("NCAM1", "SELL", "GZMK"),
  `CD56dim natural killer cell` = c("NCAM1", "FCGR3A", "PRF1", "GZMB"),
  `Eosinophil` = c("IL5RA", "CCR3", "EPX", "PRG2"),
  `Immature dendritic cell` = c("CD1A", "CD1C", "FCER1A"),
  `Macrophage` = c("CD68", "CD163", "MRC1", "MSR1"),
  `Mast cell` = c("TPSAB1", "CPA3", "HDC", "MS4A2"),
  `MDSC` = c("ITGAM", "CD33", "FUT4", "IL4R"),
  `Monocyte` = c("CD14", "FCGR3A", "CSF1R", "CCR2"),
  `Natural killer cell` = c("KLRB1", "NCR1", "NKG7", "GNLY"),
  `Natural killer T cell` = c("CD3E", "NCAM1", "KLRB1", "ZBTB16"),
  `Neutrophil` = c("CEACAM8", "FCGR3B", "CSF3R", "S100A12", "CXCL8"),
  `Plasmacytoid dendritic cell` = c("CLEC4C", "NRP1", "IL3RA", "IRF7")
)

print("正在计算 ssGSEA...")
ssgsea_param <- ssgseaParam(as.matrix(expr_symbol), charoentong_28)
ssgsea_res <- gsva(ssgsea_param)

ssgsea_df <- as.data.frame(t(ssgsea_res))
ssgsea_df$Patient_ID <- rownames(ssgsea_df)
plot_data <- merge(ssgsea_df, risk_info, by = "Patient_ID")
plot_data_long <- plot_data %>% pivot_longer(cols = names(charoentong_28), names_to = "CellType", values_to = "Score")

# ================= Figure 6A: 28 种免疫细胞全景图 =================
print("正在绘制 Figure 6A 并导出高分辨率 TIFF 与 PDF...")

# 1. 生成带有散点 (jitter) 的图形对象
p_immune <- ggboxplot(plot_data_long, x = "CellType", y = "Score", color = "RiskGroup", 
                      palette = my_pal,
                      add = "jitter", add.params = list(alpha = 0.2, size = 0.5)) +
  stat_compare_means(aes(group = RiskGroup), label = "p.signif", method = "wilcox.test", hide.ns = TRUE) +
  theme_bw() + journal_font + # 🌟 注入字体
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        legend.position = "top") +
  labs(x = "Immune Cell Types", y = "ssGSEA Enrichment Score", title = "TME Immune Cell Infiltration Profile")

# 2. 保存为 PDF (备用)
pdf("Figure6A_28_Immune_Cells.pdf", width = 14, height = 7)
print(p_immune)
dev.off()

# 🌟 3. 保存为超高分辨率 TIFF (解决 AI 卡顿的最佳方案)
# dpi = 600 保证印刷级清晰度，compression = "lzw" 保证文件不会过大
ggsave("Figure6A_28_Immune_Cells.tiff", plot = p_immune, width = 14, height = 7, 
       dpi = 600, compression = "lzw")

print("==== Figure 6A 导出完毕！请直接将生成的 .tif 文件拖入 AI 中排版 ====")


# ================= Figure 6B: 核心免疫细胞相关性散点图 (独立保存) =================
print("正在生成相关性散点图...")
all_immune_cells <- names(charoentong_28)

for (cell in all_immune_cells) {
  tmp_data <- plot_data[, c("Patient_ID", "RiskScore", cell)]
  
  p <- ggscatter(tmp_data, x = "RiskScore", y = cell, 
                 color = "#abdaec", size = 1.5, alpha = 0.8,
                 add = "reg.line", add.params = list(color = "#f6c0ce", fill = "lightgray"), conf.int = TRUE,
                 cor.coef = TRUE, cor.method = "spearman",
                 xlab = "Risk Score", ylab = paste0(cell, " (ssGSEA)"),
                 title = paste0("Correlation: RiskScore vs ", cell)) +
    journal_font # 🌟 注入字体
  
  safe_cell_name <- gsub("[ /]", "_", cell)
  file_name <- paste0("Figure6B_Cor_", safe_cell_name, ".pdf")
  
  pdf(file_name, width = 5, height = 5)
  print(p)
  graphics.off()
}


# ================= Figure 6C: 核心免疫检查点表达差异 =================
icb_genes <- c("PDCD1", "CD274", "CTLA4")
icb_exist <- intersect(icb_genes, rownames(expr_symbol))

icb_expr <- as.data.frame(t(expr_symbol[icb_exist, ]))
icb_expr$Patient_ID <- rownames(icb_expr)
data_6c <- merge(icb_expr, risk_info[, c("Patient_ID", "RiskGroup")], by = "Patient_ID")
data_6c_long <- data_6c %>% pivot_longer(cols = all_of(icb_exist), names_to = "Gene", values_to = "Expression")

pdf("Figure6C_Immune_Checkpoints.pdf", width = 6, height = 5)
p_icb <- ggboxplot(data_6c_long, x = "Gene", y = "Expression", color = "RiskGroup",
                   palette = my_pal,
                   add = "jitter", add.params = list(alpha = 0.3, size = 0.5)) +
  stat_compare_means(aes(group = RiskGroup), label = "p.signif", method = "wilcox.test") +
  theme_bw() + journal_font + # 🌟 注入字体
  labs(x = "Immune Checkpoints", y = "Expression Level (RNA-seq)", title = "Expression of Key Immune Checkpoints")
print(p_icb)
dev.off()


# ================= Figure 6D: 药物敏感性 IC50 预测 (固定种子) =================
print("开始预测 GDSC2 药物敏感性...")
GDSC2_Expr <- readRDS("../raw/drug/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res <- readRDS("../raw/drug/GDSC2_Res.rds")
candidate_drugs <- c("paclitaxel", "docetaxel", "doxorubicin", "epirubicin", "cisplatin", "5-fluorouracil", 
                     "gemcitabine", "cyclophosphamide", "tamoxifen", "fulvestrant", "letrozole", 
                     "anastrozole", "exemestane", "palbociclib", "ribociclib", "abemaciclib", "everolimus", "alpelisib")
available_cols <- grep(paste(candidate_drugs, collapse="|"), colnames(GDSC2_Res), ignore.case=TRUE, value=TRUE)

X_test <- as.data.frame(t(expr_symbol))
common_genes <- intersect(rownames(GDSC2_Expr), colnames(X_test))
X_test_matched <- as.matrix(X_test[, common_genes])
safe_genes <- make.names(common_genes, unique = TRUE)
colnames(X_test_matched) <- safe_genes

for (drug_target in available_cols) {
  y_all <- GDSC2_Res[, drug_target]
  valid_cells <- !is.na(y_all)
  y_train <- y_all[valid_cells]
  
  X_train_full <- GDSC2_Expr[, valid_cells, drop = FALSE]
  X_train <- t(X_train_full[common_genes, ])
  colnames(X_train) <- safe_genes
  
  set.seed(42)
  cv_ridge <- cv.glmnet(as.matrix(X_train), as.numeric(y_train), alpha = 0, nfolds = 5)
  pred_tmp <- as.numeric(predict(cv_ridge, s = "lambda.min", newx = X_test_matched))
  
  plot_df_tmp <- data.frame(Patient_ID = rownames(X_test_matched), IC50 = pred_tmp)
  plot_df_tmp$Patient_ID <- substr(plot_df_tmp$Patient_ID, 1, 12)
  final_tmp <- merge(plot_df_tmp, risk_info, by = "Patient_ID")
  
  if(nrow(final_tmp) > 0) {
    p_val <- wilcox.test(IC50 ~ RiskGroup, data = final_tmp)$p.value
    if(p_val < 0.05) {
      p_plot <- ggboxplot(final_tmp, x = "RiskGroup", y = "IC50",
                          color = "RiskGroup", palette = my_pal,
                          add = "jitter", add.params = list(alpha = 0.3, size = 1)) +
        stat_compare_means(comparisons = list(c("Low", "High")), method = "wilcox.test", label = "p.signif") +
        theme_bw() + journal_font + # 🌟 注入字体
        labs(x = "Risk Group", y = "Estimated IC50", title = paste(drug_target, "Sensitivity"))
      
      safe_drug_name <- gsub("[ /()]", "_", drug_target)
      pdf(paste0("Figure6D_", safe_drug_name, "_IC50.pdf"), width = 5, height = 5)
      print(p_plot)
      graphics.off()
    }
  }
}


# ================= Figure 6E: ESTIMATE 免疫评分 =================
print("正在计算 ESTIMATE 评分...")
expr_matrix <- as.matrix(expr_symbol)
gct_file <- "Breast_expr_input.gct"
con <- file(gct_file, open = "w")
writeLines("#1.2", con)
writeLines(paste(nrow(expr_matrix), ncol(expr_matrix), sep = "\t"), con)
writeLines(paste(c("Name", "Description", colnames(expr_matrix)), collapse = "\t"), con)
expr_df <- data.frame(Name = rownames(expr_matrix), Description = "na", expr_matrix, check.names = FALSE)
write.table(expr_df, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
close(con)

estimateScore(input.ds = gct_file, output.ds = "../raw/drug/Breast_estimate_score.gct", platform = "illumina")

score_raw <- read.table("../raw/drug/Breast_estimate_score.gct", skip = 2, header = TRUE, sep = "\t")
if("Description" %in% colnames(score_raw)) score_raw$Description <- NULL
rownames(score_raw) <- score_raw$NAME
score_raw$NAME <- NULL

score_res <- as.data.frame(t(score_raw))
score_res$Patient_ID <- substr(gsub("\\.", "-", rownames(score_res)), 1, 12)
score_res$ImmuneScore <- as.numeric(as.character(score_res$ImmuneScore))

final_score_data <- merge(score_res, risk_info, by = "Patient_ID")

pdf("Figure6E_ImmuneScore.pdf", width = 5, height = 5)
p_immune_est <- ggboxplot(final_score_data, x = "RiskGroup", y = "ImmuneScore",
                          color = "RiskGroup", palette = my_pal,
                          add = "jitter", add.params = list(alpha = 0.3, size = 1)) +
  stat_compare_means(comparisons = list(c("Low", "High")), method = "wilcox.test", label = "p.signif") +
  theme_bw() + journal_font + # 🌟 注入字体
  labs(x = "Risk Group", y = "ESTIMATE Immune Score", title = "Immune Infiltration Score (ESTIMATE)")
print(p_immune_est)
dev.off()

print("==== Figure 6 全套独立 PDF 图表生成完毕，代码严密，字体全刊统一！ ====")