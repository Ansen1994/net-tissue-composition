library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(stringr)
library(glmnet)
library(survival)
library(randomForest)
library(caret)
library(e1071)
library(ggplot2)
library(ggvenn) 

# ==============================================================================
# 第一部分：数据载入与 ER+ 核心样本提取
# ==============================================================================
# 1. 载入原始数据
load("brca.gdc_2024.rda")
full_exp <- assay(expquery2)
clin_info <- as.data.frame(colData(expquery2))
gene_info <- as.data.frame(rowData(expquery2))

# 2. 核心修复：使用 paper_ 前缀的真实文献级生存数据
clin_fixed <- clin_info %>%
  mutate(
    OS = ifelse(paper_vital_status == "Dead", 1, 0),
    # 死亡取 paper_days_to_death，存活取 paper_days_to_last_followup
    OS.time = ifelse(paper_vital_status == "Dead", 
                     suppressWarnings(as.numeric(as.character(paper_days_to_death))), 
                     suppressWarnings(as.numeric(as.character(paper_days_to_last_followup)))),
    Patient_ID = str_sub(sample_submitter_id, 1, 12)
  ) %>%
  filter(!is.na(OS.time) & OS.time > 0) %>%
  distinct(Patient_ID, .keep_all = TRUE)

print(paste("读取文献级生存数据后，BRCA 总有效样本数：", nrow(clin_fixed)))

# 3. 双通道提取 ER 阳性患者 (IHC 阳性 或 PAM50 LumA/B)
if("breast_carcinoma_estrogen_receptor_status" %in% colnames(clin_fixed)) {
  clin_er <- clin_fixed %>%
    filter(breast_carcinoma_estrogen_receptor_status == "Positive" | 
             paper_BRCA_Subtype_PAM50 %in% c("LumA", "LumB"))
} else {
  clin_er <- clin_fixed %>%
    filter(paper_BRCA_Subtype_PAM50 %in% c("LumA", "LumB"))
}
rownames(clin_er) <- clin_er$Patient_ID

# 4. 读取 133 个特征基因表达矩阵并对齐列名
heat_mat_old <- read.csv("Table_S2_Target_NETs_Matrix.csv", row.names = 1, check.names = FALSE)
target_ensembl <- rownames(gene_info)[gene_info$gene_name %in% rownames(heat_mat_old)]

full_exp_133 <- log2(full_exp[target_ensembl, ] + 1)
rownames(full_exp_133) <- gene_info$gene_name[match(rownames(full_exp_133), rownames(gene_info))]
colnames(full_exp_133) <- str_sub(gsub("\\.", "-", colnames(full_exp_133)), 1, 12)

# 5. 生成最终机器学习矩阵
shared_samples <- intersect(colnames(full_exp_133), rownames(clin_er))
ml_data <- as.data.frame(t(full_exp_133[, shared_samples]))
ml_data$OS <- clin_er[shared_samples, "OS"]
ml_data$OS.time <- clin_er[shared_samples, "OS.time"]

write.csv(ml_data, file = "Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv", row.names = TRUE)
print(paste("终极对齐完毕！输入机器学习的真实 ER+ 样本数：", nrow(ml_data)))


# ==============================================================================
# 第二部分：机器学习特征筛选与排版级图表绘制 (Figure 2A - 2E)
# ==============================================================================

# 提取 6 色标准色值
macaque_colors <- c("#f9eba7", "#fbd6be", "#f6c0ce", "#cac0e3", "#abdaec", "#97d1a0")

# 载入用于建模的数据
ml_data <- read.csv("Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv", row.names = 1)
x_mat <- as.matrix(ml_data[, !(colnames(ml_data) %in% c("OS", "OS.time"))])
surv_obj <- Surv(time = ml_data$OS.time, event = ml_data$OS)
os_factor <- as.factor(ifelse(ml_data$OS == 1, "Dead", "Alive"))

# ================= 算法一：Lasso 回归 (Fig 2A & 2B) =================
set.seed(2026) 
cv_lasso <- cv.glmnet(x_mat, surv_obj, family = "cox", alpha = 1, nfolds = 10)
fit_lasso <- glmnet(x_mat, surv_obj, family = "cox", alpha = 1)

lasso_coef <- coef(cv_lasso, s = "lambda.min")
lasso_genes <- rownames(lasso_coef)[which(lasso_coef != 0)]

# --- Fig 2A: Lasso Coefficient Trajectory ---
pdf("Figure2A_Lasso_Trajectory.pdf", width = 6, height = 5)
par(family = "sans") # 🌟 字体升级：注入无衬线字体
x_range <- range(-log(fit_lasso$lambda))
plot(fit_lasso, xvar = "lambda", label = FALSE, las = 1, 
     xlim = c(x_range[2], x_range[1]), 
     col = rep(macaque_colors, length.out = ncol(x_mat)), 
     lwd = 2)
title(main = "Lasso Coefficient Trajectory", line = 2.5, font.main = 2, family = "sans")
dev.off()

# --- Fig 2B: Lasso Cross-Validation ---
pdf("Figure2B_Lasso_CV_Fixed.pdf", width = 6, height = 5)
par(family = "sans") # 🌟 字体升级：注入无衬线字体
plot(cv_lasso, sign.lambda = 1)
points(log(cv_lasso$lambda), cv_lasso$cvm, pch = 20, col = "#f6c0ce", cex = 1)
title(main = "Lasso Cross-Validation", line = 2.5, font.main = 2, family = "sans")
dev.off()


# ================= 算法二：随机森林 (Fig 2C & 2D) =================
set.seed(2026)
rf_model <- randomForest(x = x_mat, y = os_factor, ntree = 500, importance = TRUE)

rf_imp <- importance(rf_model)
rf_imp_df <- data.frame(Gene = rownames(rf_imp), Gini = rf_imp[, "MeanDecreaseGini"])
rf_imp_df <- rf_imp_df[order(-rf_imp_df$Gini), ]
rf_genes <- rf_imp_df$Gene[rf_imp_df$Gini > mean(rf_imp_df$Gini)]

# --- Fig 2C: Random Forest Error Rate ---
pdf("Figure2C_RF_Error.pdf", width = 6, height = 5)
par(family = "sans") # 🌟 字体升级：注入无衬线字体
plot(rf_model, main = "Random Forest Error Rate", 
     col = c("#cac0e3", "#f6c0ce", "#abdaec"), lwd = 2.5)
legend("topright", legend = colnames(rf_model$err.rate), 
       col = c("#cac0e3", "#f6c0ce", "#abdaec"), lty = 1:3, lwd = 2.5, bty = "n")
dev.off()

# --- Fig 2D: Random Forest Variable Importance ---
library(patchwork)

# 1. 提取并格式化重要性数据
imp_df <- as.data.frame(importance(rf_model)) %>% 
  rownames_to_column("Gene")
hub_genes_highlight <- c("PLG", "CD79A", "MMP7", "ANGPT1")

# 🌟 回归顶刊标准：严格只取前 20 名，保持画板紧凑精致
needed_n <- 20 

# 2. 准备 Accuracy 数据
top_acc <- imp_df %>% 
  arrange(desc(MeanDecreaseAccuracy)) %>% 
  head(needed_n) %>%
  mutate(
    Gene = factor(Gene, levels = rev(Gene)),
    is_hub = Gene %in% hub_genes_highlight
  )

# 3. 准备 Gini 数据
top_gini <- imp_df %>% 
  arrange(desc(MeanDecreaseGini)) %>% 
  head(needed_n) %>%
  mutate(
    Gene = factor(Gene, levels = rev(Gene)),
    is_hub = Gene %in% hub_genes_highlight
  )

# 4. 构建统一的棒棒糖绘图函数
plot_lollipop <- function(data, x_var, x_label, title) {
  
  y_levels <- levels(data$Gene)
  text_colors <- ifelse(y_levels %in% hub_genes_highlight, "#d69c73", "black")
  text_faces <- ifelse(y_levels %in% hub_genes_highlight, "bold.italic", "plain")
  
  ggplot(data, aes_string(x = x_var, y = "Gene")) +
    geom_segment(aes_string(x = 0, xend = x_var, y = "Gene", yend = "Gene"), color = "grey85", linewidth = 1) +
    geom_point(aes_string(fill = x_var), shape = 21, 
               color = ifelse(data$is_hub, "#d69c73", "white"), 
               stroke = ifelse(data$is_hub, 1.8, 0.5), 
               size = 4.5) +
    scale_fill_gradient(low = "#abdaec", high = "#f6c0ce") +
    theme_minimal() +
    labs(x = x_label, y = NULL, title = title) +
    theme(
      text = element_text(family = "sans", color = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(family = "sans", face = "bold", size = 13, hjust = 0.5),
      axis.title.x = element_text(family = "sans", face = "bold", size = 11, margin = ggplot2::margin(t = 10)),
      # 注: ggplot2 新版本对矢量化文本有警告，但依然可完美渲染，安全忽略即可
      axis.text.y = element_text(family = "sans", size = 10, face = text_faces, color = text_colors),
      axis.text.x = element_text(family = "sans", size = 10, color = "black"),
      legend.position = "none"
    )
}

# 5. 生成左右两张图
# 忽略这里的 "Vectorized input..." 警告，它是新版 ggplot2 的特性，不影响最终出图
suppressWarnings(p_acc <- plot_lollipop(top_acc, "MeanDecreaseAccuracy", "MeanDecreaseAccuracy", "Accuracy Importance"))
suppressWarnings(p_gini <- plot_lollipop(top_gini, "MeanDecreaseGini", "MeanDecreaseGini", "Gini Importance"))

# 6. 使用 patchwork 拼接图表 (🌟 终极防报错语法：强制指定 ggplot2 命名空间)
p_combined_rf <- p_acc + p_gini + 
  plot_annotation(
    title = "Random Forest Variable Importance",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(family = "sans", face = "bold", size = 16, hjust = 0.5))
  )

# 动态调整画板高度：Top 20 对应 6 即可，保持紧凑美观
dynamic_height <- 6

# 导出 PDF 与高分辨率 TIFF
pdf("Figure2D_RF_Importance_Lollipop.pdf", width = 9, height = dynamic_height)
print(p_combined_rf)
dev.off()

ggsave("Figure2D_RF_Importance_Lollipop.tif", plot = p_combined_rf, width = 9, height = dynamic_height, dpi = 600, compression = "lzw")


# ================= 核心重合区：韦恩图 (Fig 2E) =================
hub_genes <- intersect(lasso_genes, rf_genes)

# --- Fig 2E: 双韦恩图 ---
venn_list_dual <- list("Lasso" = lasso_genes, 
                       "Random Forest" = rf_genes)

p_venn_dual <- ggvenn(
  venn_list_dual, 
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
  labs(title = "Identification of Hub Signature") +
  # 🌟 字体升级：强制应用 sans 并修复 margin 冲突
  theme(
    text = element_text(family = "sans", color = "black"),
    plot.title = element_text(family = "sans", hjust = 0.5, face = "bold", size = 15, margin = ggplot2::margin(b = 20))
  )

ggsave("Figure2E_Dual_Venn.pdf", plot = p_venn_dual, width = 6.5, height = 5)

# ================= 输出特征表格 =================
max_len <- max(length(lasso_genes), length(rf_genes))
lasso_col <- c(lasso_genes, rep(NA, max_len - length(lasso_genes)))
rf_col <- c(rf_genes, rep(NA, max_len - length(rf_genes)))
hub_col <- c(hub_genes, rep(NA, max_len - length(hub_genes)))

table_s4 <- data.frame(
  LASSO_Features = lasso_col, 
  RandomForest_Features = rf_col, 
  Hub_Signature = hub_col
)

write.csv(table_s4, "Table_S4_Machine_Learning_Features.csv", row.names = FALSE, na = "")
print(paste("最终锁定的核心基因为：", paste(hub_genes, collapse = ", ")))