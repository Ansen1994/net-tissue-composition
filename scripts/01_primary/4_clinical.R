library(GEOquery)
library(dplyr)
library(tibble)
library(survival)
library(survminer)
library(timeROC)
library(pheatmap)
library(ggplot2)
library(rms)

# ==============================================================================
# 第一部分：GEO 数据读取、Z-score 标准化与正确的临床终点清洗
# ==============================================================================
print("正在读取表达矩阵...")
gset <- getGEO(filename = "GSE20685_series_matrix.txt.gz", getGPL = FALSE)
exp_raw <- exprs(gset)
clin_raw <- pData(gset)

print("正在读取探针注释...")
gpl_anno <- read.delim("GPL570-55999.txt", sep = "\t", comment.char = "#", 
                       stringsAsFactors = FALSE, check.names = FALSE)

print("正在进行探针基因映射...")
probe2symbol <- gpl_anno %>%
  dplyr::select(ID, `Gene Symbol`) %>%
  dplyr::filter(`Gene Symbol` != "") %>%
  dplyr::mutate(symbol = sapply(strsplit(as.character(`Gene Symbol`), " /// "), "[", 1)) %>%
  dplyr::select(ID, symbol)

exp_annotated <- exp_raw %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::inner_join(probe2symbol, by = "ID") %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), max)) %>% 
  tibble::column_to_rownames("symbol")

# 核心：4 个 Hub 基因
hub_genes <- c("PLG", "CD79A", "MMP7", "ANGPT1")

available_genes <- intersect(hub_genes, rownames(exp_annotated))
geo_exp_hub <- as.data.frame(t(exp_annotated[available_genes, ]))

# Z-score 标准化（消除测序与芯片平台批次效应）
geo_exp_hub_scaled <- as.data.frame(scale(geo_exp_hub))

print("正在清洗临床随访数据（已修正为正确的 OS 总生存终点）...")
geo_clin <- clin_raw %>%
  dplyr::select(geo_accession, 
                OS = `event_death:ch1`,                    
                OS.time = `follow_up_duration (years):ch1`, 
                Age = `age at diagnosis:ch1`,
                T_stage = `t_stage:ch1`,
                N_stage = `n_stage:ch1`) %>%
  dplyr::mutate(
    OS = as.numeric(OS),
    OS.time = as.numeric(OS.time) * 365,                    
    Age = as.numeric(Age),
    # 🌟 逻辑安全锁：强制转换为因子，防止 rms::cph 报错
    T_stage = as.factor(T_stage),
    N_stage = as.factor(N_stage)
  )

# 合并标准化后的表达矩阵与临床数据
geo_final <- merge(geo_clin, geo_exp_hub_scaled, by.x = "geo_accession", by.y = "row.names")
geo_final <- na.omit(geo_final) 


# ==============================================================================
# 第二部分：计算 Risk Score 与截点优化、联合临床建模
# ==============================================================================
geo_data <- geo_final

# 计算 Risk Score (采用 TCGA 训练集的最精确底层 Cox 回归系数)
geo_data$RiskScore <- (0.1102 * geo_data$PLG) + (0.2295 * geo_data$ANGPT1) - (0.0854 * geo_data$CD79A) - (0.0851 * geo_data$MMP7)

# 使用 surv_cutpoint 寻找最佳截点
res.cut <- surv_cutpoint(geo_data, time = "OS.time", event = "OS", variables = "RiskScore")
res.cat <- surv_categorize(res.cut)
geo_data$RiskGroup <- factor(res.cat$RiskScore, levels = c("low", "high"), labels = c("Low", "High"))

# 按照 RiskScore 从小到大排序并赋予 PatientID（用于 4A, 4B, 4C 排序对齐）
geo_data <- geo_data[order(geo_data$RiskScore), ]
geo_data$PatientID <- 1:nrow(geo_data)

# 联合临床特征构建多因素模型，提取 CombinedScore
dd <- datadist(geo_data)
options(datadist = 'dd')
cox_combined <- cph(Surv(OS.time, OS) ~ RiskScore + Age + T_stage + N_stage, 
                    data = geo_data, x = TRUE, y = TRUE, surv = TRUE)
geo_data$CombinedScore <- predict(cox_combined, type = "lp")


# ==============================================================================
# 第三部分：Figure 4 (4A - 4E) 全套图表绘制 (🌟 全面注入顶刊 Sans 字体)
# ==============================================================================

# 定义全局 Sans 字体主题 (针对 ggplot2 图形)
journal_font <- theme(
  text = element_text(family = "sans", color = "black"),
  plot.title = element_text(family = "sans", hjust = 0.5, face = "bold", size = 14),
  axis.title = element_text(family = "sans", face = "bold"),
  axis.text = element_text(family = "sans", color = "black"),
  legend.text = element_text(family = "sans"),
  legend.title = element_text(family = "sans", face = "bold")
)

# --- Fig 4A: 风险评分散点图 ---
pdf("Figure4A_GEO_Risk_Score.pdf", width = 6, height = 4)
p_4a <- ggplot(geo_data, aes(x = PatientID, y = RiskScore)) +
  geom_point(aes(color = RiskGroup), size = 1.5, alpha = 0.8) +
  geom_vline(xintercept = sum(geo_data$RiskGroup == "Low"), lty = 2, color = "grey50") +
  scale_color_manual(values = c("Low" = "#abdaec", "High" = "#f6c0ce")) +
  theme_classic() + journal_font + # 🌟 注入字体
  labs(x = "GEO Patients (increasing risk score)", y = "Risk Score", title = "GEO Validation: Risk Score")
print(p_4a)
dev.off()

# --- Fig 4B: 生存状态散点图 ---
pdf("Figure4B_GEO_Survival_Status.pdf", width = 6, height = 4)
p_4b <- ggplot(geo_data, aes(x = PatientID, y = OS.time / 365)) + 
  geom_point(aes(color = as.factor(OS)), size = 1.5, alpha = 0.8) +
  geom_vline(xintercept = sum(geo_data$RiskGroup == "Low"), lty = 2, color = "grey50") +
  scale_color_manual(values = c("0" = "#abdaec", "1" = "#f6c0ce"), labels = c("Alive", "Dead")) +
  theme_classic() + journal_font + # 🌟 注入字体
  labs(x = "GEO Patients", y = "Survival Time (Years)", color = "Status", title = "GEO Validation: Survival Status")
print(p_4b)
dev.off()

# --- Fig 4C: 4个核心基因表达热图 ---
heat_data_geo <- t(geo_data[, hub_genes])
colnames(heat_data_geo) <- geo_data$PatientID
annotation_col_geo <- data.frame(Risk = geo_data$RiskGroup)
rownames(annotation_col_geo) <- colnames(heat_data_geo)
ann_colors_geo <- list(Risk = c(Low = "#abdaec", High = "#f6c0ce"))

pdf("Figure4C_GEO_Risk_Heatmap.pdf", width = 7, height = 4)
pheatmap(heat_data_geo, scale = "row", cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, annotation_col = annotation_col_geo, annotation_colors = ann_colors_geo,
         color = colorRampPalette(c("#5c92b8", "white", "#d69c73"))(100),
         breaks = seq(-2, 2, length.out = 101),
         fontfamily = "sans", # 🌟 注入字体 (pheatmap 专属参数)
         main = "GEO Validation: Expression Heatmap")
dev.off()

# --- Fig 4D: KM 生存曲线 ---
fit_km_geo <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = geo_data)
p_km_geo_perfect <- ggsurvplot(
  fit_km_geo, data = geo_data, pval = TRUE, 
  pval.coord = c(200, 0.15), conf.int = TRUE, conf.int.alpha = 0.18, 
  risk.table = TRUE, size = 0.7, palette = c("#abdaec", "#f6c0ce"),
  legend.title = "Risk", legend.labs = c("Low", "High"), 
  xlab = "Time (Days)", ylab = "Survival probability", title = NULL, 
  ggtheme = theme_classic() + journal_font # 🌟 注入字体
)

pdf("Figure4D_GEO_KM_Survival_CorrectedOS.pdf", width = 6, height = 6, onefile = FALSE)
print(p_km_geo_perfect)
dev.off() 

# --- Fig 4E: 临床-分子联合模型时间依赖型 ROC 曲线 ---
ROC_rt_geo_combined <- timeROC(
  T = geo_data$OS.time, delta = geo_data$OS, 
  marker = geo_data$CombinedScore, cause = 1, 
  weighting = "marginal", times = c(365, 365*3, 365*5), ROC = TRUE
) 

pdf("Figure4E_GEO_Time_ROC_ClinicalCombined.pdf", width = 6, height = 6)
par(family = "sans") # 🌟 注入字体 (Base R 专属设置)
plot(ROC_rt_geo_combined, time = 365, col = "#f6c0ce", title = FALSE, lwd = 2.5)
plot(ROC_rt_geo_combined, time = 365*3, col = "#abdaec", add = TRUE, lwd = 2.5) 
plot(ROC_rt_geo_combined, time = 365*5, col = "#97d1a0", add = TRUE, lwd = 2.5) 
legend("bottomright", 
       c(paste0("1-Year AUC = ", sprintf("%.3f", ROC_rt_geo_combined$AUC[1])),
         paste0("3-Year AUC = ", sprintf("%.3f", ROC_rt_geo_combined$AUC[2])), 
         paste0("5-Year AUC = ", sprintf("%.3f", ROC_rt_geo_combined$AUC[3]))), 
       col = c("#f6c0ce", "#abdaec", "#97d1a0"), lty = 1, lwd = 2.5, bty = "n")
title(main = "GEO Validation: Clinical-Molecular ROC", font.main = 2, family = "sans") # 🌟 注入字体
dev.off()

print("✅ Figure 4 (4A-4E) 全套高定验证图表已全部生成完毕，且全局锁定为顶刊排版字体！")