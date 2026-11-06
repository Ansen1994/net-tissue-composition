library(GEOquery)
library(dplyr)
library(tibble)
library(survival)
library(survminer)
library(timeROC)
library(pheatmap)
library(ggplot2)

# ==============================================================================
# 第一部分：GEO 数据的特征提取与清洗
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
if(length(available_genes) < length(hub_genes)) {
  warning("部分基因在GEO数据集中未找到！")
}

geo_exp_hub <- as.data.frame(t(exp_annotated[available_genes, ]))

print("正在清洗临床随访数据...")
geo_clin <- clin_raw %>%
  dplyr::select(geo_accession, 
                OS = `event_death:ch1`, 
                OS.time = `time_to_metastasis (years):ch1`) %>%
  dplyr::mutate(OS = as.numeric(OS), 
                OS.time = as.numeric(OS.time) * 365) # 转换为天数

geo_final <- merge(geo_clin, geo_exp_hub, by.x = "geo_accession", by.y = "row.names")
geo_final <- na.omit(geo_final) 

write.csv(geo_final, "GEO_Validation_Matrix.csv", row.names = FALSE)
print(paste("大功告成！纯本地外部验证矩阵已生成，包含有效患者数量：", nrow(geo_final)))


# ==============================================================================
# 第二部分：风险评分计算与图表绘制 (Figure 4A-E 顶刊马卡龙色系版)
# ==============================================================================
geo_data <- read.csv("GEO_Validation_Matrix.csv")

# 计算 Risk Score (精确底层 Cox 回归系数)
geo_data$RiskScore <- (0.1102 * geo_data$PLG) + (0.2295 * geo_data$ANGPT1) - (0.0854 * geo_data$CD79A) - (0.0851 * geo_data$MMP7)

# 使用 surv_cutpoint 寻找最佳截点
res.cut <- surv_cutpoint(geo_data, time = "OS.time", event = "OS", variables = "RiskScore")
res.cat <- surv_categorize(res.cut)
geo_data$RiskGroup <- factor(res.cat$RiskScore, levels = c("low", "high"), labels = c("Low", "High"))

# 重新排序并赋予编号
geo_data <- geo_data[order(geo_data$RiskScore), ]
geo_data$PatientID <- 1:nrow(geo_data)

write.csv(geo_data, "Table_S6_GEO_Validation_Result.csv", row.names = FALSE)


# -------- 绘图部分 (马卡龙高定配色) --------

# Fig 4A: GEO Risk Distribution (马卡龙蓝/粉)
pdf("Figure4A_GEO_Risk_Score.pdf", width = 6, height = 4)
p_risk_geo <- ggplot(geo_data, aes(x = PatientID, y = RiskScore)) +
  geom_point(aes(color = RiskGroup), size = 1.5, alpha = 0.8) +
  geom_vline(xintercept = sum(geo_data$RiskGroup == "Low"), lty = 2, color = "grey50") +
  scale_color_manual(values = c("Low" = "#abdaec", "High" = "#f6c0ce")) +
  theme_classic() +
  labs(x = "GEO Patients (increasing risk score)", y = "Risk Score", title = "GEO Validation: Risk Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_risk_geo)
dev.off()

# Fig 4B: GEO Survival Status Plot (马卡龙蓝/粉)
pdf("Figure4B_GEO_Survival_Status.pdf", width = 6, height = 4)
p_status_geo <- ggplot(geo_data, aes(x = PatientID, y = OS.time / 365)) + 
  geom_point(aes(color = as.factor(OS)), size = 1.5, alpha = 0.8) +
  geom_vline(xintercept = sum(geo_data$RiskGroup == "Low"), lty = 2, color = "grey50") +
  scale_color_manual(values = c("0" = "#abdaec", "1" = "#f6c0ce"), labels = c("Alive", "Dead")) +
  theme_classic() +
  labs(x = "GEO Patients", y = "Survival Time (Years)", color = "Status", title = "GEO Validation: Survival Status") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_status_geo)
dev.off()

# Fig 4C: GEO Expression Heatmap (高对比度加深渐变 + 马卡龙注释条)
heat_data_geo <- t(geo_data[, hub_genes])
annotation_col_geo <- data.frame(Risk = geo_data$RiskGroup)
rownames(annotation_col_geo) <- colnames(heat_data_geo)
ann_colors_geo <- list(Risk = c(Low = "#abdaec", High = "#fbd6be"))

pdf("Figure4C_GEO_Risk_Heatmap.pdf", width = 7, height = 4)
pheatmap(heat_data_geo, scale = "row", cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, annotation_col = annotation_col_geo, annotation_colors = ann_colors_geo,
         color = colorRampPalette(c("#5c92b8", "white", "#d69c73"))(100),
         breaks = seq(-2, 2, length.out = 101),
         main = "GEO Validation: Expression Heatmap")
dev.off()

# Fig 4D: GEO Kaplan-Meier Survival (纤细线条 + 马卡龙蓝粉)
# =====================================================================
library(survival)
library(survminer)

fit_km_geo <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = geo_data)

p_km_geo_perfect <- ggsurvplot(
  fit_km_geo, 
  data = geo_data, 
  pval = TRUE,                       
  pval.coord = c(200, 0.15),         # 调整 P 值在图中的坐标位置
  conf.int = TRUE,                   
  conf.int.alpha = 0.18,             # 通透的置信区间阴影
  risk.table = TRUE,                 
  size = 0.7,                        # 纤细精致的线条
  
  # 颜色严格与全篇同步：低风险马卡龙蓝，高风险马卡龙粉
  palette = c("#abdaec", "#f6c0ce"),
  
  legend.title = "Risk", 
  legend.labs = c("Low", "High"),    # 核心修改：精简图例文字，去掉 RiskGroup= 前缀
  xlab = "Time (Days)",
  ylab = "Survival probability",
  title = NULL,                      # 核心修改：去除冗余的顶部大标题，保持和 3E 一致
  ggtheme = theme_classic()
)

# 导出完美版 GEO KM 曲线
pdf("Figure4D_GEO_KM_Survival_Perfect.pdf", width = 6, height = 6, onefile = FALSE)
print(p_km_geo_perfect)
dev.off()

ggsave("Figure4D_GEO_KM_Survival_Perfect.tif", plot = print(p_km_geo_perfect), width = 6, height = 6, dpi = 600, compression = "lzw")

print("4D 验证集生存曲线已精修完毕，风格与 3E 完全一致！")

# Fig 4E: GEO Time-dependent ROC (粉、蓝、绿 三色高对比)
ROC_rt_geo <- timeROC(T = geo_data$OS.time, delta = geo_data$OS, 
                      marker = geo_data$RiskScore, cause = 1, 
                      weighting = "marginal", 
                      times = c(365, 365*3, 365*5), ROC = TRUE) 

pdf("Figure4E_GEO_Time_ROC.pdf", width = 6, height = 6)
plot(ROC_rt_geo, time = 365, col = "#f6c0ce", title = FALSE, lwd = 2.5)
plot(ROC_rt_geo, time = 365*3, col = "#abdaec", add = TRUE, lwd = 2.5) 
plot(ROC_rt_geo, time = 365*5, col = "#97d1a0", add = TRUE, lwd = 2.5) 
legend("bottomright", 
       c(paste0("1-Year AUC = ", sprintf("%.3f", ROC_rt_geo$AUC[1])),
         paste0("3-Year AUC = ", sprintf("%.3f", ROC_rt_geo$AUC[2])), 
         paste0("5-Year AUC = ", sprintf("%.3f", ROC_rt_geo$AUC[3]))), 
       col = c("#f6c0ce", "#abdaec", "#97d1a0"), lty = 1, lwd = 2.5, bty = "n")
title(main = "GEO Validation: Time-dependent ROC", font.main = 2)
dev.off()

print("✅ Figure 4 外部验证系列图表全部生成完毕，配色已完美契合全篇马卡龙风格！")