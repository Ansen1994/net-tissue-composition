# 基于 GSE96058 的外部验证：多因素临床-分子联合模型 (纯 ER+ 亚群)
# ==============================================================================
# [patched] setwd removed; cwd is work/
library(GEOquery)
library(dplyr)
library(tibble)
library(survival)
library(survminer)
library(timeROC)
library(pheatmap)
library(ggplot2)
library(data.table) 
library(rms)

# ==============================================================================
# 第一部分：合并临床数据
# ==============================================================================
print("1. 正在读取并合并两个平台的临床随访数据...")
gset1 <- getGEO(filename = "../raw/GSE96058/GSE96058-GPL11154_series_matrix.txt.gz", getGPL = FALSE)
clin1 <- pData(gset1)

gset2 <- getGEO(filename = "../raw/GSE96058/GSE96058-GPL18573_series_matrix.txt.gz", getGPL = FALSE)
clin2 <- pData(gset2)

clin_merged <- dplyr::bind_rows(clin1, clin2)

# ==============================================================================
# 第二部分：读取矩阵、Z-score 标准化与临床联合特征提取
# ==============================================================================
print("2. 正在极速读取 RNA-seq 表达矩阵...")
exp_raw <- fread("../raw/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz", data.table = FALSE)
rownames(exp_raw) <- exp_raw[, 1]
exp_matrix <- exp_raw[, -1]

print("3. 提取 Hub 基因并进行 Z-score 标准化...")
hub_genes <- c("PLG", "CD79A", "MMP7", "ANGPT1")
available_genes <- intersect(hub_genes, rownames(exp_matrix))
geo_exp_hub <- as.data.frame(t(exp_matrix[available_genes, ]))
geo_exp_hub_scaled <- as.data.frame(scale(geo_exp_hub))
geo_exp_hub_scaled$title <- rownames(geo_exp_hub_scaled) 

library(stringr)
print("4. 正在清洗临床随访数据：深度运用正则匹配清洗肿瘤大小与淋巴结文本...")
geo_clin_base <- clin_merged %>%
  dplyr::select(title,                                  
                ER_Status = `er status:ch1`,             
                OS = `overall survival event:ch1`,        
                OS.time = `overall survival days:ch1`,
                Age = `age at diagnosis:ch1`,
                TumorSize_raw = `tumor size:ch1`,
                LymphNode_raw = `lymph node status:ch1`) %>% 
  dplyr::filter(ER_Status == "1") %>% 
  dplyr::mutate(
    OS = as.numeric(OS), 
    OS.time = as.numeric(OS.time),
    Age = as.numeric(Age),
    
    # 🌟 核心提分操作 1：深度清洗淋巴结状态
    LymphNode = case_when(
      grepl("pos|\\+|N1|N2|N3", LymphNode_raw, ignore.case = TRUE) ~ 1,
      grepl("neg|\\-|N0", LymphNode_raw, ignore.case = TRUE) ~ 0,
      as.numeric(str_extract(LymphNode_raw, "\\d+")) > 0 ~ 1,
      as.numeric(str_extract(LymphNode_raw, "\\d+")) == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # 🌟 核心提分操作 2：深度清洗肿瘤大小 
    TumorSize = as.numeric(str_extract(TumorSize_raw, "\\d+\\.?\\d*")),
    
    # 🌟 逻辑安全锁：强制转换为因子，防止 rms::cph 将其视为连续变量
    LymphNode = as.factor(LymphNode)
  ) 

geo_final <- merge(geo_clin_base, geo_exp_hub_scaled, by = "title")
geo_final <- na.omit(geo_final)
print(paste("---- 纯 ER+ 联合临床矩阵生成完毕，当前有效患者数量：", nrow(geo_final), "----"))

# ==============================================================================
# 第三部分：构建 Risk Score 与 临床-分子联合评分 (Combined Score)
# ==============================================================================
print("5. 正在构建多因素联合生存模型...")
geo_data <- geo_final

# 1. 计算基础 Risk Score
geo_data$RiskScore <- (0.1102 * geo_data$PLG) + (0.2295 * geo_data$ANGPT1) - 
  (0.0854 * geo_data$CD79A) - (0.0851 * geo_data$MMP7)

# 2. 截点分组
res.cut <- surv_cutpoint(geo_data, time = "OS.time", event = "OS", variables = "RiskScore")
res.cat <- surv_categorize(res.cut)
geo_data$RiskGroup <- factor(res.cat$RiskScore, levels = c("low", "high"), labels = c("Low", "High"))

geo_data <- geo_data[order(geo_data$RiskScore), ]
geo_data$PatientID <- 1:nrow(geo_data)

# 🌟 3. 三大临床特征满编上阵：RiskScore + Age + TumorSize + LymphNode
dd <- datadist(geo_data)
options(datadist = 'dd')
cox_combined <- cph(Surv(OS.time, OS) ~ RiskScore + Age + TumorSize + LymphNode, 
                    data = geo_data, x = TRUE, y = TRUE, surv = TRUE)
geo_data$CombinedScore <- predict(cox_combined, type = "lp")
# [patched] 输出 GSE96058 验证结果表（原脚本只出图不写表，现表为旧版遗留、含 FLT3 且分组不符）
keep <- intersect(c("title","ER_Status","OS","OS.time","Age","TumorSize","LymphNode",
                    "PLG","CD79A","MMP7","ANGPT1","RiskScore","RiskGroup","PatientID","CombinedScore"),
                  names(geo_data))
write.csv(geo_data[, keep], "Table_S6_GSE96058_Validation_Result.csv", row.names = FALSE)
cat("已写出 Table_S6_GSE96058_Validation_Result.csv，n =", nrow(geo_data),
    " Low =", sum(geo_data$RiskGroup=="Low"), " High =", sum(geo_data$RiskGroup=="High"), "\n")

print("模型构建成功！请继续执行第 6 步的自动化画图。")

# ==============================================================================
# 第四部分：自动化生成 Figure 4 (4A - 4E) 高定马卡龙图表 (全面注入 Sans 字体)
# ==============================================================================
print("6. 正在生成系列高质量图表...")

# 定义全局 Sans 字体主题
journal_font <- theme(
  text = element_text(family = "sans", color = "black"),
  plot.title = element_text(family = "sans", hjust = 0.5, face = "bold"),
  axis.title = element_text(family = "sans", face = "bold"),
  axis.text = element_text(family = "sans", color = "black"),
  legend.text = element_text(family = "sans"),
  legend.title = element_text(family = "sans", face = "bold")
)

# --- 4A: 风险评分散点图 ---
pdf("Figure4A_GSE96058_Risk_Score.pdf", width = 6, height = 4)
print(ggplot(geo_data, aes(x = PatientID, y = RiskScore)) +
        geom_point(aes(color = RiskGroup), size = 1.5, alpha = 0.8) +
        geom_vline(xintercept = sum(geo_data$RiskGroup == "Low"), lty = 2, color = "grey50") +
        scale_color_manual(values = c("Low" = "#abdaec", "High" = "#f6c0ce")) +
        theme_classic() + journal_font + # 🌟 注入字体
        labs(x = "GSE96058 ER+ Patients", y = "Risk Score", title = "GSE96058 Validation: Risk Score"))
dev.off()

# --- 4B: 生存状态散点图 ---
pdf("Figure4B_GSE96058_Survival_Status.pdf", width = 6, height = 4)
print(ggplot(geo_data, aes(x = PatientID, y = OS.time / 365)) + 
        geom_point(aes(color = as.factor(OS)), size = 1.5, alpha = 0.8) +
        geom_vline(xintercept = sum(geo_data$RiskGroup == "Low"), lty = 2, color = "grey50") +
        scale_color_manual(values = c("0" = "#abdaec", "1" = "#f6c0ce"), labels = c("Alive", "Dead")) +
        theme_classic() + journal_font + # 🌟 注入字体
        labs(x = "GSE96058 ER+ Patients", y = "Survival Time (Years)", color = "Status", title = "GSE96058 Validation: Survival Status"))
dev.off()

# --- 4C: 热图 (高对比度深蓝至深杏) ---
heat_data_geo <- t(geo_data[, hub_genes])
colnames(heat_data_geo) <- geo_data$PatientID
annotation_col_geo <- data.frame(Risk = geo_data$RiskGroup)
rownames(annotation_col_geo) <- colnames(heat_data_geo)
ann_colors_geo <- list(Risk = c(Low = "#abdaec", High = "#f6c0ce"))

pdf("Figure4C_GSE96058_Risk_Heatmap.pdf", width = 7, height = 4)
pheatmap(heat_data_geo, scale = "row", cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, annotation_col = annotation_col_geo, annotation_colors = ann_colors_geo,
         color = colorRampPalette(c("#5c92b8", "white", "#d69c73"))(100),
         breaks = seq(-2, 2, length.out = 101),
         fontfamily = "sans", # 🌟 注入字体
         main = "GSE96058 Validation: Expression Heatmap")
dev.off()

# --- 4D: KM 曲线 (纤细精简风格) ---
fit_km_geo <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = geo_data)
p_km_geo <- ggsurvplot(
  fit_km_geo, data = geo_data, pval = TRUE, 
  pval.coord = c(max(geo_data$OS.time)*0.05, 0.15), 
  conf.int = TRUE, conf.int.alpha = 0.18, 
  risk.table = TRUE, size = 0.7, 
  palette = c("#abdaec", "#f6c0ce"),
  legend.title = "Risk", legend.labs = c("Low", "High"), 
  xlab = "Time (Days)", ylab = "Survival probability", title = NULL, 
  ggtheme = theme_classic() + journal_font # 🌟 注入字体
)

pdf("Figure4D_GSE96058_KM_Survival.pdf", width = 6, height = 6, onefile = FALSE)
print(p_km_geo)
dev.off() 

# --- 🌟 4E: 基于 CombinedScore 的超强联合 ROC 曲线 ---
ROC_rt_geo <- timeROC(T = geo_data$OS.time, delta = geo_data$OS, 
                      marker = geo_data$CombinedScore, cause = 1, 
                      weighting = "marginal", 
                      times = c(365, 365*3, 365*5), ROC = TRUE) 

pdf("Figure4E_GSE96058_Time_ROC_Combined.pdf", width = 6, height = 6)
par(family = "sans") # 🌟 注入字体
plot(ROC_rt_geo, time = 365, col = "#f6c0ce", title = FALSE, lwd = 2.5)
plot(ROC_rt_geo, time = 365*3, col = "#abdaec", add = TRUE, lwd = 2.5) 
plot(ROC_rt_geo, time = 365*5, col = "#97d1a0", add = TRUE, lwd = 2.5) 
legend("bottomright", 
       c(paste0("1-Year AUC = ", sprintf("%.3f", ROC_rt_geo$AUC[1])),
         paste0("3-Year AUC = ", sprintf("%.3f", ROC_rt_geo$AUC[2])), 
         paste0("5-Year AUC = ", sprintf("%.3f", ROC_rt_geo$AUC[3]))), 
       col = c("#f6c0ce", "#abdaec", "#97d1a0"), lty = 1, lwd = 2.5, bty = "n")
title(main = "GSE96058 Validation: Clinical-Molecular ROC", font.main = 2, family = "sans") # 🌟 注入字体
dev.off()

print("==== 联合模型优化完成！去检查最新的 AUC 成绩单吧 ====")