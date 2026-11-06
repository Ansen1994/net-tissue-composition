# ==============================================================================
# 基于 METABRIC 的外部验证：多因素临床-分子联合模型 (纯 ER+ 亚群)
# ==============================================================================
# [patched] setwd removed; cwd is work/
library(dplyr)
library(tibble)
library(survival)
library(survminer)
library(timeROC)
library(pheatmap)
library(ggplot2)
library(rms)
library(stringr)

# ==============================================================================
# 第一部分：读取数据与 Z-score 标准化
# ==============================================================================
print("1. 正在读取 METABRIC 表达矩阵与临床数据...")
# 跳过 cBioPortal 文件的头部注释行
clin_raw <- read.delim("../raw/METABRIC/data_clinical_patient.txt", comment.char = "#", stringsAsFactors = FALSE)
exp_raw <- read.delim("../raw/METABRIC/data_mrna_illumina_microarray.txt", stringsAsFactors = FALSE) 

print("2. 正在提取 Hub 基因并进行 Z-score 标准化...")
hub_genes <- c("PLG", "CD79A", "MMP7", "ANGPT1")

geo_exp_hub <- exp_raw %>%
  dplyr::filter(Hugo_Symbol %in% hub_genes) %>%
  dplyr::select(-2) %>% # 剔除 Entrez_Gene_Id
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), max)) %>%
  tibble::column_to_rownames("Hugo_Symbol") %>%
  t() %>%
  as.data.frame()

# 🌟 核心优化 1：Z-score 标准化
geo_exp_hub_scaled <- as.data.frame(scale(geo_exp_hub))
geo_exp_hub_scaled$PATIENT_ID <- gsub("\\.", "-", rownames(geo_exp_hub_scaled))


# ==============================================================================
# 第二部分：清洗临床数据与提取三大核心特征
# ==============================================================================
# ==============================================================================
# 第三部分：清洗临床数据与提取可用核心特征
# ==============================================================================
print("3. 正在清洗临床随访数据并筛选 ER 阳性样本...")

geo_clin_base <- clin_raw %>%
  dplyr::select(PATIENT_ID, 
                ER_Status = ER_IHC,             
                OS_Status_Raw = OS_STATUS,      
                OS.time = OS_MONTHS,
                Age = AGE_AT_DIAGNOSIS,
                LymphNode_raw = LYMPH_NODES_EXAMINED_POSITIVE) %>% 
  dplyr::filter(grepl("pos|positive", ER_Status, ignore.case = TRUE)) %>%  
  dplyr::mutate(
    OS = ifelse(grepl("DECEASED", OS_Status_Raw, ignore.case = TRUE), 1, 0),
    OS.time = as.numeric(OS.time) * 30, # 月转为天
    Age = as.numeric(Age),
    LymphNode = ifelse(as.numeric(LymphNode_raw) > 0, 1, 0) # 将淋巴结数目转化为二分类：0为阴性，1为阳性
  )

geo_final <- merge(geo_clin_base, geo_exp_hub_scaled, by = "PATIENT_ID")
geo_final <- na.omit(geo_final)

print(paste("---- 纯 ER+ 联合临床矩阵生成完毕，当前包含患者数量：", nrow(geo_final), "----"))


# ==============================================================================
# 第四部分：计算 Risk Score 与 临床-分子联合评分 (Combined Score)
# ==============================================================================
print("4. 正在构建多因素联合生存模型...")
geo_data <- geo_final

geo_data$RiskScore <- (0.1102 * geo_data$PLG) + (0.2295 * geo_data$ANGPT1) - 
  (0.0854 * geo_data$CD79A) - (0.0851 * geo_data$MMP7)

res.cut <- surv_cutpoint(geo_data, time = "OS.time", event = "OS", variables = "RiskScore")
res.cat <- surv_categorize(res.cut)
geo_data$RiskGroup <- factor(res.cat$RiskScore, levels = c("low", "high"), labels = c("Low", "High"))

geo_data <- geo_data[order(geo_data$RiskScore), ]
geo_data$PatientID <- 1:nrow(geo_data)

write.csv(geo_data, "Table_S6_METABRIC_Validation_Result.csv", row.names = FALSE)

# 🌟 构建多因素临床-分子联合模型（Age + LymphNode）
dd <- datadist(geo_data)
options(datadist = 'dd')
cox_combined <- cph(Surv(OS.time, OS) ~ RiskScore + Age + LymphNode, 
                    data = geo_data, x = TRUE, y = TRUE, surv = TRUE)
geo_data$CombinedScore <- predict(cox_combined, type = "lp")


# ==============================================================================
# 第五部分：生成图表 (马卡龙配色及无伪影优化)
# ==============================================================================
print("5. 正在生成系列高质量图表...")

# --- 4A: 风险评分散点图 ---
pdf("Figure4A_METABRIC_Risk_Score.pdf", width = 6, height = 4)
print(ggplot(geo_data, aes(x = PatientID, y = RiskScore)) +
        geom_point(aes(color = RiskGroup), size = 1.5, alpha = 0.8) +
        geom_vline(xintercept = sum(geo_data$RiskGroup == "Low"), lty = 2, color = "grey50") +
        scale_color_manual(values = c("Low" = "#abdaec", "High" = "#f6c0ce")) +
        theme_classic() +
        labs(x = "METABRIC ER+ Patients", y = "Risk Score", title = "METABRIC Validation: Risk Score") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")))
dev.off()

# --- 4B: 生存状态散点图 ---
pdf("Figure4B_METABRIC_Survival_Status.pdf", width = 6, height = 4)
print(ggplot(geo_data, aes(x = PatientID, y = OS.time / 365)) + 
        geom_point(aes(color = as.factor(OS)), size = 1.5, alpha = 0.8) +
        geom_vline(xintercept = sum(geo_data$RiskGroup == "Low"), lty = 2, color = "grey50") +
        scale_color_manual(values = c("0" = "#abdaec", "1" = "#f6c0ce"), labels = c("Alive", "Dead")) +
        theme_classic() +
        labs(x = "METABRIC ER+ Patients", y = "Survival Time (Years)", color = "Status", title = "METABRIC Validation: Survival Status") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")))
dev.off()

# --- 4C: 热图 ---
heat_data_geo <- t(geo_data[, c("PLG", "CD79A", "MMP7", "ANGPT1")])
colnames(heat_data_geo) <- geo_data$PatientID
annotation_col_geo <- data.frame(Risk = geo_data$RiskGroup)
rownames(annotation_col_geo) <- colnames(heat_data_geo)
ann_colors_geo <- list(Risk = c(Low = "#abdaec", High = "#f6c0ce"))

pdf("Figure4C_METABRIC_Risk_Heatmap.pdf", width = 7, height = 4)
pheatmap(heat_data_geo, scale = "row", cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, annotation_col = annotation_col_geo, annotation_colors = ann_colors_geo,
         color = colorRampPalette(c("#5c92b8", "white", "#d69c73"))(100),
         breaks = seq(-2, 2, length.out = 101),
         main = "METABRIC Validation: Expression Heatmap")
dev.off()

# --- 4D: KM 曲线 ---
fit_km_geo <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = geo_data)
p_km_geo <- ggsurvplot(
  fit_km_geo, data = geo_data, pval = TRUE, 
  pval.coord = c(200, 0.15), 
  conf.int = TRUE, conf.int.alpha = 0.18, 
  risk.table = TRUE, size = 0.8, 
  censor = FALSE,                        
  xlim = c(0, 6000),                     
  break.time.by = 1500,
  palette = c("#abdaec", "#f6c0ce"),
  legend.title = "Risk", legend.labs = c("Low", "High"), 
  xlab = "Time (Days)", ylab = "Survival probability", title = NULL, ggtheme = theme_classic()
)
pdf("Figure4D_METABRIC_KM_Survival.pdf", width = 6, height = 6, onefile = FALSE)
print(p_km_geo)
dev.off() 

# --- 4E: 联合临床特征的 TimeROC 曲线 ---
ROC_rt_geo <- timeROC(T = geo_data$OS.time, delta = geo_data$OS, 
                      marker = geo_data$CombinedScore, cause = 1,  
                      weighting = "marginal", 
                      times = c(365, 365*3, 365*5), ROC = TRUE) 

pdf("Figure4E_METABRIC_Time_ROC_Combined.pdf", width = 6, height = 6)
plot(ROC_rt_geo, time = 365, col = "#f6c0ce", title = FALSE, lwd = 2.5)
plot(ROC_rt_geo, time = 365*3, col = "#abdaec", add = TRUE, lwd = 2.5) 
plot(ROC_rt_geo, time = 365*5, col = "#97d1a0", add = TRUE, lwd = 2.5) 
legend("bottomright", 
       c(paste0("1-Year AUC = ", sprintf("%.3f", ROC_rt_geo$AUC[1])),
         paste0("3-Year AUC = ", sprintf("%.3f", ROC_rt_geo$AUC[2])), 
         paste0("5-Year AUC = ", sprintf("%.3f", ROC_rt_geo$AUC[3]))), 
       col = c("#f6c0ce", "#abdaec", "#97d1a0"), lty = 1, lwd = 2.5, bty = "n")
title(main = "METABRIC Validation: Clinical-Molecular ROC", font.main = 2)
dev.off()

print("==== METABRIC 全套高定验证图表运行完毕！ ====")