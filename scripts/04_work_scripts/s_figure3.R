# =====================================================================
# Figure 3: 单/多因素 Cox 回归与预后模型评估 (全套顶刊 Sans 字体版)
# =====================================================================
library(survival)
library(survminer)
library(timeROC)
library(pheatmap)
library(dplyr)
library(ggplot2)

graphics.off() 

# ================= 1. 数据准备 =================
rt <- read.csv("Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv", row.names = 1)
hub_genes <- c("PLG", "CD79A", "MMP7", "ANGPT1")

# 提取生存数据与核心靶标表达量
model_data <- rt[, c("OS.time", "OS", hub_genes)]

# ================= 2. 单因素 Cox 回归计算 (🌟 核心修复缺失逻辑) =================
print("正在计算单因素 Cox 回归...")
univ_cox_list <- lapply(hub_genes, function(gene) {
  formula_uni <- as.formula(paste("Surv(OS.time, OS) ~", gene))
  fit <- coxph(formula_uni, data = model_data)
  summary_fit <- summary(fit)
  data.frame(
    Gene = gene,
    HR = summary_fit$coefficients[1, "exp(coef)"],
    Lower_95CI = summary_fit$conf.int[1, "lower .95"],
    Upper_95CI = summary_fit$conf.int[1, "upper .95"],
    P_value = summary_fit$coefficients[1, "Pr(>|z|)"]
  )
})
univ_cox_df <- do.call(rbind, univ_cox_list)

# ================= 3. 多因素 Cox 回归计算 =================
print("正在计算多因素 Cox 回归...")
formula_str <- paste("Surv(OS.time, OS) ~", paste(hub_genes, collapse = " + "))
cox_model <- coxph(as.formula(formula_str), data = model_data)

cox_summary <- summary(cox_model)
cox_results <- data.frame(
  Gene = rownames(cox_summary$coefficients),
  Coef = cox_summary$coefficients[, "coef"],            
  HR = cox_summary$coefficients[, "exp(coef)"],        
  P_value = cox_summary$coefficients[, "Pr(>|z|)"]     
)
write.csv(cox_results, "Table_S5_Cox_Coefficients.csv", row.names = FALSE)


# =====================================================================
# 绘图部分
# =====================================================================

# ------------------------------------------------------------------------------
# 3A: 单因素 Cox 森林图 (Univariate Cox Forest Plot)
# ------------------------------------------------------------------------------
plot_uni_df <- univ_cox_df[nrow(univ_cox_df):1, ]
plot_uni_df$Y <- 1:nrow(plot_uni_df)
plot_uni_df$HR_text <- sprintf("%.2f (%.2f - %.2f)", plot_uni_df$HR, plot_uni_df$Lower_95CI, plot_uni_df$Upper_95CI)
plot_uni_df$P_text <- ifelse(plot_uni_df$P_value < 0.001, "< 0.001 ***",
                             ifelse(plot_uni_df$P_value < 0.01, sprintf("%.3f **", plot_uni_df$P_value),
                                    ifelse(plot_uni_df$P_value < 0.05, sprintf("%.3f *", plot_uni_df$P_value),
                                           sprintf("%.3f", plot_uni_df$P_value))))
plot_uni_df$ColorGroup <- ifelse(plot_uni_df$HR > 1, "Risk", "Protective")

p_uni_forest <- ggplot(plot_uni_df) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.8) +
  geom_pointrange(aes(x = HR, y = Y, xmin = Lower_95CI, xmax = Upper_95CI, color = ColorGroup),
                  shape = 15, size = 1.3, linewidth = 1.3) +
  scale_color_manual(values = c("Risk" = "#d69c73", "Protective" = "#5c92b8")) +
  # 🌟 字体升级：geom_text 必须显式声明 family = "sans"
  geom_text(aes(x = 0.5, y = Y, label = Gene, color = ColorGroup), family = "sans", fontface = "bold", hjust = 0, size = 4.5) +
  geom_text(aes(x = 1.6, y = Y, label = HR_text), family = "sans", color = "black", hjust = 0, size = 4) +
  geom_text(aes(x = 2.3, y = Y, label = P_text), family = "sans", color = "black", fontface = "italic", hjust = 0, size = 4) +
  scale_x_log10(limits = c(0.4, 2.8)) +
  scale_y_continuous(limits = c(0.5, nrow(plot_uni_df) + 0.5)) +
  theme_minimal() +
  theme(
    text = element_text(family = "sans", color = "black"), # 全局字体
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(family = "sans", size = 12, face = "bold", color = "black"),
    axis.text.x = element_text(family = "sans", size = 11, color = "black"),
    legend.position = "none",
    plot.title = element_text(family = "sans", hjust = 0.5, face = "bold", size = 16, margin = ggplot2::margin(b = 20))
  ) +
  labs(title = "Univariate Cox Regression Analysis", x = "Hazard Ratio (95% CI)")

ggsave("Figure_Univariate_Cox_Forest.pdf", plot = p_uni_forest, width = 7.5, height = 4)


# ------------------------------------------------------------------------------
# 3B: 多因素 Cox 森林图 (Multivariate Cox Forest Plot)
# ------------------------------------------------------------------------------
cox_df <- data.frame(
  Gene = rownames(cox_summary$coefficients),
  HR = cox_summary$coefficients[, "exp(coef)"],
  Lower = cox_summary$conf.int[, "lower .95"],
  Upper = cox_summary$conf.int[, "upper .95"],
  P_value = cox_summary$coefficients[, "Pr(>|z|)"]
)
cox_df <- cox_df[nrow(cox_df):1, ]
cox_df$Y <- 1:nrow(cox_df)

n_patients <- nrow(model_data)
cox_df$N_text <- paste0("(N=", n_patients, ")")
cox_df$HR_text <- sprintf("%.2f\n(%.2f - %.2f)", cox_df$HR, cox_df$Lower, cox_df$Upper)
cox_df$P_text <- ifelse(cox_df$P_value < 0.001, "< 0.001 ***",
                        ifelse(cox_df$P_value < 0.01, sprintf("%.3f **", cox_df$P_value),
                               ifelse(cox_df$P_value < 0.05, sprintf("%.3f *", cox_df$P_value),
                                      sprintf("%.3f", cox_df$P_value))))
events <- cox_model$nevent
global_p <- summary(cox_model)$logtest["pvalue"]
bottom_text <- sprintf("# Events: %d; Global p-value (Log-Rank): %s", events, signif(global_p, 5))

cox_df$ColorGroup <- ifelse(cox_df$HR > 1, "Risk", "Protective")

p_forest_multi <- ggplot(cox_df) +
  geom_rect(data = data.frame(ymin = c(0.5, 2.5), ymax = c(1.5, 3.5)),
            aes(xmin = 0.28, xmax = 2.2, ymin = ymin, ymax = ymax), fill = "grey90", color = NA) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_pointrange(aes(x = HR, y = Y, xmin = Lower, xmax = Upper, color = ColorGroup), 
                  shape = 15, size = 1.3, linewidth = 1.3) +
  scale_color_manual(values = c("Risk" = "#d69c73", "Protective" = "#5c92b8")) +
  # 🌟 字体升级：强制应用 Sans
  geom_text(aes(x = 0.3, y = Y, label = Gene, color = ColorGroup), family = "sans", fontface = "bold", hjust = 0, size = 4.5) +
  geom_text(aes(x = 0.5, y = Y, label = N_text), family = "sans", color = "black", fontface = "italic", hjust = 0.5, size = 4) +
  geom_text(aes(x = 0.7, y = Y, label = HR_text), family = "sans", color = "black", hjust = 0.5, size = 4, lineheight = 0.8) +
  geom_text(aes(x = 2.0, y = Y, label = P_text), family = "sans", color = "black", fontface = "italic", hjust = 1, size = 4) +
  annotate("text", x = 0.3, y = 0.1, label = bottom_text, family = "sans", hjust = 0, fontface = "italic", size = 3.5) +
  scale_x_log10(breaks = c(0.9, 1.0, 1.1, 1.2, 1.3, 1.4), limits = c(0.28, 2.2)) +
  scale_y_continuous(limits = c(0, 4.5)) +
  theme_minimal() +
  theme(
    text = element_text(family = "sans", color = "black"),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(family = "sans", size = 11, color = "black"),
    legend.position = "none",
    plot.title = element_text(family = "sans", hjust = 0.5, face = "bold", size = 16, margin = ggplot2::margin(b = 20))
  ) +
  labs(title = "Multivariate Cox Regression")

ggsave("Figure3B_Cox_Forest_Multi.pdf", plot = p_forest_multi, width = 8, height = 4.5)

# ------------------------------------------------------------------------------
# 风险评分分组 (Risk Score)
# ------------------------------------------------------------------------------
model_data$RiskScore <- predict(cox_model, type = "risk") 
res.cut <- median(model_data$RiskScore)
model_data$RiskGroup <- factor(ifelse(model_data$RiskScore > res.cut, "High", "Low"), levels = c("Low", "High"))
model_data <- model_data[order(model_data$RiskScore), ]
model_data$PatientID <- 1:nrow(model_data)
write.csv(model_data, "Table_S6_Patient_Risk_Score_Result.csv", row.names = TRUE)

# ------------------------------------------------------------------------------
# 3C: 风险评分分布曲线
# ------------------------------------------------------------------------------
pdf("Figure3C_Risk_Score_Curve.pdf", width = 6, height = 4)
p_curve <- ggplot(model_data, aes(x = PatientID, y = RiskScore)) +
  geom_point(aes(color = RiskGroup), size = 1.5, alpha = 0.8) +
  geom_vline(xintercept = sum(model_data$RiskGroup == "Low"), lty = 2, color = "grey50") +
  scale_color_manual(values = c("Low" = "#abdaec", "High" = "#f6c0ce")) +
  theme_classic() +
  # 🌟 字体升级
  theme(
    text = element_text(family = "sans", color = "black"),
    plot.title = element_text(family = "sans", hjust = 0.5, face = "bold")
  ) +
  labs(x = "Patients (increasing risk score)", y = "Risk Score", title = "Risk Score Distribution")
print(p_curve)
dev.off()

# ------------------------------------------------------------------------------
# 3D: 生存状态散点图
# ------------------------------------------------------------------------------
pdf("Figure3D_Survival_Status.pdf", width = 6, height = 4)
p_status <- ggplot(model_data, aes(x = PatientID, y = OS.time / 365)) + 
  geom_point(aes(color = as.factor(OS)), size = 1.5, alpha = 0.8) +
  geom_vline(xintercept = sum(model_data$RiskGroup == "Low"), lty = 2, color = "grey50") +
  scale_color_manual(values = c("0" = "#abdaec", "1" = "#f6c0ce"), labels = c("Alive", "Dead")) +
  theme_classic() +
  # 🌟 字体升级
  theme(
    text = element_text(family = "sans", color = "black"),
    plot.title = element_text(family = "sans", hjust = 0.5, face = "bold")
  ) +
  labs(x = "Patients (increasing risk score)", y = "Survival Time (Years)", color = "Status", title = "Survival Status")
print(p_status)
dev.off()

# ------------------------------------------------------------------------------
# 3E: 核心基因表达热图
# ------------------------------------------------------------------------------
heat_data <- t(model_data[, hub_genes])
annotation_col <- data.frame(Risk = model_data$RiskGroup)
rownames(annotation_col) <- colnames(heat_data)
ann_colors <- list(Risk = c(Low = "#abdaec", High = "#fbd6be"))

pdf("Figure3E_Risk_Heatmap_Fixed.pdf", width = 7, height = 4)
pheatmap(heat_data, scale = "row", cluster_cols = FALSE, cluster_rows = FALSE,
         show_colnames = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors,
         color = colorRampPalette(c("#5c92b8", "white", "#d69c73"))(100),
         breaks = seq(-2, 2, length.out = 101),
         fontfamily = "sans", # 🌟 字体升级：pheatmap 专用参数
         main = "Signature Expression Heatmap")
dev.off()

# ------------------------------------------------------------------------------
# 3F: Kaplan-Meier 生存曲线
# ------------------------------------------------------------------------------
fit_km <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = model_data)

p_km_thin <- ggsurvplot(
  fit_km, data = model_data, pval = TRUE, pval.coord = c(200, 0.25), 
  conf.int = TRUE, conf.int.alpha = 0.18, risk.table = TRUE, size = 0.7,  
  palette = c("#abdaec", "#f6c0ce"),
  legend.title = "Risk", legend.labs = c("Low", "High"),
  xlab = "Time (Days)", ylab = "Survival probability",
  # 🌟 字体升级：注入 ggplot2 标准字体主题
  ggtheme = theme_classic() + theme(text = element_text(family = "sans", color = "black"))
)

pdf("Figure3F_KM_Survival_Thin.pdf", width = 6, height = 6, onefile = FALSE)
print(p_km_thin)
dev.off()

# ------------------------------------------------------------------------------
# 3G: 时间依赖性 ROC 曲线 (Base R 绘图)
# ------------------------------------------------------------------------------
ROC_rt <- timeROC(T = model_data$OS.time, delta = model_data$OS, 
                  marker = model_data$RiskScore, cause = 1, 
                  weighting = "marginal", 
                  times = c(365, 365*3, 365*5), ROC = TRUE)

pdf("Figure3G_Time_ROC.pdf", width = 6, height = 6)
par(family = "sans") # 🌟 字体升级：Base R 绘图系统强行指定字体
plot(ROC_rt, time = 365, col = "#f6c0ce", title = FALSE, lwd = 2.5)
plot(ROC_rt, time = 365*3, col = "#abdaec", add = TRUE, lwd = 2.5)
plot(ROC_rt, time = 365*5, col = "#97d1a0", add = TRUE, lwd = 2.5)
legend("bottomright", 
       c(paste0("1-Year AUC = ", sprintf("%.3f", ROC_rt$AUC[1])),
         paste0("3-Year AUC = ", sprintf("%.3f", ROC_rt$AUC[2])),
         paste0("5-Year AUC = ", sprintf("%.3f", ROC_rt$AUC[3]))),
       col = c("#f6c0ce", "#abdaec", "#97d1a0"), lty = 1, lwd = 2.5, bty = "n")
title(main = "Time-dependent ROC", font.main = 2, family = "sans")
dev.off()

print("==== Figure 3 全套图表生成完毕，逻辑已修复，顶刊字体已全面覆盖！ ====")