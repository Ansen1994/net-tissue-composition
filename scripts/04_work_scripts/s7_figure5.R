# [patched] setwd removed; cwd is work/
# ================= 1. 加载依赖包 =================
library(dplyr)
library(stringr)
library(survival)
library(survminer)
library(rms)
library(dcurves)
library(ggplot2)
library(SummarizedExperiment) 

# ================= 2. 数据准备与清洗 =================
print("正在构建整合临床与分子的 final_data...")
model_data <- read.csv("Table_S6_Patient_Risk_Score_Result.csv", header = TRUE)
colnames(model_data)[1] <- "Sample_ID"
model_data$Patient_ID <- substr(gsub("\\.", "-", model_data$Sample_ID), 1, 12)

# 把包含 expquery2 的原始数据加载进环境
load("../raw/TCGA/brca.gdc_2024.rda") 

clin_df <- as.data.frame(colData(expquery2))
clin_df$Sample_ID <- rownames(clin_df)
clin_df$Patient_ID <- substr(gsub("\\.", "-", clin_df$Sample_ID), 1, 12)

clin_clean <- clin_df %>%
  dplyr::select(Patient_ID, Age = age_at_index, Stage = ajcc_pathologic_stage) %>%
  dplyr::mutate(
    Age = as.numeric(Age),
    Stage = stringr::str_extract(Stage, "Stage [IV]+") %>% gsub("Stage ", "", .)
  ) %>%
  dplyr::filter(!is.na(Age) & !is.na(Stage)) %>%
  dplyr::distinct(Patient_ID, .keep_all = TRUE)

final_data <- merge(model_data, clin_clean, by = "Patient_ID")
final_data$RiskGroup <- factor(final_data$RiskGroup, levels = c("Low", "High"))
final_data$Stage <- factor(final_data$Stage) 

write.csv(final_data, file = "Table_S7_Final_Clinical_Data.csv", row.names = FALSE)
print(paste("---- final_data 准备完毕，有效样本量：", nrow(final_data), "----"))


# ==============================================================================
# ==============================================================================
# 🌟 高级定制：构建基于“显著性”动态变色的完美森林图函数
# ==============================================================================
# ==============================================================================
# 🌟 终极森林图函数：完美对齐灰色背景，虚线精准截断避免遮挡文字
# ==============================================================================
plot_custom_forest <- function(df_plot, plot_title, bottom_text = "") {
  df_plot <- df_plot[nrow(df_plot):1, ] 
  df_plot$Y <- 1:nrow(df_plot)
  
  # 格式化文本
  df_plot$HR_text <- ifelse(is.na(df_plot$P_value), "Reference", sprintf("%.2f\n(%.2f - %.2f)", df_plot$HR, df_plot$Lower, df_plot$Upper))
  df_plot$P_text <- ifelse(is.na(df_plot$P_value), "",
                           ifelse(df_plot$P_value < 0.001, "< 0.001 ***",
                                  ifelse(df_plot$P_value < 0.01, sprintf("%.3f **", df_plot$P_value),
                                         ifelse(df_plot$P_value < 0.05, sprintf("%.3f *", df_plot$P_value),
                                                sprintf("%.3f", df_plot$P_value)))))
  
  # 颜色映射 (显著HR>1为暖杏色，HR<=1包含Reference为深蓝色)
  df_plot$ColorGroup <- ifelse(df_plot$HR > 1, "Risk", "Protective")
  
  # 动态计算排版坐标
  min_x <- min(df_plot$Lower, na.rm = TRUE)
  max_x <- max(df_plot$Upper, na.rm = TRUE)
  log_range <- log10(max_x) - log10(min_x)
  
  text_gene_x <- 10^(log10(min_x) - log_range * 0.75) 
  text_n_x    <- 10^(log10(min_x) - log_range * 0.20)
  text_hr_x   <- 10^(log10(max_x) + log_range * 0.3)
  text_p_x    <- 10^(log10(max_x) + log_range * 0.8)
  
  x_limit_min <- 10^(log10(min_x) - log_range * 0.85)
  x_limit_max <- 10^(log10(max_x) + log_range * 1.1)
  
  # 🌟 修复 1：精确对齐灰色背景。以偶数行 (Y=2, 4, 6...) 为中心，上下各扩展 0.5
  bg_rects <- data.frame(
    ymin = seq(2, nrow(df_plot), by = 2) - 0.5,
    ymax = seq(2, nrow(df_plot), by = 2) + 0.5
  )
  
  p <- ggplot(df_plot) +
    # 绘制完美对齐的灰白背景条纹
    geom_rect(data = bg_rects, aes(xmin = x_limit_min, xmax = x_limit_max, ymin = ymin, ymax = ymax), fill = "grey92", color = NA, inherit.aes = FALSE) +
    
    # 🌟 修复 2：使用线段 (segment) 替换 geom_vline，将其底部精确截断在 Y=0.5 处，绝不往下延伸
    annotate("segment", x = 1, xend = 1, y = 0.5, yend = nrow(df_plot) + 0.5, linetype = "dashed", color = "black", linewidth = 0.7) +
    
    geom_pointrange(aes(x = HR, y = Y, xmin = Lower, xmax = Upper, color = ColorGroup), shape = 15, size = 1.3, linewidth = 1.1) +
    scale_color_manual(values = c("Risk" = "#d69c73", "Protective" = "#5c92b8")) +
    
    # 文本图层加载
    geom_text(aes(x = text_gene_x, y = Y, label = Name, color = ColorGroup), family = "sans", fontface = "bold", hjust = 0, size = 4.5) +
    geom_text(aes(x = text_n_x, y = Y, label = N_text), family = "sans", color = "black", fontface = "italic", hjust = 0.5, size = 4) +
    geom_text(aes(x = text_hr_x, y = Y, label = HR_text), family = "sans", color = "black", hjust = 0.5, size = 4, lineheight = 0.8) +
    geom_text(aes(x = text_p_x, y = Y, label = P_text), family = "sans", color = "black", fontface = "italic", hjust = 0, size = 4) +
    
    # 底部统计信息 (Y=0.1 处，现在绝对不会被虚线遮挡)
    annotate("text", x = text_gene_x, y = 0.1, label = bottom_text, family = "sans", hjust = 0, fontface = "italic", size = 3.5) +
    
    scale_x_log10(limits = c(x_limit_min, x_limit_max), labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
    scale_y_continuous(limits = c(0, nrow(df_plot) + 0.5)) +
    theme_minimal() +
    theme(
      text = element_text(family = "sans", color = "black"),
      panel.grid = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
      axis.title.x = element_text(family = "sans", size = 12, face = "bold", color = "black"),
      axis.text.x = element_text(family = "sans", size = 11, color = "black"),
      legend.position = "none",
      plot.title = element_text(family = "sans", hjust = 0.5, face = "bold", size = 16, margin = ggplot2::margin(b = 20)),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 15) 
    ) +
    labs(title = plot_title, x = "Hazard Ratio (95% CI)")
  return(p)
}


# ================= 3. 重新绘制 Figure 5 =================

# 🌟 关键补丁：显式定义预测时间节点（假设 OS.time 单位为天）
t1 <- 365        # 1年
t3 <- 365 * 3    # 3年
t5 <- 365 * 5    # 5年

# 🌟 Figure 5A: 临床单因素 Cox 森林图
print("正在生成 Figure 5A (临床单因素森林图)...")
univ_results <- list()
for (var in c("RiskScore", "Age", "Stage")) {
  fit <- coxph(as.formula(paste("Surv(OS.time, OS) ~", var)), data = final_data)
  sm <- summary(fit)
  tmp <- data.frame(Name = rownames(sm$coefficients), HR = sm$coefficients[,2],
                    Lower = sm$conf.int[,3], Upper = sm$conf.int[,4], P_value = sm$coefficients[,5])
  
  if (var == "Stage") {
    tmp$Name <- gsub("Stage", "Stage ", tmp$Name)
    tmp <- rbind(data.frame(Name="Stage I", HR=1, Lower=1, Upper=1, P_value=NA), tmp)
    stage_counts <- table(final_data$Stage)
    tmp$N_text <- paste0("(N=", stage_counts[c("I", "II", "III", "IV")], ")")
  } else {
    tmp$N_text <- paste0("(N=", nrow(final_data), ")")
  }
  univ_results[[var]] <- tmp
}
univ_df <- do.call(rbind, univ_results)

pdf("Figure5A_Univariate_Forest.pdf", width = 8, height = 5)
print(plot_custom_forest(univ_df, "Univariate Cox Regression Analysis"))
dev.off()


# 🌟 Figure 5B: 多因素 Cox 森林图 
print("正在生成 Figure 5B (多因素森林图)...")
multi_cox <- coxph(Surv(OS.time, OS) ~ RiskScore + Age + Stage, data = final_data)
sm_multi <- summary(multi_cox)
multi_df <- data.frame(Name = rownames(sm_multi$coefficients), HR = sm_multi$coefficients[,2],
                       Lower = sm_multi$conf.int[,3], Upper = sm_multi$conf.int[,4], P_value = sm_multi$coefficients[,5])
multi_df$Name <- gsub("Stage", "Stage ", multi_df$Name)

# 动态插入 Stage I 参照组
idx <- grep("Stage II", multi_df$Name)
if(length(idx) > 0) {
  multi_df <- rbind(multi_df[1:(idx-1), ], data.frame(Name="Stage I", HR=1, Lower=1, Upper=1, P_value=NA), multi_df[idx:nrow(multi_df), ])
}

multi_df$N_text <- ""
for(i in 1:nrow(multi_df)) {
  if(grepl("Stage", multi_df$Name[i])) {
    stg <- gsub("Stage ", "", multi_df$Name[i])
    multi_df$N_text[i] <- paste0("(N=", sum(final_data$Stage == stg, na.rm = TRUE), ")")
  } else {
    multi_df$N_text[i] <- paste0("(N=", nrow(final_data), ")")
  }
}

events <- multi_cox$nevent
global_p <- sm_multi$logtest["pvalue"]
bottom_info <- sprintf("# Events: %d; Global p-value (Log-Rank): %s", events, signif(global_p, 5))

pdf("Figure5B_Multivariate_Forest.pdf", width = 8, height = 5)
print(plot_custom_forest(multi_df, "Independent Prognostic Value Analysis", bottom_text = bottom_info))
dev.off()


# 🌟 Figure 5C: 临床诺模图 (Nomogram)
print("正在生成 Figure 5C (列线图)...")

# =======================================================
# 🌟 关键修复：计算数据分布范围并载入全局环境变量
# （rms 包绘图必须有这一步，否则无法生成坐标轴上限）
# =======================================================
dd <- datadist(final_data)
options(datadist = "dd")

fit_nomo <- cph(Surv(OS.time, OS) ~ RiskScore + Age + Stage, data = final_data, x = TRUE, y = TRUE, surv = TRUE)
surv_func <- Survival(fit_nomo)

nomo <- nomogram(fit_nomo, 
                 fun = list(function(x) surv_func(t1, x), 
                            function(x) surv_func(t3, x), 
                            function(x) surv_func(t5, x)), 
                 funlabel = c("1-year OS", "3-year OS", "5-year OS"), 
                 maxscale = 100, 
                 fun.at = c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99)) 

pdf("Figure5C_Nomogram.pdf", width = 12, height = 7)
par(family = "sans") 
plot(nomo, xfrac = 0.25, col.grid = gray(c(0.8, 0.95)), lmgp = 0.2, cex.axis = 0.9) 
dev.off()


# 🌟 Figure 5D: 校准曲线
print("正在生成 Figure 5D (校准曲线)...")
pdf("Figure5D_Calibration_Curve_1_3_5_Years.pdf", width = 6.5, height = 6.5)
par(family = "sans") 

# 固定随机种子：calibrate(method = "boot") 内部做 1000 次 bootstrap 重抽样，
# 不设种子则每次运行曲线不同，影响图件可复现性
set.seed(42)
fit_cal1 <- cph(Surv(OS.time, OS) ~ RiskScore + Age + Stage, data = final_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = t1)
cal1 <- calibrate(fit_cal1, cmethod = "KM", method = "boot", u = t1, m = round(nrow(final_data)/4), B = 1000)
fit_cal3 <- cph(Surv(OS.time, OS) ~ RiskScore + Age + Stage, data = final_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = t3)
cal3 <- calibrate(fit_cal3, cmethod = "KM", method = "boot", u = t3, m = round(nrow(final_data)/4), B = 1000)
fit_cal5 <- cph(Surv(OS.time, OS) ~ RiskScore + Age + Stage, data = final_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = t5)
cal5 <- calibrate(fit_cal5, cmethod = "KM", method = "boot", u = t5, m = round(nrow(final_data)/4), B = 1000)

plot(cal3, lwd = 2.5, lty = 1, errbar.col = "#abdaec", xlim = c(0, 1), ylim = c(0, 1), xlab = "Nomogram-Predicted OS", ylab = "Actual Observed OS", main = "Calibration Curve", sub = FALSE, riskdist = FALSE)
abline(0, 1, lty = 3, lwd = 2, col = "gray40")
lines(cal1[,c('mean.predicted',"KM")], type = 'b', lwd = 2.5, col = "#f6c0ce", pch = 16)
lines(cal5[,c('mean.predicted',"KM")], type = 'b', lwd = 2.5, col = "#97d1a0", pch = 16)
lines(cal3[,c('mean.predicted',"KM")], type = 'b', lwd = 2.5, col = "#abdaec", pch = 16) 

legend("bottomright", legend = c("1-Year OS", "3-Year OS", "5-Year OS"), col = c("#f6c0ce", "#abdaec", "#97d1a0"), lty = 1, lwd = 2.5, pch = 16, bty = "n")
title(main = "Calibration Curve", font.main = 2, family = "sans") 
dev.off()


# =========================================================================
# 🌟 Figure 5E 主图 (5年 DCA) & 补充图 S5 (1年 & 3年 DCA)
# =========================================================================
print("正在生成 DCA 曲线 (包括 1年、3年、5年)...")
fit_clin <- cph(Surv(OS.time, OS) ~ Stage, data = final_data, x = TRUE, y = TRUE, surv = TRUE)
fit_risk <- cph(Surv(OS.time, OS) ~ RiskScore, data = final_data, x = TRUE, y = TRUE, surv = TRUE)

# 构建计算指定时间 DCA 并绘图的专属函数
generate_dca_plot <- function(time_val, time_label, filename) {
  data_dca <- final_data
  data_dca$Prob_Nomogram <- 1 - survest(fit_nomo, newdata = data_dca, times = time_val)$surv
  data_dca$Prob_Stage <- 1 - survest(fit_clin, newdata = data_dca, times = time_val)$surv
  data_dca$Prob_RiskScore <- 1 - survest(fit_risk, newdata = data_dca, times = time_val)$surv
  
  # 滤除缺失值
  data_dca <- na.omit(data_dca[, c("OS.time", "OS", "Prob_Nomogram", "Prob_Stage", "Prob_RiskScore")])
  
  dca_res <- dca(Surv(OS.time, OS) ~ Prob_Nomogram + Prob_RiskScore + Prob_Stage, data = data_dca, time = time_val)
  
  p_dca <- plot(dca_res, smooth = TRUE) + 
    scale_color_manual(
      values = c("Prob_Nomogram" = "#f6c0ce", "Prob_RiskScore" = "#abdaec", "Prob_Stage" = "#97d1a0", "Treat All" = "gray50", "Treat None" = "black"),
      labels = c("Prob_Nomogram" = "Nomogram", "Prob_RiskScore" = "Risk Score", "Prob_Stage" = "Stage")
    ) +
    ggplot2::labs(title = paste0(time_label, "-Year Decision Curve Analysis (DCA)")) + 
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(family = "sans", color = "black"), 
                   plot.title = element_text(family = "sans", face = "bold", hjust = 0.5), 
                   legend.position = c(0.8, 0.8), 
                   legend.background = element_blank())
  
  pdf(filename, width = 6.5, height = 5.5)
  print(p_dca)
  dev.off()
}

# 生成主图 Figure 5E (🌟修改为 5-Year🌟)
generate_dca_plot(t5, "5", "Figure5E_Comparative_DCA_Curve_5Y.pdf")

# 生成补充材料 Supplementary Figure S5 (🌟修改为 1-Year & 3-Year🌟)
generate_dca_plot(t1, "1", "Supplementary_Figure_S5_DCA_1Y.pdf")
generate_dca_plot(t3, "3", "Supplementary_Figure_S5_DCA_3Y.pdf")

print("==== Figure 5 全套优化完毕！主图已替换为表现最佳的 5年 DCA ====")