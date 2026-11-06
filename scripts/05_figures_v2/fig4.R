# ==============================================================================
# Figure 4 D–F
#   D 四种标准化方案 × 四队列 C-index 热图
#   E 增量价值（临床 / 签名 / 合并）+ ΔC + LRT p
#   F 签名本身的时间依赖 ROC（1/3/5 年）
# 输入：Table_S9_Scaling_Scheme_Cindex.csv、Table_S10a/b、Table_S6_*.csv
# ==============================================================================
source("figs_v2/style.R")

s9 <- read.csv("Table_S9_Scaling_Scheme_Cindex.csv")
s10 <- read.csv("Table_S10a_Incremental_Value.csv")

s9$scheme <- factor(s9$scheme, levels = c("S0","S1","S2","S3"))
s9$cohort <- factor(s9$cohort, levels = rev(c("TCGA","GSE20685","METABRIC","GSE96058")))
s9$lab <- sprintf("%.3f", s9$C_index)

# ---------------- A–C 三队列 KM（surv_cutpoint 分组，与表内 RiskGroup 一致）----------------
km_files <- c(GSE20685 = "Table_S6_GEO_Validation_Result.csv",
              METABRIC = "Table_S6_METABRIC_Validation_Result.csv",
              GSE96058 = "Table_S6_GSE96058_Validation_Result.csv")

make_km <- function(nm, show_legend = FALSE) {
  d <- read.csv(km_files[[nm]], check.names = FALSE)
  d <- d[!is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0, ]
  d <- d[complete.cases(d[, c("OS", "OS.time", "RiskScore", "RiskGroup")]), ]
  d$RiskGroup <- factor(d$RiskGroup, levels = c("Low", "High"))
  nL <- sum(d$RiskGroup == "Low"); nH <- sum(d$RiskGroup == "High")

  fit <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = d)
  kd <- data.frame(time = fit$time / 365.25, surv = fit$surv,
                   grp = rep(c("Low", "High"), fit$strata))
  kd <- rbind(data.frame(time = 0, surv = 1, grp = c("Low", "High")), kd)
  kd <- kd[order(kd$grp, kd$time), ]

  lr <- survdiff(Surv(OS.time, OS) ~ RiskGroup, data = d)
  pv <- 1 - pchisq(lr$chisq, df = 1)
  cx <- summary(coxph(Surv(OS.time, OS) ~ RiskGroup, data = d))
  hr <- sprintf("HR = %.2f (%.2f\u2013%.2f)", cx$conf.int[1], cx$conf.int[3], cx$conf.int[4])
  pl <- if (pv < 1e-4) "Log-rank p < 0.0001" else sprintf("Log-rank p = %.4f", pv)

  ggplot(kd, aes(time, surv, colour = grp)) +
    geom_step(linewidth = 0.5) +
    scale_colour_manual(values = c(Low = "#4C9BC4", High = "#D6604D"),
                        labels = c(sprintf("Low (n = %s)", format(nL, big.mark = ",")),
                                   sprintf("High (n = %s)", format(nH, big.mark = ","))),
                        name = NULL) +
    scale_x_continuous(limits = c(0, max(kd$time) * 1.02),
                       breaks = pretty(c(0, max(kd$time)), 4),
                       expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25),
                       expand = c(0, 0), labels = function(x) sprintf("%.2f", x)) +
    labs(x = "Time (years)", y = "Overall survival", title = nm) +
    annotate("text", x = max(kd$time) * 0.04, y = 0.145, hjust = 0, size = 2.3, label = pl) +
    annotate("text", x = max(kd$time) * 0.04, y = 0.055, hjust = 0, size = 2.3, label = hr) +
    theme_pub(7.5) +
    theme(legend.position = if (show_legend) c(0.80, 0.86) else "none",
          legend.key.size = unit(2.6, "mm"),
          legend.background = element_rect(fill = "white", colour = NA))
}

pA <- tag_panel(make_km("GSE20685", TRUE), "A")
pB <- tag_panel(make_km("METABRIC"), "B")
pC <- tag_panel(make_km("GSE96058"), "C")

# ---------------- D ----------------
pD <- ggplot(s9, aes(scheme, cohort, fill = C_index)) +
  geom_tile(colour = "white", linewidth = 1.0) +
  geom_text(aes(label = lab), size = 2.7, colour = "black") +
  scale_fill_gradientn(colours = div_cols, limits = c(0.45, 0.72),
                       breaks = c(0.50, 0.55, 0.60, 0.65, 0.70),
                       name = "C-index") +
  scale_y_discrete(labels = c("GSE96058","METABRIC","GSE20685","TCGA\n(training)")) +
  labs(x = NULL, y = NULL,
       title = "Discrimination under four scaling schemes",
       subtitle = "S0 raw scale \u00b7 S1 z-score \u00b7 S2 percentile\nS3 centred, training SD \u00b7 chance = 0.500") +
  theme_pub(8) +
  theme(axis.text.y = element_text(size = 6.8, lineheight = 0.9),
        axis.text.x = element_text(size = 7.2),
        legend.position = "top", legend.direction = "horizontal",
        legend.key.width = unit(11, "mm"),
        legend.key.height = unit(2.4, "mm"),
        legend.margin = margin(0, 0, -3, 0),
        legend.title = element_text(size = 6.5))

# ---------------- E ----------------
e <- s10 %>%
  mutate(cohort = factor(cohort, levels = rev(c("GSE20685","METABRIC","GSE96058")))) %>%
  select(cohort, C_clinical, C_signature, C_combined) %>%
  pivot_longer(-cohort, names_to = "model", values_to = "C") %>%
  mutate(model = factor(model,
                        levels = c("C_clinical","C_combined","C_signature"),
                        labels = c("Clinical","Clinical + score","Score only")))

lab_e <- s10 %>%
  mutate(cohort = factor(cohort, levels = rev(c("GSE20685","METABRIC","GSE96058"))),
         d_txt = ifelse(abs(delta_C) < 5e-4, "0.000", sprintf("%+.3f", delta_C)),
         txt = sprintf("\u0394C = %s\nLRT p = %.2f", d_txt, LRT_p))

pE <- ggplot(e, aes(C, cohort, colour = model, shape = model)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_linerange(aes(xmin = 0.5, xmax = C), position = position_dodge(width = 0.32),
                 linewidth = 0.45) +
  geom_point(position = position_dodge(width = 0.32), size = 2.0) +
  geom_text(data = lab_e, aes(x = 0.85, y = cohort, label = txt),
            inherit.aes = FALSE, size = 2.35, lineheight = 0.95, hjust = 0) +
  scale_colour_manual(values = c("Clinical" = "grey55",
                                 "Clinical + score" = "#1A1A1A",
                                 "Score only" = "#D55E00"), name = NULL) +
  scale_shape_manual(values = c("Clinical" = 1, "Clinical + score" = 16,
                                "Score only" = 17), name = NULL) +
  scale_x_continuous(limits = c(0.475, 1.20), breaks = seq(0.5, 0.8, 0.1),
                     expand = expansion(mult = c(0.005, 0))) +
  labs(x = NULL, y = NULL,
       title = "Increment over age and stage",
       subtitle = "Clinical = age and stage; dashed line = chance") +
  theme_pub(8) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box.margin = margin(t = -2, b = 0),
        legend.text = element_text(size = 6.2),
        legend.key.size = unit(2.6, "mm"),
        legend.spacing.x = unit(1.2, "mm"))

# ---------------- F ----------------
files <- c(GSE20685 = "Table_S6_GEO_Validation_Result.csv",
           METABRIC = "Table_S6_METABRIC_Validation_Result.csv",
           GSE96058 = "Table_S6_GSE96058_Validation_Result.csv")

# METABRIC 验证表不保存 CombinedScore；按 s5_metabric.R 的同一模型
# （RiskScore + Age + LymphNode）在队列内重拟合。该重拟合可精确复现
# Table S10b 的 C-index 0.657 与 AUC 0.601 / 0.668 / 0.668。
comb_marker <- function(d) {
  if ("CombinedScore" %in% names(d)) return(d$CombinedScore)
  fit <- rms::cph(Surv(OS.time, OS) ~ RiskScore + Age + LymphNode,
                  data = d, x = TRUE, y = TRUE)
  as.numeric(predict(fit, type = "lp"))
}

roc_all <- list()
for (nm in names(files)) {
  d <- read.csv(files[[nm]], check.names = FALSE)
  d <- d[!is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0, ]
  req <- c("OS", "OS.time", "RiskScore")
  if (nm == "METABRIC") req <- c(req, "Age", "LymphNode")
  d <- d[complete.cases(d[, req]), ]
  d$CombinedScore <- comb_marker(d)
  d <- d[!is.na(d$CombinedScore), ]
  # 两条曲线用同一批病人、同一批事件，差别只在预测子本身：
  # 实线是签名，虚线是队列内用签名 + 年龄 + 分期重拟合的线性预测子。
  for (mk in c("Score only", "Score + age and stage")) {
    m <- if (mk == "Score only") d$RiskScore else d$CombinedScore
    # 只画 3 年和 5 年。三个队列的 1 年事件数分别是 5、14、14，AUC 在
    # 0.46–0.65 之间摆动完全是抽样噪音，画出来既不能支持也不能反驳任何结论，
    # 反而把 3 年那组真正有力的对照挤到角落。1 年的数值完整保留在 Table S10b。
    r <- timeROC(T = d$OS.time, delta = d$OS, marker = m, cause = 1,
                 weighting = "marginal", times = c(1095, 1825), iid = FALSE)
    for (i in 1:2) {
      roc_all[[length(roc_all) + 1]] <- data.frame(
        cohort = nm, marker = mk, horizon = c("3 years","5 years")[i],
        FPR = as.numeric(r$FP[, i]), TPR = as.numeric(r$TP[, i]),
        AUC = as.numeric(r$AUC)[i])
    }
  }
}
roc <- do.call(rbind, roc_all)
roc$cohort <- factor(roc$cohort, levels = c("GSE20685","METABRIC","GSE96058"))
roc$horizon <- factor(roc$horizon, levels = c("3 years","5 years"))
roc$marker <- factor(roc$marker, levels = c("Score only","Score + age and stage"))

auc_lab <- roc %>% distinct(cohort, horizon, marker, AUC) %>%
  pivot_wider(names_from = marker, values_from = AUC) %>%
  arrange(cohort, horizon) %>%
  group_by(cohort) %>%
  summarise(lab = paste(sprintf("%.2f \u2192 %.2f", `Score only`, `Score + age and stage`),
                        collapse = "   \u00b7   "), .groups = "drop")

# 每个分面上方直接给出三个时间点的 AUC 对照（签名 → 重拟合的临床联合模型）。
# 早先把这串数字放进面板右下角，会被处于对角线以下的曲线穿过，无法阅读。
strip_lab <- setNames(
  sprintf("%s\n%s", auc_lab$cohort, auc_lab$lab),
  auc_lab$cohort)

pF <- ggplot(roc, aes(FPR, TPR, colour = marker, linetype = horizon,
                      group = interaction(marker, horizon))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.3, colour = "grey60") +
  geom_step(linewidth = 0.45) +
  facet_wrap(~ cohort, nrow = 1, labeller = as_labeller(strip_lab)) +
  scale_colour_manual(values = c("Score only" = "#4D4D4D",
                                 "Score + age and stage" = "#D55E00"), name = NULL) +
  scale_linetype_manual(values = c("3 years" = "solid", "5 years" = "22"),
                        name = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(linewidth = 0.5)),
         linetype = guide_legend(order = 2, override.aes = list(colour = "grey20",
                                                                linewidth = 0.5))) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0","0.5","1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0","0.5","1")) +
  coord_equal() +
  labs(x = "1 - specificity", y = "Sensitivity",
       title = "Time-dependent ROC in the validation cohorts",
       subtitle = "Strip: AUC at 3 and 5 years, score \u2192 score + age and stage. One-year values are in Table S10b.") +
  theme_pub(8) +
  theme(strip.text = element_text(size = 5.5, lineheight = 1.0, hjust = 0.5),
        strip.background = element_rect(fill = "grey95", colour = NA),
        panel.spacing = unit(3.5, "mm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.box.just = "center",
        legend.key.width = unit(7, "mm"),
        legend.box.margin = margin(t = -1.5, b = 0),
        legend.spacing.x = unit(2, "mm"),
        legend.text = element_text(size = 6.0))

# ---------------- 拼合 ----------------
pD <- tag_panel(pD, "D"); pE <- tag_panel(pE, "E"); pF <- tag_panel(pF, "F")

fig6 <- (pD + pE + plot_layout(widths = c(1.02, 1))) / pF +
  plot_layout(heights = c(1, 1.12))
save_fig(fig6, "Figure4_DEF", W2, 5.5)

fig <- (pA | pB | pC) / (pD | pE) / pF +
  plot_layout(heights = c(1, 1.05, 1.20))
save_fig(fig, "Figure4_ABCDEF", W2, 190 * MM)

