# ==============================================================================
# Figure 5  A–D
#   A 22 基因参照 NET 评分在四队列的 C-index
#   B 参照评分每 SD 的 HR（森林）
#   C 参照评分 vs top-50/100/200 共表达替代
#   D 评分 × 治疗交互（森林）
# 输入：j1_results.csv、j2_results.csv
# ==============================================================================
source("figs_v2/style.R")

j1 <- read.csv("j1_results.csv")
j2 <- read.csv("j2_results.csv")

j1$cohort <- factor(j1$cohort, levels = c("TCGA","GSE20685","METABRIC","GSE96058"))
ref <- subset(j1, variable == "NETs_score")
ref$tag <- factor(ifelse(ref$cohort == "TCGA", "TCGA\n(in-sample)", as.character(ref$cohort)),
                  levels = c("GSE20685","METABRIC","GSE96058","TCGA\n(in-sample)"))

# ---------------- A ----------------
pA <- ggplot(ref, aes(tag, C_index)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_segment(aes(xend = tag, y = 0.5, yend = C_index, colour = cohort), linewidth = 0.6) +
  geom_point(aes(colour = cohort), size = 2.4) +
  geom_text(aes(label = sprintf("%.3f", C_index)), vjust = -1.1, size = 2.5) +
  scale_colour_manual(values = COH, guide = "none") +
  scale_y_continuous(limits = c(0.42, 0.61), breaks = seq(0.45, 0.60, 0.05)) +
  labs(x = NULL, y = "C-index (reference NET score)",
       title = "Reference NET score is non-prognostic",
       subtitle = "22 neutrophil / NETosis genes; dashed line = chance") +
  theme_pub(8) +
  theme(axis.text.x = element_text(size = 6.8, lineheight = 0.9))

# ---------------- B ----------------
hr <- ref
hr$tag <- factor(hr$tag, levels = rev(levels(ref$tag)))
pB <- ggplot(hr, aes(HR_perSD, tag, colour = cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_errorbarh(aes(xmin = HR_lo, xmax = HR_hi), height = 0.22, linewidth = 0.5) +
  geom_point(size = 2.2) +
  geom_text(aes(x = 1.42, label = sprintf("%.2f (p=%.2f)", HR_perSD, HR_p)),
            size = 2.2, hjust = 0) +
  scale_colour_manual(values = COH, guide = "none") +
  scale_x_continuous(limits = c(0.72, 2.12), breaks = c(0.8, 1.0, 1.2),
                     trans = "log", expand = expansion(mult = c(0.02, 0))) +
  labs(x = "HR per SD (95% CI)", y = NULL,
       title = "Effect size per SD",
       subtitle = "HR per SD with 95% CI; p from Wald test") +
  theme_pub(8) +
  theme(axis.text.y = element_text(size = 6.6, lineheight = 0.9))

# ---------------- C ----------------
sets <- subset(j1, variable %in% c("NETs_score","top50_coexp","top100_coexp","top200_coexp"))
sets$variable <- factor(sets$variable,
                        levels = c("NETs_score","top50_coexp","top100_coexp","top200_coexp"),
                        labels = c("Reference","top 50","top 100","top 200"))
set_cols <- c(Reference = "#1A1A1A", "top 50" = "#0072B2",
              "top 100" = "#009E73", "top 200" = "#D55E00")

pC <- ggplot(sets, aes(cohort, C_index, colour = variable, shape = variable)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_point(position = position_dodge(width = 0.62), size = 2.1) +
  scale_colour_manual(values = set_cols, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
  scale_y_continuous(limits = c(0.40, 0.60), breaks = seq(0.40, 0.60, 0.05)) +
  labs(x = NULL, y = "C-index",
       title = "Co-expression surrogates do not recover signal",
       subtitle = "genes selected in TCGA by correlation with the reference score") +
  theme_pub(8) +
  theme(axis.text.x = element_text(size = 6.8, angle = 20, hjust = 1),
        legend.position = "top", legend.direction = "horizontal",
        legend.key.size = unit(2.8, "mm"),
        legend.margin = margin(0, 0, -3, 0))

# ---------------- D ----------------
j2$cohort_tx <- factor(c("METABRIC\nchemo","METABRIC\nendocrine",
                         "GSE96058\nchemo","GSE96058\nendocrine"),
                       levels = rev(c("METABRIC\nchemo","METABRIC\nendocrine",
                                      "GSE96058\nchemo","GSE96058\nendocrine")))
d_est <- rbind(
  data.frame(block = j2$cohort_tx, est = "Interaction",
             HR = j2$HR_interaction, lo = j2$lo, hi = j2$hi, p = j2$p_interaction),
  data.frame(block = j2$cohort_tx, est = "Treated",
             HR = j2$HR_in_treated, lo = NA, hi = NA, p = j2$p_in_treated),
  data.frame(block = j2$cohort_tx, est = "Untreated",
             HR = j2$HR_in_untreated, lo = NA, hi = NA, p = j2$p_in_untreated))
d_est$est <- factor(d_est$est, levels = c("Interaction","Treated","Untreated"))

lab_d <- d_est %>% filter(est == "Interaction") %>%
  mutate(txt = sprintf("%.2f (%.2f\u2013%.2f)\np %s", HR, lo, hi,
                       ifelse(p < 0.0001, "< 0.0001", paste("=", sprintf("%.3f", p)))))

pD <- ggplot(d_est, aes(HR, block, colour = est)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.24, linewidth = 0.5,
                position = position_dodge(width = 0.55), na.rm = TRUE) +
  geom_point(size = 2.0, position = position_dodge(width = 0.55)) +
  geom_text(data = lab_d, aes(x = 2.55, y = block, label = txt),
            inherit.aes = FALSE, size = 2.05, lineheight = 0.95, hjust = 0) +
  scale_colour_manual(values = c(Interaction = "#B2182B", Treated = "#0072B2",
                                 Untreated = "grey45"), name = NULL) +
  scale_x_continuous(limits = c(0.58, 4.9), trans = "log",
                     breaks = c(0.7, 1.0, 1.5, 2.0),
                     expand = expansion(mult = c(0.02, 0))) +
  labs(x = "Hazard ratio (log scale)", y = NULL,
       title = "Score by treatment: no consistent predictive signal",
       subtitle = "interaction HR (95% CI) annotated; circles, within-stratum HR") +
  theme_pub(8) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.key.size = unit(2.8, "mm"),
        legend.margin = margin(0, 0, -3, 0),
        axis.text.y = element_text(size = 6.6, lineheight = 0.9))

# ---------------- 拼合 ----------------
pA <- tag_panel(pA, "A"); pB <- tag_panel(pB, "B")
pC <- tag_panel(pC, "C"); pD <- tag_panel(pD, "D")

fig <- (pA + pB + plot_layout(widths = c(1.25, 1))) /
  (pC + pD + plot_layout(widths = c(1, 1.28))) +
  plot_layout(heights = c(1, 1.05))

save_fig(fig, "Figure5_ABCD", W2, 6.2)
