# ==============================================================================
# Figure 3F  训练集判别力：样本内 vs 五折交叉验证
# 输入：diag_training_optimism.R 的输出（0.690 / 0.638）
# ==============================================================================
source("figs_v2/style.R")

d <- data.frame(x = c(1, 2), C = c(0.690, 0.638))

p <- ggplot(d, aes(x, C)) +
  geom_segment(aes(x = 1, xend = 1.86, y = 0.690, yend = 0.643),
               arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
               linewidth = 0.4, colour = "grey35") +
  geom_point(size = 3.0, colour = c("#0072B2", "#D55E00")) +
  geom_text(aes(label = sprintf("%.3f", C)), vjust = -1.3, size = 2.6) +
  annotate("text", x = 1.52, y = 0.713, label = "\u0394 = \u22120.052",
           size = 2.4, colour = "grey25") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("Apparent\n(in-sample)", "5-fold\ncross-validated"),
                     limits = c(0.55, 2.45)) +
  scale_y_continuous(limits = c(0.605, 0.725), breaks = seq(0.62, 0.70, 0.04)) +
  labs(x = NULL, y = "C-index",
       title = "Optimism in the training cohort") +
  theme_pub(8) +
  theme(axis.text.x = element_text(size = 6.6, lineheight = 0.9))

save_fig(p, "Figure3F_Optimism", 62 * MM, 60 * MM)
