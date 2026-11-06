## 17_targeted_figs.R
## 可视化 16_targeted_residual.R 的结果（图内一律英文，规避本机 R 中文编码问题）
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(ggrepel) })
setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))
dir.create("fig", showWarnings = FALSE)

res <- read.csv("out/targeted_residual.csv", stringsAsFactors = FALSE)
ver <- read.csv("out/targeted_verdict.csv", stringsAsFactors = FALSE)
modes <- c("counts", "log2tpm", "zscore", "qn", "rank")
res$mode <- factor(res$mode, levels = modes)
## 按类别与稳健性排序
ord <- ver$drug[order(ver$category, ver$verdict, -abs(ver$median_d))]
res$drug <- factor(res$drug, levels = rev(ord))

lab <- function(x) ifelse(is.na(x), "", ifelse(x < 0.05, "*", ""))
res$star <- lab(res$p_res11)

## ---- 图 A：扣除共同轴（11 药留一法）后的效应量热图 ----
p1 <- ggplot(res, aes(x = mode, y = drug, fill = d_res11)) +
  geom_tile(color = "grey85", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f%s", d_res11, star)), size = 3.1, color = "grey15") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-0.65, 0.65),
                       name = "Cohen's d\n(LOO residual)", oob = scales::squish) +
  scale_x_discrete(position = "top") +
  labs(title = "Drug-specific signal after removing the common axis",
       subtitle = "11 target drugs x 5 input preprocessings; LOO axis = mean of the other 10 drugs",
       x = "Input preprocessing", y = NULL,
       caption = "* p < 0.05 (Wilcoxon on LOO residual). Positive = high-risk relatively resistant.") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 9.5),
        panel.grid = element_blank(),
        plot.caption = element_text(size = 8.5, color = "grey35"),
        legend.position = "right")
ggsave("fig/FigS_targeted_residual.png", p1, width = 9.2, height = 6.6, dpi = 300)
ggsave("fig/FigS_targeted_residual.pdf", p1, width = 9.2, height = 6.6)

## ---- 图 B：11 药轴 vs 29 药轴（log2tpm / qn），看轴规模是否改变判定 ----
r29 <- res %>% filter(!is.na(p_res29)) %>%
  select(mode, drug, d_res11, d_res29, p_res11, p_res29)
axis_cmp <- cor(r29$d_res11, r29$d_res29)
p2 <- ggplot(r29, aes(x = d_res11, y = d_res29)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_abline(slope = 1, intercept = 0, color = "grey75") +
  geom_point(aes(color = mode, shape = p_res29 < 0.05), size = 2.6, alpha = .85) +
  geom_text_repel(aes(label = drug), size = 2.5, color = "grey25",
                  max.overlaps = 20, box.padding = 0.35, min.segment.length = 0.1) +
  coord_equal(xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7)) +
  scale_color_manual(values = c(log2tpm = "#1f77b4", qn = "#d62728")) +
  labs(title = "Common-axis size does not change the verdict",
       subtitle = sprintf("LOO-11 axis vs LOO-29 axis; Pearson r = %.3f (n = %d)",
                          axis_cmp, nrow(r29)),
       x = "Cohen's d  (residual vs 11-drug axis)",
       y = "Cohen's d  (residual vs 29-drug axis)",
       color = "Preprocessing", shape = "p<0.05 (29-axis)") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave("fig/FigS_axis_comparison.png", p2, width = 7.6, height = 6.6, dpi = 300)
ggsave("fig/FigS_axis_comparison.pdf", p2, width = 7.6, height = 6.6)

cat("Figures written. LOO-11 vs LOO-29 Pearson r =", round(axis_cmp, 3), "\n")
options(width = 200); print(ver, row.names = FALSE)
