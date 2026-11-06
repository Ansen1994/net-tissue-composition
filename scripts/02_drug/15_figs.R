## 15_figs.R
## 重画敏感性热图（英文标签，避免本机 LC_CTYPE=C 的中文乱码）
suppressPackageStartupMessages(library(ggplot2))
setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))

res <- read.csv("out/sensitivity_matrix.csv", stringsAsFactors = FALSE)
ord <- read.csv("out/sensitivity_verdict.csv", stringsAsFactors = FALSE)
res$mode <- factor(res$mode, levels = c("counts", "log2tpm", "zscore", "qn", "rank"))
res$drug <- factor(res$drug, levels = rev(ord$drug))
res$star <- ifelse(res$p < 0.05, "*", "")

p1 <- ggplot(res, aes(mode, drug, fill = cohen_d)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = star), size = 5.5, vjust = 0.78, color = "grey15") +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c", midpoint = 0,
                       name = "Cohen's d") +
  labs(x = "Input preprocessing", y = NULL,
       title = "Effect size of predicted drug response vs. input preprocessing",
       subtitle = "* p < 0.05 (Wilcoxon). Blue = high-risk sensitive, red = high-risk resistant") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, color = "grey30"))
ggsave("fig/FigS_sensitivity_heatmap.pdf", p1, width = 7.6, height = 5.4)
ggsave("fig/FigS_sensitivity_heatmap.png", p1, width = 7.6, height = 5.4, dpi = 300)

## 扣除共同轴后的效应量（留一法）
r2 <- read.csv("out/residual_axis.csv", stringsAsFactors = FALSE)
r2$mode <- factor(r2$mode, levels = c("counts", "log2tpm", "zscore", "qn", "rank"))
r2$drug <- factor(r2$drug, levels = rev(c("5-Fluorouracil", "Docetaxel", "Fulvestrant", "Epirubicin")))
r2$star <- ifelse(r2$p_resid < 0.05, "*", "")

p2 <- ggplot(r2, aes(mode, drug, fill = d_resid)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = star), size = 6, vjust = 0.78, color = "grey15") +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c", midpoint = 0,
                       name = "Cohen's d\n(residual)") +
  labs(x = "Input preprocessing", y = NULL,
       title = "Effect size after removing the shared cross-drug axis (leave-one-out)",
       subtitle = "* p < 0.05. Signal surviving in all five pipelines indicates drug-specific effect") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"))
ggsave("fig/FigS_residual_axis.pdf", p2, width = 7.2, height = 3.6)
ggsave("fig/FigS_residual_axis.png", p2, width = 7.2, height = 3.6, dpi = 300)
cat("done\n")
