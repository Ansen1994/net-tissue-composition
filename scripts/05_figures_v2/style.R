# ==============================================================================
# 统一绘图风格（v2 稿件的全部新图共用）
# ==============================================================================
suppressMessages({
  library(ggplot2); library(patchwork); library(scales)
  library(survival); library(timeROC); library(dplyr); library(tidyr)
})

MM <- 1 / 25.4                      # mm -> inch
W2 <- 170 * MM                      # 双栏宽 170 mm
W1 <- 85 * MM                       # 单栏宽 85 mm

# 队列配色（Okabe-Ito，色盲友好）
COH <- c(TCGA = "#0072B2", GSE20685 = "#009E73",
         METABRIC = "#D55E00", GSE96058 = "#CC79A7")
COH_LAB <- c(TCGA = "TCGA (n = 762)", GSE20685 = "GSE20685 (n = 327)",
             METABRIC = "METABRIC (n = 1,445)", GSE96058 = "GSE96058 (n = 2,840)")

# 模块配色
MOD <- c(Stroma = "#8C6BB1", Stroma_ECM = "#8C6BB1", Prolif = "#D7301F",
         Proliferation = "#D7301F", Purity_epi = "#2166AC", Epithelial = "#2166AC")
IMM <- c(Bcell_noCD79A = "#1B7837", CD8_T = "#762A83", Treg = "#E08214",
         Cytotoxic = "#35978F", Macrophage = "#8C510A", Neutrophil = "#01665E")

# 发散色（相关系数）
DIV <- c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B")
div_cols <- colorRampPalette(DIV)(101)

theme_pub <- function(bs = 8) {
  theme_classic(base_size = bs, base_family = "sans") +
    theme(
      axis.line   = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.35, colour = "black"),
      axis.text   = element_text(colour = "black", size = bs - 0.5),
      axis.title  = element_text(colour = "black", size = bs),
      plot.title  = element_text(face = "bold", size = bs + 0.5, hjust = 0,
                                 margin = ggplot2::margin(b = 3)),
      plot.subtitle = element_text(size = bs - 1, colour = "grey25",
                                   margin = ggplot2::margin(b = 3)),
      strip.background = element_blank(),
      strip.text  = element_text(face = "bold", size = bs - 0.3),
      legend.title = element_text(size = bs - 0.5),
      legend.text  = element_text(size = bs - 1),
      legend.key.size = unit(3.0, "mm"),
      legend.background = element_blank(),
      panel.grid  = element_blank(),
      plot.tag    = element_text(face = "bold", size = bs + 3),
      plot.margin = ggplot2::margin(4, 6, 4, 4)
    )
}

# 面板字母（用 labs(tag=)，避免 patchwork 的 & 运算符兼容问题）
tag_panel <- function(p, lab) p + labs(tag = lab)

save_fig <- function(p, name, w, h) {
  dir.create("figures_v2", showWarnings = FALSE)
  ggsave(file.path("figures_v2", paste0(name, ".pdf")), p,
         width = w, height = h, units = "in", device = cairo_pdf)
  cat("saved:", name, sprintf("(%.1f x %.1f mm)\n", w * 25.4, h * 25.4))
}
