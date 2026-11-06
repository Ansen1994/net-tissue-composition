## 04_figs.R —— 药敏结果可视化
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))
res <- read.csv("out/drug_results.csv", stringsAsFactors = FALSE)
pred <- read.csv("out/pred_matrix.csv", row.names = 1, check.names = FALSE)
clin <- readRDS("tcga_log2tpm.rds")$clin

cat_order <- c("Endocrine","Chemo","HER2","PI3K-AKT-mTOR","PARP","CDK4-6","CDK4-6-control")
cat_lab <- c("Endocrine therapy","Chemotherapy","Anti-HER2","PI3K/AKT/mTOR","PARP inhibitor",
             "CDK4/6 inhibitor","CDK4/6 (pan-CDK ref.)")
names(cat_lab) <- cat_order
res$category <- factor(res$category, levels = cat_order, labels = cat_lab)
res <- res[order(res$category, res$diff), ]
res$drug <- factor(res$drug, levels = unique(res$drug))

## ---- 1. 森林图：High vs Low 预测 ln(IC50) 差值 ----
res$sig <- ifelse(res$fdr < 0.05, "FDR < 0.05",
           ifelse(res$p_wilcox < 0.05, "nominal p < 0.05", "ns"))
p1 <- ggplot(res, aes(x = diff, y = drug, color = sig)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.5) +
  geom_point(size = 2.2) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("FDR < 0.05" = "#c0392b",
                                "nominal p < 0.05" = "#e59866",
                                "ns" = "grey60")) +
  labs(x = expression(Delta * " predicted ln(IC50)  (High-risk - Low-risk)"),
       y = NULL, color = NULL,
       title = "Predicted drug sensitivity by risk group") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
        strip.background = element_rect(fill = "grey94"),
        legend.position = "top")
ggsave("fig/forest_all_drugs.pdf", p1, width = 8.5, height = 9)
ggsave("fig/forest_all_drugs.png", p1, width = 8.5, height = 9, dpi = 300)

## ---- 2. 各类别代表药物箱线图 ----
sel <- c("Fulvestrant","Tamoxifen","Epirubicin","Docetaxel","Paclitaxel",
         "Lapatinib","Alpelisib","Ipatasertib","Olaparib","Palbociclib")
sel <- sel[sel %in% colnames(pred)]
long <- pred[, sel, drop = FALSE]
long$Patient_ID <- rownames(pred)
long <- merge(long, clin[, c("Patient_ID","RiskGroup","RiskScore")], by = "Patient_ID")
lg <- pivot_longer(long, all_of(sel), names_to = "drug", values_to = "z")

pl <- lapply(sel, function(d) {
  sub <- lg[lg$drug == d, ]
  r <- res[res$drug == d, ]
  lab <- sprintf("p = %.2g\nFDR = %.2g", r$p_wilcox[1], r$fdr[1])
  ggplot(sub, aes(RiskGroup, z, fill = RiskGroup)) +
    geom_boxplot(outlier.size = 0.4, width = 0.6) +
    scale_fill_manual(values = c(Low = "#abdaec", High = "#f6c0ce")) +
    annotate("text", x = 1.5, y = max(sub$z) * 1.02, label = lab, size = 2.7) +
    labs(title = d, x = NULL, y = "Predicted sensitivity (z)") +
    theme_bw(base_size = 9) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))
})
p2 <- wrap_plots(pl, ncol = 5)
ggsave("fig/boxplots_key_drugs.pdf", p2, width = 15, height = 6)
ggsave("fig/boxplots_key_drugs.png", p2, width = 15, height = 6, dpi = 300)

## ---- 3. 类别层面汇总 ----
catsum <- res %>% group_by(category) %>%
  summarise(n = n(),
            n_sig = sum(fdr < 0.05),
            mean_d = round(mean(cohen_d), 2),
            med_p = signif(median(p_wilcox), 2), .groups = "drop")
write.csv(catsum, "out/category_summary.csv", row.names = FALSE)
print(catsum)

p3 <- ggplot(res, aes(x = cohen_d, y = category, color = fdr < 0.05)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_jitter(width = 0, height = 0.15, size = 2, alpha = 0.85) +
  scale_color_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "grey65"),
                     labels = c("FDR >= 0.05", "FDR < 0.05")) +
  labs(x = "Cohen's d (High vs Low risk)", y = NULL, color = NULL,
       title = "Effect size distribution by drug class") +
  theme_bw(base_size = 10) + theme(legend.position = "top")
ggsave("fig/category_effects.pdf", p3, width = 7, height = 5)

## ---- 4. 免疫微环境 ----
imm <- read.csv("out/immune_proxy.csv", stringsAsFactors = FALSE)
keep <- c("T-cell inflamed GEP (18-gene)","T cell infiltration","Cytolytic activity (GZMA/PRF1)",
          "IFN-gamma response","ImmuneScore","CD274","CD8A")
im <- imm[imm$variable %in% keep, ]
im$variable <- factor(im$variable, levels = im$variable[order(im$diff)])
p4 <- ggplot(im, aes(x = diff, y = variable, fill = diff < 0)) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = 0, color = "grey40") +
  scale_fill_manual(values = c(`TRUE` = "#abdaec", `FALSE` = "#f6c0ce"), guide = "none") +
  labs(x = "High-risk minus Low-risk (ssGSEA / score units)", y = NULL,
       title = "Immune microenvironment proxies") +
  theme_bw(base_size = 10)
ggsave("fig/immune_proxy.pdf", p4, width = 7.5, height = 4.2)

cat("图表输出完成\n")
