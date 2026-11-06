## 08_class_score.R —— 方案/类别层面合成评分
setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))
pred <- read.csv("out/pred_matrix.csv", row.names = 1, check.names = FALSE)
clin <- readRDS("tcga_log2tpm.rds")$clin
res  <- read.csv("out/drug_results.csv", stringsAsFactors = FALSE)

cls <- list(
  `Anthracycline (Epirubicin, Mitoxantrone)` = c("Epirubicin","Mitoxantrone"),
  `Taxane (Paclitaxel, Docetaxel)`            = c("Paclitaxel","Docetaxel"),
  `Antimetabolite (5-FU, Gemcitabine)`        = c("5-Fluorouracil","Gemcitabine"),
  `Vinca (Vinorelbine, Vinblastine)`          = c("Vinorelbine","Vinblastine"),
  `Platinum (Cisplatin)`                      = c("Cisplatin"),
  `Endocrine (Tamoxifen, Fulvestrant)`        = c("Tamoxifen","Fulvestrant"),
  `Anti-HER2 TKI (Lapatinib, Afatinib, Sapitinib)` = c("Lapatinib","Afatinib","Sapitinib"),
  `PI3K/AKT/mTOR (n=8)`                       = c("Alpelisib","Buparlisib","Pictilisib","Taselisib",
                                                  "Ipatasertib","MK-2206","Rapamycin","AZD2014"),
  `PARP (Olaparib, Talazoparib, Niraparib)`   = c("Olaparib","Talazoparib","Niraparib"),
  `CDK4/6 (Palbociclib, Ribociclib)`          = c("Palbociclib","Ribociclib")
)

cohen_d <- function(a, b) {
  sp <- sqrt(((length(a)-1)*var(a) + (length(b)-1)*var(b))/(length(a)+length(b)-2))
  (mean(a)-mean(b))/sp
}

out <- do.call(rbind, lapply(names(cls), function(k) {
  ds <- intersect(cls[[k]], colnames(pred))
  s  <- rowMeans(pred[, ds, drop = FALSE])
  hi <- s[clin$RiskGroup == "High"]; lo <- s[clin$RiskGroup == "Low"]
  set.seed(42)
  bd <- replicate(1000, mean(sample(hi, length(hi), TRUE)) - mean(sample(lo, length(lo), TRUE)))
  w  <- wilcox.test(hi, lo)
  data.frame(regimen = k, n_drugs = length(ds),
             diff = round(mean(hi) - mean(lo), 3),
             ci_lo = round(quantile(bd, .025), 3), ci_hi = round(quantile(bd, .975), 3),
             p = w$p.value, d = round(cohen_d(hi, lo), 3), stringsAsFactors = FALSE)
}))
out$fdr <- p.adjust(out$p, "BH")
out <- out[order(out$diff), ]
print(out, row.names = FALSE)
write.csv(out, "out/regimen_scores.csv", row.names = FALSE)

library(ggplot2)
out$regimen <- factor(out$regimen, levels = out$regimen)
out$sig <- ifelse(out$fdr < 0.05, "FDR < 0.05", ifelse(out$p < 0.05, "nominal p < 0.05", "ns"))
p <- ggplot(out, aes(diff, regimen, color = sig)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = .25, linewidth = .5) +
  geom_point(size = 2.6) +
  scale_color_manual(values = c("FDR < 0.05" = "#c0392b", "nominal p < 0.05" = "#e59866", "ns" = "grey60")) +
  labs(x = "Delta predicted sensitivity score (High - Low risk)", y = NULL, color = NULL,
       title = "Regimen-level predicted sensitivity by risk group") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("fig/regimen_forest.pdf", p, width = 8.6, height = 4.6)
ggsave("fig/regimen_forest.png", p, width = 8.6, height = 4.6, dpi = 300)
cat("完成\n")
