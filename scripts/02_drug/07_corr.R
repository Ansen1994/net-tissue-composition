## 07_corr.R —— 药物间预测相关性 + 高危组内部异质性
setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))
pred <- read.csv("out/pred_matrix.csv", row.names = 1, check.names = FALSE)
clin <- readRDS("tcga_log2tpm.rds")$clin
res  <- read.csv("out/drug_results.csv", stringsAsFactors = FALSE)

cat("预测矩阵:", dim(pred), "\n")

## 1. 药物两两相关（Spearman）
cm <- cor(pred, method = "spearman")
write.csv(round(cm, 3), "out/drug_correlation.csv")

## 2. 全局共同轴：每位患者跨药物平均 z
glob <- rowMeans(pred)
cat("\n全局平均耐药分数 vs RiskScore: rho =",
    round(cor(glob, clin$RiskScore, method = "spearman"), 3),
    " p =", signif(cor.test(glob, clin$RiskScore, method="spearman")$p.value, 3), "\n")
cat("全局平均 vs 分组: ", round(tapply(glob, clin$RiskGroup, mean), 3), "\n")

## 3. 三个显著药物的相互相关
sig3 <- c("Sapitinib","Docetaxel","Fulvestrant")
cat("\n三个显著药物间 Spearman:\n"); print(round(cm[sig3, sig3], 3))

## 4. 去全局轴后的偏相关：残差化再检验
resid_mat <- apply(pred, 2, function(x) residuals(lm(x ~ glob)))
cm_res <- cor(resid_mat, method = "spearman")
cat("\n去全局轴后 Docetaxel vs Sapitinib:",
    round(cm_res["Docetaxel","Sapitinib"], 3), "\n")

## 5. 高危组内：多西他赛预测耐药的四分位分布
hi <- clin$RiskGroup == "High"
q <- quantile(pred$Docetaxel, c(.75, .9))
cat("\n多西他赛预测值 P75 =", round(q[1],2), " P90 =", round(q[2],2), "\n")
cat("高危组中处于全体前 25% 耐药者:",
    sum(pred$Docetaxel[hi] > q[1]), "/", sum(hi), "\n")
cat("低危组中处于全体前 25% 耐药者:",
    sum(pred$Docetaxel[!hi] > q[1]), "/", sum(!hi), "\n")

## 6. 与分期的关系
cat("\n多西他赛预测值按分期:\n"); print(round(tapply(pred$Docetaxel, clin$Stage, mean), 3))
cat("Sapitinib 按分期:\n"); print(round(tapply(pred$Sapitinib, clin$Stage, mean), 3))
cat("Fulvestrant 按分期:\n"); print(round(tapply(pred$Fulvestrant, clin$Stage, mean), 3))

## Kruskal-Wallis
for (d in sig3) cat(d, "分期 KW p =", signif(kruskal.test(pred[[d]] ~ clin$Stage)$p.value, 3), "\n")
