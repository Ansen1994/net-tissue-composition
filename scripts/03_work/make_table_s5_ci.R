# ==============================================================================
# 重出 Table S5：多因素 Cox 系数 + 标准误 + 95% 置信区间
#
# 原来的 Table_S5_Cox_Coefficients.csv 只有 Gene / Coef / HR / P_value，
# 而 Figure 3A 印的是 HR (95% CI)，图表对不上。这里用与 s_figure3.R 完全相同的
# 输入与公式重拟合一次模型，把 summary(coxph)$conf.int 的上下限写进表里。
#
# 不手算 CI：coxph 的 conf.int 来自 Wald 区间 exp(coef ± 1.96*se)，
# 自己按 p 反推 SE 会在 p 极小时放大误差，直接用模型对象最稳。
#
# 跑完后与备份逐列比对 Coef / HR / P_value，确认是同一个模型而非重新拟合出偏差。
# ==============================================================================
suppressMessages(library(survival))

rt  <- read.csv("Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv", row.names = 1)
hub <- c("PLG", "CD79A", "MMP7", "ANGPT1")

md  <- rt[, c("OS.time", "OS", hub)]
fit <- coxph(as.formula(paste("Surv(OS.time, OS) ~", paste(hub, collapse = " + "))),
             data = md)
s   <- summary(fit)

out <- data.frame(
  Gene    = rownames(s$coefficients),
  Coef    = s$coefficients[, "coef"],
  SE      = s$coefficients[, "se(coef)"],
  HR      = s$conf.int[, "exp(coef)"],
  CI_low  = s$conf.int[, "lower .95"],
  CI_high = s$conf.int[, "upper .95"],
  z       = s$coefficients[, "z"],
  P_value = s$coefficients[, "Pr(>|z|)"],
  row.names = NULL,
  stringsAsFactors = FALSE
)
out <- out[match(hub, out$Gene), ]

write.csv(out, "Table_S5_Cox_Coefficients.csv", row.names = FALSE)

cat(sprintf("n = %d, events = %d\n", fit$n, fit$nevent))
print(out[, c("Gene", "HR", "CI_low", "CI_high", "P_value")], digits = 6)

# ---- 与旧表比对：Coef / HR / P_value 必须一致，否则说明模型被改动过 ----
old <- read.csv("/tmp/S5_backup.csv", stringsAsFactors = FALSE)
chk <- merge(old, out, by = "Gene", suffixes = c(".old", ".new"))
for (v in c("Coef", "HR", "P_value")) {
  d <- max(abs(chk[[paste0(v, ".old")]] - chk[[paste0(v, ".new")]]))
  cat(sprintf("max|%s_old - %s_new| = %.3e\n", v, v, d))
  if (d > 1e-10) stop("模型结果与旧表不一致，停止")
}
cat("校验通过：同一模型，仅新增 SE / CI / z 三列\n")
