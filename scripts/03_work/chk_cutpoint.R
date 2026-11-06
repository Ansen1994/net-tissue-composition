## 分组方式敏感性：median vs surv_cutpoint
## 目的：核实方法 2.4 描述与代码是否一致，并检验结论对分组方式是否敏感
suppressPackageStartupMessages({
  library(survival); library(survminer)
})

files <- c(
  Training_TCGA      = "Table_S6_Patient_Risk_Score_Result.csv",
  GSE20685           = "Table_S6_GEO_Validation_Result.csv",
  GSE96058           = "Table_S6_GSE96058_Validation_Result.csv",
  METABRIC           = "Table_S6_METABRIC_Validation_Result.csv"
)

out <- data.frame()
for (nm in names(files)) {
  d <- read.csv(files[[nm]], stringsAsFactors = FALSE)
  d$OS <- as.numeric(d$OS); d$OS.time <- as.numeric(d$OS.time)
  d <- d[complete.cases(d[, c("OS", "OS.time", "RiskScore")]) & d$OS.time > 0, ]
  if (nrow(d) < 30) next

  # 方案 A：中位数
  gA <- ifelse(d$RiskScore > median(d$RiskScore), "High", "Low")
  # 方案 B：maximally selected rank statistic
  cp <- tryCatch(surv_cutpoint(d, time = "OS.time", event = "OS",
                               variables = "RiskScore", minprop = 0.1),
                 error = function(e) NULL)
  if (is.null(cp)) { cutB <- median(d$RiskScore) } else { cutB <- cp$cutpoint$cutpoint[1] }
  gB <- ifelse(d$RiskScore > cutB, "High", "Low")

  for (k in c("A", "B")) {
    g <- if (k == "A") gA else gB
    cutv <- if (k == "A") median(d$RiskScore) else cutB
    d$g <- factor(g, levels = c("Low", "High"))
    fit <- survdiff(Surv(OS.time, OS) ~ g, data = d)
    p <- 1 - pchisq(fit$chisq, df = length(fit$n) - 1)
    cx <- coxph(Surv(OS.time, OS) ~ g, data = d)
    s <- summary(cx)
    out <- rbind(out, data.frame(
      Cohort = nm, Method = ifelse(k == "A", "median", "surv_cutpoint"),
      N = nrow(d), Events = sum(d$OS),
      Cut = round(cutv, 4),
      nHigh = sum(g == "High"), nLow = sum(g == "Low"),
      HR = round(s$conf.int[1], 3),
      Lower = round(s$conf.int[3], 3), Upper = round(s$conf.int[4], 3),
      P_logrank = signif(p, 3)
    ))
  }
}
print(out, row.names = FALSE)
write.csv(out, "Table_S16_Cutpoint_Sensitivity.csv", row.names = FALSE)
