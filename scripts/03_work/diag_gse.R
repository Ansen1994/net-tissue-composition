suppressMessages({library(survival); library(survminer)})
d <- read.csv("Table_S6_GSE96058_Validation_Result.csv", stringsAsFactors = FALSE)
for (mp in c(0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5)) {
  r <- try(surv_cutpoint(d, time="OS.time", event="OS", variables="RiskScore", minprop=mp), silent=TRUE)
  if (inherits(r,"try-error")) { cat(sprintf("minprop=%.2f 失败\n", mp)); next }
  cat_cut <- surv_categorize(r)
  tb <- table(cat_cut$RiskScore)
  cat(sprintf("minprop=%.2f  low=%d high=%d  截点=%.4f\n", mp, tb["low"], tb["high"], r$cutpoint$cutpoint))
}
