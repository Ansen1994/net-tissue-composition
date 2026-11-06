suppressMessages({library(survival); library(timeROC)})
ci <- function(lp, d) as.numeric(concordance(Surv(d$OS.time,d$OS) ~ lp, reverse=TRUE)$concordance)
rep1 <- function(f, nm, clin_formula, clin_vars) {
  d <- read.csv(f, check.names=FALSE)
  d <- d[!is.na(d$OS)&!is.na(d$OS.time)&d$OS.time>0,]
  for (v in clin_vars) if (is.character(d[[v]])) d[[v]] <- factor(d[[v]])
  d <- d[complete.cases(d[, c("OS","OS.time","RiskScore",clin_vars)]),]
  m1 <- coxph(as.formula(paste("Surv(OS.time,OS) ~", clin_formula)), data=d)
  m2 <- coxph(Surv(OS.time,OS) ~ RiskScore, data=d)
  m3 <- coxph(as.formula(paste("Surv(OS.time,OS) ~ RiskScore +", clin_formula)), data=d)
  cat(sprintf("\n== %s (n=%d, 事件=%d) ==\n", nm, nrow(d), sum(d$OS)))
  cat(sprintf("  临床单独      C=%.3f\n", ci(predict(m1,type='lp'), d)))
  cat(sprintf("  签名单独      C=%.3f\n", ci(predict(m2,type='lp'), d)))
  cat(sprintf("  临床+签名     C=%.3f  (Δ vs 临床 = %+.3f)\n", ci(predict(m3,type='lp'), d),
              ci(predict(m3,type='lp'), d) - ci(predict(m1,type='lp'), d)))
  lrt <- anova(m1, m3, test="LRT")
  cat(sprintf("  似然比检验 签名增量 p = %.3g\n", lrt$`Pr(>|Chi|)`[2]))
}
rep1("Table_S6_GEO_Validation_Result.csv","GSE20685","Age + T_stage + N_stage", c("T_stage","N_stage"))
rep1("Table_S6_METABRIC_Validation_Result.csv","METABRIC","Age + LymphNode", c("LymphNode"))
rep1("Table_S6_GSE96058_Validation_Result.csv","GSE96058","Age + TumorSize + LymphNode", c("LymphNode"))
