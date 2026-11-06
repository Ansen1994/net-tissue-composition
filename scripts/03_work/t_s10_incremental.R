# ==============================================================================
# Table S10：增量价值（临床 / 签名 / 临床+签名）与按变量的时间依赖 AUC
# ==============================================================================
suppressMessages({library(survival); library(timeROC)})
ci <- function(lp, d) as.numeric(concordance(Surv(d$OS.time, d$OS) ~ lp, reverse = TRUE)$concordance)

files <- c(GSE20685 = "Table_S6_GEO_Validation_Result.csv",
           METABRIC = "Table_S6_METABRIC_Validation_Result.csv",
           GSE96058 = "Table_S6_GSE96058_Validation_Result.csv")
clin <- list(GSE20685 = list(f = "Age + T_stage + N_stage", v = c("T_stage","N_stage")),
             METABRIC = list(f = "Age + LymphNode",          v = c("LymphNode")),
             GSE96058 = list(f = "Age + TumorSize + LymphNode", v = c("LymphNode")))

inc <- list(); auc <- list()
for (nm in names(files)) {
  d <- read.csv(files[[nm]], check.names = FALSE)
  d <- d[!is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0, ]
  cv <- clin[[nm]]$v
  for (v in cv) if (is.character(d[[v]])) d[[v]] <- factor(d[[v]])
  d <- d[complete.cases(d[, c("OS","OS.time","RiskScore", cv)]), ]
  m1 <- coxph(as.formula(paste("Surv(OS.time,OS) ~", clin[[nm]]$f)), data = d)
  m2 <- coxph(Surv(OS.time, OS) ~ RiskScore, data = d)
  m3 <- coxph(as.formula(paste("Surv(OS.time,OS) ~ RiskScore +", clin[[nm]]$f)), data = d)
  c1 <- ci(predict(m1, type = "lp"), d); c2 <- ci(predict(m2, type = "lp"), d)
  c3 <- ci(predict(m3, type = "lp"), d)
  lrt <- anova(m1, m3, test = "LRT")
  inc[[nm]] <- data.frame(cohort = nm, n = nrow(d), events = sum(d$OS),
                          C_clinical = c1, C_signature = c2, C_combined = c3,
                          delta_C = c3 - c1, LRT_p = lrt$`Pr(>|Chi|)`[2])
  for (v in intersect(c("RiskScore","CombinedScore"), names(d))) {
    r <- tryCatch(timeROC(T = d$OS.time, delta = d$OS, marker = d[[v]], cause = 1,
                          weighting = "marginal", times = c(365,1095,1825), iid = FALSE),
                  error = function(e) NULL)
    auc[[length(auc)+1]] <- data.frame(cohort = nm, variable = v,
      C_index = as.numeric(concordance(Surv(d$OS.time, d$OS) ~ d[[v]], reverse = TRUE)$concordance),
      AUC_1y = if (is.null(r)) NA else as.numeric(r$AUC)[1],
      AUC_3y = if (is.null(r)) NA else as.numeric(r$AUC)[2],
      AUC_5y = if (is.null(r)) NA else as.numeric(r$AUC)[3])
  }
}
inc <- do.call(rbind, inc); auc <- do.call(rbind, auc)
inc[, 4:7] <- lapply(inc[, 4:7], function(x) round(x, 3)); inc$LRT_p <- signif(inc$LRT_p, 3)
auc[, 3:6] <- lapply(auc[, 3:6], function(x) round(x, 3))
write.csv(inc, "Table_S10a_Incremental_Value.csv", row.names = FALSE)
write.csv(auc, "Table_S10b_AUC_by_Variable.csv", row.names = FALSE)
cat("\n===== 增量价值 =====\n"); print(inc, row.names = FALSE)
cat("\n===== 按变量 AUC =====\n"); print(auc, row.names = FALSE)
