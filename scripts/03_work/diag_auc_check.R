# 用图对应的结果表核对：C-index 与 timeROC AUC 是否自洽
suppressMessages({library(survival); library(timeROC)})

chk <- function(f, nm) {
  d <- read.csv(f, check.names = FALSE)
  d <- d[!is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$RiskScore), ]
  cat("\n=====", nm, " n =", nrow(d), " 事件 =", sum(d$OS), "=====\n")
  ci <- as.numeric(concordance(Surv(d$OS.time, d$OS) ~ d$RiskScore, reverse = TRUE)$concordance)
  cat("  C-index =", round(ci, 3), "\n")
  for (tt in c(365, 365*3, 365*5)) {
    ev <- sum(d$OS == 1 & d$OS.time <= tt)
    at <- sum(d$OS.time > tt)
    cat(sprintf("  %4.0f 天：该时点前事件 %4d 例，存活至该时点 %4d 例", tt, ev, at))
    r <- tryCatch(timeROC(T = d$OS.time, delta = d$OS, marker = d$RiskScore,
                          cause = 1, times = tt, iid = FALSE), error = function(e) NULL)
    if (!is.null(r)) cat("  | timeROC AUC =", round(as.numeric(r$AUC[2]), 3))
    cat("\n")
  }
}

chk("Table_S6_GEO_Validation_Result.csv", "GSE20685")
chk("Table_S6_METABRIC_Validation_Result.csv", "METABRIC")
chk("Table_S6_GSE96058_Validation_Result.csv", "GSE96058")
