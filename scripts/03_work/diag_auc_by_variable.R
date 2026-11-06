suppressMessages({library(survival); library(timeROC)})
chk <- function(f, nm) {
  d <- read.csv(f, check.names = FALSE)
  d <- d[!is.na(d$OS)&!is.na(d$OS.time)&d$OS.time>0, ]
  cat("\n=====", nm, "n =", nrow(d), "=====\n")
  for (v in c("RiskScore","CombinedScore")) {
    if (!v %in% names(d)) { cat("  缺", v, "\n"); next }
    ci <- as.numeric(concordance(Surv(d$OS.time,d$OS) ~ d[[v]], reverse=TRUE)$concordance)
    r <- tryCatch(timeROC(T=d$OS.time, delta=d$OS, marker=d[[v]], cause=1,
                          weighting="marginal", times=c(365,1095,1825), iid=FALSE), error=function(e) NULL)
    cat(sprintf("  %-14s C-index %.3f | AUC 1/3/5y = %s\n", v, ci,
        if(is.null(r)) "ERR" else paste(round(as.numeric(r$AUC),3), collapse=" / ")))
  }
}
chk("Table_S6_GEO_Validation_Result.csv","GSE20685")
chk("Table_S6_METABRIC_Validation_Result.csv","METABRIC")
chk("Table_S6_GSE96058_Validation_Result.csv","GSE96058")
