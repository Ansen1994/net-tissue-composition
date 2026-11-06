suppressMessages({library(survival); library(timeROC)})
d <- read.csv("Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv", row.names=1, check.names=FALSE)
hub <- c("PLG","CD79A","MMP7","ANGPT1")
d <- d[complete.cases(d[,c("OS.time","OS",hub)]) & d$OS.time>0, ]
cat("TCGA n =", nrow(d), "事件 =", sum(d$OS), "\n")
m <- coxph(Surv(OS.time,OS) ~ PLG+CD79A+MMP7+ANGPT1, data=d)
d$rs <- predict(m, type="risk")
cat("样本内 C-index =", round(as.numeric(concordance(Surv(d$OS.time,d$OS)~d$rs, reverse=TRUE)$concordance),3), "\n")
r <- timeROC(T=d$OS.time, delta=d$OS, marker=d$rs, cause=1, weighting="marginal", times=c(365,1095,1825), iid=FALSE)
cat("样本内 AUC 1/3/5y =", paste(round(as.numeric(r$AUC),3), collapse=" / "), "\n")
set.seed(42)
K <- 5; folds <- sample(rep(1:K, length.out=nrow(d)))
pred <- rep(NA, nrow(d))
for (k in 1:K) {
  tr <- d[folds!=k, ]; te <- d[folds==k, ]
  f <- coxph(Surv(OS.time,OS) ~ PLG+CD79A+MMP7+ANGPT1, data=tr)
  pred[folds==k] <- predict(f, newdata=te, type="risk")
}
cat("5 折交叉验证 C-index =", round(as.numeric(concordance(Surv(d$OS.time,d$OS)~pred, reverse=TRUE)$concordance),3), "\n")
