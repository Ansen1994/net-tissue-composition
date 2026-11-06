# 尺度一致性检验：验证队列的 risk score 若改用原始 log2 表达（与训练集一致），结论是否改变
suppressMessages({ library(data.table) })
e <- fread("../raw/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",
           data.table = FALSE)
rownames(e) <- e[[1]]; e <- as.matrix(e[, -1])
fd <- read.csv("Table_S6_GSE96058_Validation_Result.csv", stringsAsFactors = FALSE)
s <- intersect(colnames(e), fd$title)
cat("样本:", length(s), "\n")

hub <- c("PLG","CD79A","MMP7","ANGPT1")
coef <- c(PLG = 0.1102, ANGPT1 = 0.2295, CD79A = -0.0854, MMP7 = -0.0851)

raw4 <- t(e[hub, s])
z4   <- as.data.frame(scale(raw4))
rs_raw <- as.numeric(as.matrix(raw4) %*% coef[hub])
rs_z   <- as.numeric(as.matrix(z4)   %*% coef[hub])

prolif <- c("ANLN","BIRC5","BUB1","CCNB1","CDC20","CDC25B","CEP55","CKS2","EXO1",
            "KPNA2","MELK","MKI67","MYBL2","NUF2","ORC6L","PTTG1","RRM2","TYMS","UBE2C")
pg <- intersect(prolif, rownames(e))
pm <- e[pg, s]; pz <- t(scale(t(pm))); pz[is.na(pz)] <- 0
prolif_score <- colMeans(pz)

cat("\n=== GSE96058 ===\n")
cat(sprintf("4 基因原始值范围: PLG %.2f-%.2f | CD79A %.2f-%.2f | MMP7 %.2f-%.2f | ANGPT1 %.2f-%.2f\n",
            min(raw4[,1]), max(raw4[,1]), min(raw4[,2]), max(raw4[,2]),
            min(raw4[,3]), max(raw4[,3]), min(raw4[,4]), max(raw4[,4])))
cat(sprintf("risk score (原始值): 范围 %.2f ~ %.2f, 均值 %.2f\n", min(rs_raw), max(rs_raw), mean(rs_raw)))
cat(sprintf("risk score (z-score): 范围 %.2f ~ %.2f, 均值 %.2f\n", min(rs_z), max(rs_z), mean(rs_z)))
cat(sprintf("两种算法下 risk score 的 Spearman 相关: %.3f\n",
            cor(rs_raw, rs_z, method = "spearman")))
cat(sprintf("与增殖模块相关 (原始值算法): rho = %.3f\n",
            cor(rs_raw, prolif_score, method = "spearman")))
cat(sprintf("与增殖模块相关 (z-score 算法): rho = %.3f\n",
            cor(rs_z,   prolif_score, method = "spearman")))
cat(sprintf("已发表表里的 RiskScore 与原始值算法相关: %.3f\n",
            cor(fd$RiskScore[match(s, fd$title)], rs_raw, method = "spearman")))

suppressMessages({library(survival); library(survminer)})
d <- data.frame(title=s, rs_raw=rs_raw, rs_z=rs_z, prolif=prolif_score)
d <- merge(d, fd[,c("title","OS","OS.time")], by="title")
d <- d[!is.na(d$OS)&!is.na(d$OS.time)&d$OS.time>0,]
cat("\n=== 中位数分组 KM (n=", nrow(d), ", 事件=", sum(d$OS), ") ===\n")
for (v in c("rs_raw","rs_z")) {
  d$G <- as.numeric(d[[v]] > median(d[[v]]))
  f <- summary(coxph(Surv(OS.time,OS) ~ G, data=d))
  lr <- survdiff(Surv(OS.time,OS) ~ G, data=d)
  cat(sprintf("  %-7s HR = %.2f (%.2f-%.2f), Cox p = %.4f, log-rank p = %.4f\n",
      v, f$conf.int[1], f$conf.int[3], f$conf.int[4], f$coefficients[1,5],
      1-pchisq(lr$chisq,1)))
}
# surv_cutpoint 最优截点
r <- surv_cutpoint(d, time="OS.time", event="OS", variables="rs_raw")
cat(sprintf("  surv_cutpoint 截点: %.3f, 分组 %d / %d\n", r$cutpoint$cutpoint,
            sum(d$rs_raw <= r$cutpoint$cutpoint), sum(d$rs_raw > r$cutpoint$cutpoint)))
d$G2 <- as.numeric(d$rs_raw > r$cutpoint$cutpoint)
f2 <- summary(coxph(Surv(OS.time,OS) ~ G2, data=d))
lr2 <- survdiff(Surv(OS.time,OS) ~ G2, data=d)
cat(sprintf("  最优截点分组: HR = %.2f (%.2f-%.2f), log-rank p = %.4f\n",
            f2$conf.int[1], f2$conf.int[3], f2$conf.int[4], 1-pchisq(lr2$chisq,1)))
