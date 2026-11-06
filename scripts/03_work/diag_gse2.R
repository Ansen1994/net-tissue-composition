suppressMessages({library(data.table); library(GEOquery); library(survival); library(survminer); library(dplyr)})
cat("读取表达矩阵...\n")
exp_raw <- fread("../raw/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz", data.table=FALSE)
rownames(exp_raw) <- exp_raw[,1]; exp_matrix <- exp_raw[,-1]
hub <- c("PLG","CD79A","MMP7","ANGPT1")
sc_all <- as.data.frame(t(exp_matrix[intersect(hub, rownames(exp_matrix)), ]))
sc_all <- as.data.frame(scale(sc_all))
sc_all$RiskScore <- 0.1102*sc_all$PLG + 0.2295*sc_all$ANGPT1 - 0.0854*sc_all$CD79A - 0.0851*sc_all$MMP7
cat("全部样本 n =", nrow(sc_all), "\n")
cat("全队列中位数 =", round(median(sc_all$RiskScore),4), " 按中位数分组: low", sum(sc_all$RiskScore <= median(sc_all$RiskScore)), "high", sum(sc_all$RiskScore > median(sc_all$RiskScore)), "\n")

g1 <- getGEO(filename="../raw/GSE96058/GSE96058-GPL11154_series_matrix.txt.gz", getGPL=FALSE)
g2 <- getGEO(filename="../raw/GSE96058/GSE96058-GPL18573_series_matrix.txt.gz", getGPL=FALSE)
c1 <- as.data.frame(pData(g1)); c2 <- as.data.frame(pData(g2))
clin <- rbind(c1[,intersect(names(c1),names(c2))], c2[,intersect(names(c1),names(c2))])
clin <- clin[!duplicated(clin$title),]
clin$ER <- clin[["er status:ch1"]]
cat("临床样本", nrow(clin), " ER+ =", sum(clin$ER %in% c("positive","Positive","1"), na.rm=TRUE), "\n")
m <- merge(sc_all, clin[,c("title","ER","overall survival event:ch1","overall survival days:ch1")], by.x="row.names", by.y="title")
cat("合并后", nrow(m), "\n")

m$OS <- as.numeric(m[["overall survival event:ch1"]])
m$OS.time <- as.numeric(m[["overall survival days:ch1"]])
mm <- m[!is.na(m$OS) & !is.na(m$OS.time),]
cat("全队列有随访 n =", nrow(mm), "\n")
r_all <- surv_cutpoint(mm, time="OS.time", event="OS", variables="RiskScore")
cat("全队列 surv_cutpoint 截点 =", round(r_all$cutpoint$cutpoint,4),
    " low", sum(mm$RiskScore <= r_all$cutpoint$cutpoint), "high", sum(mm$RiskScore > r_all$cutpoint$cutpoint), "\n")
er <- mm[mm$ER %in% c("positive","Positive","1"),]
cat("ER+ 有随访 n =", nrow(er), "\n")
cat("用全队列截点划分 ER+ : low", sum(er$RiskScore <= r_all$cutpoint$cutpoint), "high", sum(er$RiskScore > r_all$cutpoint$cutpoint), "\n")
cat("ER+ 内 median 划分: low", sum(er$RiskScore <= median(er$RiskScore)), "high", sum(er$RiskScore > median(er$RiskScore)), "\n")
r_er <- surv_cutpoint(er, time="OS.time", event="OS", variables="RiskScore")
cat("ER+ 内 surv_cutpoint: low", sum(er$RiskScore <= r_er$cutpoint$cutpoint), "high", sum(er$RiskScore > r_er$cutpoint$cutpoint), "\n")
