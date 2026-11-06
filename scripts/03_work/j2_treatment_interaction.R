# ==============================================================================
# 判定 2：这个评分能不能预测「谁从治疗中获益」（predictive），
#          区别于它已经失败的「预测预后」（prognostic）
#
# 关键：检验交互项 RiskScore × Treatment。
#       主效应为零 ≠ 交互项为零，两者数学上独立。
# 队列：METABRIC（有 CHEMOTHERAPY / HORMONE_THERAPY）
#       GSE96058（有 chemo treated / endocrine treated）
# 输出：j2_results.csv
# ==============================================================================
suppressMessages({ library(data.table); library(survival) })

# ------------------------------------------------------------------ 工具
inter_test <- function(d, treat_col, label, covariates = c("Age")) {
  d <- d[!is.na(d[[treat_col]]) & !is.na(d$OS.time) & d$OS.time > 0, ]
  d$T <- as.numeric(d[[treat_col]])
  d$s <- as.numeric(scale(d$RiskScore))
  covs <- intersect(covariates, names(d))
  covs <- covs[sapply(d[covs], function(x) sum(!is.na(x)) > 0.9 * nrow(d))]
  fml <- as.formula(paste("Surv(OS.time, OS) ~ s * T +", paste(covs, collapse = " + ")))
  fit <- coxph(fml, data = d)
  sm <- summary(fit)
  ix <- grep("^s:T|^T:s", rownames(sm$coefficients))
  if (!length(ix)) return(NULL)

  # 分层 HR：治疗组内与非治疗组内，risk 每 SD 的 HR
  hr_in <- function(sub) {
    if (nrow(sub) < 60 || sum(sub$OS) < 15) return(c(NA, NA, NA, NA))
    f <- coxph(Surv(OS.time, OS) ~ s, data = sub)
    z <- summary(f)
    c(z$conf.int[1], z$conf.int[3], z$conf.int[4], z$coefficients[1,5])
  }
  h1 <- hr_in(d[d$T == 1, ]); h0 <- hr_in(d[d$T == 0, ])

  cat(sprintf("\n---- %s ----\n", label))
  cat(sprintf("  n=%d  治疗组 %d (%.0f%%)  事件 %d\n",
              nrow(d), sum(d$T == 1), 100*mean(d$T == 1), sum(d$OS)))
  cat(sprintf("  risk 主效应 (T=0 基准): HR/SD = %.2f (%.2f-%.2f), p = %.3f\n",
              sm$conf.int[1,1], sm$conf.int[1,3], sm$conf.int[1,4], sm$coefficients[1,5]))
  cat(sprintf("  治疗主效应:            HR = %.2f (%.2f-%.2f), p = %.3f\n",
              sm$conf.int["T",1], sm$conf.int["T",3], sm$conf.int["T",4],
              sm$coefficients["T",5]))
  cat(sprintf("  >>> 交互项 s:T:        HR = %.2f (%.2f-%.2f), p = %.4f  %s\n",
              sm$conf.int[ix[1],1], sm$conf.int[ix[1],3], sm$conf.int[ix[1],4],
              sm$coefficients[ix[1],5],
              ifelse(sm$coefficients[ix[1],5] < 0.05, "<<< 显著", "")))
  cat(sprintf("  治疗组内 risk HR/SD = %.2f (%.2f-%.2f) p=%.3f\n", h1[1],h1[2],h1[3],h1[4]))
  cat(sprintf("  未治组内 risk HR/SD = %.2f (%.2f-%.2f) p=%.3f\n", h0[1],h0[2],h0[3],h0[4]))

  data.frame(cohort = label, n = nrow(d), n_treated = sum(d$T == 1),
             events = sum(d$OS),
             HR_main = sm$conf.int[1,1], p_main = sm$coefficients[1,5],
             HR_treat = sm$conf.int["T",1], p_treat = sm$coefficients["T",5],
             HR_interaction = sm$conf.int[ix[1],1],
             lo = sm$conf.int[ix[1],3], hi = sm$conf.int[ix[1],4],
             p_interaction = sm$coefficients[ix[1],5],
             HR_in_treated = h1[1], p_in_treated = h1[4],
             HR_in_untreated = h0[1], p_in_untreated = h0[4],
             stringsAsFactors = FALSE)
}

res <- list()

# ==============================================================================
# METABRIC
# ==============================================================================
cat("========== METABRIC ==========\n")
clin <- read.delim("../raw/METABRIC/data_clinical_patient.txt",
                   comment.char = "#", stringsAsFactors = FALSE)
cat("CHEMOTHERAPY:\n"); print(table(clin$CHEMOTHERAPY, useNA = "ifany"))
cat("HORMONE_THERAPY:\n"); print(table(clin$HORMONE_THERAPY, useNA = "ifany"))

dmb <- read.csv("Table_S6_METABRIC_Validation_Result.csv", stringsAsFactors = FALSE)
dmb$PID <- gsub("\\.", "-", dmb$PATIENT_ID)
clin$PID <- gsub("\\.", "-", clin$PATIENT_ID)
m <- merge(dmb, clin[, c("PID","CHEMOTHERAPY","HORMONE_THERAPY","NPI")], by = "PID")
cat("合并后 n =", nrow(m), "\n")
m$OS.time <- as.numeric(m$OS.time)
m$LN <- as.numeric(m$LymphNode)
m$Chemo <- ifelse(grepl("YES", m$CHEMOTHERAPY, ignore.case = TRUE), 1,
           ifelse(grepl("NO", m$CHEMOTHERAPY, ignore.case = TRUE), 0, NA))
m$Horm <- ifelse(grepl("YES", m$HORMONE_THERAPY, ignore.case = TRUE), 1,
          ifelse(grepl("NO", m$HORMONE_THERAPY, ignore.case = TRUE), 0, NA))
cat(sprintf("Chemo 可用 %d (接受 %d)  Horm 可用 %d (接受 %d)\n",
            sum(!is.na(m$Chemo)), sum(m$Chemo == 1, na.rm=TRUE),
            sum(!is.na(m$Horm)), sum(m$Horm == 1, na.rm=TRUE)))

res[[1]] <- inter_test(m, "Chemo", "METABRIC_chemo", c("Age","LN"))
res[[2]] <- inter_test(m, "Horm",  "METABRIC_endocrine", c("Age","LN"))

# ==============================================================================
# GSE96058
# ==============================================================================
cat("\n\n========== GSE96058 ==========\n")
suppressMessages(library(GEOquery))
g1 <- getGEO(filename = "../raw/GSE96058/GSE96058-GPL11154_series_matrix.txt.gz", getGPL = FALSE)
g2 <- getGEO(filename = "../raw/GSE96058/GSE96058-GPL18573_series_matrix.txt.gz", getGPL = FALSE)
cl <- rbind(pData(g1)[, intersect(c("title","chemo treated:ch1","endocrine treated:ch1"), names(pData(g1)))],
            pData(g2)[, intersect(c("title","chemo treated:ch1","endocrine treated:ch1"), names(pData(g2)))])
names(cl) <- c("title","chemo","endo")
d96 <- read.csv("Table_S6_GSE96058_Validation_Result.csv", stringsAsFactors = FALSE)
m9 <- merge(d96, cl, by = "title")
cat("合并后 n =", nrow(m9), "\n")
m9$OS.time <- as.numeric(m9$OS.time)
m9$LN <- as.numeric(as.character(m9$LymphNode))
m9$Chemo <- as.numeric(m9$chemo)
m9$Horm  <- as.numeric(m9$endo)
cat(sprintf("Chemo 可用 %d (接受 %d)  Endo 可用 %d (接受 %d)\n",
            sum(!is.na(m9$Chemo)), sum(m9$Chemo==1, na.rm=TRUE),
            sum(!is.na(m9$Horm)), sum(m9$Horm==1, na.rm=TRUE)))

res[[3]] <- inter_test(m9, "Chemo", "GSE96058_chemo", c("Age","LN"))
res[[4]] <- inter_test(m9, "Horm",  "GSE96058_endocrine", c("Age","LN"))

R <- do.call(rbind, Filter(Negate(is.null), res))
write.csv(R, "j2_results.csv", row.names = FALSE)

cat("\n\n================ 判定 2 汇总 ================\n")
cat(sprintf("%-22s %6s %8s %10s %8s\n", "队列_治疗", "n", "交互HR", "交互p", "解读"))
for (i in seq_len(nrow(R))) {
  cat(sprintf("%-22s %6d %8.2f %10.4f %8s\n", R$cohort[i], R$n[i],
              R$HR_interaction[i], R$p_interaction[i],
              ifelse(R$p_interaction[i] < 0.05, "有信号",
              ifelse(R$p_interaction[i] < 0.15, "边缘", "无"))))
}
cat("\n判读：交互 p < 0.05 且两个队列方向一致 → 可作为 predictive biomarker 报告\n")
cat("      交互 p > 0.15 或方向相反 → 「指导用药」这条路关闭\n")
