# ==============================================================================
# 判定 2b：GSE96058 化疗交互信号是真交互，还是适应证混杂？
# 化疗是"病重者才做"，高 risk 者更可能被化疗 → 化疗组内 risk HR 自然更高
# 三项排查：基线平衡 / 加肿瘤大小调整 / 淋巴结分层
# ==============================================================================
suppressMessages({ library(survival); library(GEOquery) })

g1 <- getGEO(filename = "../raw/GSE96058/GSE96058-GPL11154_series_matrix.txt.gz", getGPL = FALSE)
g2 <- getGEO(filename = "../raw/GSE96058/GSE96058-GPL18573_series_matrix.txt.gz", getGPL = FALSE)
p1 <- pData(g1); p2 <- pData(g2)
cols <- c("title","chemo treated:ch1","endocrine treated:ch1")
cl <- rbind(p1[, cols], p2[, cols]); names(cl) <- c("title","chemo","endo")
d <- read.csv("Table_S6_GSE96058_Validation_Result.csv", stringsAsFactors = FALSE)
m <- merge(d, cl, by = "title")
m$OS.time <- as.numeric(m$OS.time); m$OS <- as.numeric(m$OS)
m$LN <- as.numeric(as.character(m$LymphNode))
m$Chemo <- as.numeric(as.character(m$chemo))
m$Age <- as.numeric(m$Age); m$TumorSize <- as.numeric(m$TumorSize)
m <- m[!is.na(m$Chemo) & !is.na(m$OS.time) & m$OS.time > 0, ]

cat("================ 1. 化疗组与未化疗组的基线差异 ================\n")
cat("（若 RiskScore 在化疗组显著更高，说明存在适应证混杂）\n\n")
for (v in c("RiskScore","Age","TumorSize","LN")) {
  a <- m[[v]][m$Chemo == 1]; b <- m[[v]][m$Chemo == 0]
  wt <- suppressWarnings(wilcox.test(a, b))
  cat(sprintf("  %-10s 化疗组 %.2f ± %.2f   未化疗 %.2f ± %.2f   p = %.2e\n",
              v, mean(a, na.rm=TRUE), sd(a, na.rm=TRUE),
              mean(b, na.rm=TRUE), sd(b, na.rm=TRUE), wt$p.value))
}
cat(sprintf("  %-10s 化疗组事件率 %.1f%%   未化疗 %.1f%%\n", "event rate",
            100*mean(m$OS[m$Chemo==1]), 100*mean(m$OS[m$Chemo==0])))

m$s <- as.numeric(scale(m$RiskScore))
run <- function(fml, lab, dd = m) {
  f <- coxph(as.formula(fml), data = dd)
  sm <- summary(f)
  ix <- grep("s:Chemo|Chemo:s", rownames(sm$coefficients))
  if (!length(ix)) return(invisible(NULL))
  cat(sprintf("  %-42s 交互 HR=%.2f (%.2f-%.2f)  p=%.4f   n=%d\n", lab,
              sm$conf.int[ix[1],1], sm$conf.int[ix[1],3], sm$conf.int[ix[1],4],
              sm$coefficients[ix[1],5], nrow(dd)))
}

cat("\n================ 2. 逐步加入协变量，看交互是否衰减 ================\n")
run("Surv(OS.time, OS) ~ s * Chemo", "仅交互，不调整")
run("Surv(OS.time, OS) ~ s * Chemo + Age", "+ 年龄")
run("Surv(OS.time, OS) ~ s * Chemo + Age + LN", "+ 年龄 + 淋巴结")
run("Surv(OS.time, OS) ~ s * Chemo + Age + LN + TumorSize", "+ 年龄 + 淋巴结 + 肿瘤大小")
cat("\n判读：若加入肿瘤大小后交互大幅衰减或变不显著 → 适应证混杂\n")

cat("\n================ 3. 淋巴结分层内的交互 ================\n")
for (lv in sort(unique(na.omit(m$LN)))) {
  sub <- m[!is.na(m$LN) & m$LN == lv, ]
  if (nrow(sub) < 200) next
  f <- coxph(Surv(OS.time, OS) ~ s * Chemo + Age, data = sub)
  sm <- summary(f)
  ix <- grep("s:Chemo|Chemo:s", rownames(sm$coefficients))
  cat(sprintf("  LN=%d (n=%d, 化疗 %d, 事件 %d): 交互 HR=%.2f p=%.3f\n",
              lv, nrow(sub), sum(sub$Chemo==1), sum(sub$OS),
              sm$conf.int[ix[1],1], sm$coefficients[ix[1],5]))
}

cat("\n================ 4. 同一模型在 METABRIC（化疗 9%）================\n")
cat("检验力不足是 METABRIC 无法验证的主因：\n")
clin <- read.delim("../raw/METABRIC/data_clinical_patient.txt", comment.char = "#", stringsAsFactors = FALSE)
dmb <- read.csv("Table_S6_METABRIC_Validation_Result.csv", stringsAsFactors = FALSE)
dmb$PID <- gsub("\\.","-",dmb$PATIENT_ID); clin$PID <- gsub("\\.","-",clin$PATIENT_ID)
mm <- merge(dmb, clin[, c("PID","CHEMOTHERAPY")], by="PID")
mm$Chemo <- ifelse(grepl("YES", mm$CHEMOTHERAPY, ignore.case=TRUE), 1,
            ifelse(grepl("NO", mm$CHEMOTHERAPY, ignore.case=TRUE), 0, NA))
mm$OS.time <- as.numeric(mm$OS.time); mm$LN <- as.numeric(mm$LymphNode)
mm <- mm[!is.na(mm$Chemo), ]
mm$s <- as.numeric(scale(mm$RiskScore))
cat(sprintf("  化疗组 n=%d，事件=%d（交互检验力极低）\n",
            sum(mm$Chemo==1), sum(mm$OS[mm$Chemo==1])))
f <- coxph(Surv(OS.time, OS) ~ s * Chemo + Age + LN, data=mm)
sm <- summary(f); ix <- grep("s:Chemo|Chemo:s", rownames(sm$coefficients))
cat(sprintf("  交互 HR=%.2f (%.2f-%.2f) p=%.3f  ← 方向与 GSE96058 相反\n",
            sm$conf.int[ix[1],1], sm$conf.int[ix[1],3], sm$conf.int[ix[1],4],
            sm$coefficients[ix[1],5]))
