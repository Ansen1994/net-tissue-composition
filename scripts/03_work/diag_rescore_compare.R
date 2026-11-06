# ============================================================
# risk score 标准化方案对比（四队列 × 四方案）
# S0 现行：训练集用原始 log2 线性预测子；验证集用队列内 z 后套同一系数
# S1 全队列基因级 z（含训练集），系数在 z 尺度上重新拟合
# S2 全队列基因级百分位（经验 CDF，队列内），系数在百分位尺度上重新拟合
# S3 队列内中心化 + 训练集 SD 缩放（保留训练集相对权重），系数沿用原始拟合
# ============================================================
suppressMessages({
  library(survival); library(data.table); library(GEOquery)
  library(dplyr); library(tibble); library(stringr)
})
hub <- c("PLG","CD79A","MMP7","ANGPT1")
options(stringsAsFactors = FALSE)

# ---------------- 1. TCGA ----------------
tc <- read.csv("Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv", row.names = 1, check.names = FALSE)
tc <- tc[complete.cases(tc[, c("OS.time","OS",hub)]), ]
tc <- tc[tc$OS.time > 0, ]
tc_dat <- data.frame(OS.time = tc$OS.time, OS = tc$OS, tc[, hub], check.names = FALSE)

# ---------------- 2. GSE20685 ----------------
g1 <- getGEO(filename = "../raw/GSE20685/GSE20685_series_matrix.txt.gz", getGPL = FALSE)
e1 <- exprs(g1); p1 <- pData(g1)
probe <- read.delim("../raw/annotation/GPL570-55999.txt", sep = "\t", comment.char = "#",
                    stringsAsFactors = FALSE, check.names = FALSE)
pcol <- grep("Gene Symbol", names(probe), value = TRUE)[1]
map <- data.frame(ID = probe[[1]], symbol = probe[[pcol]]) %>%
  filter(!is.na(symbol), symbol != "") %>%
  mutate(symbol = sapply(strsplit(as.character(symbol), " /// "), "[", 1))
keep1 <- rep(TRUE, nrow(p1))  # 与原脚本一致：GSE20685 全队列 327，未做 ER 过滤
e1 <- e1[, keep1]; p1 <- p1[keep1, ]
idx <- match(rownames(e1), map$ID)
e1 <- e1[!is.na(idx), ]; map1 <- map[idx[!is.na(idx)], ]
e1 <- rowsum(e1, map1$symbol); e1 <- e1 / as.vector(table(map1$symbol)[rownames(e1)]) # 均值合并
x1 <- as.data.frame(t(e1[intersect(hub, rownames(e1)), , drop = FALSE]))
d1 <- data.frame(OS.time = as.numeric(p1[["follow_up_duration (years):ch1"]]) * 365,
                 OS = as.numeric(p1[["event_death:ch1"]]))
d1 <- cbind(d1, x1); d1 <- d1[complete.cases(d1) & d1$OS.time > 0, ]

# ---------------- 3. METABRIC ----------------
mb <- fread("../raw/METABRIC/data_mrna_illumina_microarray.txt", data.table = FALSE)
mb <- mb[mb$Hugo_Symbol %in% hub, ]
mm <- mb %>% dplyr::select(-2) %>% dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), max)) %>%
  tibble::column_to_rownames("Hugo_Symbol")
x2 <- as.data.frame(t(as.matrix(mm))); x2$PATIENT_ID <- gsub("\\.", "-", rownames(x2))
cl2 <- fread("../raw/METABRIC/data_clinical_patient.txt", data.table = FALSE, skip = 4)
cl2 <- cl2 %>% dplyr::transmute(PATIENT_ID, ER_Status = ER_IHC, OS_Status_Raw = OS_STATUS,
                                OS.time = as.numeric(OS_MONTHS) * 30) %>%
  dplyr::filter(grepl("pos|positive", ER_Status, ignore.case = TRUE)) %>%
  dplyr::mutate(OS = ifelse(grepl("DECEASED", OS_Status_Raw, ignore.case = TRUE), 1, 0))
d2 <- merge(cl2, x2, by = "PATIENT_ID"); d2 <- d2[complete.cases(d2) & d2$OS.time > 0, ]

# ---------------- 4. GSE96058 ----------------
ge <- fread("../raw/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",
            data.table = FALSE)
rownames(ge) <- ge[, 1]; gm <- as.matrix(ge[, -1])
x3 <- as.data.frame(t(gm[intersect(hub, rownames(gm)), , drop = FALSE]))
x3$title <- rownames(x3)
s1 <- getGEO(filename = "../raw/GSE96058/GSE96058-GPL11154_series_matrix.txt.gz", getGPL = FALSE)
s2 <- getGEO(filename = "../raw/GSE96058/GSE96058-GPL18573_series_matrix.txt.gz", getGPL = FALSE)
c1 <- as.data.frame(pData(s1)); c2 <- as.data.frame(pData(s2))
cm <- rbind(c1[, intersect(names(c1), names(c2))], c2[, intersect(names(c1), names(c2))])
cm <- cm[!duplicated(cm$title), ]
cl3 <- data.frame(title = cm$title, ER_Status = cm[["er status:ch1"]],
                  OS = as.numeric(cm[["overall survival event:ch1"]]),
                  OS.time = as.numeric(cm[["overall survival days:ch1"]])) %>%
  dplyr::filter(ER_Status == "1")
d3 <- merge(cl3, x3, by = "title"); d3 <- d3[complete.cases(d3) & d3$OS.time > 0, ]

# 统一成 OS.time / OS / 四基因
mk <- function(df, nm) {
  out <- data.frame(OS.time = df$OS.time, OS = df$OS)
  for (g in hub) out[[g]] <- df[[g]]
  cat(sprintf("%-8s n = %4d  事件 = %3d\n", nm, nrow(out), sum(out$OS)))
  out
}
cat("===== 各队列入组情况 =====\n")
D <- list(TCGA = mk(tc_dat, "TCGA"), GSE20685 = mk(d1, "GSE20685"),
          METABRIC = mk(d2, "METABRIC"), GSE96058 = mk(d3, "GSE96058"))

cat("\n===== 各队列 4 基因零值比例与离散度 =====\n")
for (nm in names(D)) {
  cat("--", nm, "\n")
  for (g in hub) {
    v <- D[[nm]][[g]]
    cat(sprintf("   %-7s 零值 %5.1f%%  mean %7.3f  sd %6.3f\n",
                g, mean(v == 0) * 100, mean(v), sd(v)))
  }
}

# ---------------- 方案构造 ----------------
zscore <- function(v) (v - mean(v)) / sd(v)
pct <- function(v) rank(v, ties.method = "average") / length(v)

score_sets <- list()
# S0 现行：训练集原始尺度拟合
f_raw <- coxph(Surv(OS.time, OS) ~ PLG + CD79A + MMP7 + ANGPT1, data = D$TCGA)
b_raw <- coef(f_raw)
cat("\n===== S0 训练集原始尺度系数（per unit log2）=====\n"); print(round(b_raw, 4))

S0 <- list(); for (nm in names(D)) {
  d <- D[[nm]]
  if (nm == "TCGA") S0[[nm]] <- as.numeric(as.matrix(d[, hub]) %*% b_raw)
  else {
    dz <- as.data.frame(lapply(d[, hub], zscore))
    S0[[nm]] <- as.numeric(as.matrix(dz) %*% b_raw)
  }
}

# S1 全队列 z + 在 z 尺度重拟合
S1 <- list()
for (nm in names(D)) {
  dz <- as.data.frame(lapply(D[[nm]][, hub], zscore)); names(dz) <- hub
  S1[[nm]] <- dz
}
f1 <- coxph(Surv(D$TCGA$OS.time, D$TCGA$OS) ~ PLG + CD79A + MMP7 + ANGPT1, data = S1$TCGA)
b1 <- coef(f1); cat("\n===== S1 z 尺度系数（per SD）=====\n"); print(round(b1, 4))
S1s <- lapply(S1, function(dz) as.numeric(as.matrix(dz[, hub]) %*% b1))

# S2 全队列百分位 + 重拟合
S2 <- list()
for (nm in names(D)) {
  dp <- as.data.frame(lapply(D[[nm]][, hub], pct)); names(dp) <- hub
  S2[[nm]] <- dp
}
f2 <- coxph(Surv(D$TCGA$OS.time, D$TCGA$OS) ~ PLG + CD79A + MMP7 + ANGPT1, data = S2$TCGA)
b2 <- coef(f2); cat("\n===== S2 百分位尺度系数 =====\n"); print(round(b2, 4))
S2s <- lapply(S2, function(dp) as.numeric(as.matrix(dp[, hub]) %*% b2))

# S3 队列内中心化 + 训练集 SD 缩放，沿用原始系数
sd_tr <- sapply(D$TCGA[, hub], sd)
S3s <- list(); for (nm in names(D)) {
  d <- D[[nm]]
  dz <- as.data.frame(lapply(hub, function(g) (d[[g]] - mean(d[[g]])) / sd_tr[g]))
  names(dz) <- hub
  S3s[[nm]] <- as.numeric(as.matrix(dz) %*% b_raw)
}
cat("\n===== 训练集 SD =====\n"); print(round(sd_tr, 3))

# ---------------- 评价 ----------------
ev <- function(score, d, nm) {
  s <- score; if (anyNA(s) || sd(s, na.rm=TRUE) == 0) return(NULL)
  ci <- tryCatch(as.numeric(concordance(Surv(d$OS.time, d$OS) ~ s, reverse = TRUE)$concordance),
                 error = function(e) NA)
  sz <- as.numeric(scale(s))
  f <- summary(coxph(Surv(d$OS.time, d$OS) ~ sz))
  g <- as.numeric(s > median(s))
  lr <- survdiff(Surv(d$OS.time, d$OS) ~ g)
  fg <- summary(coxph(Surv(d$OS.time, d$OS) ~ g))
  data.frame(cohort = nm, Cindex = round(ci, 3),
             HRperSD = round(f$conf.int[1], 3),
             p_perSD = signif(f$coefficients[1, 5], 3),
             HRsplit = round(fg$conf.int[1], 3),
             p_split = signif(1 - pchisq(lr$chisq, 1), 3))
}
cat("\n\n===== 四方案对比（C-index / 每 SD 的 HR / 中位分组 HR）=====\n")
for (scheme in c("S0","S1","S2","S3")) {
  S <- switch(scheme, S0 = S0, S1 = S1s, S2 = S2s, S3 = S3s)
  cat("\n--", scheme, "--\n")
  for (nm in names(D)) {
    r <- ev(S[[nm]], D[[nm]], nm)
    if (!is.null(r)) print(r, row.names = FALSE)
  }
}
