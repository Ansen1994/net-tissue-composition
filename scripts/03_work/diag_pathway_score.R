# ============================================================
# 关键判断：132 个候选 NETs 基因的整体评分，能否跨队列复现
# 评分方式：队列内基因级 z-score 后取均值（与 4 基因签名同一构造逻辑）
# ============================================================
suppressMessages({
  library(survival); library(timeROC); library(data.table)
  library(GEOquery); library(dplyr); library(tibble)
})
pat <- read.csv("Table_S2_Target_NETs_Matrix.csv", row.names = 1, check.names = FALSE)
G132 <- rownames(pat)
cat("候选基因数 =", length(G132), "\n")

eval_score <- function(s, d, nm) {
  ok <- !is.na(s) & !is.na(d$OS) & !is.na(d$OS.time) & d$OS.time > 0
  s <- s[ok]; dd <- d[ok, ]
  ci <- as.numeric(concordance(Surv(dd$OS.time, dd$OS) ~ s, reverse = TRUE)$concordance)
  r <- tryCatch(timeROC(T = dd$OS.time, delta = dd$OS, marker = s, cause = 1,
                        weighting = "marginal", times = c(365, 1095, 1825), iid = FALSE),
                error = function(e) NULL)
  g <- as.numeric(s > median(s))
  lr <- survdiff(Surv(dd$OS.time, dd$OS) ~ g)
  fg <- summary(coxph(Surv(dd$OS.time, dd$OS) ~ g))
  sz <- as.numeric(scale(s)); f <- summary(coxph(Surv(dd$OS.time, dd$OS) ~ sz))
  cat(sprintf("%-9s n=%4d 事件=%3d | C-index %.3f | AUC %s | HR/SD %.2f (p=%.2g) | 中位分组 HR %.2f (logrank p=%.3f)\n",
              nm, nrow(dd), sum(dd$OS), ci,
              if (is.null(r)) "ERR" else paste(round(as.numeric(r$AUC), 3), collapse = "/"),
              f$conf.int[1], f$coefficients[1, 5], fg$conf.int[1], 1 - pchisq(lr$chisq, 1)))
  invisible(ci)
}

# ---- TCGA ----
tc <- read.csv("Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv", row.names = 1, check.names = FALSE)
common <- intersect(G132, colnames(tc))
cat("TCGA 命中基因 =", length(common), "\n")
m <- as.matrix(pat[common, , drop = FALSE]); rownames(m) <- common
colnames(m) <- substr(colnames(m), 1, 12)
# 用样本 ID 对齐
ids <- intersect(colnames(m), rownames(tc))
cat("TCGA 可对齐样本 =", length(ids), "| Table_S3 样本 =", nrow(tc), "\n")
if (length(ids) < 100) {
  # 列名与行名对不上时，按顺序对齐（Table_S3 由 Table_S2 转置派生）
  m2 <- as.matrix(pat[common, seq_len(min(ncol(pat), nrow(tc))), drop = FALSE])
  s_tcga <- rowMeans(scale(t(m2)), na.rm = TRUE)
  d_tcga <- tc[seq_len(length(s_tcga)), ]
} else {
  m2 <- m[, ids, drop = FALSE]
  s_tcga <- rowMeans(scale(t(m2)), na.rm = TRUE)
  d_tcga <- tc[ids, ]
}
eval_score(as.numeric(s_tcga), d_tcga, "TCGA")

# ---- GSE20685 ----
g1 <- getGEO(filename = "../raw/GSE20685/GSE20685_series_matrix.txt.gz", getGPL = FALSE)
e1 <- exprs(g1); p1 <- pData(g1)
probe <- read.delim("../raw/annotation/GPL570-55999.txt", sep = "\t", comment.char = "#",
                    stringsAsFactors = FALSE, check.names = FALSE)
map <- data.frame(ID = probe[[1]], symbol = probe[["Gene Symbol"]], stringsAsFactors = FALSE) %>%
  filter(!is.na(symbol), symbol != "") %>%
  mutate(symbol = sapply(strsplit(symbol, " /// "), "[", 1))
idx <- match(rownames(e1), map$ID)
e1 <- e1[!is.na(idx), , drop = FALSE]; map1 <- map[idx[!is.na(idx)], ]
e1 <- rowsum(e1, map1$symbol); e1 <- e1 / as.vector(table(map1$symbol)[rownames(e1)])
gg <- intersect(G132, rownames(e1))
cat("GSE20685 命中基因 =", length(gg), "\n")
s1 <- rowMeans(scale(t(e1[gg, , drop = FALSE])), na.rm = TRUE)
d1 <- data.frame(OS.time = as.numeric(p1[["follow_up_duration (years):ch1"]]) * 365,
                 OS = as.numeric(p1[["event_death:ch1"]]))
eval_score(as.numeric(s1), d1, "GSE20685")

# ---- METABRIC ----
mb <- fread("../raw/METABRIC/data_mrna_illumina_microarray.txt", data.table = FALSE)
mb <- mb[mb$Hugo_Symbol %in% G132, ]
mm <- mb %>% dplyr::select(-2) %>% dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), max)) %>%
  tibble::column_to_rownames("Hugo_Symbol")
cat("METABRIC 命中基因 =", nrow(mm), "\n")
s2 <- rowMeans(scale(t(as.matrix(mm))), na.rm = TRUE)
cl2 <- fread("../raw/METABRIC/data_clinical_patient.txt", data.table = FALSE, skip = 4)
cl2 <- cl2 %>% dplyr::transmute(PATIENT_ID = gsub("\\.", "-", PATIENT_ID), ER_Status = ER_IHC,
                                OS_Status_Raw = OS_STATUS, OS.time = as.numeric(OS_MONTHS) * 30) %>%
  dplyr::filter(grepl("pos|positive", ER_Status, ignore.case = TRUE)) %>%
  dplyr::mutate(OS = ifelse(grepl("DECEASED", OS_Status_Raw, ignore.case = TRUE), 1, 0))
ds2 <- data.frame(PATIENT_ID = gsub("\\.", "-", names(s2)), s = as.numeric(s2))
d2 <- merge(cl2, ds2, by = "PATIENT_ID")
eval_score(d2$s, d2, "METABRIC")

# ---- GSE96058 ----
ge <- fread("../raw/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",
            data.table = FALSE)
rownames(ge) <- ge[, 1]; gm <- as.matrix(ge[, -1])
gg3 <- intersect(G132, rownames(gm))
cat("GSE96058 命中基因 =", length(gg3), "\n")
s3 <- rowMeans(scale(t(gm[gg3, , drop = FALSE])), na.rm = TRUE)
s1g <- getGEO(filename = "../raw/GSE96058/GSE96058-GPL11154_series_matrix.txt.gz", getGPL = FALSE)
s2g <- getGEO(filename = "../raw/GSE96058/GSE96058-GPL18573_series_matrix.txt.gz", getGPL = FALSE)
c1 <- as.data.frame(pData(s1g)); c2 <- as.data.frame(pData(s2g))
cm <- rbind(c1[, intersect(names(c1), names(c2))], c2[, intersect(names(c1), names(c2))])
cm <- cm[!duplicated(cm$title), ]
cl3 <- data.frame(title = cm$title, ER_Status = cm[["er status:ch1"]],
                  OS = as.numeric(cm[["overall survival event:ch1"]]),
                  OS.time = as.numeric(cm[["overall survival days:ch1"]])) %>%
  dplyr::filter(ER_Status == "1")
ds3 <- data.frame(title = colnames(gm), s = as.numeric(s3))
d3 <- merge(cl3, ds3, by = "title")
eval_score(d3$s, d3, "GSE96058")
