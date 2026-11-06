# ==============================================================================
# Table S9：四种跨平台标准化方案 × 四队列的判别力（导出 csv）
# 逻辑与 diag_rescore_compare.R 完全一致，仅把结果写成表供作图与补充材料使用
# ==============================================================================
suppressMessages({
  library(survival); library(data.table); library(GEOquery)
  library(dplyr); library(tibble); library(stringr)
})
hub <- c("PLG","CD79A","MMP7","ANGPT1")
options(stringsAsFactors = FALSE)

tc <- read.csv("Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv", row.names = 1, check.names = FALSE)
tc <- tc[complete.cases(tc[, c("OS.time","OS",hub)]), ]
tc <- tc[tc$OS.time > 0, ]
tc_dat <- data.frame(OS.time = tc$OS.time, OS = tc$OS, tc[, hub], check.names = FALSE)

g1 <- getGEO(filename = "../raw/GSE20685/GSE20685_series_matrix.txt.gz", getGPL = FALSE)
e1 <- exprs(g1); p1 <- pData(g1)
probe <- read.delim("../raw/annotation/GPL570-55999.txt", sep = "\t", comment.char = "#",
                    stringsAsFactors = FALSE, check.names = FALSE)
pcol <- grep("Gene Symbol", names(probe), value = TRUE)[1]
map <- data.frame(ID = probe[[1]], symbol = probe[[pcol]]) %>%
  filter(!is.na(symbol), symbol != "") %>%
  mutate(symbol = sapply(strsplit(as.character(symbol), " /// "), "[", 1))
idx <- match(rownames(e1), map$ID)
e1 <- e1[!is.na(idx), ]; map1 <- map[idx[!is.na(idx)], ]
e1 <- rowsum(e1, map1$symbol); e1 <- e1 / as.vector(table(map1$symbol)[rownames(e1)])
x1 <- as.data.frame(t(e1[intersect(hub, rownames(e1)), , drop = FALSE]))
d1 <- data.frame(OS.time = as.numeric(p1[["follow_up_duration (years):ch1"]]) * 365,
                 OS = as.numeric(p1[["event_death:ch1"]]))
d1 <- cbind(d1, x1); d1 <- d1[complete.cases(d1) & d1$OS.time > 0, ]

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

mk <- function(df, nm) {
  out <- data.frame(OS.time = df$OS.time, OS = df$OS)
  for (g in hub) out[[g]] <- df[[g]]
  out
}
D <- list(TCGA = mk(tc_dat, "TCGA"), GSE20685 = mk(d1, "GSE20685"),
          METABRIC = mk(d2, "METABRIC"), GSE96058 = mk(d3, "GSE96058"))

zscore <- function(v) (v - mean(v)) / sd(v)
pct <- function(v) rank(v, ties.method = "average") / length(v)

f_raw <- coxph(Surv(OS.time, OS) ~ PLG + CD79A + MMP7 + ANGPT1, data = D$TCGA)
b_raw <- coef(f_raw)
S0 <- list(); for (nm in names(D)) {
  d <- D[[nm]]
  if (nm == "TCGA") S0[[nm]] <- as.numeric(as.matrix(d[, hub]) %*% b_raw)
  else { dz <- as.data.frame(lapply(d[, hub], zscore)); S0[[nm]] <- as.numeric(as.matrix(dz) %*% b_raw) }
}
S1 <- list(); for (nm in names(D)) { dz <- as.data.frame(lapply(D[[nm]][, hub], zscore)); names(dz) <- hub; S1[[nm]] <- dz }
b1 <- coef(coxph(Surv(D$TCGA$OS.time, D$TCGA$OS) ~ PLG + CD79A + MMP7 + ANGPT1, data = S1$TCGA))
S1s <- lapply(S1, function(dz) as.numeric(as.matrix(dz[, hub]) %*% b1))
S2 <- list(); for (nm in names(D)) { dp <- as.data.frame(lapply(D[[nm]][, hub], pct)); names(dp) <- hub; S2[[nm]] <- dp }
b2 <- coef(coxph(Surv(D$TCGA$OS.time, D$TCGA$OS) ~ PLG + CD79A + MMP7 + ANGPT1, data = S2$TCGA))
S2s <- lapply(S2, function(dp) as.numeric(as.matrix(dp[, hub]) %*% b2))
sd_tr <- sapply(D$TCGA[, hub], sd)
S3s <- list(); for (nm in names(D)) {
  d <- D[[nm]]
  dz <- as.data.frame(lapply(hub, function(g) (d[[g]] - mean(d[[g]])) / sd_tr[g])); names(dz) <- hub
  S3s[[nm]] <- as.numeric(as.matrix(dz) %*% b_raw)
}

ev <- function(score, d, nm, scheme) {
  s <- score; if (anyNA(s) || sd(s, na.rm = TRUE) == 0) return(NULL)
  ci <- as.numeric(concordance(Surv(d$OS.time, d$OS) ~ s, reverse = TRUE)$concordance)
  sz <- as.numeric(scale(s))
  f <- summary(coxph(Surv(d$OS.time, d$OS) ~ sz))
  data.frame(scheme = scheme, cohort = nm, n = nrow(d), events = sum(d$OS),
             C_index = ci, HR_perSD = f$conf.int[1], p_perSD = f$coefficients[1, 5])
}
rows <- list()
for (scheme in c("S0","S1","S2","S3")) {
  S <- switch(scheme, S0 = S0, S1 = S1s, S2 = S2s, S3 = S3s)
  for (nm in names(D)) {
    r <- ev(S[[nm]], D[[nm]], nm, scheme)
    if (!is.null(r)) rows[[length(rows)+1]] <- r
  }
}
out <- do.call(rbind, rows)
out$C_index <- round(out$C_index, 3)
out$HR_perSD <- round(out$HR_perSD, 3)
write.csv(out, "Table_S9_Scaling_Scheme_Cindex.csv", row.names = FALSE)

cat("\n===== S9 写入完成 =====\n")
print(out, row.names = FALSE)
cat("\nC-index 范围（外部队列）:",
    paste(range(out$C_index[out$cohort != "TCGA"]), collapse = " ~ "), "\n")
