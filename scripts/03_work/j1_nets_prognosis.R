# ==============================================================================
# 判定 1：NETs 生物学本身有没有预后价值？
#   这决定「换基因集能不能提升判别力」这条路值不值得走
#
# 设计：
#   (a) NETs 参考评分（中性粒 + NETosis 执行集）的直接预后能力 —— 四队列
#   (b) 在 TCGA 中取「与 NETs 评分最相关」的 top 50/100/200 基因，
#       构成绩效评分，再到 METABRIC / GSE96058 外部验证
#       —— 这一步直接回答用户想法：换选法能否突破天花板
# 输出：j1_results.csv
# ==============================================================================
suppressMessages({
  library(data.table); library(survival)
  library(SummarizedExperiment); library(org.Hs.eg.db)
})

nets_core <- c("PADI4","ELANE","MPO","GSDMD","HMGB1","NCF1","CYBB","AZU1","CTSG",
               "DNASE1L3","FCGR3B","CSF3R","CXCR2","S100A8","S100A9","FPR1","MMP9",
               "CAMP","LCN2","ITGAM","ITGB2","TLR4","CAMP")

cidx <- function(time, event, score) {
  ok <- is.finite(time) & is.finite(event) & is.finite(score) & time > 0
  if (sum(ok) < 50 || sum(event[ok]) < 10) return(NA)
  as.numeric(concordance(Surv(time[ok], event[ok]) ~ score[ok], reverse = TRUE)$concordance)
}

hr_sd <- function(d, sc) {
  d$s <- as.numeric(scale(d[[sc]]))
  f <- coxph(Surv(OS.time, OS) ~ s, data = d)
  s <- summary(f)
  c(HR = s$conf.int[1], lo = s$conf.int[3], hi = s$conf.int[4], p = s$coefficients[1,5])
}

km_p <- function(d, sc) {
  g <- ifelse(d[[sc]] > median(d[[sc]], na.rm = TRUE), 1, 0)
  sd <- survdiff(Surv(OS.time, OS) ~ g, data = d)
  pchisq(sd$chisq, length(sd$n) - 1, lower.tail = FALSE)
}

zscore <- function(mat) { z <- t(scale(t(mat))); z[is.na(z)] <- 0; z }

score_of <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 5) return(NULL)
  colMeans(zscore(mat[g, , drop = FALSE]))
}

# ==============================================================================
# 1. TCGA：算 NETs 评分，选 top-N 共表达基因
# ==============================================================================
cat("========== TCGA ==========\n")
load("../raw/TCGA/brca.gdc_2024.rda")
expr <- assay(expquery2)
colnames(expr) <- substr(gsub("\\.", "-", colnames(expr)), 1, 12)
ens <- gsub("\\..*", "", rownames(expr))
sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(expr) <- sym
expr <- expr[!is.na(rownames(expr)), ]
expr <- expr[!duplicated(rownames(expr)), ]
expr <- as.matrix(expr)
cat("  基因数:", nrow(expr), "\n")

fd <- read.csv("Table_S7_Final_Clinical_Data.csv", stringsAsFactors = FALSE)
s <- intersect(colnames(expr), fd$Patient_ID)
cat("  匹配样本:", length(s), "\n")
E <- expr[, s, drop = FALSE]
rm(expr, expquery2); gc()

nets_sc <- score_of(E, nets_core)
cat("  NETs 可用基因:", length(intersect(nets_core, rownames(E))), "/", length(unique(nets_core)), "\n")

d_tcga <- data.frame(Patient_ID = s, NETs = as.numeric(nets_sc), stringsAsFactors = FALSE)
d_tcga <- merge(d_tcga, fd[, c("Patient_ID","OS","OS.time","RiskScore")], by = "Patient_ID")
d_tcga <- d_tcga[!is.na(d_tcga$OS.time) & d_tcga$OS.time > 0, ]
cat("  分析集 n =", nrow(d_tcga), " 事件 =", sum(d_tcga$OS), "\n\n")

# 与 NETs 评分的相关性（全基因）
keep_g <- rownames(E)[rowMeans(E[, s], na.rm = TRUE) > 1]
cat("  参与筛选的基因:", length(keep_g), "\n")
rhos <- apply(E[keep_g, s], 1, function(x)
  suppressWarnings(cor(x, nets_sc, method = "spearman")))
rhos <- rhos[is.finite(rhos)]
ord <- order(rhos, decreasing = TRUE)
cat(sprintf("  相关分布: 最高 %.3f, 中位 %.3f, 最低 %.3f\n",
            rhos[ord[1]], median(rhos), rhos[ord[length(ord)]]))
cat("  NETs core 基因自身排名（前 5）:\n")
print(round(head(sort(rhos[nets_core], decreasing = TRUE), 5), 3))

topN <- list()
for (N in c(50, 100, 200)) topN[[as.character(N)]] <- names(rhos)[ord[1:N]]
saveRDS(topN, "j1_topN_genes.rds")
saveRDS(d_tcga, "j1_tcga_clin.rds")

# ==============================================================================
# 2. 各队列：NETs 评分 + top-N 评分
# ==============================================================================
res <- list()
add <- function(cohort, var, d, sc) {
  ci <- cidx(d$OS.time, d$OS, d[[sc]])
  h  <- hr_sd(d, sc)
  p  <- km_p(d, sc)
  cat(sprintf("  %-14s C=%5.3f  HR/SD=%.2f (%.2f-%.2f) p=%.3g  KM p=%.2e\n",
              var, ci, h[1], h[2], h[3], h[4], p))
  res[[length(res)+1]] <<- data.frame(cohort = cohort, variable = var,
    C_index = ci, HR_perSD = h[1], HR_lo = h[2], HR_hi = h[3],
    HR_p = h[4], KM_logrank_p = p, n = nrow(d), events = sum(d$OS),
    stringsAsFactors = FALSE)
}

cat("\n========== 各队列判定 ==========\n")
cat("---- TCGA (训练) ----\n")
d_tcga$NETs <- as.numeric(d_tcga$NETs)
add("TCGA", "NETs_score", d_tcga, "NETs")
add("TCGA", "risk_4gene", d_tcga, "RiskScore")

for (N in c(50, 100, 200)) {
  sc <- score_of(E, topN[[as.character(N)]])
  d_tcga[[paste0("top", N)]] <- as.numeric(sc)
  add("TCGA", paste0("top", N, "_coexp"), d_tcga, paste0("top", N))
}
need <- unique(c(nets_core, unlist(topN)))

# ---- METABRIC ----
cat("\n---- METABRIC ----\n")
mb <- fread("../raw/METABRIC/data_mrna_illumina_microarray.txt", data.table = FALSE, showProgress = FALSE)
mb <- mb[mb$Hugo_Symbol %in% need, ]
mb <- mb[!duplicated(mb$Hugo_Symbol), ]
rownames(mb) <- mb$Hugo_Symbol
mb <- as.matrix(mb[, -c(1,2), drop = FALSE])
colnames(mb) <- gsub("\\.", "-", colnames(mb))
dmb <- read.csv("Table_S6_METABRIC_Validation_Result.csv", stringsAsFactors = FALSE)
dmb$PID <- gsub("\\.", "-", dmb$PATIENT_ID)
cm <- intersect(colnames(mb), dmb$PID)
dmb <- dmb[match(cm, dmb$PID), ]
dmb$OS.time <- as.numeric(dmb$OS.time); dmb$OS <- as.numeric(dmb$OS)
dmb$NETs <- as.numeric(score_of(mb[, cm, drop = FALSE], nets_core))
for (N in c(50,100,200)) dmb[[paste0("top",N)]] <- as.numeric(score_of(mb[, cm, drop=FALSE], topN[[as.character(N)]]))
add("METABRIC", "NETs_score", dmb, "NETs")
add("METABRIC", "risk_4gene", dmb, "RiskScore")
for (N in c(50,100,200)) add("METABRIC", paste0("top",N,"_coexp"), dmb, paste0("top",N))
rm(mb); gc()

# ---- GSE96058 ----
cat("\n---- GSE96058 ----\n")
s96 <- fread("../raw/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",
             data.table = FALSE, showProgress = FALSE)
rn <- s96[[1]]; s96 <- s96[,-1, drop = FALSE]; rownames(s96) <- rn
s96 <- as.matrix(s96[rownames(s96) %in% need, , drop = FALSE])
rm(rn); gc()
d96 <- read.csv("Table_S6_GSE96058_Validation_Result.csv", stringsAsFactors = FALSE)
c9 <- intersect(colnames(s96), d96$title)
d96 <- d96[match(c9, d96$title), ]
d96$NETs <- as.numeric(score_of(s96[, c9, drop=FALSE], nets_core))
for (N in c(50,100,200)) d96[[paste0("top",N)]] <- as.numeric(score_of(s96[, c9, drop=FALSE], topN[[as.character(N)]]))
add("GSE96058", "NETs_score", d96, "NETs")
add("GSE96058", "risk_4gene", d96, "RiskScore")
for (N in c(50,100,200)) add("GSE96058", paste0("top",N,"_coexp"), d96, paste0("top",N))
rm(s96); gc()

# ---- GSE20685 ----
cat("\n---- GSE20685 ----\n")
suppressMessages(library(GEOquery))
gset <- getGEO(filename = "../raw/GSE20685/GSE20685_series_matrix.txt.gz", getGPL = FALSE)
e20 <- exprs(gset)
gpl <- read.delim("../raw/annotation/GPL570-55999.txt", sep = "\t", comment.char = "#",
                  stringsAsFactors = FALSE, check.names = FALSE)
p2s <- data.frame(ID = gpl$ID,
                  symbol = sapply(strsplit(as.character(gpl$`Gene Symbol`), " /// "), "[", 1),
                  stringsAsFactors = FALSE)
p2s <- p2s[!is.na(p2s$symbol) & p2s$symbol != "" & p2s$symbol %in% need, ]
idx <- match(rownames(e20), p2s$ID); ridx <- which(!is.na(idx))
e20 <- e20[ridx, , drop = FALSE]; sy <- p2s$symbol[idx[ridx]]
o <- order(sy, -rowMeans(e20, na.rm = TRUE)); e20 <- e20[o,]; sy <- sy[o]
kp <- !duplicated(sy); e20 <- e20[kp, , drop = FALSE]; rownames(e20) <- sy[kp]
d20 <- read.csv("Table_S6_GEO_Validation_Result.csv", stringsAsFactors = FALSE)
c2 <- intersect(colnames(e20), d20$geo_accession)
d20 <- d20[match(c2, d20$geo_accession), ]
d20$NETs <- as.numeric(score_of(e20[, c2, drop=FALSE], nets_core))
for (N in c(50,100,200)) d20[[paste0("top",N)]] <- as.numeric(score_of(e20[, c2, drop=FALSE], topN[[as.character(N)]]))
add("GSE20685", "NETs_score", d20, "NETs")
add("GSE20685", "risk_4gene", d20, "RiskScore")
for (N in c(50,100,200)) add("GSE20685", paste0("top",N,"_coexp"), d20, paste0("top",N))

R <- do.call(rbind, res)
write.csv(R, "j1_results.csv", row.names = FALSE)

# ==============================================================================
# 3. 汇总
# ==============================================================================
cat("\n\n================ 判定 1 汇总：C-index ================\n")
vars <- c("NETs_score","top50_coexp","top100_coexp","top200_coexp","risk_4gene")
cat(sprintf("%-16s %8s %10s %10s %10s\n", "变量","TCGA","GSE20685","METABRIC","GSE96058"))
for (v in vars) {
  row <- sapply(c("TCGA","GSE20685","METABRIC","GSE96058"), function(cc) {
    z <- R[R$cohort == cc & R$variable == v, "C_index"]; if (length(z)) round(z,3) else NA })
  cat(sprintf("%-16s %8s %10s %10s %10s\n", v, row[1], row[2], row[3], row[4]))
}
cat("\n判读标准：若 top-N 四项外部 C-index 仍全部 < 0.58，则换基因集无法突破天花板。\n")
