# ==============================================================================
# 免疫轴外部复现 —— 第一步：四队列统一模块评分构建
# 关键设计：B 细胞模块分「含 CD79A」与「不含 CD79A」两版
#   因为 risk score 本身含 CD79A（系数 -0.0854），含 CD79A 的版本会构造性相关
# 输出：ia_scores_<cohort>.rds（每队列一个，含 RiskScore + 各模块评分 + 关键基因）
# ==============================================================================
suppressMessages({
  library(data.table)
  library(survival)
})

# ---- 统一基因集（刻意把 CD79A 单列，用于拆解循环性）----
modules <- list(
  Bcell_noCD79A = c("MS4A1","CD19","CD79B","TNFRSF17","MZB1","JCHAIN",
                    "IGKC","IGHM","IGHD","PAX5","BANK1","FCRL1"),
  Bcell_withCD79A = c("CD79A","MS4A1","CD19","CD79B","TNFRSF17","MZB1",
                      "JCHAIN","IGKC","IGHM","IGHD","PAX5","BANK1","FCRL1"),
  CD8_T      = c("CD8A","CD8B","GZMK","GZMA","PRF1","NKG7","CXCL9","CXCL10","CD3D","CD3E"),
  Treg       = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF4","CCR8","IL10"),
  Cytotoxic  = c("GZMB","GNLY","KLRD1","PRF1","NKG7","GZMH"),
  Macrophage = c("CD68","CD163","MRC1","CSF1R","MSR1","C1QA","C1QB","CD14"),
  Mast       = c("TPSAB1","CPA3","MS4A2","KIT"),
  Neutrophil = c("FCGR3B","CSF3R","CXCR2","S100A8","S100A9","FPR1"),
  Stroma     = c("COL1A1","COL1A2","COL3A1","FAP","ACTA2","PDGFRB","LUM","TAGLN"),
  Prolif     = c("MKI67","AURKA","BIRC5","CCNB1","MYBL2","ANLN","UBE2C","TOP2A","CCNB2","TYMS"),
  Purity_epi = c("EPCAM","KRT8","KRT18","KRT19","CDH1","KRT7","ELF3"),
  Hub        = c("PLG","CD79A","MMP7","ANGPT1")
)
allg <- unique(unlist(modules))

# 模块评分：基因队列内 z-score 后对样本取均值（行=基因，列=样本）
module_score <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 3) return(NULL)
  m <- mat[g, , drop = FALSE]
  z <- t(scale(t(m)))            # 每个基因跨样本 z
  z[is.na(z)] <- 0
  sc <- colMeans(z)              # 每个样本的模块均分
  attr(sc, "n") <- length(g)
  sc
}

build_cohort <- function(mat, risk_df, id_col, cohort) {
  # mat: 基因×样本（已 log2 尺度）；risk_df 含 id_col 与 RiskScore
  ids <- risk_df[[id_col]]
  common <- intersect(colnames(mat), ids)
  cat(sprintf("[%s] 表达样本 %d，评分表 %d，交集 %d\n",
              cohort, ncol(mat), nrow(risk_df), length(common)))
  if (length(common) < 50) { cat("  样本过少，跳过\n"); return(NULL) }
  mat <- mat[, common, drop = FALSE]
  out <- data.frame(sample = common,
                    RiskScore = risk_df$RiskScore[match(common, ids)],
                    stringsAsFactors = FALSE)
  for (nm in names(modules)) {
    sc <- module_score(mat, modules[[nm]])
    if (is.null(sc)) { cat(sprintf("  %-16s 可用基因不足，跳过\n", nm)); next }
    out[[nm]] <- as.numeric(sc)
    cat(sprintf("  %-16s 使用 %d/%d 个基因\n", nm, attr(sc, "n"), length(modules[[nm]])))
  }
  # 单基因表达也留着，便于后续拆解
  for (g in c("CD79A","PLG","MMP7","ANGPT1","PTPRC","CD8A","FOXP3")) {
    if (g %in% rownames(mat)) out[[paste0("g_", g)]] <- as.numeric(mat[g, common])
  }
  saveRDS(out, sprintf("ia_scores_%s.rds", cohort))
  cat(sprintf("  -> ia_scores_%s.rds  (%d 样本)\n\n", cohort, nrow(out)))
  invisible(out)
}

# ==============================================================================
# 1. TCGA-BRCA（训练队列）
# ==============================================================================
cat("========== 1. TCGA-BRCA ==========\n")
suppressMessages(library(SummarizedExperiment)); suppressMessages(library(org.Hs.eg.db))
load("../raw/TCGA/brca.gdc_2024.rda")
expr <- assay(expquery2)
colnames(expr) <- substr(gsub("\\.", "-", colnames(expr)), 1, 12)
ens <- gsub("\\..*", "", rownames(expr))
sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
rownames(expr) <- sym
expr <- expr[!is.na(rownames(expr)), ]
expr <- expr[!duplicated(rownames(expr)), ]
expr <- as.matrix(expr)
# 只保留需要的基因，释放内存
expr <- expr[intersect(allg, rownames(expr)), , drop = FALSE]
cat("  TCGA 提取基因数:", nrow(expr), "\n")
tcga_fd <- read.csv("Table_S7_Final_Clinical_Data.csv", stringsAsFactors = FALSE)
build_cohort(expr, tcga_fd, "Patient_ID", "TCGA")
rm(expr, expquery2); gc()

# ==============================================================================
# 2. GSE20685（芯片，GPL570）
# ==============================================================================
cat("========== 2. GSE20685 ==========\n")
suppressMessages(library(GEOquery))
gset <- getGEO(filename = "../raw/GSE20685/GSE20685_series_matrix.txt.gz", getGPL = FALSE)
exp_raw <- exprs(gset)
gpl <- read.delim("../raw/annotation/GPL570-55999.txt", sep = "\t", comment.char = "#",
                  stringsAsFactors = FALSE, check.names = FALSE)
p2s <- data.frame(ID = gpl$ID,
                  symbol = sapply(strsplit(as.character(gpl$`Gene Symbol`), " /// "), "[", 1),
                  stringsAsFactors = FALSE)
# 注意：symbol 含 NA，不能用 `%in% & != ""` 组合筛选（会产生 NA 行索引）
p2s <- p2s[!is.na(p2s$symbol) & p2s$symbol != "" & p2s$symbol %in% allg, ]
idx <- match(rownames(exp_raw), p2s$ID)   # 长度 = 探针数，值为 p2s 行号
ridx <- which(!is.na(idx))
exp20685 <- exp_raw[ridx, , drop = FALSE]
sym_sub  <- p2s$symbol[idx[ridx]]          # 必须用 idx[ridx] 取值，不能直接掩码 p2s
# 同一基因多探针：保留平均表达最高的探针（避免 avereps 的 NA 分组问题）
ord <- order(sym_sub, -rowMeans(exp20685, na.rm = TRUE))
exp20685 <- exp20685[ord, , drop = FALSE]
sym_sub  <- sym_sub[ord]
keep <- !duplicated(sym_sub)
exp20685 <- exp20685[keep, , drop = FALSE]
rownames(exp20685) <- sym_sub[keep]
cat("  GSE20685 基因数:", nrow(exp20685), "\n")
d20685 <- read.csv("Table_S6_GEO_Validation_Result.csv", stringsAsFactors = FALSE)
build_cohort(exp20685, d20685, "geo_accession", "GSE20685")
rm(exp_raw, gset, exp20685); gc()

# ==============================================================================
# 3. METABRIC（芯片，Illumina HT-12 v3）
# ==============================================================================
cat("========== 3. METABRIC ==========\n")
mb <- fread("../raw/METABRIC/data_mrna_illumina_microarray.txt",
            data.table = FALSE, showProgress = FALSE)
cat("  METABRIC 原始维度:", dim(mb), "\n")
mb2 <- mb[mb$Hugo_Symbol %in% allg, ]
mb2 <- mb2[!duplicated(mb2$Hugo_Symbol), ]
rownames(mb2) <- mb2$Hugo_Symbol
mb2 <- as.matrix(mb2[, -c(1, 2), drop = FALSE])
rm(mb); gc()
cat("  METABRIC 基因数:", nrow(mb2), "\n")
dmb <- read.csv("Table_S6_METABRIC_Validation_Result.csv", stringsAsFactors = FALSE)
dmb$PID2 <- gsub("\\.", "-", dmb$PATIENT_ID)
colnames(mb2) <- gsub("\\.", "-", colnames(mb2))
build_cohort(mb2, dmb, "PID2", "METABRIC")
rm(mb2); gc()

# ==============================================================================
# 4. GSE96058（RNA-seq，3273 样本）
# ==============================================================================
cat("========== 4. GSE96058 ==========\n")
s96 <- fread("../raw/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",
             data.table = FALSE, showProgress = FALSE)
cat("  GSE96058 原始维度:", dim(s96), "\n")
rn <- s96[[1]]
s96 <- s96[, -1, drop = FALSE]
rownames(s96) <- rn
s96 <- s96[rownames(s96) %in% allg, , drop = FALSE]
s96 <- as.matrix(s96)
rm(rn); gc()
cat("  GSE96058 基因数:", nrow(s96), "\n")
d96 <- read.csv("Table_S6_GSE96058_Validation_Result.csv", stringsAsFactors = FALSE)
build_cohort(s96, d96, "title", "GSE96058")
rm(s96); gc()

cat("\n== 全部队列评分构建完成 ==\n")
