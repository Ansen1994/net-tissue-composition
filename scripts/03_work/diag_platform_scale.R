# ============================================================
# 诊断：四个队列中 4 个 hub 基因的绝对水平与离散度是否可比
# 目的：决定 risk score 的跨队列标准化方案
# ============================================================
suppressMessages({
  library(data.table); library(GEOquery); library(dplyr); library(tibble)
})
hub <- c("PLG","CD79A","MMP7","ANGPT1")

# ---------- TCGA 训练集 ----------
rt <- read.csv("Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv", row.names = 1, check.names = FALSE)
cat("===== TCGA 训练集 (n =", nrow(rt), ") =====\n")
for (g in hub) cat(sprintf("  %-7s mean=%7.3f  sd=%6.3f  median=%7.3f\n",
                           g, mean(rt[[g]]), sd(rt[[g]]), median(rt[[g]])))

stat <- function(x, lab) {
  cat("=====", lab, "(n =", nrow(x), ") =====\n")
  for (g in hub) {
    v <- x[[g]]
    cat(sprintf("  %-7s mean=%7.3f  sd=%6.3f  median=%7.3f  range=[%.2f, %.2f]\n",
                g, mean(v, na.rm=TRUE), sd(v, na.rm=TRUE), median(v, na.rm=TRUE),
                min(v, na.rm=TRUE), max(v, na.rm=TRUE)))
  }
}

# ---------- GSE20685 ----------
g <- getGEO(filename="../raw/GSE20685/GSE20685_series_matrix.txt.gz", getGPL=FALSE)
e <- exprs(g); p <- pData(g)
# 只取 ER+ 子集
er_col <- grep("^er_status", names(p), value=TRUE, ignore.case=TRUE)
if (length(er_col)) {
  keep <- p[[er_col[1]]] %in% c("positive","Positive","1","ER+")
  cat("\n[GSE20685] ER 字段:", er_col[1], "| ER+ n =", sum(keep, na.rm=TRUE), "/", nrow(p), "\n")
  e <- e[, keep]
}
x <- as.data.frame(t(e[intersect(hub, rownames(e)), , drop=FALSE]))
stat(x, "GSE20685 ER+")

# ---------- METABRIC ----------
m <- fread("../raw/METABRIC/data_mrna_illumina_microarray.txt", data.table=FALSE)
m <- m[m$Hugo_Symbol %in% hub, ]
mm <- m %>% dplyr::select(-2) %>% dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(dplyr::across(dplyr::where(is.numeric), max)) %>%
  tibble::column_to_rownames("Hugo_Symbol")
x2 <- as.data.frame(t(as.matrix(mm)))
cat("\n[METABRIC] 全队列 n =", nrow(x2), "（此处仅看平台尺度，未过滤 ER）\n")
for (g2 in intersect(hub, colnames(x2))) {
  v <- x2[[g2]]
  cat(sprintf("  %-7s mean=%7.3f  sd=%6.3f  median=%7.3f  range=[%.2f, %.2f]\n",
              g2, mean(v, na.rm=TRUE), sd(v, na.rm=TRUE), median(v, na.rm=TRUE),
              min(v, na.rm=TRUE), max(v, na.rm=TRUE)))
}

# ---------- GSE96058 ----------
g3 <- fread("../raw/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",
            data.table=FALSE)
rownames(g3) <- g3[,1]; gm <- as.matrix(g3[,-1])
x3 <- as.data.frame(t(gm[intersect(hub, rownames(gm)), , drop=FALSE]))
cat("\n[GSE96058] 全队列 n =", nrow(x3), "（仅看平台尺度，未过滤 ER）\n")
for (g2 in intersect(hub, colnames(x3))) {
  v <- x3[[g2]]
  cat(sprintf("  %-7s mean=%7.3f  sd=%6.3f  median=%7.3f  range=[%.2f, %.2f]\n",
              g2, mean(v, na.rm=TRUE), sd(v, na.rm=TRUE), median(v, na.rm=TRUE),
              min(v, na.rm=TRUE), max(v, na.rm=TRUE)))
}

cat("\n===== 结论指标：各队列 SD / 训练集 SD =====\n")
