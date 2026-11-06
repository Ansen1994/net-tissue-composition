## 03_immune.R
## 目的：免疫治疗在 GDSC 细胞系筛选中无对应药物（ICI 均为抗体，不在小分子库内），
##       改用公认的免疫微环境代理指标评估免疫治疗潜在获益
## 同时输出内分泌/CDK4-6/PI3K 通路靶点表达，供治疗选择判断

suppressPackageStartupMessages({
  library(GSVA)
  library(estimate)
})

OUT <- Sys.getenv("WB_OUTPUT_DIR", unset = "out")
setwd(OUT)

tcga <- readRDS("tcga_log2tpm.rds")
expr <- tcga$expr; clin <- tcga$clin
cat("维度:", dim(expr), "\n")

## ============ 1. 基因集 ============
sets <- list(
  `T-cell inflamed GEP (18-gene)` = c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG",
                                     "CD274","CCL5","CXCL11","CXCR6","LAG3","NKG7",
                                     "CD8A","CMKLR1","PSMB10","HLA-E","CD276","TIGIT"),
  `Cytolytic activity (GZMA/PRF1)` = c("GZMA","PRF1"),
  `T cell infiltration` = c("CD3D","CD3E","CD3G","CD2","CD8A","CD8B"),
  `IFN-gamma response` = c("IFNG","STAT1","IRF1","CXCL9","CXCL10","GBP1","GBP2","IDO1","HLA-DRA"),
  `Checkpoint module` = c("PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT","IDO1","BTLA","VSIR"),
  `MHC class II` = c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","CIITA")
)
sets <- lapply(sets, function(g) intersect(g, rownames(expr)))
cat("基因集命中数:\n"); print(lengths(sets))

ss <- gsva(ssgseaParam(as.matrix(expr), sets), verbose = FALSE)
ss <- t(ss)
cat("ssGSEA 结果:", dim(ss), "\n")

## ============ 2. 检查点单基因 ============
ck <- c("CD274","PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","IDO1","CD8A")
ck <- intersect(ck, rownames(expr))
ck_mat <- t(expr[ck, , drop = FALSE])

## ============ 3. 治疗靶点表达 ============
tg <- c("ESR1","PGR","AR","GATA3","FOXA1","ERBB2","EGFR","MKI67","CCND1",
        "CDK4","CDK6","RB1","PIK3CA","PTEN","AKT1","MTOR","BRCA1","BRCA2","MMP7","PLG","ANGPT1","CD79A")
tg <- intersect(tg, rownames(expr))
tg_mat <- t(expr[tg, , drop = FALSE])

## ============ 4. ESTIMATE（基于 log2 TPM，修正原代码用 counts 的问题）============
cat("== 运行 ESTIMATE ==\n")
gct <- "out/tcga_log2tpm.gct"
con <- file(gct, "w")
writeLines("#1.2", con)
writeLines(paste(nrow(expr), ncol(expr), sep = "\t"), con)
writeLines(paste(c("Name", "Description", colnames(expr)), collapse = "\t"), con)
write.table(data.frame(Name = rownames(expr), Description = "na", expr,
                       check.names = FALSE),
            con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
close(con)

est_ok <- tryCatch({
  estimateScore(gct, "out/tcga_estimate.gct", platform = "illumina")
  TRUE
}, error = function(e) { cat("ESTIMATE 失败:", conditionMessage(e), "\n"); FALSE })

est <- NULL
if (est_ok) {
  s <- read.table("out/tcga_estimate.gct", skip = 2, header = TRUE, sep = "\t",
                  check.names = FALSE)
  rownames(s) <- s$NAME; s$NAME <- NULL; s$Description <- NULL
  est <- as.data.frame(t(s))
  for (cc in intersect(c("ImmuneScore","StromalScore","ESTIMATEScore","TumorPurity"), colnames(est)))
    est[[cc]] <- as.numeric(est[[cc]])
  if (!"TumorPurity" %in% colnames(est) && all(c("ESTIMATEScore") %in% colnames(est)))
    est$TumorPurity <- cos(0.6049872018 + 0.0001467884 * est$ESTIMATEScore)
  cat("ESTIMATE 完成:", dim(est), "\n")
  print(summary(est[, intersect(c("ImmuneScore","StromalScore","TumorPurity"), colnames(est))]))
}

## ============ 5. 分组比较 ============
allv <- cbind(as.data.frame(ss), ck_mat, est)
allv$Patient_ID <- rownames(allv)
df <- merge(allv, clin[, c("Patient_ID","RiskScore","RiskGroup","Stage","Age")], by = "Patient_ID")
cat("合并样本:", nrow(df), "\n")

vars <- setdiff(colnames(allv), "Patient_ID")
out <- do.call(rbind, lapply(vars, function(v) {
  x <- df[[v]]; g <- df$RiskGroup
  hi <- x[g == "High"]; lo <- x[g == "Low"]
  w <- wilcox.test(hi, lo)
  r <- suppressWarnings(cor.test(x, df$RiskScore, method = "spearman"))
  data.frame(variable = v,
             mean_Low = round(mean(lo, na.rm = TRUE), 3),
             mean_High = round(mean(hi, na.rm = TRUE), 3),
             diff = round(mean(hi, na.rm = TRUE) - mean(lo, na.rm = TRUE), 3),
             p = w$p.value, rho_risk = round(r$estimate, 3),
             p_rho = r$p.value, stringsAsFactors = FALSE)
}))
out$fdr <- p.adjust(out$p, "BH")
out$group <- ifelse(out$variable %in% names(sets), "ssGSEA signature",
             ifelse(out$variable %in% ck, "Checkpoint gene", "ESTIMATE"))
out <- out[order(out$p), ]
write.csv(out, "out/immune_proxy.csv", row.names = FALSE)
print(out[, c("variable","group","diff","p","fdr","rho_risk")], row.names = FALSE)

## 治疗靶点表达
tg_df <- merge(data.frame(tg_mat, Patient_ID = rownames(tg_mat), check.names = FALSE),
               clin[, c("Patient_ID","RiskScore","RiskGroup")], by = "Patient_ID")
tgo <- do.call(rbind, lapply(colnames(tg_mat), function(v) {
  x <- tg_df[[v]]; g <- tg_df$RiskGroup
  w <- wilcox.test(x[g == "High"], x[g == "Low"])
  data.frame(gene = v, mean_Low = round(mean(x[g=="Low"]),3), mean_High = round(mean(x[g=="High"]),3),
             p = w$p.value, stringsAsFactors = FALSE)
}))
tgo$fdr <- p.adjust(tgo$p, "BH")
tgo <- tgo[order(tgo$p), ]
write.csv(tgo, "out/target_expression.csv", row.names = FALSE)
cat("\n== 治疗靶点表达 ==\n"); print(tgo, row.names = FALSE)

saveRDS(list(df = df, sets = names(sets)), "out/immune_df.rds")
cat("== 完成 ==\n")
