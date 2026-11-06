# ==============================================================================
# 诊断：GSE96058 为何独独不复现免疫轴？
# 排查方向：模块内部一致性（数据质量） vs 权重实际构成（评分本身）
# ==============================================================================
cohorts <- c("TCGA","GSE20685","METABRIC","GSE96058")
D <- list(); for (c1 in cohorts) D[[c1]] <- readRDS(sprintf("ia_scores_%s.rds", c1))

cat("================ 一、评分实际由哪个基因驱动 ================\n")
cat("（risk 与各 hub 基因表达的 Spearman，四队列对比）\n\n")
cat(sprintf("%-10s %10s %10s %10s %10s   %s\n", "队列","PLG","CD79A","MMP7","ANGPT1","主导基因"))
for (c1 in cohorts) {
  d <- D[[c1]]
  rr <- sapply(c("g_PLG","g_CD79A","g_MMP7","g_ANGPT1"), function(g)
    suppressWarnings(cor(d$RiskScore, d[[g]], method="spearman")))
  cat(sprintf("%-10s %10.3f %10.3f %10.3f %10.3f   %s\n",
              c1, rr[1], rr[2], rr[3], rr[4], names(which.max(abs(rr)))))
}

cat("\n================ 二、hub 基因间的相互相关（Spearman）================\n")
for (c1 in cohorts) {
  d <- D[[c1]]
  m <- cor(d[, c("g_PLG","g_CD79A","g_MMP7","g_ANGPT1")], method="spearman")
  cat(sprintf("---- %s ----\n", c1))
  print(round(m, 2))
  cat(sprintf("  平均 |r| = %.3f\n\n", mean(abs(m[upper.tri(m)]))))
}

cat("\n================ 三、B 细胞模块内部一致性 ================\n")
cat("（若某队列内部一致性低，说明是平台/数据问题，不是生物学差异）\n\n")
# 用原始基因表达（rds 里只存了 CD79A/PTPRC/CD8A/FOXP3 单基因，需回读矩阵）
# 这里改用间接指标：模块评分 vs 单基因表达的一致性
cat(sprintf("%-10s %14s %14s %14s\n","队列","B模块 vs CD8A","CD8_T vs CD8A","CD8_T vs FOXP3"))
for (c1 in cohorts) {
  d <- D[[c1]]
  f <- function(a,b) if (a %in% names(d) && b %in% names(d))
    round(suppressWarnings(cor(d[[a]], d[[b]], method="spearman")),3) else NA
  cat(sprintf("%-10s %14s %14s %14s\n", c1,
              f("Bcell_noCD79A","g_CD8A"), f("CD8_T","g_CD8A"), f("CD8_T","g_FOXP3")))
}

cat("\n================ 四、risk 与免疫相关的方向，按 hub 基因分层检验 ================\n")
cat("若高 risk 组的免疫低下完全由 CD79A 造成，则剔除 CD79A 后应消失\n\n")
for (c1 in cohorts) {
  d <- D[[c1]]
  d$grp <- ifelse(d$RiskScore > median(d$RiskScore, na.rm=TRUE), "High","Low")
  # 只在 CD79A 高表达（B 细胞存在）的样本内比较
  q <- quantile(d$g_CD79A, 0.5, na.rm=TRUE)
  sub <- d[d$g_CD79A >= q, ]
  a <- median(sub$CD8_T[sub$grp=="High"], na.rm=TRUE)
  b <- median(sub$CD8_T[sub$grp=="Low"], na.rm=TRUE)
  p <- suppressWarnings(wilcox.test(sub$CD8_T[sub$grp=="High"],
                                    sub$CD8_T[sub$grp=="Low"]))$p.value
  cat(sprintf("%-10s CD79A 上半区 (n=%d): CD8_T  High %.3f vs Low %.3f  Δ=%.3f  p=%.2g\n",
              c1, nrow(sub), a, b, a-b, p))
}

cat("\n================ 五、GSE96058 表达值分布（判断是否可比）================\n")
for (c1 in cohorts) {
  d <- D[[c1]]
  cat(sprintf("%-10s PTPRC %.2f  CD8A %.2f  CD79A %.2f  FOXP3 %.2f  EPCAM %.2f (均值)\n",
              c1, mean(d$g_PTPRC), mean(d$g_CD8A), mean(d$g_CD79A),
              mean(d$g_FOXP3), mean(d$Purity_epi)))
}
