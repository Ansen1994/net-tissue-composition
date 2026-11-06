# ==============================================================================
# 免疫轴信号溯源：把 risk score 拆成「CD79A 部分」与「其余三基因部分」
# 目的：免疫轴到底是独立的生物学关联，还是评分内源变量的算术必然
# ==============================================================================
cohorts <- c("TCGA","GSE20685","METABRIC","GSE96058")
D <- list(); for (c1 in cohorts) D[[c1]] <- readRDS(sprintf("ia_scores_%s.rds", c1))

mods <- c("Bcell_noCD79A","CD8_T","Treg","Cytotoxic","Macrophage","Neutrophil","Stroma")
out <- list()

cat("================ 评分拆解（队列内 z 化后套原系数）================\n\n")
for (c1 in cohorts) {
  d <- D[[c1]]
  z <- scale(d[, c("g_PLG","g_CD79A","g_MMP7","g_ANGPT1")])
  d$risk_full    <- 0.1102*z[,1] + 0.2295*z[,4] - 0.0854*z[,2] - 0.0851*z[,3]
  d$risk_noCD79A <- 0.1102*z[,1] + 0.2295*z[,4] - 0.0851*z[,3]
  d$risk_CD79A   <- -0.0854*z[,2]
  # 核对：拆解后与原 RiskScore 的秩相关应接近 1
  cat(sprintf("%-10s cor(risk_full, 原RiskScore) = %.3f\n", c1,
              cor(d$risk_full, d$RiskScore, method="spearman")))

  cat(sprintf("  %-14s %18s %18s %14s\n", "免疫模块", "全评分", "去CD79A", "仅CD79A"))
  for (m in mods) {
    if (!m %in% names(d)) next
    f <- function(x) round(suppressWarnings(cor(x, d[[m]], method="spearman")), 3)
    cat(sprintf("  %-14s %18.3f %18.3f %14.3f\n", m,
                f(d$risk_full), f(d$risk_noCD79A), f(d$risk_CD79A)))
    out[[length(out)+1]] <- data.frame(cohort = c1, module = m,
      rho_full = f(d$risk_full), rho_noCD79A = f(d$risk_noCD79A),
      rho_CD79Aonly = f(d$risk_CD79A), stringsAsFactors = FALSE)
  }
  cat("\n")
}
O <- do.call(rbind, out)
write.csv(O, "Table_ImmuneAxis_ScoreDecomposition.csv", row.names = FALSE)

cat("\n================ 汇总：去 CD79A 后免疫关联保留比例 ================\n")
for (m in mods) {
  s <- O[O$module == m, ]
  if (nrow(s) < 4) next
  cat(sprintf("%-14s 全评分均值 %+.3f | 去CD79A均值 %+.3f | 仅CD79A均值 %+.3f | 负相关队列数(全/去CD79A) %d/%d\n",
              m, mean(s$rho_full), mean(s$rho_noCD79A), mean(s$rho_CD79Aonly),
              sum(s$rho_full < 0), sum(s$rho_noCD79A < 0)))
}
