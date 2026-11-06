# ==============================================================================
# 免疫轴外部复现 —— 第二步：统计检验与循环性拆解
# 核心问题：risk 与免疫的相关，有多少是真实的跨队列信号，
#           有多少是「评分里含 CD79A」这一构造事实的同义反复
# ==============================================================================

modules_all <- c("Bcell_noCD79A","Bcell_withCD79A","CD8_T","Treg","Cytotoxic",
                 "Macrophage","Mast","Neutrophil","Stroma","Prolif","Purity_epi")

cohorts <- c("TCGA","GSE20685","METABRIC","GSE96058")
D <- list()
for (c1 in cohorts) {
  f <- sprintf("ia_scores_%s.rds", c1)
  if (file.exists(f)) D[[c1]] <- readRDS(f) else cat("缺:", f, "\n")
}
cat("载入队列:", paste(names(D), collapse=", "), "\n")
for (c1 in names(D)) cat(sprintf("  %-10s n = %d\n", c1, nrow(D[[c1]])))

sp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 30) return(c(rho = NA, p = NA, n = sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

# 秩残差偏相关（消除单调非线性）
pcor <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  if (sum(ok) < 30) return(c(rho = NA, p = NA, n = sum(ok)))
  rx <- resid(lm(rank(x[ok]) ~ rank(z[ok])))
  ry <- resid(lm(rank(y[ok]) ~ rank(z[ok])))
  ct <- suppressWarnings(cor.test(rx, ry, method = "pearson"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

# ------------------------------------------------------------------ 主表
res <- list()
for (c1 in names(D)) {
  d <- D[[c1]]
  for (m in modules_all) {
    if (!m %in% names(d)) next
    r <- sp(d$RiskScore, d[[m]])
    res[[length(res)+1]] <- data.frame(cohort = c1, analysis = "primary",
                                       variable = m, rho = r[1], p = r[2], n = r[3],
                                       stringsAsFactors = FALSE)
  }
  # 单基因层面
  for (g in grep("^g_", names(d), value = TRUE)) {
    r <- sp(d$RiskScore, d[[g]])
    res[[length(res)+1]] <- data.frame(cohort = c1, analysis = "single_gene",
                                       variable = g, rho = r[1], p = r[2], n = r[3],
                                       stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, res)

# 队列内 BH 校正
R$p_BH <- NA
for (c1 in unique(R$cohort)) {
  for (a in unique(R$analysis)) {
    i <- which(R$cohort == c1 & R$analysis == a)
    R$p_BH[i] <- p.adjust(R$p[i], method = "BH")
  }
}

# ------------------------------------------------------------------ 循环性拆解
cat("\n\n================ 一、B 细胞相关性：真实信号 vs 构造性同义反复 ================\n")
cat(sprintf("%-10s %14s %14s %14s\n", "队列",
            "含CD79A模块", "不含CD79A模块", "控制CD79A后"))
decomp <- list()
for (c1 in names(D)) {
  d <- D[[c1]]
  if (!all(c("Bcell_withCD79A","Bcell_noCD79A","g_CD79A") %in% names(d))) next
  a <- sp(d$RiskScore, d$Bcell_withCD79A)
  b <- sp(d$RiskScore, d$Bcell_noCD79A)
  cc <- pcor(d$RiskScore, d$Bcell_noCD79A, d$g_CD79A)
  cat(sprintf("%-10s %6.3f (p=%.1e) %6.3f (p=%.2f) %6.3f (p=%.2f)   n=%d\n",
              c1, a[1], a[2], b[1], b[2], cc[1], cc[2], a[3]))
  decomp[[c1]] <- data.frame(cohort = c1,
    rho_withCD79A = a[1], p_withCD79A = a[2],
    rho_noCD79A   = b[1], p_noCD79A   = b[2],
    rho_noCD79A_ctrlCD79A = cc[1], p_noCD79A_ctrlCD79A = cc[2], n = a[3])
}
DEC <- do.call(rbind, decomp)

# ------------------------------------------------------------------ 混杂校正
cat("\n\n================ 二、控制纯度与增殖后的免疫相关性 ================\n")
cat("（Purity_epi = 上皮/肿瘤细胞模块，作为肿瘤纯度的代理变量）\n\n")
adj <- list()
for (c1 in names(D)) {
  d <- D[[c1]]
  cat(sprintf("---- %s (n=%d) ----\n", c1, nrow(d)))
  cat(sprintf("  %-16s %18s %22s %22s\n", "免疫模块", "原始 rho",
              "控制纯度后", "控制纯度+增殖后"))
  for (m in c("Bcell_noCD79A","CD8_T","Treg","Cytotoxic","Macrophage","Neutrophil")) {
    if (!m %in% names(d)) next
    raw <- sp(d$RiskScore, d[[m]])
    a1 <- pcor(d$RiskScore, d[[m]], d$Purity_epi)
    # 双协变量：先分别对两个协变量回归取残差
    ok <- is.finite(d$RiskScore) & is.finite(d[[m]]) & is.finite(d$Purity_epi) & is.finite(d$Prolif)
    rx <- resid(lm(rank(d$RiskScore[ok]) ~ rank(d$Purity_epi[ok]) + rank(d$Prolif[ok])))
    ry <- resid(lm(rank(d[[m]][ok]) ~ rank(d$Purity_epi[ok]) + rank(d$Prolif[ok])))
    a2r <- suppressWarnings(cor.test(rx, ry))
    a2 <- c(rho = unname(a2r$estimate), p = a2r$p.value)
    cat(sprintf("  %-16s %7.3f (p=%.1e) %7.3f (p=%.2f)  %7.3f (p=%.2f)\n",
                m, raw[1], raw[2], a1[1], a1[2], a2[1], a2[2]))
    adj[[length(adj)+1]] <- data.frame(cohort = c1, module = m,
      rho_raw = raw[1], p_raw = raw[2],
      rho_adjPurity = a1[1], p_adjPurity = a1[2],
      rho_adjPurityProlif = a2[1], p_adjPurityProlif = a2[2], n = raw[3])
  }
  cat("\n")
}
ADJ <- do.call(rbind, adj)

# ------------------------------------------------------------------ 跨队列一致性
cat("\n\n================ 三、跨队列一致性汇总 ================\n")
cat(sprintf("%-16s %s\n", "模块", "各队列 rho（括号内 p）"))
cons <- list()
for (m in modules_all) {
  sub <- R[R$analysis == "primary" & R$variable == m, ]
  if (nrow(sub) < 2) next
  sub <- sub[match(intersect(cohorts, sub$cohort), sub$cohort), ]
  txt <- paste(sprintf("%s:%.2f(%.2g)", sub$cohort, sub$rho, sub$p), collapse = "  ")
  same <- length(unique(sign(sub$rho[!is.na(sub$rho)]))) == 1
  # Stouffer 合并（方向一致时）
  z <- qnorm(sub$p/2) * sign(sub$rho)
  zc <- sum(z) / sqrt(length(z))
  cat(sprintf("%-16s %s   [符号%s, 合并 p=%.2g]\n", m, txt,
              ifelse(same, "一致", "不一致"), 2*pnorm(-abs(zc))))
  cons[[m]] <- data.frame(module = m, n_cohort = nrow(sub),
    sign_consistent = same, n_neg = sum(sub$rho < 0, na.rm = TRUE),
    rho_min = min(sub$rho, na.rm = TRUE), rho_max = max(sub$rho, na.rm = TRUE),
    combined_p = 2*pnorm(-abs(zc)), stringsAsFactors = FALSE)
}
CONS <- do.call(rbind, cons)

# ------------------------------------------------------------------ 高低风险组
cat("\n\n================ 四、高/低风险组免疫模块差异（组内中位数分组）================\n")
grp <- list()
for (c1 in names(D)) {
  d <- D[[c1]]
  d$grp <- ifelse(d$RiskScore > median(d$RiskScore, na.rm = TRUE), "High", "Low")
  cat(sprintf("---- %s ----\n", c1))
  for (m in c("Bcell_noCD79A","CD8_T","Treg","Cytotoxic","Macrophage","Neutrophil","Purity_epi")) {
    if (!m %in% names(d)) next
    x <- d[[m]][d$grp == "High"]; y <- d[[m]][d$grp == "Low"]
    wt <- suppressWarnings(wilcox.test(x, y))
    dbar <- (median(x, na.rm=TRUE) - median(y, na.rm=TRUE)) /
            sd(d[[m]], na.rm = TRUE)
    cat(sprintf("  %-14s High %.3f vs Low %.3f   Δ(std) = %+.2f   p = %.2g\n",
                m, median(x, na.rm=TRUE), median(y, na.rm=TRUE), dbar, wt$p.value))
    grp[[length(grp)+1]] <- data.frame(cohort = c1, module = m,
      median_High = median(x, na.rm=TRUE), median_Low = median(y, na.rm=TRUE),
      delta_std = dbar, p = wt$p.value)
  }
  cat("\n")
}
GRP <- do.call(rbind, grp)

# ------------------------------------------------------------------ 导出
write.csv(R,    "Table_ImmuneAxis_AllCorrelations.csv", row.names = FALSE)
write.csv(DEC,  "Table_ImmuneAxis_BcellDecomposition.csv", row.names = FALSE)
write.csv(ADJ,  "Table_ImmuneAxis_Adjusted.csv", row.names = FALSE)
write.csv(CONS, "Table_ImmuneAxis_Consistency.csv", row.names = FALSE)
write.csv(GRP,  "Table_ImmuneAxis_GroupDiff.csv", row.names = FALSE)
cat("\n已导出 5 张表：Table_ImmuneAxis_*.csv\n")
saveRDS(list(R=R, DEC=DEC, ADJ=ADJ, CONS=CONS, GRP=GRP), "ia_results.rds")
