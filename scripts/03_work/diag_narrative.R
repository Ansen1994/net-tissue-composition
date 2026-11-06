# ==============================================================================
# 叙事可行性验证：这个 risk score 到底是不是"独立于增殖的轴"
# 目的：在决定是否转向"非增殖免疫轴"之前，用现有数据给出判定
# 输出：相关性 / 偏相关 / 偏C-index / 分层KM / 随机基因集对照
# ==============================================================================
suppressMessages({
  library(SummarizedExperiment)
  library(org.Hs.eg.db)
  library(survival)
})

load("../raw/TCGA/brca.gdc_2024.rda")
expr <- assay(expquery2)
colnames(expr) <- substr(gsub("\\.", "-", colnames(expr)), 1, 12)
ens <- gsub("\\..*", "", rownames(expr))
sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(expr) <- sym
expr <- expr[!is.na(rownames(expr)), ]
expr <- expr[!duplicated(rownames(expr)), ]
expr <- as.matrix(expr)

fd <- read.csv("Table_S7_Final_Clinical_Data.csv", stringsAsFactors = FALSE)
s <- intersect(colnames(expr), fd$Patient_ID)
cat("匹配样本:", length(s), "\n")

# ---- 基因集 ----
prolif <- c("ANLN","BIRC5","BUB1","CCNB1","CDC20","CDC25B","CEP55","CKS2","EXO1",
            "KPNA2","MELK","MKI67","MYBL2","NUF2","ORC6L","PTTG1","RRM2","TYMS","UBE2C")
prolif5 <- c("MKI67","AURKA","BIRC5","CCNB1","MYBL2")
immune <- c("CD8A","CD8B","GZMA","GZMB","PRF1","CXCL9","CXCL10","PDCD1","CTLA4",
            "FOXP3","MS4A1","CD79A","CD79B","TNFRSF17","IGHM")
stroma <- c("COL1A1","COL1A2","COL3A1","COL5A1","FAP","ACTA2","TAGLN","PDGFRB","THY1","LUM")

zscore_score <- function(genes, samples) {
  g <- intersect(genes, rownames(expr))
  if (length(g) < 3) return(rep(NA, length(samples)))
  m <- expr[g, samples, drop = FALSE]
  z <- t(scale(t(log2(m + 1))))
  z[is.na(z)] <- 0
  colMeans(z)
}

df <- data.frame(Patient_ID = s,
                 Prolif   = zscore_score(prolif, s),
                 Prolif5  = zscore_score(prolif5, s),
                 Immune   = zscore_score(immune, s),
                 Stroma   = zscore_score(stroma, s),
                 check.names = FALSE)
df <- merge(df, fd[, c("Patient_ID","RiskScore","Age","Stage","OS","OS.time","RiskGroup")],
            by = "Patient_ID")
df <- df[!is.na(df$OS) & !is.na(df$OS.time) & df$OS.time > 0, ]
cat("纳入分析:", nrow(df), " 事件:", sum(df$OS == 1), "\n\n")

sp <- function(a, b) {
  ct <- suppressWarnings(cor.test(a, b, method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value)
}

cat("=== 1. risk score 与各模块的相关性 (Spearman) ===\n")
for (v in c("Prolif","Prolif5","Immune","Stroma")) {
  r <- sp(df$RiskScore, df[[v]])
  cat(sprintf("  %-8s rho = %6.3f   p = %.3g\n", v, r[1], r[2]))
}

# ---- 偏相关：控制增殖后 risk 与免疫的关系 ----
pcor <- function(x, y, z) {
  rx <- resid(lm(x ~ z)); ry <- resid(lm(y ~ z))
  cor(rx, ry, method = "spearman")
}
cat("\n=== 2. 偏相关（控制增殖模块）===\n")
cat(sprintf("  risk vs Immune | Prolif = %6.3f\n", pcor(df$RiskScore, df$Immune, df$Prolif)))
cat(sprintf("  risk vs Stroma | Prolif = %6.3f\n", pcor(df$RiskScore, df$Stroma, df$Prolif)))

# ---- 随机基因集对照：无关是特异的还是普遍的 ----
set.seed(42)
allg <- rownames(expr)[rowMeans(expr[, s]) > 1]
nullrho <- numeric(1000)
for (i in 1:1000) {
  g <- sample(allg, length(prolif))
  sc <- zscore_score(g, s)
  d2 <- data.frame(Patient_ID = s, sc = sc)
  d2 <- merge(d2, df[, c("Patient_ID","RiskScore")], by = "Patient_ID")
  nullrho[i] <- suppressWarnings(cor(d2$sc, d2$RiskScore, method = "spearman"))
}
obs <- sp(df$RiskScore, df$Prolif)[1]
cat("\n=== 3. 随机 19 基因集对照 (1000 次) ===\n")
cat(sprintf("  随机集 rho 分布: 中位 %.3f, 95%% 区间 [%.3f, %.3f]\n",
            median(nullrho), quantile(nullrho, .025), quantile(nullrho, .975)))
cat(sprintf("  实际增殖集 rho = %.3f, 百分位 = %.1f%%\n",
            obs, mean(nullrho < obs) * 100))

# ---- 多因素 Cox：risk 与增殖是否互相独立 ----
df$Stage2 <- factor(ifelse(grepl("IV", df$Stage), "IV",
                    ifelse(grepl("III", df$Stage), "III",
                    ifelse(grepl("II", df$Stage), "II", "I"))), levels = c("I","II","III","IV"))
df$ProlifGrp <- ifelse(df$Prolif > median(df$Prolif), "High", "Low")

cat("\n=== 4. 多因素 Cox（risk + 增殖 + 年龄 + 分期）===\n")
fit <- coxph(Surv(OS.time, OS) ~ RiskScore + Prolif + Age + Stage2, data = df)
print(summary(fit)$conf.int[, c(1,3,4)], digits = 3)
cat("  p 值:\n"); print(round(summary(fit)$coefficients[, "Pr(>|z|)"], 4))

# ---- 偏 C-index：增殖之上，risk 还能加多少 ----
cidx <- function(f, d) {
  if (inherits(f, "try-error")) return(NA)
  lp <- try(predict(f, newdata = d, type = "lp"), silent = TRUE)
  if (inherits(lp, "try-error")) return(NA)
  as.numeric(concordance(Surv(d$OS.time, d$OS) ~ lp, reverse = TRUE)$concordance)
}
c_prolif <- cidx(coxph(Surv(OS.time, OS) ~ Prolif, data = df), df)
c_risk   <- cidx(coxph(Surv(OS.time, OS) ~ RiskScore, data = df), df)
c_both   <- cidx(coxph(Surv(OS.time, OS) ~ Prolif + RiskScore, data = df), df)
c_clin   <- cidx(coxph(Surv(OS.time, OS) ~ Age + Stage2, data = df), df)
c_clin_r <- cidx(coxph(Surv(OS.time, OS) ~ Age + Stage2 + Prolif, data = df), df)
c_all    <- cidx(coxph(Surv(OS.time, OS) ~ Age + Stage2 + Prolif + RiskScore, data = df), df)
cat("\n=== 5. C-index 对比 ===\n")
cat(sprintf("  增殖单独            %.4f\n", c_prolif))
cat(sprintf("  risk 单独           %.4f\n", c_risk))
cat(sprintf("  增殖 + risk         %.4f   (Δ vs 增殖 = %+.4f)\n", c_both, c_both - c_prolif))
cat(sprintf("  临床单独            %.4f\n", c_clin))
cat(sprintf("  临床 + 增殖         %.4f\n", c_clin_r))
cat(sprintf("  临床 + 增殖 + risk  %.4f   (Δ vs 临床+增殖 = %+.4f)\n", c_all, c_all - c_clin_r))

# ---- 分层：增殖高/低亚组内，risk 还有没有分层能力 ----
cat("\n=== 6. 按增殖分层的 risk 组 HR ===\n")
for (lv in c("Low","High")) {
  sub <- df[df$ProlifGrp == lv, ]
  sub$RGnum <- as.numeric(sub$RiskScore > median(sub$RiskScore))
  f <- coxph(Surv(OS.time, OS) ~ RGnum, data = sub)
  sm <- summary(f)
  cat(sprintf("  增殖%s组 (n=%d, 事件=%d): HR = %.2f (%.2f-%.2f), p = %.4f\n",
              lv, nrow(sub), sum(sub$OS), sm$conf.int[1], sm$conf.int[3],
              sm$conf.int[4], sm$coefficients[1, 5]))
}
df$RGnum <- as.numeric(df$RiskScore > median(df$RiskScore))
lr <- survdiff(Surv(OS.time, OS) ~ RGnum, data = df)
cat(sprintf("  全体 log-rank p = %.3g\n", 1 - pchisq(lr$chisq, length(lr$n) - 1)))

# ---- 与 PAM50 近似分型的关系 ----
cat("\n=== 7. 与受体状态的关系 ===\n")
for (g in c("ESR1","PGR","ERBB2")) {
  if (g %in% rownames(expr)) {
    r <- sp(df$RiskScore, expr[g, df$Patient_ID])
    cat(sprintf("  %-6s rho = %6.3f  p = %.3g\n", g, r[1], r[2]))
  }
}

saveRDS(df, "feasibility_narrative.rds")
cat("\n完成，结果对象已存 feasibility_narrative.rds\n")
