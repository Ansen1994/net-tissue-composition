# ==============================================================================
# Table S15b —— GSE96058 化疗与非化疗患者的基线比较
#
# Results 3.5 引用了三组基线数字（年龄 55.2 vs 67.3 岁、肿瘤大小 22.0 vs 18.3 mm、
# 淋巴结阳性 57% vs 26%），此前没有对应表格。本脚本按 j2_treatment_interaction.R
# 相同的合并方式重建该表，供正文引用。
#
# 输入：Table_S6_GSE96058_Validation_Result.csv
#       ../raw/GSE96058/GSE96058-GPL11154_series_matrix.txt.gz
#       ../raw/GSE96058/GSE96058-GPL18573_series_matrix.txt.gz
# 输出：Table_S15b_Baseline_by_chemo.csv
# ==============================================================================
suppressMessages({ library(GEOquery) })

d96 <- read.csv("Table_S6_GSE96058_Validation_Result.csv", stringsAsFactors = FALSE)

g1 <- getGEO(filename = "../raw/GSE96058/GSE96058-GPL11154_series_matrix.txt.gz", getGPL = FALSE)
g2 <- getGEO(filename = "../raw/GSE96058/GSE96058-GPL18573_series_matrix.txt.gz", getGPL = FALSE)
cols <- c("title", "chemo treated:ch1")
cl <- rbind(pData(g1)[, intersect(cols, names(pData(g1)))],
            pData(g2)[, intersect(cols, names(pData(g2)))])
names(cl) <- c("title", "chemo")

m <- merge(d96, cl, by = "title")
m$Chemo <- as.numeric(m$chemo)
m$LN    <- as.numeric(as.character(m$LymphNode))
m$Age   <- as.numeric(m$Age)
m$Size  <- as.numeric(m$TumorSize)

cat("合并后 n =", nrow(m), " 化疗字段可用 =", sum(!is.na(m$Chemo)),
    " 接受化疗 =", sum(m$Chemo == 1, na.rm = TRUE), "\n")

tr <- m[which(m$Chemo == 1), ]   # which() 必要：NA 会被 m[NA,] 复制成整行 NA
un <- m[which(m$Chemo == 0), ]

fmt <- function(x, dig = 1) sprintf(paste0("%.", dig, "f"), mean(x, na.rm = TRUE))
sdf <- function(x, dig = 1) sprintf(paste0("%.", dig, "f"), sd(x, na.rm = TRUE))

rows <- list()
add <- function(var, unit, f, test, dig = 1, pct = FALSE) {
  a <- f(tr[[var]]); b <- f(un[[var]])
  p <- tryCatch(test(tr[[var]], un[[var]]), error = function(e) NA_real_)
  rows[[length(rows) + 1]] <<- data.frame(
    variable = var, unit = unit,
    chemo_mean = a, chemo_sd = sdf(tr[[var]], dig),
    no_chemo_mean = b, no_chemo_sd = sdf(un[[var]], dig),
    p = signif(p, 3), n_chemo = sum(!is.na(tr[[var]])),
    n_no_chemo = sum(!is.na(un[[var]])), stringsAsFactors = FALSE)
}

add("Age", "years", fmt, function(a, b) wilcox.test(a, b)$p.value)

# 肿瘤大小：0 值视为缺失（GSE96058 用 0 表示未记录）
m$Size[m$Size == 0] <- NA
tr$Size <- m$Size[which(m$Chemo == 1)]; un$Size <- m$Size[which(m$Chemo == 0)]
add("Size", "mm", fmt, function(a, b) wilcox.test(a, b)$p.value)

# 淋巴结阳性率
add("LN", "proportion positive", function(x) sprintf("%.1f", 100 * mean(x, na.rm = TRUE)),
    function(a, b) chisq.test(table(c(a, b), c(rep(1, length(a)), rep(0, length(b)))))$p.value,
    dig = 1)

out <- do.call(rbind, rows)
# 比例行的离散度不报告 SD，留空以免被误读为百分比的标准差
out$chemo_sd[out$variable == "LN"] <- ""
out$no_chemo_sd[out$variable == "LN"] <- ""
out$n_total <- nrow(m)
write.csv(out, "Table_S15b_Baseline_by_chemo.csv", row.names = FALSE)

cat("\n==== Table S15b ====\n")
print(out[, c("variable","unit","chemo_mean","chemo_sd","no_chemo_mean","no_chemo_sd","p")],
      row.names = FALSE)
cat("\n化疗 n =", nrow(tr), " 未化疗 n =", nrow(un), "\n")
