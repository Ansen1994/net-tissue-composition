# ==============================================================================
# Figure 3 全图重画（A–G），统一 v2 风格
# 输入：Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv（n = 762, 99 events）
# A 多因素 Cox 森林 · B 风险评分排序 · C 生存状态 · D 四基因热图
# E KM（median 分组）· F 乐观度 · G 时间依赖 ROC
# ==============================================================================
source("figs_v2/style.R")

d0 <- read.csv("Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv",
               row.names = 1, check.names = FALSE)
genes <- c("PLG", "CD79A", "MMP7", "ANGPT1")
d <- d0[, c("OS.time", "OS", genes)]
d <- d[complete.cases(d) & d$OS.time > 0, ]

mf <- coxph(Surv(OS.time, OS) ~ PLG + CD79A + MMP7 + ANGPT1, data = d)
d$risk <- as.numeric(predict(mf, type = "risk"))
d$lp   <- as.numeric(predict(mf, type = "lp"))
d$grp  <- factor(ifelse(d$risk > median(d$risk), "High", "Low"),
                 levels = c("Low", "High"))
n_low  <- sum(d$grp == "Low")
d <- d[order(d$risk), ]
d$idx <- seq_len(nrow(d))
cat("n =", nrow(d), " events =", sum(d$OS), " nLow =", n_low, "\n")

GRP <- c(Low = "#4C9BC4", High = "#D6604D")

# ---------------------------------------------------------------- A 森林图
s  <- summary(mf)
co <- data.frame(gene = rownames(s$conf.int),
                 HR = s$conf.int[, 1], lo = s$conf.int[, 3], hi = s$conf.int[, 4],
                 p  = s$coefficients[, 5], stringsAsFactors = FALSE)
co <- co[order(co$HR), ]
co$y <- seq_len(nrow(co))
co$p_lab <- ifelse(co$p < 0.001, "<0.001", sprintf("%.3f", co$p))
# 显著性标记：原来只印 p 值，读者要自己比对 0.05，这里按惯例补星号
co$stars <- ifelse(co$p < 0.001, "***",
            ifelse(co$p < 0.01,  "**",
            ifelse(co$p < 0.05,  "*", "")))
co$p_lab <- ifelse(co$stars == "", co$p_lab, paste0(co$p_lab, " ", co$stars))
co$hr_lab <- sprintf("%.2f (%.2f\u2013%.2f)", co$HR, co$lo, co$hi)

pA_txt <- ggplot(co, aes(y = y)) +
  geom_text(aes(x = 0, label = gene), hjust = 0, size = 2.7, fontface = "italic") +
  geom_text(aes(x = 0.95, label = hr_lab), hjust = 0, size = 2.5) +
  geom_text(aes(x = 2.55, label = p_lab), hjust = 0, size = 2.5) +
  annotate("text", x = 0.95, y = max(co$y) + 0.85, label = "HR (95% CI)",
           hjust = 0, size = 2.3, colour = "grey30") +
  annotate("text", x = 2.55, y = max(co$y) + 0.85, label = "p",
           hjust = 0, size = 2.3, colour = "grey30") +
  scale_x_continuous(limits = c(-0.05, 3.55)) +
  scale_y_continuous(limits = c(0.4, max(co$y) + 1.1)) +
  labs(title = "Multivariate Cox model, training cohort") +
  theme_void(base_size = 7.5) +
  theme(plot.title = element_text(face = "bold", size = 8, hjust = 0,
                                  margin = margin(b = 3)),
        plot.margin = margin(4, 0, 4, 4))

pA_f <- ggplot(co, aes(y = y)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45,
                 colour = "grey25") +
  geom_point(aes(x = HR), shape = 15, size = 2.4, colour = "#0072B2") +
  scale_x_log10(breaks = c(0.9, 1.0, 1.1, 1.2, 1.3, 1.4),
                limits = c(0.82, 1.45),
                labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(limits = c(0.4, max(co$y) + 1.1)) +
  labs(x = "Hazard ratio (log scale)", y = NULL) +
  theme_pub(7.5) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = margin(4, 4, 4, 2))

pA <- (pA_txt | pA_f) + plot_layout(widths = c(0.44, 0.56))
pA <- wrap_elements(full = pA) + labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = 10.5),
        plot.tag.position = c(0, 1.055))

# ---------------------------------------------------------------- B 风险评分
pB <- ggplot(d, aes(idx, risk)) +
  geom_point(aes(colour = grp), size = 0.35, stroke = 0) +
  geom_vline(xintercept = n_low + 0.5, linetype = 2, linewidth = 0.3,
             colour = "grey35") +
  scale_colour_manual(values = GRP, guide = "none") +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(x = "Patients ranked by risk score", y = "Risk score",
       title = "Risk-score distribution") +
  theme_pub(7.5) +
  annotate("text", x = 8, y = max(d$risk), label = "Low", hjust = 0,
           vjust = 1.2, size = 2.3, colour = GRP["Low"]) +
  annotate("text", x = nrow(d) - 8, y = max(d$risk), label = "High", hjust = 1,
           vjust = 1.2, size = 2.3, colour = GRP["High"])

# ---------------------------------------------------------------- C 生存状态
pC <- ggplot(d, aes(idx, OS.time / 365.25)) +
  geom_point(aes(colour = factor(OS), shape = factor(OS)), size = 0.45,
             stroke = 0.15) +
  geom_vline(xintercept = n_low + 0.5, linetype = 2, linewidth = 0.3,
             colour = "grey35") +
  scale_colour_manual(values = c("0" = "#4C9BC4", "1" = "#D6604D"),
                      labels = c("Alive", "Dead"), name = NULL) +
  scale_shape_manual(values = c("0" = 1, "1" = 16),
                     labels = c("Alive", "Dead"), name = NULL) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(x = "Patients ranked by risk score", y = "Follow-up (years)",
       title = "Survival status") +
  theme_pub(7.5) +
  theme(legend.position = c(0.85, 0.74),
        legend.key.size = unit(2.4, "mm"),
        legend.background = element_rect(fill = "white", colour = NA))

# ---------------------------------------------------------------- D 热图
hm <- do.call(rbind, lapply(genes, function(g) {
  v <- as.numeric(d[[g]])
  data.frame(gene = g, idx = d$idx, z = (v - mean(v)) / sd(v))
}))
hm$gene <- factor(hm$gene, levels = rev(c("ANGPT1", "PLG", "MMP7", "CD79A")))

pD <- ggplot(hm, aes(idx, gene, fill = z)) +
  geom_raster() +
  geom_vline(xintercept = n_low + 0.5, linetype = 2, linewidth = 0.3,
             colour = "black") +
  scale_fill_gradientn(colours = div_cols, limits = c(-2.2, 2.2),
                       oob = squish, name = "z-score",
                       breaks = c(-2, 0, 2)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(x = "Patients ranked by risk score", y = NULL,
       title = "Hub-gene expression") +
  theme_pub(7.5) +
  theme(axis.text.y = element_text(face = "italic", size = 6.6),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        legend.position = "right",
        legend.key.width = unit(2.0, "mm"),
        legend.key.height = unit(5.0, "mm"),
        legend.title = element_text(size = 6))

# ---------------------------------------------------------------- E KM
fit <- survfit(Surv(OS.time, OS) ~ grp, data = d)
kd  <- data.frame(time = fit$time / 365.25,
                  surv = fit$surv,
                  n.censor = fit$n.censor,
                  grp = rep(c("Low", "High"), fit$strata))
kd  <- kd[order(kd$grp, kd$time), ]
# 每个分组补齐起点
kd  <- rbind(data.frame(time = 0, surv = 1, n.censor = 0,
                        grp = c("Low", "High")), kd)
lr  <- survdiff(Surv(OS.time, OS) ~ grp, data = d)
pLR <- 1 - pchisq(lr$chisq, df = 1)
pLR_lab <- if (pLR < 1e-4) "Log-rank p < 0.0001" else sprintf("Log-rank p = %.4f", pLR)
hr  <- summary(coxph(Surv(OS.time, OS) ~ grp, data = d))
hrT <- sprintf("HR = %.2f (%.2f\u2013%.2f)", hr$conf.int[1], hr$conf.int[3], hr$conf.int[4])

pE <- ggplot(kd, aes(time, surv, colour = grp)) +
  geom_step(linewidth = 0.55) +
  scale_colour_manual(values = GRP,
                      labels = c(sprintf("Low (n = %d)", n_low),
                                 sprintf("High (n = %d)", nrow(d) - n_low)),
                      name = NULL) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, 2),
                     expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25),
                     expand = c(0, 0), labels = function(x) sprintf("%.2f", x)) +
  labs(x = "Time (years)", y = "Overall survival",
       title = "Kaplan\u2013Meier, median split") +
  annotate("text", x = 0.25, y = 0.10, hjust = 0, size = 2.4, label = pLR_lab) +
  annotate("text", x = 0.25, y = 0.03, hjust = 0, size = 2.4, label = hrT) +
  theme_pub(7.5) +
  theme(legend.position = c(0.70, 0.60),
        legend.key.size = unit(2.6, "mm"),
        legend.background = element_rect(fill = "white", colour = NA))

# ---------------------------------------------------------------- F 乐观度
df <- data.frame(x = c(1, 2), C = c(0.690, 0.638))
pF <- ggplot(df, aes(x, C)) +
  geom_segment(aes(x = 1, xend = 1.86, y = 0.690, yend = 0.645),
               arrow = arrow(length = unit(1.4, "mm"), type = "closed"),
               linewidth = 0.35, colour = "grey35") +
  geom_point(size = 2.6, colour = c("#0072B2", "#D55E00")) +
  geom_text(aes(label = sprintf("%.3f", C)), vjust = -1.25, size = 2.4) +
  annotate("text", x = 1.50, y = 0.712, label = "\u0394 = \u22120.052",
           size = 2.3, colour = "grey25") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("Apparent", "5-fold CV"),
                     limits = c(0.58, 2.42)) +
  scale_y_continuous(limits = c(0.605, 0.725), breaks = seq(0.62, 0.70, 0.04)) +
  labs(x = NULL, y = "C-index", title = "Optimism") +
  theme_pub(7.5) +
  theme(axis.text.x = element_text(size = 6.2))

# ---------------------------------------------------------------- G ROC
tp <- c(365, 1095, 1825)
ro <- timeROC(T = d$OS.time, delta = d$OS, marker = d$lp, cause = 1,
              times = tp, iid = FALSE)
lab_t <- c("1 year", "3 years", "5 years")
# 每个时间窗内实际发生的死亡数。1 年只有 8 例，AUC 0.817 是在这 8 个事件上
# 算出来的，不确定度极大；不标出来读者会把它误读成"模型判别力很好"。
ev <- sapply(tp, function(t) sum(d$OS == 1 & d$OS.time <= t))
roc <- do.call(rbind, lapply(seq_along(tp), function(i) {
  data.frame(FP = ro$FP[, i], TP = ro$TP[, i],
             t = sprintf("%s \u00b7 AUC %.3f \u00b7 %d events",
                         lab_t[i], ro$AUC[i], ev[i]))
}))
roc$t <- factor(roc$t, levels = unique(roc$t))
ROC_COL <- c("#0072B2", "#009E73", "#D55E00")

pG <- ggplot(roc, aes(FP, TP, colour = t)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3,
              colour = "grey55") +
  geom_step(linewidth = 0.55) +
  scale_colour_manual(values = ROC_COL, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = c(0.01, 0)) +
  labs(x = "1 \u2212 specificity", y = "Sensitivity",
       title = "Time-dependent ROC (in-sample)") +
  theme_pub(7.5) +
  theme(legend.position = c(0.985, 0.045),
        legend.justification = c(1, 0),
        legend.key.size = unit(2.4, "mm"),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.text = element_text(size = 5.5))

# ---------------------------------------------------------------- 拼版
fig <- (pA) /
       (tag_panel(pB, "B") | tag_panel(pC, "C")) /
       (tag_panel(pD, "D") | tag_panel(pE, "E")) /
       (tag_panel(pF, "F") | tag_panel(pG, "G")) +
  plot_layout(heights = c(1.02, 1, 1.12, 1))
# 首行的 A 字母原本距纸张上缘只有 0.4 pt，制版裁切可能切到字形；
# 顶部加一个约 2 mm 的空白行（patchwork 对 plot.margin 不生效）
fig <- (plot_spacer() / fig) + plot_layout(heights = c(0.003, 1))

save_fig(fig, "Figure3_ABCDEFG", 170 * MM, 205 * MM)
