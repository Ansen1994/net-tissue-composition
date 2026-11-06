# ==============================================================================
# Figure 2 — 双机器学习把候选池收窄到 4 个 hub 基因（v2 重排版）
# A LASSO 系数路径 | B LASSO 十折交叉验证 | C RF 重要性 | D 交集
# 数字与 Table S4 一致：LASSO 11、RF 58、交集 4
# ==============================================================================
source("figs_v2/style.R")
suppressMessages({ library(glmnet); library(randomForest) })

BLU <- "#0072B2"; GRN <- "#009E73"; ORG <- "#D55E00"
MAG <- "#CC79A7"; SKY <- "#56B4E9"; YEL <- "#E69F00"; GRY <- "#B4B4B4"

HUB <- c("PLG", "CD79A", "MMP7", "ANGPT1")

ml <- read.csv("Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv", row.names = 1, check.names = FALSE)
X <- as.matrix(ml[, setdiff(colnames(ml), c("OS", "OS.time"))])
surv <- Surv(ml$OS.time, ml$OS)

# ---------------------------------------------------------------- 拟合
set.seed(2026)
fit <- glmnet(X, surv, family = "cox", alpha = 1)
set.seed(2026)
cv  <- cv.glmnet(X, surv, family = "cox", alpha = 1, nfolds = 10)

# 两条参考线（A、B 共用，保证纵向对齐）
lmin <- log(cv$lambda.min)
l1se <- log(cv$lambda.1se)

# ---------------------------------------------------------------- A 路径
beta <- as.matrix(fit$beta)
traj <- data.frame(
  loglam = rep(log(fit$lambda), each = nrow(beta)),
  coef   = as.vector(beta),
  gene   = rep(rownames(beta), times = ncol(beta))
)
traj$hub <- traj$gene %in% HUB
traj$gene <- factor(traj$gene, levels = rownames(beta))

hub_cols <- setNames(c(ORG, BLU, GRN, MAG), HUB)

# 顶部副轴用：log(lambda) 递增排列对应的非零系数个数
.o <- order(log(fit$lambda))
xs_df <- data.frame(x = log(fit$lambda)[.o], y = fit$df[.o])
BRK <- c(-9, -7, -5)

pA <- ggplot() +
  geom_line(data = traj[!traj$hub, ], aes(loglam, coef, group = gene),
            colour = "#D6D6D6", linewidth = 0.22) +
  geom_line(data = traj[traj$hub, ], aes(loglam, coef, group = gene, colour = gene),
            linewidth = 0.62) +
  scale_colour_manual(values = hub_cols, name = NULL) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  geom_vline(xintercept = l1se, linetype = 2, linewidth = 0.35, colour = GRY) +
  geom_vline(xintercept = lmin, linetype = 2, linewidth = 0.35, colour = BLU) +
  annotate("text", x = lmin + 0.10, y = 1.245, label = "11", hjust = 0,
           size = 1.95, fontface = "bold", colour = BLU) +
  annotate("segment", x = BRK, xend = BRK, y = 1.155, yend = 1.215,
           linewidth = 0.25, colour = "grey35") +
  annotate("text", x = BRK, y = 1.265,
           label = as.character(round(approx(xs_df$x, xs_df$y, xout = BRK, rule = 2)$y)),
           size = 1.85, colour = "grey25") +
  annotate("text", x = -10.42, y = 1.265, label = "Non-zero", hjust = 0,
           size = 1.85, colour = "grey35") +
  annotate("text", x = -10.42, y = 1.175, label = "features", hjust = 0,
           size = 1.85, colour = "grey35") +
  scale_x_continuous(limits = c(-10.45, -2.95), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.05, 1.42), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  labs(x = expression(log(lambda)), y = "Coefficient",
       title = "LASSO coefficient paths") +
  theme_pub(7.2) +
  theme(legend.position = c(0.185, 0.135),
        legend.key.size = unit(2.2, "mm"),
        legend.text = element_text(size = 5.8, face = "italic"),
        legend.spacing.y = unit(0.4, "mm"),
        legend.background = element_rect(fill = "white", colour = NA)) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

# ---------------------------------------------------------------- B CV
cvd <- data.frame(loglam = log(cv$lambda), cvm = cv$cvm,
                  lo = cv$cvlo, hi = cv$cvup, nz = cv$nzero)

pB <- ggplot(cvd, aes(loglam, cvm)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.28, colour = "#BEBEBE") +
  geom_point(aes(colour = loglam >= lmin), size = 0.85, stroke = 0) +
  scale_colour_manual(values = c("FALSE" = GRY, "TRUE" = ORG), guide = "none") +
  geom_vline(xintercept = lmin, linetype = 2, linewidth = 0.35, colour = BLU) +
  geom_vline(xintercept = l1se, linetype = 2, linewidth = 0.35, colour = GRY) +
  annotate("text", x = l1se, y = max(cvd$hi) * 0.965, label = "lambda[1*se]",
           parse = TRUE, hjust = 1.05, size = 2.0, colour = "grey40") +
  annotate("text", x = lmin, y = max(cvd$hi) * 0.855, label = "lambda[min]",
           parse = TRUE, hjust = 1.05, size = 2.0, colour = BLU) +
  scale_x_continuous(limits = c(-10.45, -2.95), expand = c(0, 0)) +
  labs(x = expression(log(lambda)), y = "Partial likelihood deviance",
       title = "Ten-fold cross-validation") +
  theme_pub(7.2)

# ---------------------------------------------------------------- C RF
set.seed(2026)
rf <- randomForest(x = X, y = as.factor(ifelse(ml$OS == 1, "Dead", "Alive")),
                   ntree = 500, importance = TRUE)
imp <- as.data.frame(importance(rf))
imp$gene <- rownames(imp)

# 已选特征集取自 Table S4（原始分析的 58 个），在本环境重算其重要性数值
s4 <- read.csv("Table_S4_Machine_Learning_Features.csv", stringsAsFactors = FALSE)
rf_sel <- s4$RandomForest_Features[s4$RandomForest_Features != ""]
lasso_sel <- s4$LASSO_Features[s4$LASSO_Features != ""]
hub_sel <- s4$Hub_Signature[s4$Hub_Signature != ""]
imp_sel <- imp[imp$gene %in% rf_sel, ]

NTOP <- 20
ord <- imp_sel[order(-imp_sel$MeanDecreaseGini), ]
top <- ord[seq_len(NTOP), ]
## 排名仍由 Gini 决定，但四个 hub 基因无论排到第几都保留 —— 否则 ANGPT1
## （accuracy 排第 4、Gini 排第 38）会被整行删掉，读者看到的只剩 MMP7。
sel <- ord[ord$gene %in% c(top$gene, HUB), ]
rank_hub <- sapply(HUB, function(g) rank(-imp_sel$MeanDecreaseGini)[match(g, imp_sel$gene)])

imp_long <- rbind(
  data.frame(gene = sel$gene, metric = "Mean decrease in accuracy", value = sel$MeanDecreaseAccuracy),
  data.frame(gene = sel$gene, metric = "Mean decrease in Gini",     value = sel$MeanDecreaseGini)
)
imp_long$metric <- factor(imp_long$metric, levels = c("Mean decrease in accuracy", "Mean decrease in Gini"))
imp_long$gene <- factor(imp_long$gene, levels = rev(sel$gene))
imp_long$hub <- imp_long$gene %in% HUB

mk_imp <- function(d, xlab, keep_y) {
  g <- ggplot(d, aes(value, gene)) +
    geom_segment(aes(x = 0, xend = value, yend = gene), colour = "grey88", linewidth = 0.3) +
    geom_point(aes(fill = hub, size = hub), shape = 21, colour = "grey35", stroke = 0.3) +
    scale_fill_manual(values = c("FALSE" = SKY, "TRUE" = ORG), guide = "none") +
    scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 2.4), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = xlab, y = NULL) +
    theme_pub(7.2) +
    theme(panel.grid.major.x = element_line(linewidth = 0.18, colour = "grey93"))
  if (keep_y) g else g + theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(), axis.line.y = element_blank())
}
pC1 <- mk_imp(imp_long[imp_long$metric == levels(imp_long$metric)[1], ],
              "Mean decrease in accuracy", TRUE)
pC2 <- mk_imp(imp_long[imp_long$metric == levels(imp_long$metric)[2], ],
              "Mean decrease in Gini", FALSE)
pC1 <- pC1 + theme(axis.text.y = element_text(size = 5.9, face = "italic"))
pC1 <- pC1 + theme(plot.margin = ggplot2::margin(7, 4, 4, 4))
pC <- (pC1 | pC2) + plot_layout(widths = c(1.06, 1))
pC <- wrap_elements(full = pC) + labs(tag = "C") +
  theme(plot.tag = element_text(face = "bold", size = 11),
        plot.tag.position = c(0, 1.0))

# ---------------------------------------------------------------- D 韦恩
circ <- function(cx, cy, r, n = 240) {
  th <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(th), y = cy + r * sin(th))
}
c1 <- circ(-0.48, 0, 0.90); c1$set <- "LASSO"
c2 <- circ( 0.48, 0, 0.90); c2$set <- "Random forest"

lasso_g <- lasso_sel; rf_g <- rf_sel
nL <- length(lasso_g); nR <- length(rf_g); nI <- length(intersect(lasso_g, rf_g))

pD <- ggplot() +
  geom_polygon(data = c1, aes(x, y, fill = set), alpha = 0.55, colour = BLU, linewidth = 0.4) +
  geom_polygon(data = c2, aes(x, y, fill = set), alpha = 0.55, colour = ORG, linewidth = 0.4) +
  scale_fill_manual(values = c("LASSO" = SKY, "Random forest" = YEL), guide = "none") +
  annotate("text", x = -1.02, y = 0, label = nL - nI, size = 3.0) +
  annotate("text", x =  0.00, y = 0, label = nI, size = 3.4, fontface = "bold") +
  annotate("text", x =  1.00, y = 0, label = nR - nI, size = 3.0) +
  annotate("text", x = -0.58, y = 1.14, label = sprintf("LASSO\n%d", nL),
           size = 2.2, colour = BLU, lineheight = 0.95, vjust = 0, fontface = "bold") +
  annotate("text", x =  0.58, y = 1.14, label = sprintf("Random forest\n%d", nR),
           size = 2.2, colour = "#B87A00", lineheight = 0.95, vjust = 0, fontface = "bold") +
  coord_equal(xlim = c(-1.52, 1.52), ylim = c(-0.96, 1.48)) +
  labs(title = "Consensus features") +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 8.5, hjust = 0,
                                  margin = ggplot2::margin(0, 0, 2, 0)),
        plot.margin = ggplot2::margin(2, 2, 2, 2))

# ---------------------------------------------------------------- 拼版
row1 <- (tag_panel(pA, "A") | tag_panel(pB, "B")) + plot_layout(widths = c(1, 1))
row2 <- pC
row3 <- (plot_spacer() | pD | plot_spacer()) + plot_layout(widths = c(0.27, 0.46, 0.27))
row3 <- wrap_elements(full = row3) + labs(tag = "D") +
  theme(plot.tag = element_text(face = "bold", size = 11),
        plot.tag.position = c(0, 1.0))

final <- row1 / row2 / row3 + plot_layout(heights = c(0.93, 1.36, 1.02))
save_fig(final, "Figure2_ABCD", W2, 7.55)

cat("\nLASSO:", nL, " RF:", nR, " 交集:", nI, " ->", paste(sort(intersect(lasso_g, rf_g)), collapse = ","), "\n")
cat("Hub rank within 58 selected:", paste(names(rank_hub), rank_hub, sep = "=", collapse = " "), "\n")
cat("S4 hub:", paste(hub_sel, collapse = ","), "\n")
