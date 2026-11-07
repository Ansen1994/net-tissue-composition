# ==============================================================================
# Table S12h - ESTIMATE score comparison (convergent, expression-based check)
#
# Run from 05_数据与结果/work/ :  Rscript s12h_estimate_check.R
# All paths are relative so that the script also works under LC_CTYPE=C on
# Windows, where R cannot open multibyte directory names.
#
# Inputs
#   final_data.csv                         686 TCGA-BRCA ER+ tumours, risk score
#   ../raw/drug/Breast_estimate_score.gct  ESTIMATE scores for the same 686 tumours
#                                          (estimate R package, Yoshihara 2013 [31])
#   j4_module_scores.csv                   nine generalisation modules (Methods 2.9)
#
# Output
#   Table_S12h_ESTIMATE.csv                consumed by build_supplementary_tables.py
#
# Interpretation: ESTIMATE scores are inferred from the SAME bulk expression
# matrix as the risk score. They are an algorithmically independent readout of
# stromal/immune admixture, NOT an independent histological gold standard.
# ==============================================================================

options(stringsAsFactors = FALSE)

fd <- read.csv("final_data.csv", check.names = FALSE)
gct <- readLines("../raw/drug/Breast_estimate_score.gct")
hdr <- strsplit(gct[3], "\t")[[1]]
samples <- gsub("\\.", "-", hdr[-(1:2)])
mat <- do.call(rbind, lapply(gct[-(1:3)], function(z) {
  as.numeric(strsplit(z, "\t")[[1]][-(1:2)])
}))
rownames(mat) <- vapply(gct[-(1:3)], function(z) strsplit(z, "\t")[[1]][1], "")
colnames(mat) <- samples
est <- as.data.frame(t(mat))
est$Patient_ID <- rownames(est)

j4 <- read.csv("j4_module_scores.csv", check.names = FALSE)
j4$Patient_ID <- substr(j4[[1]], 1, 12)
mods <- c("Stroma_ECM", "Epithelial", "T_cells", "Myeloid",
          "Proliferation", "B_cells")
j4s <- j4[, c("Patient_ID", mods)]
# 11 barcodes recur because some patients contributed several aliquots to the
# 1,111-tumour generalisation matrix. Collapse to one row per patient (mean of
# the aliquot z-scores) before linking to the 686-patient risk-score table.
j4s <- aggregate(j4s[, mods], by = list(Patient_ID = j4s$Patient_ID), FUN = mean)

d <- merge(fd[, c("Patient_ID", "RiskScore")], est, by = "Patient_ID")
d <- merge(d, j4s, by = "Patient_ID")
cat("Matched tumours:", nrow(d), "\n")
stopifnot(nrow(d) == 686)

sp <- function(x, y) {
  ct <- suppressWarnings(cor.test(d[[x]], d[[y]], method = "spearman",
                                  exact = TRUE))
  p <- ct$p.value
  # R's exact Spearman p underflows to 0 for |rho| > 0.7 at n = 686; report the
  # asymptotic t-approximation in those three rows (flagged in the table note).
  if (!is.finite(p) || p == 0) {
    p <- suppressWarnings(cor.test(d[[x]], d[[y]], method = "spearman",
                                   exact = FALSE))$p.value
  }
  c(rho = unname(ct$estimate), p = p, n = sum(is.finite(d[[x]]) &
                                              is.finite(d[[y]])))
}

# Rank-residual partial correlation, same construction as ia2_stats.R
pcor <- function(x, y, z) {
  ok <- is.finite(d[[x]]) & is.finite(d[[y]]) & is.finite(d[[z]])
  rx <- resid(lm(rank(d[[x]][ok]) ~ rank(d[[z]][ok])))
  ry <- resid(lm(rank(d[[y]][ok]) ~ rank(d[[z]][ok])))
  ct <- cor.test(rx, ry)
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

rows <- list()
add <- function(analysis, comparison, r, note) {
  rows[[length(rows) + 1L]] <<- data.frame(
    analysis = analysis, comparison = comparison,
    rho = round(r[["rho"]], 6), p = signif(r[["p"]], 6),
    n = as.integer(r[["n"]]), note = note, stringsAsFactors = FALSE)
}

# --- 1. primary: risk score vs ESTIMATE scores -------------------------------
for (k in c("StromalScore", "ImmuneScore", "ESTIMATEScore")) {
  add("primary", paste("RiskScore ~", k), sp("RiskScore", k),
      "Risk score vs ESTIMATE score; BH family = 7 primary/construct rows")
}

# --- 2. construct validity: do ESTIMATE scores measure composition here? -----
pairs <- list(c("StromalScore", "Stroma_ECM"),
              c("ImmuneScore", "T_cells"),
              c("ImmuneScore", "Myeloid"),
              c("ESTIMATEScore", "Epithelial"))
for (pr in pairs) {
  add("construct", paste(pr[1], "~", pr[2]), sp(pr[1], pr[2]),
      "Construct check against the nine-module panel of Methods 2.9 / Table S13e")
}

# --- 3. partial: does controlling for stromal content remove the immune link? -
for (k in c("B_cells", "T_cells")) {
  raw <- sp("RiskScore", k)
  pr <- pcor("RiskScore", k, "StromalScore")
  add("partial", paste("RiskScore ~", k, "| StromalScore"), pr,
      sprintf("Rank-residual partial correlation; raw rho = %.3f; sensitivity only, excluded from BH",
              raw[["rho"]]))
}

out <- do.call(rbind, rows)

# Benjamini-Hochberg within the 7 primary + construct tests only
bh <- rep(NA_real_, nrow(out))
idx <- out$analysis %in% c("primary", "construct")
bh[idx] <- p.adjust(out$p[idx], method = "BH")
out$p_BH <- ifelse(is.na(bh), "", signif(bh, 6))
out <- out[, c("analysis", "comparison", "rho", "p", "n", "p_BH", "note")]

write.csv(out, "Table_S12h_ESTIMATE.csv", row.names = FALSE, na = "")
cat("\n")
print(out[, 1:6], row.names = FALSE)
cat("\n[written] Table_S12h_ESTIMATE.csv\n")
