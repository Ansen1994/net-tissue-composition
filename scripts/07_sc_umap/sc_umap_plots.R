## =============================================================================
## sc_umap_plots.R - draw Supplementary Figure S4 from the REUSABLE OUTPUTS.
## Reads the csv.gz tables only; never re-reads the 2.1 GB matrices, never
## recomputes UMAP. Safe to run from another project with WB_SC_UMAP_DIR set.
##
## Output: Supplementary_Figure_S3.pdf (2 pages, A4 portrait) in WB_SC_UMAP_DIR.
## =============================================================================

args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
here <- dirname(normalizePath(this_file))
source(file.path(here, "sc_umap_config.R"))
source(file.path(here, "sc_umap_load.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

audit_dir <- Sys.getenv("WB_SC_AUDIT_DIR", unset = "")
COL9 <- c("Cancer Epithelial" = "#D55E00", "Normal Epithelial" = "#E69F00",
          "T-cells" = "#0072B2", "B-cells" = "#CC79A7",
          "Plasmablasts" = "#7B2D8E", "Myeloid" = "#009E73",
          "CAFs" = "#56B4E9", "PVL" = "#8C6D1F", "Endothelial" = "#999999")
COL3 <- c("ER+" = "#0072B2", "HER2+" = "#D55E00", "TNBC" = "#009E73")

emb <- load_sc_umap(dir = WB_SC_UMAP_DIR, validate = FALSE)
con <- gzfile(file.path(WB_SC_UMAP_DIR, "sc_umap_features.csv.gz"), "rt")
feat <- utils::read.csv(con, check.names = FALSE, stringsAsFactors = FALSE); close(con)
audit <- NULL; marker <- NULL
if (nzchar(audit_dir)) {
  p <- file.path(audit_dir, "j3_annotation_audit.csv")
  q <- file.path(audit_dir, "j3_marker_verification.csv")
  if (file.exists(p)) audit <- utils::read.csv(p, stringsAsFactors = FALSE)
  if (file.exists(q)) marker <- utils::read.csv(q, stringsAsFactors = FALSE)
}

base_umap <- function(d) ggplot(d, aes(UMAP_1, UMAP_2)) +
  coord_equal() +
  theme_void() +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 7, colour = "grey35", hjust = 0),
        legend.text = element_text(size = 6.5), legend.title = element_blank(),
        legend.key.size = grid::unit(0.28, "cm"))

## ---------------------------------------------------------------------------
## Display sub-sample. The PDF is a vector device, so every plotted point costs
## file size; 100,064 points in each of two panels alone exceeded 20 MB. The
## embedding shape is fully resolved at 30,000 points, and the sub-sample is
## stated in the legend. Percentages are always recomputed on the full table.
set.seed(WB_SEED_UMAP)
n_show <- min(30000L, nrow(emb))
emb_show <- emb[sample.int(nrow(emb), n_show), , drop = FALSE]

pA <- base_umap(emb_show) +
  geom_point(aes(colour = major), size = 0.06, alpha = 0.75) +
  scale_colour_manual(values = COL9) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1),
                               nrow = 3, byrow = TRUE)) +
  labs(title = "A  UMAP, nine published major classes",
       subtitle = "Labels reused unchanged, no re-clustering") +
  theme(legend.position = "bottom")

pB <- base_umap(emb_show) +
  geom_point(aes(colour = subtype), size = 0.06, alpha = 0.7) +
  scale_colour_manual(values = COL3) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 1)) +
  labs(title = "B  Same embedding, clinical subtype",
       subtitle = "11 ER+ / 5 HER2+ / 10 TNBC tumours (GEO clinical_subtype)") +
  theme(legend.position = "bottom")

pC <- ggplot() +
  geom_point(data = emb_show, aes(UMAP_1, UMAP_2), colour = "grey88", size = 0.05) +
  geom_point(data = emb_show[emb_show$subtype == "ER+", ],
             aes(UMAP_1, UMAP_2), colour = COL3["ER+"], size = 0.06, alpha = 0.8) +
  coord_equal() + theme_void() +
  labs(title = "C  ER+ subset used for hub-gene attribution",
       subtitle = "38,241 cells from 11 ER+ tumours; all other cells greyed out") +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 7, colour = "grey35", hjust = 0),
        plot.margin = margin(4, 4, 4, 12))

if (!is.null(audit)) {
  d <- audit[audit$level == "major", ]
  d$label <- factor(d$label, levels = d$label[order(d$n)])
  pD <- ggplot(d, aes(n, label, fill = label)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.12, size = 2.4) +
    scale_fill_manual(values = COL9, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "D  Cell counts by major class", x = "cells", y = NULL,
         subtitle = "This annotation contains no neutrophil class") +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey35"),
          axis.text.y = element_text(size = 7))
} else pD <- ggplot() + theme_void()

if (!is.null(marker)) {
  # MKI67 is a proliferation marker, not a lineage marker: at the nine-major
  # level it peaks in Myeloid (cycling myeloid cells), so it is excluded from
  # this lineage check and remains only in the source audit table.
  marker <- marker[marker$marker != "MKI67", , drop = FALSE]
  marker$label <- factor(marker$marker, levels = marker$marker[order(marker$mean_expr)])
  pE <- ggplot(marker, aes(mean_expr, label, fill = top_major)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = top_major), hjust = -0.1, size = 2.2) +
    scale_fill_manual(values = COL9, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.6))) +
    labs(title = "E  Lineage-marker verification", x = "CP10K-log mean", y = NULL,
         subtitle = "Each marker peaks in its expected compartment (Table S11a)") +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey35"),
          axis.text.y = element_text(size = 7))
} else pE <- ggplot() + theme_void()

page1 <- ((pA | pB) / (pC | (pD / pE))) +
  plot_layout(heights = c(0.60, 1)) +
  plot_annotation(
    title = "Supplementary Figure S3  Recomputed UMAP of GSE176078",
    subtitle = "Display only: no re-clustering, no re-annotation; coordinates are not those of Wu et al. [32].")

## ------------------------------------------------------------ page 2: features
## All 12 genes in the released feature table are displayed in a 4 x 3 grid on
## A4 portrait; order follows WB_SC_FEATURE_GENES (four hub, then eight NET /
## granule genes). A single colour scale is shared by all 12 panels so that the
## panels can be compared with each other; the upper limit is set by the
## highest-expressing transcript in the set.
COLS_EXPR <- c("#4575B4", "#FFFFBF", "#D73027")
fgenes <- intersect(WB_SC_FEATURE_GENES, names(feat))
E <- merge(emb[, c("cell_id", "UMAP_1", "UMAP_2")], feat, by = "cell_id", sort = FALSE)
vmax_all <- max(vapply(fgenes, function(g) max(E[[g]], na.rm = TRUE), numeric(1)))
set.seed(WB_SEED_UMAP)
plist <- lapply(fgenes, function(g) {
  v <- E[[g]]
  bg <- which(v <= 0)
  if (length(bg) > 2500) bg <- sample(bg, 2500)
  fg <- which(v > 0)
  ## Down-sample for display only; the printed percentage uses all fg cells.
  fg_show <- fg
  if (length(fg_show) > 5000) fg_show <- sample(fg_show, 5000)
  d0 <- E[bg, c("UMAP_1", "UMAP_2")]
  d1 <- data.frame(UMAP_1 = E$UMAP_1[fg_show], UMAP_2 = E$UMAP_2[fg_show],
                   val = v[fg_show])
  n_pct <- sprintf("%.2f%% detected", 100 * length(fg) / nrow(E))
  ggplot() +
    geom_point(data = d0, aes(UMAP_1, UMAP_2), colour = "grey90", size = 0.04) +
    geom_point(data = d1, aes(UMAP_1, UMAP_2, colour = val), size = 0.16, alpha = 0.85) +
    scale_colour_gradientn(colours = COLS_EXPR,
                           limits = c(0, vmax_all), oob = scales::squish,
                           guide = "none") +
    coord_equal() + theme_void() +
    labs(title = g, subtitle = n_pct) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 6.5, colour = "grey35", hjust = 0.5),
          plot.background = element_rect(colour = "grey75", fill = NA, linewidth = 0.3))
})

## Shared colour bar for the 12 panels on page 2.
cbar <- ggplot(data.frame(x = seq(0, vmax_all, length.out = 200), y = 1), aes(x, y, fill = x)) +
  geom_tile() +
  scale_fill_gradientn(colours = COLS_EXPR, limits = c(0, vmax_all)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.008, 0.008)), breaks = scales::pretty_breaks(6)) +
  labs(x = "log1p counts per 10,000 UMIs (one scale for all 12 panels)", y = NULL, fill = NULL) +
  theme_minimal(base_size = 6.5) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), legend.position = "none",
        plot.margin = margin(1, 46, 3, 46))

page2 <- (patchwork::wrap_plots(plist, ncol = 3) / cbar) +
  plot_layout(heights = c(1, 0.035)) +
  plot_annotation(
    title = "Supplementary Figure S3 (continued)  Gene expression on the recomputed UMAP",
    subtitle = paste0(
      "Four hub genes and eight NET effector / granule transcripts, all 12 genes in the released feature table.\n",
      "One shared colour scale: a pale panel means low expression, not a different scale. Percentages are over\n",
      "all 100,064 cells, not the ER+ subset of Figure 6D; points are down-sampled for display."))

out <- file.path(WB_SC_UMAP_DIR, "Supplementary_Figure_S3.pdf")
grDevices::cairo_pdf(out, width = 8.27, height = 11.69, onefile = TRUE)
print(page1)
print(page2)
grDevices::dev.off()
cat("[plots] wrote", out, "\n")
