## =============================================================================
## sc_umap_config.R - central configuration for the C-tier UMAP recomputation.
##
## Source this from sc_umap_pipeline.R / sc_umap_plots.R, or override any value
## with an environment variable. There are no machine-specific paths here.
##
## WB_SC_RAW_DIR    directory with one sub-folder per tumour (metadata.csv,
##                  count_matrix_genes.tsv, count_matrix_barcodes.tsv,
##                  count_matrix_sparse.mtx)
## WB_SC_UMAP_DIR   output directory for all reusable artifacts
## =============================================================================

WB_SC_RAW_DIR  <- Sys.getenv("WB_SC_RAW_DIR",  unset = "data/raw/scRNA/sc")
WB_SC_UMAP_DIR <- Sys.getenv("WB_SC_UMAP_DIR", unset = "out/sc_umap")

WB_SEED_UMAP       <- as.integer(Sys.getenv("WB_SEED_UMAP",      unset = "20261106"))
WB_SC_HVG_N        <- as.integer(Sys.getenv("WB_SC_HVG_N",       unset = "2000"))
WB_SC_NPCS         <- as.integer(Sys.getenv("WB_SC_NPCS",        unset = "30"))
WB_SC_NNEIGHBORS   <- as.integer(Sys.getenv("WB_SC_NNEIGHBORS",  unset = "30"))
WB_SC_MINDIST      <- as.numeric(Sys.getenv("WB_SC_MINDIST",     unset = "0.3"))
WB_SC_SAMPLE_LIMIT <- as.integer(Sys.getenv("WB_SC_SAMPLE_LIMIT", unset = "0"))
WB_SC_SKIP_HARMONY <- Sys.getenv("WB_SC_SKIP_HARMONY", unset = "0") == "1"

WB_SC_FEATURE_GENES <- c("PLG", "CD79A", "MMP7", "ANGPT1",
                         "PADI4", "ELANE", "MPO", "GSDMD",
                         "FCGR3B", "CXCR2", "S100A8", "S100A9")
WB_SC_REDUCTION         <- "harmony"
WB_SC_PIPELINE_VERSION  <- "1.0.0"
WB_SC_CACHE_KEY_VERSION <- 1L

dir.create(WB_SC_UMAP_DIR, showWarnings = FALSE, recursive = TRUE)
