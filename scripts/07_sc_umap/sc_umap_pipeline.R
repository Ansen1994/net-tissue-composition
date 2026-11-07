## =============================================================================
## sc_umap_pipeline.R - recompute a display UMAP embedding for GSE176078.
##
## IMPORTANT SCIENTIFIC SCOPE
##   This script recomputes an embedding from the deposited count matrices for
##   VISUALISATION ONLY. It does NOT re-cluster cells and does NOT re-annotate
##   them: cell types remain the curated labels published by Wu et al. (2021).
##   Coordinates are therefore NOT the original atlas coordinates.
##
## Input : 26 tumour folders under WB_SC_RAW_DIR (GEO GSE176078 processed files).
## Output: sc_umap_embeddings.csv.gz, sc_umap_features.csv.gz,
##         sc_umap_object.rds, sc_umap_manifest.json,
##         sessionInfo_sc_umap.txt  (all in WB_SC_UMAP_DIR).
##
## Caching: a SHA-256 cache key is derived from the parameter set, package
##   versions and an input-file fingerprint. If the manifest and all outputs
##   are present and the key matches, the script exits without recomputing.
##   Set WB_SC_FORCE=1 to force a rebuild.
## =============================================================================

args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
here <- dirname(normalizePath(this_file))
source(file.path(here, "sc_umap_config.R"))

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(harmony)
  library(jsonlite)
  library(digest)
})

msg <- function(...) cat(sprintf(...), "\n", sep = "")

## ---------------------------------------------------------------- input scan
dirs <- list.dirs(WB_SC_RAW_DIR, recursive = FALSE, full.names = TRUE)
dirs <- dirs[file.exists(file.path(dirs, "metadata.csv")) &
             file.exists(file.path(dirs, "count_matrix_sparse.mtx"))]
dirs <- sort(dirs)
if (WB_SC_SAMPLE_LIMIT > 0) dirs <- head(dirs, WB_SC_SAMPLE_LIMIT)
if (!length(dirs)) stop("No tumour folders found under ", WB_SC_RAW_DIR)
samples <- basename(dirs)
msg("[info] %d tumour folders, first=%s last=%s", length(dirs), samples[1], tail(samples, 1))

fingerprint_files <- unlist(lapply(dirs, function(d) {
  f <- file.path(d, c("metadata.csv", "count_matrix_genes.tsv",
                      "count_matrix_barcodes.tsv", "count_matrix_sparse.mtx"))
  fi <- file.info(f)
  paste(basename(d), basename(f), fi$size, as.integer(fi$mtime), sep = ":")
}))
input_fp <- digest::digest(paste(fingerprint_files, collapse = "|"), algo = "sha256")

pkg_ver <- function(p) as.character(packageVersion(p))
versions <- list(Seurat = pkg_ver("Seurat"), SeuratObject = pkg_ver("SeuratObject"),
                 harmony = pkg_ver("harmony"), uwot = pkg_ver("uwot"),
                 Matrix = pkg_ver("Matrix"), R = as.character(getRversion()))

params <- list(seed = WB_SEED_UMAP, hvg_n = WB_SC_HVG_N, npcs = WB_SC_NPCS,
               n_neighbors = WB_SC_NNEIGHBORS, min_dist = WB_SC_MINDIST,
               skip_harmony = WB_SC_SKIP_HARMONY, sample_limit = WB_SC_SAMPLE_LIMIT,
               reduction = WB_SC_REDUCTION, pipeline_version = WB_SC_PIPELINE_VERSION,
               cache_key_version = WB_SC_CACHE_KEY_VERSION)

cache_key <- digest::digest(list(params = params, versions = versions,
                                 input_fingerprint = input_fp,
                                 genes = WB_SC_FEATURE_GENES), algo = "sha256")

out_names <- c("sc_umap_embeddings.csv.gz", "sc_umap_features.csv.gz",
               "sc_umap_object.rds", "sc_umap_manifest.json",
               "sessionInfo_sc_umap.txt")
manifest_path <- file.path(WB_SC_UMAP_DIR, "sc_umap_manifest.json")
all_present <- all(file.exists(file.path(WB_SC_UMAP_DIR, out_names)))
if (all_present && file.exists(manifest_path) &&
    identical(jsonlite::fromJSON(manifest_path)$cache_key, cache_key) &&
    Sys.getenv("WB_SC_FORCE", "0") != "1") {
  msg("[cache hit] key=%s ; all 5 outputs present, nothing recomputed.", substr(cache_key, 1, 12))
  quit(save = "no", status = 0)
}
msg("[info] cache miss/force, key=%s", substr(cache_key, 1, 12))

## ------------------------------------------------------------ read 10x matrices
mats <- list(); metas <- list()
for (d in dirs) {
  sid <- basename(d)
  m <- Matrix::readMM(file.path(d, "count_matrix_sparse.mtx"))
  m <- as(m, "CsparseMatrix")
  gn <- readLines(file.path(d, "count_matrix_genes.tsv"), warn = FALSE)
  bc <- readLines(file.path(d, "count_matrix_barcodes.tsv"), warn = FALSE)
  if (nrow(m) != length(gn)) stop("gene dimension mismatch in ", sid)
  if (ncol(m) != length(bc)) stop("barcode dimension mismatch in ", sid)
  rownames(m) <- make.unique(gn)
  colnames(m) <- bc
  md <- read.csv(file.path(d, "metadata.csv"), row.names = 1,
                 stringsAsFactors = FALSE, check.names = FALSE)
  idx <- match(colnames(m), rownames(md))
  if (anyNA(idx)) stop("barcode/metadata mismatch in ", sid, ": ",
                       sum(is.na(idx)), " unmatched cells")
  md <- md[idx, , drop = FALSE]
  md$sample <- sid
  mats[[sid]] <- m
  metas[[sid]] <- md
  msg("  read %-10s cells=%-6d genes=%d nnz=%d", sid, ncol(m), nrow(m), length(m@x))
}

## ------------------------------------------------------------- merge + Seurat
cell_ids <- unlist(lapply(mats, colnames), use.names = FALSE)
counts <- do.call(cbind, mats)
colnames(counts) <- cell_ids
rm(mats); invisible(gc(FALSE))
mr <- do.call(rbind, metas)
rm(metas); invisible(gc(FALSE))
stopifnot(ncol(counts) == nrow(mr), length(cell_ids) == ncol(counts))
meta <- data.frame(cell_id = cell_ids, sample = mr$sample,
                   subtype = mr$subtype, major = mr$celltype_major,
                   minor = mr$celltype_minor, subset = mr$celltype_subset,
                   nCount_RNA = mr$nCount_RNA, nFeature_RNA = mr$nFeature_RNA,
                   percent_mito = mr$percent.mito, stringsAsFactors = FALSE)
rownames(meta) <- cell_ids
stopifnot(identical(colnames(counts), rownames(meta)))
msg("[info] merged matrix %d genes x %d cells, nnz=%d",
    nrow(counts), ncol(counts), length(counts@x))

obj <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta,
                                        min.cells = 0, min.features = 0,
                                        project = "GSE176078_UMAP")
stopifnot(ncol(obj) == nrow(meta))
rm(counts); invisible(gc(FALSE))

obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize",
                             scale.factor = 1e4, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst",
                                    nfeatures = WB_SC_HVG_N, verbose = FALSE)
obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(obj),
                         verbose = FALSE)
obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(obj),
                      npcs = WB_SC_NPCS, verbose = FALSE)
if (!WB_SC_SKIP_HARMONY) {
  obj <- harmony::RunHarmony(obj, group.by.vars = "sample",
                             reduction.use = "pca", verbose = FALSE)
  red <- "harmony"
} else {
  red <- "pca"
}
obj <- Seurat::RunUMAP(obj, reduction = red, dims = seq_len(WB_SC_NPCS),
                       seed.use = WB_SEED_UMAP, n.neighbors = WB_SC_NNEIGHBORS,
                       min.dist = WB_SC_MINDIST, verbose = FALSE)
msg("[info] UMAP done on reduction '%s'", red)

## ------------------------------------------------------------------- outputs
emb <- as.data.frame(Seurat::Embeddings(obj, "umap"))
names(emb) <- c("UMAP_1", "UMAP_2")
embeddings <- data.frame(cell_id = colnames(obj),
                         sample = obj$sample, subtype = obj$subtype,
                         major = obj$major, minor = obj$minor,
                         subset = obj$subset,
                         nCount_RNA = obj$nCount_RNA,
                         nFeature_RNA = obj$nFeature_RNA,
                         percent_mito = obj$percent_mito,
                         UMAP_1 = emb$UMAP_1, UMAP_2 = emb$UMAP_2,
                         stringsAsFactors = FALSE)

feat_genes <- intersect(WB_SC_FEATURE_GENES, rownames(obj))
dat <- Seurat::GetAssayData(obj, layer = "data", assay = "RNA")
feat <- as.matrix(t(as.matrix(dat[feat_genes, , drop = FALSE])))
features <- data.frame(cell_id = rownames(feat), feat,
                       check.names = FALSE, stringsAsFactors = FALSE)

con <- gzfile(file.path(WB_SC_UMAP_DIR, "sc_umap_embeddings.csv.gz"), "wt")
utils::write.csv(embeddings, con, row.names = FALSE); close(con)
con <- gzfile(file.path(WB_SC_UMAP_DIR, "sc_umap_features.csv.gz"), "wt")
utils::write.csv(features, con, row.names = FALSE); close(con)

light <- list(embeddings = as.matrix(emb), meta = meta,
              feature_expr = feat, hvg = Seurat::VariableFeatures(obj),
              pca_loadings = Seurat::Loadings(obj, "pca"),
              reduction_used = red,
              note = "Lightweight R object: reductions/meta/selected genes only; NOT a full Seurat object.")
saveRDS(light, file.path(WB_SC_UMAP_DIR, "sc_umap_object.rds"), compress = TRUE)
writeLines(utils::capture.output(utils::sessionInfo()),
           file.path(WB_SC_UMAP_DIR, "sessionInfo_sc_umap.txt"))

outputs <- lapply(out_names, function(f) {
  p <- file.path(WB_SC_UMAP_DIR, f)
  list(file = f, bytes = if (file.exists(p)) as.numeric(file.size(p)) else NA_real_)
})
names(outputs) <- out_names

manifest <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  dataset = "GSE176078 (Wu et al. Nat Genet 2021;53:1334-1347)",
  scope = "Display UMAP recomputed by the present authors; curated cell-type labels reused unchanged; NO de novo clustering or re-annotation; coordinates are NOT the original atlas coordinates.",
  n_samples = length(dirs), n_cells = ncol(obj), n_genes = nrow(obj),
  samples = samples, params = params, package_versions = versions,
  input_fingerprint_method = "sha256 over basename:file:size:mtime of the 4 input files per tumour (not a content hash)",
  input_fingerprint = input_fp, cache_key = cache_key,
  feature_genes = feat_genes, outputs = outputs,
  validation = list(major_classes = length(unique(embeddings$major)),
                    minor_classes = length(unique(embeddings$minor)),
                    subtypes = as.list(table(embeddings$subtype)),
                    neutrophil_labels = sum(grepl("neutro|granulo|pmn",
                                                  unique(c(embeddings$major, embeddings$minor)),
                                                  ignore.case = TRUE))))
writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, digits = 6),
           manifest_path)
msg("[done] %d cells written; manifest cache_key=%s", nrow(embeddings), substr(cache_key, 1, 12))
