## =============================================================================
## sc_umap_load.R - load-only helper. NEVER recomputes anything.
##
##   source("sc_umap_load.R")
##   emb <- load_sc_umap(subtype = "ER+")              # 38,241 rows
##   emb <- load_sc_umap(with_features = TRUE)         # + gene expression columns
##
## Set WB_SC_UMAP_DIR, or pass dir=, to point at this folder from another
## project. Returns a data.frame; the manifest is attached as attr(,"manifest").
## =============================================================================

load_sc_umap <- function(dir = Sys.getenv("WB_SC_UMAP_DIR", unset = "out/sc_umap"),
                         subtype = NULL, with_features = FALSE, validate = TRUE) {
  ## R under the C locale cannot pass a non-ASCII absolute path to gzfile().
  ## setwd() accepts it, so read via relative filenames and restore the cwd.
  old_wd <- getwd()
  if (inherits(try(setwd(dir), silent = TRUE), "try-error"))
    stop("Cannot access sc_umap directory: ", dir,
         "\nSet WB_SC_UMAP_DIR or pass dir= to the folder holding the outputs.")
  on.exit(setwd(old_wd), add = TRUE)
  mp <- "sc_umap_manifest.json"
  ep <- "sc_umap_embeddings.csv.gz"
  fp <- "sc_umap_features.csv.gz"
  if (!file.exists(ep)) stop("Embeddings not found in: ", dir,
                             "\nSet WB_SC_UMAP_DIR or pass dir= to the sc_umap folder.")
  con <- gzfile(ep, "rt"); emb <- utils::read.csv(con, stringsAsFactors = FALSE,
                                                 check.names = FALSE); close(con)

  man <- NULL
  if (file.exists(mp)) {
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      man <- jsonlite::fromJSON(mp)
    } else {
      txt <- paste(readLines(mp, warn = FALSE), collapse = "\n")
      man <- list(cache_key = sub(".*\"cache_key\"\\s*:\\s*\"([^\"]+)\".*", "\\1", txt))
      warning("jsonlite not installed: manifest parsed as plain text (cache_key only).")
    }
  }

  if (validate) {
    need <- c("cell_id", "sample", "subtype", "major", "minor",
              "UMAP_1", "UMAP_2")
    miss <- setdiff(need, names(emb))
    if (length(miss)) stop("Malformed embeddings table, missing: ", paste(miss, collapse = ", "))
    if (!is.null(man) && !is.null(man$n_cells) && nrow(emb) != man$n_cells)
      warning("Row count ", nrow(emb), " differs from manifest n_cells ", man$n_cells)
    cat(sprintf("[load_sc_umap] %d cells, %d samples, %d major classes; cache_key=%s\n",
                nrow(emb), length(unique(emb$sample)), length(unique(emb$major)),
                if (!is.null(man$cache_key)) substr(man$cache_key, 1, 12) else "n/a"))
  }

  if (!is.null(subtype)) emb <- emb[emb$subtype %in% subtype, , drop = FALSE]
  if (with_features) {
    if (!file.exists(fp)) stop("Feature table not found: ", fp)
    con <- gzfile(fp, "rt"); fe <- utils::read.csv(con, stringsAsFactors = FALSE,
                                                   check.names = FALSE); close(con)
    emb <- merge(emb, fe, by = "cell_id", all.x = TRUE, sort = FALSE)
  }
  attr(emb, "manifest") <- man
  emb
}
