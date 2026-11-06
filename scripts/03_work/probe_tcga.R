load("../raw/TCGA/brca.gdc_2024.rda")
objs <- ls()
cat("OBJ:", paste(objs, collapse = ","), "\n")
for (o in objs) {
  x <- get(o)
  if (is.matrix(x) || is.data.frame(x)) {
    cat("\n--", o, "dim:", paste(dim(x), collapse = "x"), "\n")
    cat("ROW:", paste(head(rownames(x), 3), collapse = " | "), "\n")
    cat("COL:", paste(head(colnames(x), 3), collapse = " | "), "\n")
    v <- as.matrix(x[seq_len(min(80, nrow(x))), seq_len(min(30, ncol(x))), drop = FALSE])
    cat("RANGE:", paste(round(range(v, na.rm = TRUE), 2), collapse = " ~ "), "\n")
  } else if (is.character(x) || is.factor(x)) {
    cat("\n--", o, "length:", length(x), " sample:", paste(head(x, 3), collapse = " | "), "\n")
  }
}
