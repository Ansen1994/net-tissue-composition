suppressMessages({ library(SummarizedExperiment) })
load("../raw/TCGA/brca.gdc_2024.rda")
x <- expquery2
cat("class:", class(x), "\n")
cat("dim:", paste(dim(x), collapse = "x"), "\n")
cat("assays:", paste(assayNames(x), collapse = ", "), "\n")
cat("rowData cols:", paste(head(names(rowData(x)), 10), collapse = ", "), "\n")
cat("colData cols:", paste(head(names(colData(x)), 12), collapse = ", "), "\n")
cat("rowname示例:", paste(head(rownames(x), 3), collapse = " | "), "\n")
cat("colname示例:", paste(head(colnames(x), 3), collapse = " | "), "\n")
a <- assay(x, 1)
cat("assay1 class:", class(a), " range:",
    paste(round(range(as.matrix(a[seq_len(min(200, nrow(a))), seq_len(min(50, ncol(a)))]), na.rm = TRUE), 2), collapse = " ~ "), "\n")
if ("unstranded" %in% assayNames(x)) {
  b <- assay(x, "unstranded")
  cat("unstranded range:",
      paste(round(range(as.matrix(b[seq_len(min(200, nrow(b))), seq_len(min(50, ncol(b)))]), na.rm = TRUE), 2), collapse = " ~ "), "\n")
}
cat("\n--- rowData 前3行 ---\n")
print(head(as.data.frame(rowData(x)), 3))
cat("\n--- colData 前2行(前8列) ---\n")
print(head(as.data.frame(colData(x))[, 1:min(8, ncol(colData(x)))]), 2)
