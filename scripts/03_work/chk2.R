cat("ggvenn:", requireNamespace("ggvenn", quietly=TRUE), "| version:", as.character(packageVersion("ggvenn")), "\n")
cat("ggVennDiagram:", requireNamespace("ggVennDiagram", quietly=TRUE), "\n")
cat("enrichKEGG 在线:", tryCatch({clusterProfiler::enrichKEGG("TP53", organism="hsa"); "ok"}, error=function(e) paste("ERR", conditionMessage(e))), "\n")
