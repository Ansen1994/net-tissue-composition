## Central configuration for the portable code archive.
## Set WB_DATA_ROOT and WB_OUTPUT_DIR as environment variables before running,
## or edit the defaults below. No script in this archive contains a machine path.

WB_DATA_ROOT  <- Sys.getenv("WB_DATA_ROOT",  unset = "data/raw")
WB_OUTPUT_DIR <- Sys.getenv("WB_OUTPUT_DIR", unset = "out")

## Fixed random seeds used in the manuscript.
WB_SEED_FEATURE_SELECTION <- 2026L  # LASSO cv.glmnet and randomForest (Methods 2.3)
WB_SEED_CROSS_VALIDATION  <- 42L    # 5-fold CV optimism, module scores, drug panel

dir.create(WB_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
