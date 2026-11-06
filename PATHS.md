# Path portability audit

Every absolute path found in the original scripts is listed below together with its replacement in this archive. **The originals under `04_代码/` and `05_数据与结果/work/` were not modified.**

Replacements use `WB_DATA_ROOT` (raw input data) and `WB_OUTPUT_DIR` (derived tables and figures); both are defined in `00_config.R` and can be overridden by environment variables.

| File | Original line | Portable line |
|---|---|---|
| `scripts/01_primary/Supplementary Tables.R` | `# setwd("G:\\Project\\nets_breast canner_1")` | `# setwd(Sys.getenv("WB_DATA_ROOT", unset = "data/raw"))` |
| `scripts/02_drug/01_prep.R` | `G   <- "G:/Project/nets_breast canner_1"` | `G   <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")` |
| `scripts/02_drug/01_prep.R` | `OUT <- "E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel"` | `OUT <- Sys.getenv("WB_OUTPUT_DIR", unset = "out")` |
| `scripts/02_drug/02_drug_panel.R` | `G   <- "G:/Project/nets_breast canner_1"` | `G   <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")` |
| `scripts/02_drug/02_drug_panel.R` | `OUT <- "E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel"` | `OUT <- Sys.getenv("WB_OUTPUT_DIR", unset = "out")` |
| `scripts/02_drug/03_immune.R` | `OUT <- "E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel"` | `OUT <- Sys.getenv("WB_OUTPUT_DIR", unset = "out")` |
| `scripts/02_drug/04_figs.R` | `setwd("E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel")` | `setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))` |
| `scripts/02_drug/05_dupcheck.R` | `G <- "G:/Project/nets_breast canner_1"` | `G <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")` |
| `scripts/02_drug/05_dupcheck.R` | `setwd("E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel")` | `setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))` |
| `scripts/02_drug/06_fix_dup.R` | `G <- "G:/Project/nets_breast canner_1"` | `G <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")` |
| `scripts/02_drug/06_fix_dup.R` | `setwd("E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel")` | `setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))` |
| `scripts/02_drug/07_corr.R` | `setwd("E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel")` | `setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))` |
| `scripts/02_drug/08_class_score.R` | `setwd("E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel")` | `setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))` |
| `scripts/02_drug/13_full_sensitivity.R` | `G <- "G:/Project/nets_breast canner_1"` | `G <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")` |
| `scripts/02_drug/13_full_sensitivity.R` | `setwd("E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel")` | `setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))` |
| `scripts/02_drug/14_residual_axis.R` | `G <- "G:/Project/nets_breast canner_1"` | `G <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")` |
| `scripts/02_drug/14_residual_axis.R` | `setwd("E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel")` | `setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))` |
| `scripts/02_drug/15_figs.R` | `setwd("E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel")` | `setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))` |
| `scripts/02_drug/16_targeted_residual.R` | `G <- "G:/Project/nets_breast canner_1"` | `G <- Sys.getenv("WB_DATA_ROOT", unset = "data/raw")` |
| `scripts/02_drug/16_targeted_residual.R` | `setwd("E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel")` | `setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))` |
| `scripts/02_drug/17_targeted_figs.R` | `setwd("E:/Work Buddy/2026-10-31-19-25-29/projects/submission-optimization/drug_panel")` | `setwd(Sys.getenv("WB_OUTPUT_DIR", unset = "out"))` |
| `scripts/04_work_scripts/s9_supp_tables.R` | `# setwd("G:\\Project\\nets_breast canner_1")` | `# setwd(Sys.getenv("WB_DATA_ROOT", unset = "data/raw"))` |

Total replacements: **22** across **15** files.

Two commented-out `setwd()` lines that used Windows backslash paths are also listed; they are inert but were rewritten for consistency.
