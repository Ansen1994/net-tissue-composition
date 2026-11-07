# Script manifest

Generated 2026-11-07 14:47 from `04_代码/` and `05_数据与结果/work/`.

Directories are ordered by analysis stage. `图表脚本对照表.md` maps each figure and table to the script that produced it.

### `scripts/01_primary/` (11 files)

Original analysis and figure scripts (DESeq2, LASSO, randomForest, Cox, KM).

- `scripts/01_primary/4_clinical.R`
- `scripts/01_primary/4_metabric.R`
- `scripts/01_primary/Supplementary Tables.R` (paths rewritten)
- `scripts/01_primary/figure1.R`
- `scripts/01_primary/figure1_extract_degs.R`; renamed from `figure1_提取degs.R`
- `scripts/01_primary/figure2.R`
- `scripts/01_primary/figure3.R`
- `scripts/01_primary/figure4.R`
- `scripts/01_primary/figure4_GSE96058.R`
- `scripts/01_primary/figure5.R`
- `scripts/01_primary/figure6.R`

### `scripts/02_drug/` (13 files)

Exploratory drug-sensitivity analysis (GDSC2 ridge regression, Table S17).

- `scripts/02_drug/01_prep.R` (paths rewritten)
- `scripts/02_drug/02_drug_panel.R` (paths rewritten)
- `scripts/02_drug/03_immune.R` (paths rewritten)
- `scripts/02_drug/04_figs.R` (paths rewritten)
- `scripts/02_drug/05_dupcheck.R` (paths rewritten)
- `scripts/02_drug/06_fix_dup.R` (paths rewritten)
- `scripts/02_drug/07_corr.R` (paths rewritten)
- `scripts/02_drug/08_class_score.R` (paths rewritten)
- `scripts/02_drug/13_full_sensitivity.R` (paths rewritten)
- `scripts/02_drug/14_residual_axis.R` (paths rewritten)
- `scripts/02_drug/15_figs.R` (paths rewritten)
- `scripts/02_drug/16_targeted_residual.R` (paths rewritten)
- `scripts/02_drug/17_targeted_figs.R` (paths rewritten)

### `scripts/03_work/` (41 files)

Re-evaluation analyses: cross-platform scaling, reference NET score, treatment interaction, immune-axis decomposition, single-cell origin, generalisation set.

- `scripts/03_work/chk2.R`
- `scripts/03_work/chk_cutpoint.R`
- `scripts/03_work/diag_auc_by_variable.R`
- `scripts/03_work/diag_auc_check.R`
- `scripts/03_work/diag_gse.R`
- `scripts/03_work/diag_gse2.R`
- `scripts/03_work/diag_incremental_value.R`
- `scripts/03_work/diag_narrative.R`
- `scripts/03_work/diag_narrative_gse96058.R`
- `scripts/03_work/diag_nets_core.R`
- `scripts/03_work/diag_pathway_score.R`
- `scripts/03_work/diag_platform_scale.R`
- `scripts/03_work/diag_poolbias.R`
- `scripts/03_work/diag_rescore_compare.R`
- `scripts/03_work/diag_riskscore.R`
- `scripts/03_work/diag_scale_check.R`
- `scripts/03_work/diag_training_optimism.R`
- `scripts/03_work/diag_verify686.R`
- `scripts/03_work/fig1_build_enrich.R`
- `scripts/03_work/fig1d_build.R`
- `scripts/03_work/fig1d_cellorigin.R`
- `scripts/03_work/ia1_build_scores.R`
- `scripts/03_work/ia2_stats.R`
- `scripts/03_work/ia3_diagnose_gse96058.R`
- `scripts/03_work/ia4_decompose.R`
- `scripts/03_work/j1_nets_prognosis.R`
- `scripts/03_work/j2_treatment_interaction.R`
- `scripts/03_work/j2b_confounding.R`
- `scripts/03_work/j3_singlecell.R`
- `scripts/03_work/j4_generalize.R`
- `scripts/03_work/j4b_generalize.R`
- `scripts/03_work/j4c_generalize.R`
- `scripts/03_work/make_table_s13f.R`
- `scripts/03_work/make_table_s15b.R`
- `scripts/03_work/make_table_s5_ci.R`
- `scripts/03_work/make_table_s5c.R`
- `scripts/03_work/probe_tcga.R`
- `scripts/03_work/probe_tcga2.R`
- `scripts/03_work/s12h_estimate_check.R`
- `scripts/03_work/t_s10_incremental.R`
- `scripts/03_work/t_s9_scheme_table.R`

### `scripts/04_work_scripts/` (9 files)

Full re-run copies of the primary scripts, kept for the byte-for-byte verification.

- `scripts/04_work_scripts/s4_clinical.R`
- `scripts/04_work_scripts/s5_metabric.R`
- `scripts/04_work_scripts/s6_gse96058.R`
- `scripts/04_work_scripts/s7_figure5.R`
- `scripts/04_work_scripts/s8_figure6.R`
- `scripts/04_work_scripts/s9_supp_tables.R` (paths rewritten)
- `scripts/04_work_scripts/s_figure1.R`
- `scripts/04_work_scripts/s_figure2.R`
- `scripts/04_work_scripts/s_figure3.R`

### `scripts/05_figures_v2/` (9 files)

Final v2 figure code (unified Okabe-Ito palette and typography).

- `scripts/05_figures_v2/fig1.R`
- `scripts/05_figures_v2/fig2.R`
- `scripts/05_figures_v2/fig3.R`
- `scripts/05_figures_v2/fig3f.R`
- `scripts/05_figures_v2/fig4.R`
- `scripts/05_figures_v2/fig5.R`
- `scripts/05_figures_v2/fig6.R`
- `scripts/05_figures_v2/fig7.R`
- `scripts/05_figures_v2/style.R`

### `scripts/06_python/` (12 files)

Reproducibility helpers: supplementary workbook, manuscript package, reference audit, review-fix patches, TIFF conversion.

- `scripts/06_python/apply_ref_fixes.py`
- `scripts/06_python/apply_review_fixes.py`
- `scripts/06_python/build_code_archive.py`
- `scripts/06_python/build_manuscript_package.py`
- `scripts/06_python/build_supplementary_tables.py`
- `scripts/06_python/count_words.py`
- `scripts/06_python/patch.py`
- `scripts/06_python/patch_refs.py`
- `scripts/06_python/patch_supp_s17.py`
- `scripts/06_python/verify_refs.py`
- `scripts/06_python/verify_specific.py`
- `scripts/06_python/make_tiff.py` (600 dpi TIFF conversion)

### `scripts/07_sc_umap/` (5 files)

Reusable single-cell asset: recomputes and plots the display-only UMAP embedding of Supplementary Figure S3 (standalone; enters no statistical model).

- `scripts/07_sc_umap/sc_umap_config.R`
- `scripts/07_sc_umap/sc_umap_load.R`
- `scripts/07_sc_umap/sc_umap_pipeline.R`
- `scripts/07_sc_umap/sc_umap_plots.R`
- `scripts/07_sc_umap/sc_umap_load.py` (read-only loader for the cached UMAP embeddings)

## Verification

- Portable R scripts: **87** parsed with `parse()`; failures: **none**.
- Portable Python helpers: **13**.
- Absolute paths rewritten: **22** (see `PATHS.md`).

Parsing checks syntax only. A full numerical re-run additionally requires the public input data listed in `README.md`; that step was not re-executed while building this archive.

## Filename normalisation

R running under the `C` locale cannot open filenames containing multibyte characters. The following copies were renamed (contents unchanged):

- `scripts/01_primary/figure1_提取degs.R` → `scripts/01_primary/figure1_extract_degs.R`
