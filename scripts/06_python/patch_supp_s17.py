# -*- coding: utf-8 -*-
"""
Add Table S17a-S17f (exploratory drug-sensitivity results) to the consolidated
supplementary workbook, and rename the workbook S1-S16 -> S1-S17 everywhere.

Pure exact-string patching; every replacement asserts a unique match.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, "build_supplementary_tables.py")
PKG = os.path.join(HERE, "build_manuscript_package.py")

S17_SPECS = '''    ("Table S17a", "S17a", "Predicted log IC50 differences for 29 compounds (primary analysis)",
     ["../\\u836f\\u654f\\u91cd\\u7b97/drug_results.csv"], "2.10, 3.6",
     "Ridge regression (glmnet, alpha = 0, 5-fold CV) trained on GDSC2. diff = median(high-risk) - median(low-risk); fdr is Benjamini-Hochberg across all 29 compounds. Only sapitinib (FDR 9.9e-09) and fulvestrant (FDR 0.017) survive FDR < 0.05."),
    ("Table S17b", "S17b", "Drug-sensitivity predictions under five preprocessing strategies",
     ["../\\u836f\\u654f\\u91cd\\u7b97/\\u654f\\u611f\\u6027/sensitivity_matrix.csv"], "2.10, 3.6",
     "29 compounds x 5 input modes (counts, log2tpm, zscore, quantile-normalised, rank). outside_frac = fraction of predictions outside the GDSC training response range; cv_r2 = cross-validation R2. Absolute effect sizes are not comparable across modes."),
    ("Table S17c", "S17c", "Robustness of direction across the five preprocessing strategies",
     ["../\\u836f\\u654f\\u91cd\\u7b97/\\u654f\\u611f\\u6027/sensitivity_verdict.csv"], "2.10, 3.6",
     "n_sig = modes with p < 0.05; robust requires >= 4 significant modes and a consistent direction. Fulvestrant is directionally unstable before common-axis adjustment."),
    ("Table S17d", "S17d", "Leave-one-out common-axis adjustment for 11 targeted compounds",
     ["../\\u836f\\u654f\\u91cd\\u7b97/\\u654f\\u611f\\u6027/targeted_verdict.csv"], "2.10, 3.6",
     "After removing the shared prediction axis (leave-one-out over 11 targeted drugs), only sapitinib, lapatinib and fulvestrant retain a consistent direction across all five modes (ROBUST-SIGNAL); tamoxifen is ROBUST-NULL."),
    ("Table S17e", "S17e", "Per-mode residual effect sizes underlying Table S17d",
     ["../\\u836f\\u654f\\u91cd\\u7b97/\\u654f\\u611f\\u6027/targeted_residual.csv"], "2.10, 3.6",
     "d_res11 = Cohen's d after LOO-11 common-axis removal; d_res29 after LOO-29 removal. n_axis11 / n_axis29 give the number of compounds contributing to the common axis."),
    ("Table S17f", "S17f", "Sensitivity of the overall drug panel to input preprocessing",
     ["../\\u836f\\u654f\\u91cd\\u7b97/\\u654f\\u611f\\u6027/sensitivity_by_mode.csv"], "2.10, 3.6",
     "n_sig = compounds with p < 0.05 under each mode; rho_avg_risk = mean correlation between predicted IC50 and risk score; pc1_pct = variance explained by the first principal component (the common axis)."),
'''.replace("\\u836f\\u654f\\u91cd\\u7b97", "药敏重算").replace("\\u654f\\u611f\\u6027", "敏感性")

S16_TAIL = '''    ("Table S16", "S16", "Sensitivity of Kaplan-Meier results to the dichotomisation rule",
     ["Table_S16_Cutpoint_Sensitivity.csv"], "2.4, 3.4",
     "Median cut-point (pre-specified) vs maximally selected rank statistic. External log-rank significance depends on the rule."),
]
'''

PATCHES_SUPP = [
    ("Output: 03_补充材料/Supplementary_Tables_S1-S16.xlsx",
     "Output: 03_补充材料/Supplementary_Tables_S1-S17.xlsx"),
    ("(S1-S16). S1-S7 keep the legacy file prefixes; S8-S16 follow 改稿说明_v2.md §5.",
     "(S1-S17). S1-S7 keep the legacy file prefixes; S8-S16 follow 改稿说明_v2.md §5; S17a-S17f are the exploratory drug-sensitivity results."),
    ('XLSX = os.path.join(OUT_DIR, "Supplementary_Tables_S1-S16.xlsx")',
     'XLSX = os.path.join(OUT_DIR, "Supplementary_Tables_S1-S17.xlsx")'),
    (S16_TAIL, S16_TAIL[:-2] + S17_SPECS + "]\n"),
    ('"1. Table IDs follow the citations in the v2 manuscript (S1-S16). Lettered sheets are sub-panels of the same table number.",',
     '"1. Table IDs follow the citations in the v2 manuscript (S1-S17). Lettered sheets are sub-panels of the same table number.",'),
    ('"7. DEG count: 13,082 rows meet padj < 0.05 and |log2FC| > 1; they collapse to 12,984 unique HGNC symbols (98 duplicated Ensembl mappings). The manuscript reports 12,984, i.e. the symbol-level count used by intersect().",',
     '"7. DEG count: 13,082 rows meet padj < 0.05 and |log2FC| > 1; they collapse to 12,984 unique HGNC symbols (98 duplicated Ensembl mappings). The manuscript reports 12,984, i.e. the symbol-level count used by intersect().",\n'
     '        "8. Table S17a-S17f are the complete exploratory drug-sensitivity results (Methods 2.10, Results 3.6). Predictions come from a direct glmnet ridge regression, not the pRRophetic package; only two of 29 compounds survive FDR in the primary analysis and only three remain directionally stable after common-axis adjustment. These tables are supplied for completeness and are not used to support any pharmacological claim.",'),
    ('"**Exploratory analyses.** Drug-sensitivity predictions (Methods 2.10) are hypothesis-generating only and are not included in the consolidated tables.",',
     '"**Exploratory analyses.** The complete drug-sensitivity predictions behind Methods 2.10 and Results 3.6 are now supplied in Table S17a-S17f. They are hypothesis-generating only: only sapitinib and fulvestrant survive FDR correction in the primary analysis, and the direction of most compounds changes with preprocessing or disappears after removal of the shared prediction axis. No pharmacological conclusion is drawn from them.",'),
    ('"Generated %s. All supplementary tables are consolidated in `Supplementary_Tables_S1-S16.xlsx` (one worksheet per table / sub-panel)." % datetime.now().strftime("%Y-%m-%d")',
     '"Generated %s. All supplementary tables are consolidated in `Supplementary_Tables_S1-S17.xlsx` (one worksheet per table / sub-panel)." % datetime.now().strftime("%Y-%m-%d")'),
    ('"Tables are numbered S1-S16 exactly as cited in the manuscript. Lettered suffixes (for example S6a-S6d, S13a-S13f) denote sub-panels of the same table number and correspond to separate worksheets in the workbook."',
     '"Tables are numbered S1-S17 exactly as cited in the manuscript. Lettered suffixes (for example S6a-S6d, S13a-S13f, S17a-S17f) denote sub-panels of the same table number and correspond to separate worksheets in the workbook."'),
    ('for f in ["Supplementary_Tables_S1-S16.xlsx", "Supplementary_Materials_v2.md",',
     'for f in ["Supplementary_Tables_S1-S17.xlsx", "Supplementary_Materials_v2.md",'),
    ('"Deposit the analysis scripts in a public repository and mint a DOI (Zenodo) for the data-availability statement.",',
     '"Insert the public repository URL, release tag and Zenodo DOI into the [REPOSITORY URL] / [VERSION] / [DOI] placeholders in the manuscript Data availability section and the Title page Code availability section; the portable code archive itself is ready in 11_代码公开包/.",'),
]

PATCHES_PKG = [
    ('    "Supplementary_Tables_S1-S16.xlsx",',
     '    "Supplementary_Tables_S1-S17.xlsx",'),
    ("│   ├── Supplementary_Tables_S1-S16.xlsx   数据表合订（39 工作表）",
     "│   ├── Supplementary_Tables_S1-S17.xlsx   数据表合订（45 工作表）"),
    ("| 补充 | `Supplementary/Supplementary_Tables_S1-S16.xlsx` | `03_补充材料/` | 39 个工作表 = README + 38 张表（S1–S16 及子面板） |",
     "| 补充 | `Supplementary/Supplementary_Tables_S1-S17.xlsx` | `03_补充材料/` | 45 个工作表 = README + 44 张表（S1–S17 及子面板，含 S17a–f 药敏全表） |"),
    ("7 份 TIFF / 7 份 PDF 齐全、补充工作簿 39 个工作表、旧硬伤表述零残留。",
     "7 份 TIFF / 7 份 PDF 齐全、补充工作簿 45 个工作表、旧硬伤表述零残留。"),
]


def apply(path, patches):
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    print("修补 %s" % os.path.basename(path))
    for old, new in patches:
        n = text.count(old)
        if n != 1:
            print("  [FAIL] 匹配 %d 次: %s ..." % (n, old[:68].replace("\n", "\\n")))
            sys.exit(1)
        text = text.replace(old, new)
        print("  [OK]   %s" % old[:68].replace("\n", "\\n"))
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print("  [written] %s\n" % os.path.basename(path))


if __name__ == "__main__":
    apply(SUPP, PATCHES_SUPP)
    apply(PKG, PATCHES_PKG)
