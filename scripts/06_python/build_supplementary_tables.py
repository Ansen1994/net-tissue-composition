# -*- coding: utf-8 -*-
"""
Build the consolidated supplementary-table workbook for manuscript v2
(tissue-composition / measurement-reproducibility version).

Input : CSV result files in this directory (work/)
Output: 03_补充材料/Supplementary_Tables_S1-S17.xlsx
        03_补充材料/Supplementary_Materials_v2.md
        03_补充材料/Supplementary_Table_Legends.md

Table numbering follows the citations in manuscript_v2_tissue_composition.md
(S1-S17). S1-S7 keep the legacy file prefixes; S8-S16 follow 改稿说明_v2.md §5; S17a-S17f are the exploratory drug-sensitivity results.
"""
import csv
import os
import re
import sys
from datetime import datetime

from openpyxl import Workbook
from openpyxl.cell import WriteOnlyCell
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

csv.field_size_limit(50 * 1024 * 1024)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT_DIR = os.path.join(ROOT, "03_补充材料")
XLSX = os.path.join(OUT_DIR, "Supplementary_Tables_S1-S17.xlsx")
MD_ALL = os.path.join(OUT_DIR, "Supplementary_Materials_v2.md")
MD_LEG = os.path.join(OUT_DIR, "Supplementary_Table_Legends.md")

HDR_FILL = PatternFill("solid", fgColor="1F4E79")
HDR_FONT = Font(bold=True, color="FFFFFF", size=10)
NOTE_FONT = Font(bold=True, color="9C0006")


def read_csv(name):
    with open(os.path.join(HERE, name), encoding="utf-8-sig", errors="replace", newline="") as fh:
        r = csv.reader(fh)
        try:
            header = next(r)
        except StopIteration:
            return [], []
        rows = [row for row in r]
    return header, rows


def typed(v):
    """Keep numbers numeric in Excel; never let text become a formula."""
    if v is None or v == "":
        return None
    if isinstance(v, (int, float, bool)):
        return v
    s = str(v).strip()
    if s in ("NA", "NaN", "NULL", "Inf", "-Inf"):
        return s
    if re.fullmatch(r"[+-]?\d+", s):
        try:
            n = int(s)
            return n if abs(n) < 2 ** 53 else s
        except ValueError:
            return s
    if re.fullmatch(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?", s):
        try:
            return float(s)
        except ValueError:
            return s
    if s[:1] in ("=", "+", "@"):
        return "'" + s
    return s[:32000]


def clean_header(header, ncols):
    out = []
    for i, h in enumerate(header):
        h = (h or "").strip()
        out.append(h if h else "Column_%d" % (i + 1))
    while len(out) < ncols:
        out.append("Column_%d" % (len(out) + 1))
    return out


# --------------------------------------------------------------------------
# Derived / newly assembled tables
# --------------------------------------------------------------------------

def first_col(name):
    _, rows = read_csv(name)
    return [r[0] for r in rows if r and r[0].strip()]


def build_s1_candidates():
    """The 132 unique candidate symbols, annotated with ER contrast and cell origin."""
    s2_genes = first_col("Table_S2_Target_NETs_Matrix.csv")
    uniq = []
    for g in s2_genes:
        base = re.sub(r"\.\d+$", "", g)
        if base not in uniq:
            uniq.append(base)
    er = {}
    _, rows = read_csv("Table_CandidatePool_ERcontrast.csv")
    for r in rows:
        if r and r[0]:
            er[r[0]] = r
    origin = {}
    _, rows = read_csv("fig1d_cell_origin_final.csv")
    for r in rows:
        if r and r[0]:
            origin[r[0]] = r
    header = ["GeneSymbol", "log2FoldChange_ERpos_vs_ERneg", "padj",
              "SingleCell_top_celltype", "SingleCell_class", "In_S2_matrix"]
    out = []
    for g in uniq:
        e = er.get(g, ["", "", ""])
        o = origin.get(g, ["", "", "", "", "", ""])
        out.append([g, e[1] if len(e) > 1 else "", e[2] if len(e) > 2 else "",
                    o[1] if len(o) > 1 else "", o[5] if len(o) > 5 else "", "Yes"])
    return header, out


def build_s8_summary():
    """Numeric characterisation of the candidate pool quoted in Results 3.1."""
    h_all, r_all = read_csv("Table_S1_All_DEGs.csv")
    h_can, r_can = read_csv("Table_CandidatePool_ERcontrast.csv")
    ix_fc = h_all.index("log2FoldChange") if "log2FoldChange" in h_all else 2
    ix_p = h_all.index("padj") if "padj" in h_all else 6
    ix_sym = h_all.index("Real_Symbol") if "Real_Symbol" in h_all else 0

    sig_syms = []
    for r in r_all:
        if len(r) <= max(ix_fc, ix_p, ix_sym):
            continue
        try:
            if abs(float(r[ix_fc])) > 1 and float(r[ix_p]) < 0.05:
                sig_syms.append(r[ix_sym])
        except (ValueError, TypeError):
            pass

    def nums(rows, ix):
        v = []
        for r in rows:
            if len(r) <= ix:
                continue
            try:
                x = float(r[ix])
                if x == x:
                    v.append(x)
            except ValueError:
                pass
        return v

    def median(v):
        v = sorted(v)
        n = len(v)
        return (v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2) if n else ""

    all_fc = nums(r_all, ix_fc)
    can_fc = [x for x in nums(r_can, 1)]
    can_p = nums(r_can, 2)
    pool = [r[0] for r in r_can if r and r[0]]
    prolif = ["MKI67", "CCNB1", "TOP2A", "AURKA", "BIRC5"]
    rows = [
        ["All tested genes in DESeq2 output (rows, n)", len(r_all), "Table_S1_All_DEGs.csv"],
        ["DEG rows meeting padj < 0.05 and |log2FC| > 1 (n)", len(sig_syms), "computed"],
        ["Unique gene symbols among DEGs (n)", len(set(sig_syms)), "computed"],
        ["Duplicated symbols removed by intersect() (n)",
         len(sig_syms) - len(set(sig_syms)), "computed"],
        ["Candidate pool (unique symbols, n)", len(pool), "Table_CandidatePool_ERcontrast.csv"],
        ["Median |log2FC| all tested genes", round(median([abs(x) for x in all_fc]), 4), "computed"],
        ["Median |log2FC| candidate pool", round(median([abs(x) for x in can_fc]), 4), "computed"],
        ["Candidate pool with padj < 1e-10 (%)",
         round(100 * sum(1 for x in can_p if x < 1e-10) / max(len(can_p), 1), 2), "computed"],
        ["Candidate pool lower in ER+ than ER- (%)",
         round(100 * sum(1 for x in can_fc if x < 0) / max(len(can_fc), 1), 2), "computed"],
        ["ESR1 in candidate pool", "Yes" if "ESR1" in pool else "No", "computed"],
        ["ERBB2 in candidate pool", "Yes" if "ERBB2" in pool else "No", "computed"],
    ]
    for g in prolif:
        rows.append(["Proliferation gene %s in pool" % g, "Yes" if g in pool else "No", "computed"])
    return ["Property", "Value", "Source"], rows


def parse_r_list(path, varname, stop_at=None):
    s = open(os.path.join(HERE, path), encoding="utf-8").read()
    start = s.index(varname + " <- list(")
    body = s[start + len(varname + " <- list("):]
    if stop_at and stop_at in body:
        body = body[:body.index(stop_at)]
    out = {}
    for m in re.finditer(r"([A-Za-z0-9_]+)\s*=\s*c\((.*?)\)", body, re.S):
        out[m.group(1)] = re.findall(r'"([^"]+)"', m.group(2))
    return out


def build_s12_modules():
    mods = parse_r_list("ia1_build_scores.R", "modules", stop_at="\nallg <-")
    header = ["Module", "Gene", "Gene_order", "n_genes_in_module", "Note"]
    out = []
    for name, genes in mods.items():
        note = ""
        if name == "Hub":
            note = "Risk-score genes themselves; not counted among the eleven compartment modules"
        elif name.endswith("CD79A"):
            note = "B-cell module version; CD79A included" if "with" in name else "B-cell module version; CD79A excluded"
        for i, g in enumerate(genes, 1):
            out.append([name, g, i, len(genes), note])
    return header, out


def build_s13_signatures():
    sigs = parse_r_list("j4c_generalize.R", "sigs", stop_at="\n# 已知系数")
    # w_sigs is a list of named numeric vectors, e.g. BRCA_A = c(ENO1 = 0.012, ...)
    s = open(os.path.join(HERE, "j4c_generalize.R"), encoding="utf-8").read()
    wblock = s[s.index("w_sigs <- list("):s.index("neutrophil_ref <- c(")]
    weights = {}
    for m in re.finditer(r"([A-Za-z0-9_]+)\s*=\s*c\((.*?)\)", wblock, re.S):
        sig = m.group(1)
        if sig == "w_sigs":
            continue
        weights[sig] = {g: float(v) for g, v in
                        re.findall(r"([A-Z0-9]+)\s*=\s*(-?[\d.]+)", m.group(2))}
    src = {}
    _, rows = read_csv("j4_signature_purity.csv")
    for r in rows:
        if len(r) > 1 and r[0]:
            src[r[0]] = r[1]
    refmap = {"BRCA_A": "[14]", "BRCA_B": "[15]", "BRCA_C": "[16]",
              "GC_A": "[17]", "GC_B": "[18]", "CRC_A": "[19]",
              "Ours": "This study"}
    header = ["Signature", "Reference", "Source", "Gene", "Published_weight", "Weight_status"]
    out = []
    for sig, genes in sigs.items():
        w = weights.get(sig, {})
        for g in genes:
            out.append([sig, refmap.get(sig, ""), src.get(sig, ""), g, w.get(g, ""),
                        "Published Cox coefficient used" if w else "Equal weights (no coefficients published)"])
    return header, out


def build_s13_consensus():
    # fig1d_build.R defines `neut <- c(...)` (a vector, not a list) -> parse directly
    s = open(os.path.join(HERE, "fig1d_build.R"), encoding="utf-8").read()
    m = re.search(r"neut <- c\((.*?)\)", s, re.S)
    c16 = re.findall(r'"([^"]+)"', m.group(1)) if m else []
    s2 = open(os.path.join(HERE, "j4c_generalize.R"), encoding="utf-8").read()
    m2 = re.search(r"neutrophil_ref <- c\((.*?)\)", s2, re.S)
    c18 = re.findall(r'"([^"]+)"', m2.group(1)) if m2 else []
    header = ["Gene", "In_consensus_16_Methods_2_9_Figure_1D",
              "In_counted_18_j4c_signature_counts", "Status"]
    union = list(c16) + [g for g in c18 if g not in c16]
    out = []
    for g in union:
        a, b = g in c16, g in c18
        out.append([g, "Yes" if a else "No", "Yes" if b else "No",
                    "Both lists" if a and b else ("Consensus 16 only" if a else "Counted 18 only")])
    return header, out


def build_s14_reference():
    s = open(os.path.join(HERE, "j1_nets_prognosis.R"), encoding="utf-8").read()
    m = re.search(r"nets_core <- c\((.*?)\)", s, re.S)
    genes = re.findall(r'"([^"]+)"', m.group(1)) if m else []
    seen, uniq = set(), []
    for g in genes:
        if g not in seen:
            seen.add(g)
            uniq.append(g)
    return ["Order", "Gene", "Note"], [[i, g, "22 unique loci after removing duplicated CAMP"]
                                       for i, g in enumerate(uniq, 1)]


# --------------------------------------------------------------------------
# Sheet specification
# --------------------------------------------------------------------------
# (sheet, table id, title, [source csv] or builder, cited-in, note)
SPECS = [
    ("Table S1", "S1", "All differentially expressed genes, ER+ vs ER- (TCGA-BRCA)",
     ["Table_S1_All_DEGs.csv"], "3.1",
     "Column In_132_candidates flags the intersection with the 509 GeneCards NET entries: 133 rows are flagged because ELFN2 has two Ensembl mappings, corresponding to 132 unique symbols."),
    ("Table S1b", "S1b", "The 132 candidate genes (unique symbols) with ER contrast and single-cell origin",
     "build_s1_candidates", "3.1",
     "Derived. 133 rows in Table S2 collapse to 132 unique symbols (ELFN2 has two Ensembl mappings)."),
    ("Table S2", "S2", "Expression matrix of the candidate pool in TCGA-BRCA tumours (ER+ and ER-)",
     ["Table_S2_Target_NETs_Matrix.csv"], "3.1",
     "133 rows x 977 samples (780 ER+, 197 ER-). ELFN2 and ELFN2.1 are two Ensembl IDs for one symbol; modelling used the 132 unique symbols and the 762 ER+ training tumours only."),
    ("Table S3", "S3", "Machine-learning input matrix (762 ER+ tumours, 132 genes + survival)",
     ["Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv"], "2.3",
     "762 ER-positive tumours with overall survival data (99 deaths)."),
    ("Table S4", "S4", "LASSO (lambda.min), random-forest and hub gene sets",
     ["Table_S4_Machine_Learning_Features.csv"], "3.2",
     "LASSO 11 features; random forest 58; intersection 4 (PLG, CD79A, MMP7, ANGPT1)."),
    ("Table S4b", "S4b", "Four hub genes of the signature",
     ["Table_S4_Hub_Genes.csv"], "3.2", ""),
    ("Table S5", "S5", "Multivariate Cox coefficients with 95% confidence intervals used to build the risk score",
     ["Table_S5_Cox_Coefficients.csv"], "3.3",
     "Joint model of the four hub genes (n = 762, 99 events). Coef is the per-unit log-hazard coefficient; SE, z and the Wald 95% confidence intervals are taken from summary(coxph). Distinct from the univariate estimates in Table S5b; the proportional-hazards test for this model is in Table S5c."),
    ("Table S5b", "S5b", "Univariate Cox estimates for the four hub genes",
     ["Table_Univariate_Cox_Results.csv"], "3.3", ""),
    ("Table S5c", "S5c", "Proportional-hazards test for the joint model (scaled Schoenfeld residuals)",
     ["Table_S5c_PH_Test.csv"], "3.3",
     "survival::cox.zph applied to the joint model of Table S5 (n = 762, 99 events). CD79A, MMP7 and ANGPT1 depart from proportional hazards (global p = 0.004), so the hazard ratios in Table S5 and Figure 3A are time-averaged rather than constant. The risk score is a linear predictor of that model, and its discrimination (Figures 3F and 4D) does not rest on this assumption."),
    ("Table S6a", "S6a", "Risk scores and survival, TCGA-BRCA training cohort (n = 762)",
     ["Table_S6_Patient_Risk_Score_Result.csv"], "3.3", "Training cohort; dichotomised at the cohort median."),
    ("Table S6b", "S6b", "Risk scores and survival, GSE20685 (n = 327)",
     ["Table_S6_GEO_Validation_Result.csv"], "3.4", "Validation; surv_cutpoint dichotomisation."),
    ("Table S6c", "S6c", "Risk scores and survival, METABRIC ER+ (n = 1,445)",
     ["Table_S6_METABRIC_Validation_Result.csv"], "3.4", "Validation; surv_cutpoint dichotomisation."),
    ("Table S6d", "S6d", "Risk scores and survival, GSE96058 ER+ (n = 2,840)",
     ["Table_S6_GSE96058_Validation_Result.csv"], "3.4", "Validation; surv_cutpoint dichotomisation."),
    ("Table S7", "S7", "Clinical data used for the age/stage-adjusted models (n = 686)",
     ["Table_S7_Final_Clinical_Data.csv"], "2.1",
     "686 patients with both age and pathological stage recorded; not the full 762-patient training set."),
    ("Table S8", "S8", "ER+ vs ER- contrast for each of the 132 candidate genes",
     ["Table_CandidatePool_ERcontrast.csv"], "2.2, 3.1", ""),
    ("Table S8b", "S8b", "Numeric characterisation of the candidate pool",
     "build_s8_summary", "3.1",
     "Derived from Table S1 and Table S8. The 12,984 DEGs quoted in the text are unique HGNC symbols after removing 98 duplicated mappings from 13,082 significant rows."),
    ("Table S9", "S9", "C-index under four cross-platform scaling schemes in three external cohorts",
     ["Table_S9_Scaling_Scheme_Cindex.csv"], "3.4", 'S0 training coefficients on raw expression (validation cohorts z-scored); S1 z-scores in all cohorts with coefficients re-expressed per SD; S2 within-cohort percentile ranks with coefficients refitted in TCGA; S3 within-cohort centring scaled by training-cohort SDs with the original coefficients. Only S0/S1 reproduce the in-sample TCGA C-index of 0.690; S2 (0.654) and S3 (0.666) differ because the transformation changes the linear predictor. This table required expression and survival only (METABRIC n = 1,498 / 879 events; GSE96058 n = 2,935 / 264 events), whereas Tables S10 and S15 additionally required clinical covariates (1,445 / 848 and 2,840 / 247). Methods 2.1 states both complete-case definitions with their event counts.'),
    ("Table S10a", "S10a", "Incremental value of the risk score over age and stage",
     ["Table_S10a_Incremental_Value.csv"], "3.4", 'Complete cases with age, stage/nodal status, risk score and survival: n = 327 (83 events), 1,445 (848) and 2,840 (247). Likelihood-ratio test of combined model vs clinical model. C_clinical / C_combined are distinct fits from the refitted ROC marker in Table S10b.'),
    ("Table S10b", "S10b", "C-index and time-dependent AUC by predictor variable",
     ["Table_S10b_AUC_by_Variable.csv"], "3.4", 'RiskScore = signature only; CombinedScore = refitted score + age + stage. The METABRIC validation table does not store CombinedScore; Figure 4F refits it within the cohort from score, age and nodal status and reproduces this row exactly (C-index 0.657; AUC 0.601/0.668/0.668), so the row is retained. It is not used to claim incremental value. All three horizons are tabulated here, but only the 3- and 5-year curves are drawn in Figure 4F: the 1-year estimate rests on 5, 14 and 14 events, too few for a stable ROC.'),
    ("Table S11a", "S11a", "Single-cell mean expression and detection rate of the four hub genes by cell type",
     ["j3_celltype_expression.csv"], "2.6, 3.6", "GSE176078, ER+ tumours, 26 donors, 38,241 cells. Mean expression is CP10K-log (log1p counts per 10,000 full-transcriptome UMIs) and detection rate is the fraction of cells with at least one UMI. Labels were reused from the published annotation without re-clustering or dimensionality reduction; the annotation scheme has no neutrophil or granulocyte class (Methods 2.6, Supplementary Materials: single-cell annotation audit)."),
    ("Table S11b", "S11b", "Detection of NET effector transcripts in single cells",
     ["j3_nets_detection.csv"], "3.6",
     "Detection rates are given for all 100,064 cells and for the 38,241 ER+ cells; Figure 6D plots the ER+ column only. Results 3.6 reports both columns."),
    ("Table S11c", "S11c", "Cell-of-origin classification for all 132 candidate genes",
     ["fig1d_cell_origin_final.csv", "fig1d_summary.csv"], "3.1",
     "Two blocks in one sheet: per-gene assignment, then class summary."),
    ("Table S12a", "S12a", "Gene lists of the eleven compartment modules",
     "build_s12_modules", "2.5",
     "Derived from ia1_build_scores.R. Eleven module scores are counted, with the CD79A-free and CD79A-including B-cell versions counted separately; the 'Hub' term is the risk score itself and is not one of the eleven. This panel is distinct from the nine generalisation modules used in Methods 2.9 and Figure 7D (neutrophil, NETosis, stroma/ECM, endothelium, epithelial, proliferation, T cells, B cells, myeloid); 'NETosis' there is a six-gene effector module, not the 22-gene reference score of Table S14b."),
    ("Table S12b", "S12b", "Spearman correlations between risk score and module scores in four cohorts",
     ["Table_ImmuneAxis_AllCorrelations.csv"], "3.6", ""),
    ("Table S12c", "S12c", "Cross-cohort sign consistency and Fisher combined p-values",
     ["Table_ImmuneAxis_Consistency.csv"], "3.6", ""),
    ("Table S12d", "S12d", "Correlations adjusted for epithelial purity and proliferation",
     ["Table_ImmuneAxis_Adjusted.csv"], "3.6", "rho_raw is unadjusted; rho_adjPurity controls for epithelial purity; rho_adjPurityProlif additionally controls for proliferation. The three values quoted in Results 3.6 (−0.377, −0.268, −0.269) are the Bcell_noCD79A rows of rho_adjPurity in TCGA, GSE20685 and METABRIC."),
    ("Table S12e", "S12e", "B-cell module: effect of removing CD79A from the score",
     ["Table_ImmuneAxis_BcellDecomposition.csv"], "3.6", ""),
    ("Table S12f", "S12f", "Risk-score decomposition into full, CD79A-free and CD79A-only terms",
     ["Table_ImmuneAxis_ScoreDecomposition.csv"], "3.6", ""),
    ("Table S12g", "S12g", "Module score differences between high- and low-risk groups",
     ["Table_ImmuneAxis_GroupDiff.csv"], "3.6", ""),
    ("Table S13a", "S13a", "Six published NET-related signatures and the present signature",
     "build_s13_signatures", "2.9",
     'Derived from j4c_generalize.R and j4_signature_purity.csv. Six previously published signatures (37 loci) plus the present four-gene signature (4 loci) give seven signatures and 41 loci in total; 8 of 41 loci (19.5%) are neutrophil-restricted, or 8 of 37 (21.6%) among published signatures alone. The Reference column maps BRCA_A, BRCA_B, BRCA_C, GC_A, GC_B and CRC_A to manuscript citations [14], [15], [16], [17], [18] and [19]. Published Cox weights were available for BRCA_A and Ours only.'),
    ("Table S13b", "S13b", "Neutrophil-restricted / granule gene lists used in this study",
     "build_s13_consensus", "2.9",
     'Audit trail for the neutrophil gene lists. The 16-gene consensus (Methods 2.9, Figure 1D) is now used throughout; the 18-gene list formerly used by j4c_generalize.R is retained for comparison. Both give identical counts for all seven signatures (six published plus the present one: 2/2/0/2/2/0/0; 8 of 41 loci, 19.5%), so no published value or figure changes. 11 genes are shared.'),
    ("Table S13c", "S13c", "Neutrophil-restricted gene content of each signature",
     ["j4c_nets_content.csv"], "3.7", "Counts against the 16-gene consensus of Table S13b; identical to the superseded 18-gene list."),
    ("Table S13d", "S13d", "Detectability of NET effector transcripts in 1,111 TCGA tumours",
     ["j4c_nets_detectability.csv"], "3.7", ""),
    ("Table S13e", "S13e", "Signature x module correlations, raw and partial (controlling epithelial content)",
     ["j4c_signature_module_correlations.csv"], "3.7", "Genes shared between signature and module removed before correlation."),
    ("Table S13f", "S13f", "Signature purity: module correlation profile and dominant compartment",
     ["j4_signature_purity.csv"], "3.7", "Earlier j4 run; retained as supporting information."),
    ("Table S14", "S14", "Reference NET score and co-expression surrogate top-N sets in four cohorts",
     ["j1_results.csv"], "2.7", ""),
    ("Table S14b", "S14b", "The 22-gene reference NET score gene list",
     "build_s14_reference", "2.7", "Derived from j1_nets_prognosis.R; 22 unique loci after removing a duplicated CAMP entry."),
    ("Table S15", "S15", "Score x treatment interaction in METABRIC and GSE96058",
     ["j2_results.csv"], "2.8", ""),
    ("Table S16", "S16", "Sensitivity of Kaplan-Meier results to the dichotomisation rule",
     ["Table_S16_Cutpoint_Sensitivity.csv"], "2.4, 3.4",
     "Median cut-point (pre-specified) vs maximally selected rank statistic. External log-rank significance depends on the rule."),
    ("Table S17a", "S17a", "Predicted log IC50 differences for 29 compounds (primary analysis)",
     ["../药敏重算/drug_results.csv"], "2.10, 3.6",
     "Ridge regression (glmnet, alpha = 0, 5-fold CV) trained on GDSC2. diff = median(high-risk) - median(low-risk); fdr is Benjamini-Hochberg across all 29 compounds. Only sapitinib (FDR 9.9e-09) and fulvestrant (FDR 0.017) survive FDR < 0.05."),
    ("Table S17b", "S17b", "Drug-sensitivity predictions under five preprocessing strategies",
     ["../药敏重算/敏感性/sensitivity_matrix.csv"], "2.10, 3.6",
     "The 12 compounds carried into the sensitivity panel x 5 input modes (counts, log2tpm, zscore, quantile-normalised, rank). outside_frac = fraction of predictions outside the GDSC training response range; cv_r2 = cross-validation R2. Absolute effect sizes are not comparable across modes."),
    ("Table S17c", "S17c", "Robustness of direction across the five preprocessing strategies",
     ["../药敏重算/敏感性/sensitivity_verdict.csv"], "2.10, 3.6",
     "n_sig = modes with p < 0.05; robust requires >= 4 significant modes and a consistent direction. Fulvestrant is directionally unstable before common-axis adjustment."),
    ("Table S17d", "S17d", "Leave-one-out common-axis adjustment for 11 targeted compounds",
     ["../药敏重算/敏感性/targeted_verdict.csv"], "2.10, 3.6",
     "After removing the shared prediction axis (leave-one-out over 11 targeted drugs), only sapitinib, lapatinib and fulvestrant retain a consistent direction across all five modes (ROBUST-SIGNAL); tamoxifen is ROBUST-NULL."),
    ("Table S17e", "S17e", "Per-mode residual effect sizes underlying Table S17d",
     ["../药敏重算/敏感性/targeted_residual.csv"], "2.10, 3.6",
     "d_res11 = Cohen's d after LOO-11 common-axis removal; d_res29 after LOO-29 removal. n_axis11 / n_axis29 give the number of compounds contributing to the common axis."),
    ("Table S17f", "S17f", "Sensitivity of the overall drug panel to input preprocessing",
     ["../药敏重算/敏感性/sensitivity_by_mode.csv"], "2.10, 3.6",
     "n_sig = of the 12 sensitivity-panel compounds, the number with p < 0.05 under each mode; rho_avg_risk = mean correlation between predicted IC50 and risk score; pc1_pct = variance explained by the first principal component (the common axis)."),
]

BUILDERS = {
    "build_s1_candidates": build_s1_candidates,
    "build_s8_summary": build_s8_summary,
    "build_s12_modules": build_s12_modules,
    "build_s13_signatures": build_s13_signatures,
    "build_s13_consensus": build_s13_consensus,
    "build_s14_reference": build_s14_reference,
}


def resolve(spec):
    """Return (header, rows) for a sheet."""
    src = spec[3]
    if isinstance(src, str) and src.startswith("build_"):
        return BUILDERS[src]()
    header, rows = None, []
    for k, f in enumerate(src):
        h, r = read_csv(f)
        if len(src) == 1:
            header, rows = h, r
            break
        if header is None:
            header = h
            rows.extend(r)
        else:
            # stack several small tables vertically with a labelled separator row
            rows.append([])
            rows.append(["--- %s ---" % f] + [""] * max(0, len(header) - 1))
            rows.append(list(h) + [""] * max(0, len(header) - len(h)))
            rows.extend(r)
    ncols = max((len(r) for r in rows), default=0)
    header = clean_header(header or [], ncols)
    if spec[0] == "Table S1":
        cand = set(re.sub(r"\.\d+$", "", g)
                   for g in first_col("Table_S2_Target_NETs_Matrix.csv"))
        # Table_S1_All_DEGs.csv stores Ensembl IDs in GeneSymbol; the HGNC symbol
        # is in Real_Symbol, which is what intersect() was applied to.
        ix = header.index("Real_Symbol") if "Real_Symbol" in header else 0
        header.append("In_132_candidates")
        for r in rows:
            sym = r[ix] if len(r) > ix else ""
            del r[len(header) - 1:]
            r += [""] * (len(header) - 1 - len(r))
            r.append("Yes" if sym in cand else "No")
    return header, rows


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    wb = Workbook(write_only=True)

    # ---------------- README / index ----------------
    ws = wb.create_sheet("README")
    ws.freeze_panes = "A3"
    title = WriteOnlyCell(ws, value="Supplementary Tables - manuscript v2 (tissue composition / measurement reproducibility)")
    title.font = Font(bold=True, size=13, color="1F4E79")
    ws.append([title])
    ws.append(["Generated", datetime.now().strftime("%Y-%m-%d %H:%M"),
               "Source directory", "05_数据与结果/work/",
               "Manuscript", "01_手稿/_md/manuscript_v2_tissue_composition.md"])
    ws.append([])
    ws.append(["Sheet", "Table ID", "Title", "Source file(s)", "Rows (data)", "Cols", "Cited in", "Note"])
    meta = []
    for spec in SPECS:
        header, rows = resolve(spec)
        meta.append((spec, header, rows))
        ws.append([spec[0], spec[1], spec[2],
                   "; ".join(spec[3]) if isinstance(spec[3], list) else "derived (see builder script)",
                   len(rows), len(header), spec[4], spec[5]])
    ws.append([])
    ws.append(["NOTES"])
    for line in [
        "1. Table IDs follow the citations in the v2 manuscript (S1-S17). Lettered sheets are sub-panels of the same table number.",
        "2. Table S2 has 133 rows (ELFN2 / ELFN2.1 are two Ensembl IDs for one symbol); all modelling used 132 features and the manuscript reports 132 candidates.",
        "3. Table S13b retains the audit trail of two neutrophil gene lists. j4c_generalize.R now uses the same 16-gene consensus as Methods 2.9 and Figure 1D; re-running it reproduced all three result files byte-for-byte, because both lists give identical counts for every signature.",
        "4. Table S6a uses a median cut-point; S6b-d use surv_cutpoint. Table S16 quantifies the sensitivity to this choice.",
        "5. The 58 random-forest features in Table S4 come from a single run under set.seed(2026); randomForest importance rankings are not byte-identical across package versions, so this set is the reference and no alternative seed was evaluated. Figure 2 importance values were recalculated within those 58 features.",
        "6. All p-values are two-sided; no multiple-testing adjustment beyond that stated in Methods 2.11 unless a p_BH column is present.",
        "7. DEG count: 13,082 rows meet padj < 0.05 and |log2FC| > 1; they collapse to 12,984 unique HGNC symbols (98 duplicated Ensembl mappings). The manuscript reports 12,984, i.e. the symbol-level count used by intersect().",
        "8. Table S17a-S17f are the complete exploratory drug-sensitivity results (Methods 2.10, Results 3.6). Predictions come from a direct glmnet ridge regression, not the pRRophetic package; only two of 29 compounds survive FDR in the primary analysis and only three remain directionally stable after common-axis adjustment. These tables are supplied for completeness and are not used to support any pharmacological claim.",
        "9. Six published signatures (37 loci) plus the present signature (4 loci) = 41 loci in Table S13a; the Reference column maps each signature to citations [14]-[19]. The 11-module panel of Table S12a is distinct from the 9-module generalisation panel of Table S13e/Figure 7D.",
    ]:
        ws.append([line])

    # ---------------- data sheets ----------------
    for spec, header, rows in meta:
        name = spec[0]
        ws = wb.create_sheet(name[:31])
        hdr_cells = []
        for h in header:
            c = WriteOnlyCell(ws, value=h)
            c.font = HDR_FONT
            c.fill = HDR_FILL
            c.alignment = Alignment(vertical="center", wrap_text=True)
            hdr_cells.append(c)
        ws.append(hdr_cells)
        ws.freeze_panes = "A2"
        ncols = len(header)
        if ncols < 40:
            for i, h in enumerate(header, 1):
                ws.column_dimensions[get_column_letter(i)].width = min(max(11, len(str(h)) + 3), 42)
        else:
            ws.column_dimensions["A"].width = 22
        for r in rows:
            vals = [typed(x) for x in list(r)[:ncols]]
            vals += [None] * (ncols - len(vals))
            ws.append(vals)
        print("[sheet] %-18s %7d rows x %4d cols" % (name, len(rows), ncols), flush=True)

    wb.save(XLSX)
    print("[saved]", XLSX, os.path.getsize(XLSX), "bytes")
    write_markdown(meta)
    return meta


def write_markdown(meta):
    lines = []
    lines.append("# Supplementary Materials")
    lines.append("")
    lines.append("**Manuscript v2 — tissue composition and measurement reproducibility of a NET-related prognostic signature in ER-positive breast cancer**")
    lines.append("")
    lines.append("Generated %s. All supplementary tables are consolidated in `Supplementary_Tables_S1-S17.xlsx` (one worksheet per table / sub-panel)." % datetime.now().strftime("%Y-%m-%d"))
    lines.append("")
    lines.append("## Supplementary Tables")
    lines.append("")
    lines.append("Tables are numbered S1-S17 exactly as cited in the manuscript. Lettered suffixes (for example S6a-S6d, S13a-S13f, S17a-S17f) denote sub-panels of the same table number and correspond to separate worksheets in the workbook.")
    lines.append("")
    for spec, header, rows in meta:
        tid = spec[1]
        lines.append("**Table %s. %s** (worksheet %s; %d data rows; cited in Section %s)." %
                     (tid, spec[2], spec[0], len(rows), spec[4]))
        if spec[5]:
            lines.append("")
            lines.append("> %s" % spec[5])
        lines.append("")
        last_tid = tid
    lines.append("## Supplementary Figures")
    lines.append("")
    fig_leg = os.path.join(OUT_DIR, "Supplementary_Figure_Legends.md")
    if os.path.exists(fig_leg):
        txt = open(fig_leg, encoding="utf-8").read()
        txt = re.sub(r"^# Supplementary Figure Legends\s*", "", txt).strip()
        lines.append(txt)
        lines.append("")
    lines.append("## Notes on data conventions")
    lines.append("")
    for n in [
        "**Candidate-pool size.** The intersection of 12,984 ER+ vs ER- DEGs with 509 GeneCards NET entries yielded 132 unique gene symbols. In Table S1, 13,082 rows meet the stated DESeq2 threshold; these collapse to 12,984 unique HGNC symbols because 98 symbols have more than one Ensembl mapping, and R's `intersect()` operates on the de-duplicated symbol vector. The expression matrix in Table S2 contains 133 rows because ELFN2 maps to two Ensembl identifiers; all downstream modelling used 132 features.",
        "**Signature count.** Six previously published signatures (37 loci) plus the present four-gene signature give seven signatures and 41 loci; 8/41 = 19.5% is the all-seven fraction and 8/37 = 21.6% the published-only fraction. Table S13a carries the reference-number mapping to citations [14]-[19].",
        "**Neutrophil gene lists (harmonised).** Methods 2.9 and Figure 1D use a 16-gene consensus of neutrophil-restricted / granule genes (ELANE, MPO, PADI4, AZU1, CTSG, PRTN3, BPI, CEACAM8, LCN2, LTF, DEFA3, DEFA4, CD177, CXCR1, PRSS57, LPO). The signature-level counts in Table S13c were originally produced with an 18-gene `neutrophil_ref` list in `j4c_generalize.R`. The lists share 11 genes; the 7 genes unique to the 18-gene list (CAMP, CSF3R, CXCR2, FCGR3B, FPR1, S100A8, S100A9) occur in none of the six published signatures, so both lists give exactly the same counts (2/2/0/2/2/0/0; 8 of 41 loci, 19.5%). `j4c_generalize.R` has been updated to the 16-gene consensus and re-run, reproducing all three result files byte-for-byte. Table S13b preserves both lists for audit; no manuscript value and no panel of Figure 7 changes.",
        "**Dichotomisation.** Training-cohort Kaplan-Meier curves use the cohort median; external cohorts use the maximally selected rank statistic (`surv_cutpoint`) as in the original reports. Table S16 shows that external log-rank significance depends entirely on this choice.",
        "**Random-forest reproducibility.** Feature selection was run once under `set.seed(2026)` and is deterministic under the recorded package versions. The 58-feature set in Table S4 is the reference result because randomForest importance rankings are not byte-identical across package versions; no alternative seed was evaluated. Figure 2 displays importance values recalculated within these 58 features only.",
        "**Sample counts.** 762 ER-positive TCGA tumours with survival data (99 deaths) were used for model fitting; the 686 patients with complete age and stage were used for adjusted models and immune-axis analyses. Methods 2.1 reports both complete-case definitions with event counts. External cohorts: GSE20685 n = 327 (83 events) in both definitions; METABRIC ER+ n = 1,498 (879 events) for expression-plus-survival analyses and n = 1,445 (848 events) when clinical covariates were also required; GSE96058 ER+ n = 2,935 (264 events) and n = 2,840 (247 events), respectively. Table S9 uses the former; Tables S6, S10 and S15 and Figure 4E use the latter, which is the subset quoted elsewhere in the manuscript unless stated otherwise.",
        "**Single-cell annotation audit.** The single-cell analysis reused the curated labels of GSE176078 without quality-control filtering, re-clustering or re-annotation; expression was normalised to CP10K-log (log1p counts per 10,000 full-transcriptome UMIs) for gene-level attribution only. The attribution analysis therefore contains no clustering or dimensionality-reduction step. A UMAP embedding was recomputed separately from the deposited matrices for display only (Supplementary Figure S3); its coordinates are not those of Wu et al. [32] and were not used for clustering, annotation or inference. Across all 26 tumours the nine major classes contain 100,064 cells: Cancer Epithelial 24,489; T-cells 35,214; Myeloid 9,675; Endothelial 7,605; CAFs 6,573; PVL 5,423; Normal Epithelial 4,355; Plasmablasts 3,524; B-cells 3,206 (subtypes: ER+ 38,241; HER2+ 19,311; TNBC 42,512). Neither the nine major nor the 29 minor labels contain a neutrophil or granulocyte class; the myeloid minor labels are Macrophage (5,929), Monocyte (2,328), DCs (955) and Cycling_Myeloid (463), 9,675 cells in total; NK cells (1,846) are annotated as a separate lineage and are not part of the myeloid class. Lineage-marker verification (mean CP10K-log expression in ER+ cells): EPCAM -> Cancer Epithelial 0.85, CD3D -> T-cells 1.42, MS4A1 -> B-cells 2.03, LYZ -> Myeloid 2.31, PECAM1 -> Endothelial 1.96, COL1A1 -> CAFs 3.40; these six lineage markers are displayed in Figure S3E. MKI67 is a proliferation marker, not a lineage marker: at the nine-major level it peaks in Myeloid (0.03) and is excluded from Figure S3E but retained in the audit table. Per-class counts and marker values are reproducible from the released j3_singlecell.R, which writes j3_annotation_audit.csv and j3_marker_verification.csv.",
        "**Single-cell subtype denominators and reusable UMAP outputs.** Hub-gene attribution, marker verification and Figure 6D use the 11 ER+ tumours (38,241 cells); all-cell rates over 100,064 cells are reported alongside, and Figure 1D compartment assignment uses all 26 tumours without subtype filtering. The deposited subtype field is a three-category split (ER+ 11 / HER2+ 5 / TNBC 10), whereas the GEO sample characteristics record two HER2+/ER+ tumours (CID3586, CID4066; 11,487 cells) and one label discrepancy (CID3963: ER+ in GEO, TNBC in the processed metadata; 3,527 cells). The ER+ subset is therefore not an exhaustive sample of ER-positive disease. Recomputed coordinates, per-cell metadata, 12 gene values, parameters, package versions and a SHA-256 cache key are archived in `12_可复用资源/sc_umap/` with `sc_umap_load.R` / `sc_umap_load.py` loaders (see REUSE.md).",
        "**Two module panels.** Table S12a holds the eleven compartment modules of Methods 2.5 (the two B-cell versions counted separately). Figure 7D and Table S13e use a separate nine-module generalisation panel defined in `j4c_generalize.R`; its NETosis module is a six-gene effector set, not the 22-gene reference NET score of Table S14b.",
        "**Exploratory analyses.** The complete drug-sensitivity predictions behind Methods 2.10 and Results 3.6 are now supplied in Table S17a-S17f. They are hypothesis-generating only: only sapitinib and fulvestrant survive FDR correction in the primary analysis, and the direction of most compounds changes with preprocessing or disappears after removal of the shared prediction axis. No pharmacological conclusion is drawn from them.",
    ]:
        lines.append("- " + n)
    lines.append("")
    lines.append("## File inventory")
    lines.append("")
    for f in ["Supplementary_Tables_S1-S17.xlsx", "Supplementary_Materials_v2.md",
              "Supplementary_Table_Legends.md", "Supplementary_Figure_Legends.md",
              "Supplementary Figure S1.pdf", "Supplementary Figure S2.tif",
              "Supplementary Figure S3.pdf", "源文件/ (editable artwork)"]:
        lines.append("- `03_补充材料/%s`" % f)
    lines.append("")
    txt = "\n".join(lines)
    open(MD_ALL, "w", encoding="utf-8").write(txt)

    # A table-legends-only file, ready to paste into the journal's template
    leg = ["# Supplementary Table Legends", ""]
    for spec, header, rows in meta:
        leg.append("**Table %s. %s**" % (spec[1], spec[2]))
        if spec[5]:
            leg.append("")
            leg.append("%s" % spec[5])
        leg.append("")
    open(MD_LEG, "w", encoding="utf-8").write("\n".join(leg))
    print("[saved]", MD_ALL)
    print("[saved]", MD_LEG)


if __name__ == "__main__":
    main()
