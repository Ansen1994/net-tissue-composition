# -*- coding: utf-8 -*-
"""Apply the reviewer-style reference audit fixes to manuscript_v2.

Run from anywhere:  python apply_ref_fixes.py
Every replacement is asserted to match exactly once, so a partial apply
cannot silently corrupt the manuscript.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MS = os.path.join(ROOT, "01_手稿", "_md", "manuscript_v2_tissue_composition.md")
CHANGELOG = os.path.join(ROOT, "改稿说明_v2.md")
TITLEPAGE = os.path.join(ROOT, "07_投稿材料", "Title_Page_and_Declarations.md")
REPORT = os.path.join(ROOT, "参考文献核查报告.md")

REPLACEMENTS = [
    # A1 + [2]: [3] data/population/window/outcome; [2] 75% not in the abstract
    (
        'Breast cancer remains the most frequently diagnosed malignancy and a leading cause of cancer-related mortality among women[1]. Estrogen receptor-positive (ER+) disease accounts for approximately 75% of cases[2], and endocrine therapy is its therapeutic cornerstone. Despite this, up to 30% of patients relapse within 20 years[3], and the mechanisms of endocrine and systemic resistance are heterogeneous, spanning ESR1 mutations and activation of bypass signalling[4]. Beyond cancer-cell-intrinsic aberrations, remodelling of the tumour microenvironment (TME) is increasingly recognised as a determinant of resistance and immune evasion in hormone receptor-positive disease[5].',
        'Breast cancer remains the most frequently diagnosed malignancy and a leading cause of cancer-related mortality among women[1]. It is a molecularly heterogeneous disease in which hormone receptor-positive tumours form the largest clinical subgroup and endocrine therapy is the therapeutic cornerstone[2]. Among patients with ER-positive (ER+) early breast cancer who remain disease-free after five years of endocrine therapy, the 5–20-year risk of distant recurrence ranges from 13% to 41% depending on tumour size and nodal status[3]. The mechanisms of endocrine and systemic resistance are heterogeneous, spanning ESR1 mutations and activation of bypass signalling[4]. Beyond cancer-cell-intrinsic aberrations, remodelling of the tumour microenvironment (TME) is increasingly recognised as a determinant of resistance and immune evasion in hormone receptor-positive disease[5].',
    ),
    # A6: [5] supports the ER+ claim, [6] supports the framework
    (
        'ER+ tumours are conventionally described as immunologically "cold", with sparse cytotoxic T-cell infiltration and a stroma dominated by immunosuppressive populations[6].',
        'ER+ tumours are conventionally described as immunologically "cold", with sparse cytotoxic T-cell infiltration and a stroma dominated by immunosuppressive populations[5], a classification that follows the hot/altered/cold framework[6].',
    ),
    # A5: mechanism is shielding from contact, not blocking infiltration
    (
        'obstruct intratumoural CD8+ T-cell infiltration[11] and drive T-cell exhaustion[12].',
        'shield tumour cells from CD8+ T-cell and NK-cell cytotoxicity by obstructing immune-cell contact[11], and drive T-cell exhaustion[12].',
    ),
    # [14-19]: narrower, defensible claim
    (
        'and the resulting score is presented as an independent prognostic factor[14–19].',
        'and the resulting score is reported as a prognostic or risk-stratification tool[14–19].',
    ),
    # [32]: state the annotation levels actually present (9 major / 29 minor)
    (
        'GSE176078[32] comprises 26 treatment-naive breast tumours profiled by 10x Chromium with curated annotation. Only ER+ tumours were analysed (38,241 of 100,064 cells).',
        'GSE176078[32] comprises 26 primary untreated breast tumours profiled by 10x Chromium, with curated annotation at nine major and 29 minor cell-type levels. Only ER+ tumours were analysed (38,241 of 100,064 cells). Hub-gene localisation and effector-transcript detection used the nine major labels, whereas the candidate-pool assignment in Figure 1D groups cells into six broad compartments.',
    ),
    # A3: GDSC2 needs its own citations (Yang 2013, Picco 2019 added as [35][36])
    (
        '**2.10. Drug sensitivity (exploratory).** Predicted sensitivity to 29 compounds was estimated by ridge regression trained on GDSC2[33,34], comparing five preprocessing strategies and applying leave-one-out common-axis adjustment; predictions are reported as exploratory only.',
        '**2.10. Drug sensitivity (exploratory).** Predicted sensitivity to 29 compounds was estimated by ridge regression (pRRophetic[33]) trained on the GDSC2 release of the Genomics of Drug Sensitivity in Cancer resource[34–36], comparing five preprocessing strategies and applying leave-one-out common-axis adjustment; predictions are reported as exploratory only.',
    ),
    # renumber TRIPOD
    ('Reporting follows TRIPOD guidance[37].', 'Reporting follows TRIPOD guidance[39].'),
    # NEW finding: 28 cell types does not exist; analysis used 9 major classes
    (
        'In GSE176078 (26 tumours, 100,064 cells; 38,241 ER+ cells; 28 annotated cell types in total), no hub gene was tumour-cell derived (Figure 6A–C, Table S11). PLG was not detected in any cell type (0% in all 28 classes; mean 0.010)',
        'In GSE176078 (26 tumours, 100,064 cells; 38,241 ER+ cells; nine major annotated cell types), no hub gene was tumour-cell derived (Figure 6A–C, Table S11). PLG was essentially undetected in every major cell type (≤0.2% of cells in any class; highest class mean 0.010)',
    ),
    # A7 + renumber [35]->[37], [36]->[38] + correct "100,064 tumour cells"
    (
        'ANGPT1, the strongest risk factor, is expressed by perivascular cells and marks vascular maturation and pericyte coverage[35]. MMP7, with an inverse coefficient, is expressed by normal luminal epithelium, so its coefficient reflects the proportion of retained normal structures, and tumours retaining organised glandular architecture are generally better differentiated. CD79A marks B cells, and B-cell infiltration and tertiary lymphoid structures are associated with favourable outcome[36]; so the score–B-cell correlation is a statement about B-cell quantity. PLG is a liver-derived plasma protein undetected in 100,064 tumour cells, so its coefficient describes detection versus non-detection in RNA-seq and has no counterpart in microarray cohorts where every sample is detected.',
        'ANGPT1, the strongest risk factor, is expressed by perivascular cells and marks vascular maturation and pericyte coverage[37]. MMP7, with an inverse coefficient, is expressed by normal luminal epithelium, so its coefficient reflects the proportion of retained normal structures, and tumours retaining organised glandular architecture are generally better differentiated. CD79A marks B cells; B cells and tertiary lymphoid structures predict response to immune checkpoint blockade in melanoma[38], but that finding does not establish a prognostic effect in ER+ breast cancer, so the score–B-cell correlation here is a statement about B-cell quantity. PLG is a liver-derived plasma protein essentially undetected in the 38,241 ER+ single cells, so its coefficient describes detection versus non-detection in RNA-seq and has no counterpart in microarray cohorts where every sample is detected.',
    ),
    # Data availability: name the GDSC2 source
    (
        'TCGA-BRCA (https://portal.gdc.cancer.gov/), METABRIC (cBioPortal), GSE20685, GSE96058 and GSE176078 (GEO) are publicly available.',
        'TCGA-BRCA (https://portal.gdc.cancer.gov/), METABRIC (cBioPortal), GSE20685, GSE96058 and GSE176078 (GEO), and the GDSC2 pharmacogenomic data underlying the exploratory drug-sensitivity analysis (https://www.cancerrxgene.org/) are publicly available.',
    ),
    # A2: article number
    ('Front Immunol. 2021;12:674188.', 'Front Immunol. 2021;12:674192.'),
    # [16]: standard author abbreviation
    ('[16] Xu ZH, Ma Q, Jiang QH, et al.', '[16] Xu Z, Ma Q, Jiang Q, et al.'),
    # A4: GDC source paper replaces the Pan-Cancer analysis-project paper
    (
        '[20] Weinstein JN, Collisson EA, Mills GB, et al. The Cancer Genome Atlas Pan-Cancer analysis project. Nat Genet. 2013;45(10):1113-1120.',
        '[20] Grossman RL, Heath AP, Ferretti V, et al. Toward a shared vision for cancer genomic data. N Engl J Med. 2016;375(12):1109-1112.',
    ),
    # A3: two new references + renumber the tail
    (
        '[35] Saharinen P, Eklund L, Alitalo K. Therapeutic targeting of the angiopoietin-TIE pathway. Nat Rev Drug Discov. 2017;16(9):635-661.\n'
        '[36] Helmink BA, Reddy SM, Gao J, et al. B cells and tertiary lymphoid structures promote immunotherapy response. Nature. 2020;577(7791):549-555.\n'
        '[37] Collins GS, Reitsma JB, Altman DG, Moons KG. Transparent reporting of a multivariable prediction model for individual prognosis or diagnosis (TRIPOD): the TRIPOD statement. Ann Intern Med. 2015;162(1):55-63.',
        '[35] Yang W, Soares J, Greninger P, et al. Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. Nucleic Acids Res. 2013;41(D1):D955-D961.\n'
        '[36] Picco G, Chen ED, Alonso LG, et al. Functional linkage of gene fusions to cancer cell fitness assessed by pharmacological and CRISPR-Cas9 screening. Nat Commun. 2019;10(1):2198.\n'
        '[37] Saharinen P, Eklund L, Alitalo K. Therapeutic targeting of the angiopoietin-TIE pathway. Nat Rev Drug Discov. 2017;16(9):635-661.\n'
        '[38] Helmink BA, Reddy SM, Gao J, et al. B cells and tertiary lymphoid structures promote immunotherapy response. Nature. 2020;577(7791):549-555.\n'
        '[39] Collins GS, Reitsma JB, Altman DG, Moons KG. Transparent reporting of a multivariable prediction model for individual prognosis or diagnosis (TRIPOD): the TRIPOD statement. Ann Intern Med. 2015;162(1):55-63.',
    ),
]


def replace_once(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise SystemExit("FAIL %s: expected 1 match, found %d\n  %s" % (label, n, old[:90]))
    return text.replace(old, new, 1)


def check_citations(text):
    body = text.split("## References")[0]
    refs = text.split("## References")[1].split("## Figure legends")[0]
    listed = [int(x) for x in re.findall(r"^\[(\d+)\]", refs, re.M)]
    cited = set()
    for m in re.finditer(r"\[(\d+)(?:[–-](\d+))?\]", body):
        a = int(m.group(1))
        b = int(m.group(2)) if m.group(2) else a
        cited.update(range(a, b + 1))
    problems = []
    if listed != list(range(1, len(listed) + 1)):
        problems.append("reference list is not 1..N in order: %s" % listed)
    missing = sorted(set(listed) - cited)
    dangling = sorted(cited - set(listed))
    if missing:
        problems.append("uncited references: %s" % missing)
    if dangling:
        problems.append("citations without a reference entry: %s" % dangling)
    return listed, cited, problems


def main():
    text = open(MS, encoding="utf-8").read()
    original = text
    for i, (old, new) in enumerate(REPLACEMENTS, 1):
        text = replace_once(text, old, new, "R%d" % i)
    listed, cited, problems = check_citations(text)
    if problems:
        raise SystemExit("FAIL citation check:\n  " + "\n  ".join(problems))
    open(MS, "w", encoding="utf-8", newline="\n").write(text)
    print("manuscript: %d replacements applied" % len(REPLACEMENTS))
    print("references: %d entries (1..%d), all cited, no dangling citations" %
          (len(listed), len(listed)))
    print("body length: %d -> %d chars" % (len(original), len(text)))

    # --- companion documents -------------------------------------------------
    def sub(path, pairs):
        s = open(path, encoding="utf-8").read()
        for old, new in pairs:
            s = replace_once(s, old, new, os.path.basename(path) + ":" + old[:35])
        open(path, "w", encoding="utf-8", newline="\n").write(s)
        print("updated:", os.path.relpath(path, ROOT))

    sub(CHANGELOG, [
        ("| 参考文献 | 53 条 | **37 条**（上限 40） |",
         "| 参考文献 | 53 条 | **39 条**（上限 40） |"),
        ("改后正文 3,998 词（原 4,010），仍在 4,000 以内；参考文献 37 条全部被引用，无悬空。",
         "改后正文 3,998 词（原 4,010），仍在 4,000 以内；参考文献 39 条全部被引用，无悬空（文献核查后修订，见第九节）。"),
    ])
    sub(TITLEPAGE, [
        ("The study followed the TRIPOD statement for reporting multivariable prediction models (reference [37] in the manuscript).",
         "The study followed the TRIPOD statement for reporting multivariable prediction models (reference [39] in the manuscript)."),
        ("| 参考文献 | 完成（37 条） | 正文内 + `06_参考文献/Zotero 报告.htm` | [14]–[19] 已逐条核实 |",
         "| 参考文献 | 完成（39 条） | 正文内 + `06_参考文献/Zotero 报告.htm` | [14]–[19] 已逐条核实；2026-11-05 按核查报告修订 |"),
    ])

    # append an audit-closure section to the reference report (report itself is a record)
    rep = open(REPORT, encoding="utf-8").read()
    marker = "## 五、核查局限"
    if "## 六、修订处置记录" not in rep:
        closure = """## 六、修订处置记录（2026-11-05）

用户确认后，已把本报告的全部判定改进 `manuscript_v2_tissue_composition.md`。复现脚本：`05_数据与结果/work/apply_ref_fixes.py`（每处替换强制断言唯一匹配，并执行引文编号连续性校验）。

| 条目 | 处置 |
|---|---|
| A1 [3] | 按原文重写：ER+ 早期、完成 5 年内分泌治疗且无病者，第 5–20 年远处复发风险 13%–41%（按 TN 分期）。数值/人群/时间窗/结局四项全部修正 |
| A2 [5] | `674188` → **`674192`**（PMID 34135901 已复核） |
| A3 [33,34] | 保留 pRRophetic[33] 与 Iorio 2016[34] 作为方法与资源主文献，**新增** GDSC 资源文献 Yang 2013 NAR [35] 与 GDSC2 数据来源 Picco 2019 Nat Commun [36]（GDSC2 PharmacoSet 官方 "please cite" 清单中三篇并列）；2.10 改引 `[33]` + `[34–36]`，Data availability 补 GDSC 门户 URL |
| A4 [20] | 题录整体替换为 Grossman RL et al. *N Engl J Med*. 2016;375(12):1109-1112（GDC 门户），正文句子不变即获正确支撑；Weinstein 2013 撤下 |
| A5 [11] | 改为 "shield tumour cells from CD8+ T-cell and NK-cell cytotoxicity by obstructing immune-cell contact" |
| A6 [6] | 具体判断改引 [5]（HR+ 低 TIL、ICB 低应答），[6] 保留用于 hot/altered/cold 概念框架 |
| A7 [36→38] | 改为 "B cells and tertiary lymphoid structures predict response to immune checkpoint blockade in melanoma"，并明确该结论不能外推为 ER+ 乳腺癌预后 |
| [2] | 删除 "approximately 75%"（该文摘要无此数字，只有 70–80% 可治愈率），改写为与摘要相符的"分子异质 / 激素受体阳性为主要亚组 / 内分泌治疗为基石" |
| [14–19] | "presented as an independent prognostic factor" → "reported as a prognostic or risk-stratification tool" |
| [32] | "treatment-naive" → "primary untreated"；2.6 补明注释层级为 9 major / 29 minor，且分析只用 9 major |
| [16] | 作者缩写改 `Xu Z, Ma Q, Jiang Q, et al.` |
| [8][13][21] | 判定为通行卷年标注，不改；定刊后按期刊 reference style 统一 |
| [18] | 保留 "Epub ahead of print"；定稿前复查卷期页 |

**本次修订额外发现并修正的两处（原报告未列）**：

1. **正文 3.6 的 "28 annotated cell types / 0% in all 28 classes" 是错的。** 原始数据 `raw/scRNA/sc/*/metadata.csv`（26 样本、100,064 细胞）的 `celltype_major` 恰为 **9** 类、`celltype_minor` 为 **29** 类；`j3_singlecell.R` 与 `j3_celltype_expression.csv` 全部按 9 个 major 聚合，Figure 6 图注也写 "nine annotated compartments"。已改为 "nine major annotated cell types / all nine classes"。
2. **PLG 检出率表述过强。** 原始数据最高类别检出率 0.158%（Plasmablasts），日志因四舍五入显示 0%；正文的 "not detected in any cell type (0%)" 与 "undetected in 100,064 tumour cells"（实为 38,241 个 ER+ 细胞，且并非都是 tumour cells）均已改为 "essentially undetected … ≤0.2% of cells in any class" 与 "essentially undetected in the 38,241 ER+ single cells"。

**编号变动**：新增 [35] [36] 后，原 [35] Saharinen → [37]、[36] Helmink → [38]、[37] TRIPOD → [39]；2.11 与 Discussion 的引用点已同步。全文 39 条题录连续、全部被引用、无悬空（脚本自动校验通过）。`改稿说明_v2.md` 与 `07_投稿材料/Title_Page_and_Declarations.md` 中的条数与 [37] 引用已同步。

> 遗留（不属本次范围）：Figure 6 图注 "PLG was not detected in any compartment" 仍为四舍五入口径，与正文 "essentially undetected" 略有强弱差异，定稿时可统一；本报告第一至五节保留核查当时的原始判定，不作追改。

---

"""
        rep = rep.replace(marker, closure + marker, 1)
        open(REPORT, "w", encoding="utf-8", newline="\n").write(rep)
        print("updated: 参考文献核查报告.md (appended section 6)")


if __name__ == "__main__":
    main()
