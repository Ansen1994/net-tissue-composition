# -*- coding: utf-8 -*-
"""
Apply the second-round peer-review revisions to manuscript v2.

Covers: code availability, random-forest seed, drug-sensitivity data and
method wording, 10x detection caveat, cross-platform mechanism in the
abstract, PLG blood-content interpretation, and R version correction.

Each replacement asserts a unique match, so re-running after a source edit
fails loudly instead of silently producing a partially fixed manuscript.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MS = os.path.join(ROOT, "01_手稿", "_md", "manuscript_v2_tissue_composition.md")

FIXES = [
    # 1. Abstract: cross-platform heterogeneity as a mechanism of non-transferability
    (
        "Abstract / cross-platform mechanism",
        "**Conclusions:** In ER+ breast cancer these signatures index tissue composition rather than neutrophil biology, "
        "and bulk transcriptomes lack the substrate to measure NET formation. Scores of this type should be interpreted accordingly.",
        "**Conclusions:** In ER+ breast cancer these signatures index tissue composition rather than neutrophil biology, "
        "and bulk transcriptomes lack the substrate to measure NET formation. Heterogeneous detection limits across "
        "platforms \u2014 PLG was zero in 53.9% of TCGA RNA-seq samples but a near-constant in METABRIC microarrays \u2014 "
        "further limit transferability. Scores of this type should be interpreted accordingly.",
    ),
    # 2. Methods 2.3: state the fixed random-forest seed and the version caveat
    (
        "Methods 2.3 / random-forest seed",
        "random forest[26] ranked features by MeanDecreaseGini on binary survival status with 500 trees and retained "
        "those above mean importance. The signature was the intersection of both selections.",
        "random forest[26] ranked features by MeanDecreaseGini on binary survival status with 500 trees and retained "
        "those above mean importance. Both selectors were run under a fixed seed (set.seed(2026)); because Gini "
        "importance is sensitive to the randomForest build, the released 58-feature set (Table S4) is treated as "
        "definitive and the plotted importances were recomputed within that fixed set. The signature was the "
        "intersection of both selections.",
    ),
    # 3. Methods 2.10: describe the implementation actually used (glmnet, not the pRRophetic package)
    (
        "Methods 2.10 / actual ridge implementation",
        "**2.10. Drug sensitivity (exploratory).** Predicted sensitivity to 29 compounds was estimated by ridge "
        "regression (pRRophetic[33]) trained on the GDSC2 release of the Genomics of Drug Sensitivity in Cancer "
        "resource[34\u201336], comparing five preprocessing strategies and applying leave-one-out common-axis "
        "adjustment; predictions are reported as exploratory only.",
        "**2.10. Drug sensitivity (exploratory).** Predicted log IC50 for 29 compounds was estimated by ridge "
        "regression implemented directly in glmnet (alpha = 0, five-fold cross-validation) and trained on the GDSC2 "
        "expression and drug-response matrices of the Genomics of Drug Sensitivity in Cancer resource[34\u201336]. "
        "This reproduces the ridge-regression principle of pRRophetic[33] but is a separate implementation, so "
        "predicted IC50 magnitudes are not directly comparable with published pRRophetic output. Five preprocessing "
        "strategies were compared and a leave-one-out common-axis adjustment was applied; full results are given in "
        "Table S17 and all predictions are treated as exploratory.",
    ),
    # 4. Methods 2.11: R version measured on this machine
    (
        "Methods 2.11 / R version",
        "Analyses used R 4.4.3.",
        "Analyses used R 4.4.2; package versions are recorded in sessionInfo.txt accompanying the analysis code.",
    ),
    # 5. Results 3.6: front-load the 10x neutrophil-capture caveat
    (
        "Results 3.6 / 10x caveat up front",
        "The same data show why a NET score cannot be anchored in these tumours. PADI4, the enzyme executing histone "
        "citrullination, was detected in 0.3% of all cells",
        "The same data show why a NET score cannot be anchored in these tumours. Because droplet-based 10x capture "
        "systematically under-represents neutrophils (see Limitations), the following rates describe a detection "
        "floor rather than proof of absent neutrophils. PADI4, the enzyme executing histone citrullination, was "
        "detected in 0.3% of all cells",
    ),
    # 6. Results 3.6: report the actual drug-sensitivity numbers
    (
        "Results 3.6 / concrete drug predictions",
        "Predicted drug-sensitivity differences were likewise attributable to the composition axis: after common-axis "
        "adjustment the only stable differences involved anti-HER2 kinase inhibitors and fulvestrant, and their "
        "direction followed the stromal and epithelial content of the risk groups. These predictions are reported "
        "for completeness only.",
        "Predicted drug-sensitivity differences likewise followed the composition axis (Table S17). In the primary "
        "analysis only two compounds survived Benjamini\u2013Hochberg correction: high-risk tumours were predicted more "
        "resistant to sapitinib (difference in log IC50 +0.321, FDR = 9.9\u00d710\u207b\u2079, Cohen's d = 0.47) and more "
        "sensitive to fulvestrant (\u22120.109, FDR = 0.017, d = \u22120.267). After leave-one-out common-axis adjustment "
        "three compounds retained a consistent direction across all five preprocessing modes \u2014 sapitinib (median "
        "d = +0.360), lapatinib (+0.331) and fulvestrant (\u22120.338) \u2014 while the remaining targeted compounds were "
        "directionally unstable or null (Table S17). Because the differences track the stromal and epithelial content "
        "of the risk groups and are attenuated under alternative preprocessing, they are reported as an illustration "
        "of the composition axis, not as pharmacological findings.",
    ),
    # 7. Discussion: PLG as blood content / perfusion, aligned with ANGPT1
    (
        "Discussion / PLG blood content",
        "PLG is a liver-derived plasma protein essentially undetected in the 38,241 ER+ single cells, so its "
        "coefficient describes detection versus non-detection in RNA-seq and has no counterpart in microarray "
        "cohorts where every sample is detected.",
        "PLG is a liver-derived plasma protein essentially undetected in the 38,241 ER+ single cells, so its bulk "
        "RNA-seq signal represents residual blood content and vascular perfusion of the excised specimen rather than "
        "transcription by any tumour-resident cell; its coefficient therefore describes detection versus "
        "non-detection, has no counterpart in microarray cohorts where every sample is detected, and is pathologically "
        "continuous with ANGPT1 as a second readout of vascularisation.",
    ),
    # 8. Data availability: replace local-path / on-request wording with a public archive
    (
        "Data availability / public code archive",
        "TCGA-BRCA (https://portal.gdc.cancer.gov/), METABRIC (cBioPortal), GSE20685, GSE96058 and GSE176078 (GEO), "
        "and the GDSC2 pharmacogenomic data underlying the exploratory drug-sensitivity analysis "
        "(https://www.cancerrxgene.org/) are publicly available. Processed matrices, result tables and R scripts "
        "supporting the findings are available from the corresponding author on reasonable request.",
        "TCGA-BRCA (https://portal.gdc.cancer.gov/), METABRIC (cBioPortal), GSE20685, GSE96058 and GSE176078 (GEO), "
        "and the GDSC2 pharmacogenomic data underlying the exploratory drug-sensitivity analysis "
        "(https://www.cancerrxgene.org/) are publicly available. All processed result tables supporting the figures "
        "and conclusions are provided in Supplementary Tables S1\u2013S17. The analysis code, run instructions and "
        "package-version record are archived in a public versioned repository at [REPOSITORY URL] (release [VERSION]; "
        "Zenodo DOI [DOI]).",
    ),
]


def check_citations(text):
    listed = [int(n) for n in re.findall(r"^\[(\d+)\]\s", text, re.M)]
    body = text.split("## References")[0]
    cited = set()
    for m in re.finditer(r"\[(\d+)(?:[\u2013-](\d+))?\]", body):
        a, b = int(m.group(1)), int(m.group(2) or m.group(1))
        cited.update(range(a, b + 1))
    problems = []
    if listed != list(range(1, len(listed) + 1)):
        problems.append("题录编号不连续: %s" % listed)
    missing = sorted(cited - set(listed))
    unused = sorted(set(listed) - cited)
    if missing:
        problems.append("引用了不存在的编号: %s" % missing)
    if unused:
        problems.append("未被正文引用的题录: %s" % unused)
    return listed, cited, problems


def word_count(text):
    abstract = text.split("## Abstract")[1].split("## 1. Introduction")[0]
    abstract = re.sub(r"\*\*[^*]+:\*\*", "", abstract.split("**Keywords:**")[0])
    body = text.split("## 1. Introduction")[1].split("## Data availability")[0]
    return len(abstract.split()), len(body.split())


def main():
    with open(MS, encoding="utf-8") as fh:
        text = fh.read()

    print("修订 %s" % os.path.relpath(MS, ROOT))
    before_abs, before_body = word_count(text)

    for label, old, new in FIXES:
        n = text.count(old)
        if n != 1:
            print("  [FAIL] %s \u2014 匹配 %d 次（应为 1 次），已中止" % (label, n))
            sys.exit(1)
        text = text.replace(old, new)
        print("  [OK]   %s" % label)

    listed, cited, problems = check_citations(text)
    abs_w, body_w = word_count(text)
    stale = [s for s in ["R 4.4.3", "reported for completeness only",
                         "available from the corresponding author on reasonable request",
                         "ridge regression (pRRophetic[33])"]
             if s in text]

    print("\n校验:")
    print("  题录 %d 条连续, 正文引用编号 %d 个, 问题: %s" %
          (len(listed), len(cited), problems or "无"))
    print("  表引用 S17: %s" % ("Table S17" in text))
    print("  摘要词数 %d -> %d (+%d); 正文词数 %d -> %d (+%d)" %
          (before_abs, abs_w, abs_w - before_abs, before_body, body_w, body_w - before_body))
    print("  旧表述残留: %s" % (stale or "无"))

    if problems or stale:
        print("  [FAIL] 校验未通过，未写入文件")
        sys.exit(1)

    with open(MS, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print("\n[written] %s" % MS)


if __name__ == "__main__":
    main()
