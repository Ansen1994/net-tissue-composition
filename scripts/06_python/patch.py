import re, pathlib, sys
MAP = {
 "brca.gdc_2024.rda": "../raw/TCGA/brca.gdc_2024.rda",
 "GSE20685_series_matrix.txt.gz": "../raw/GSE20685/GSE20685_series_matrix.txt.gz",
 "GPL570-55999.txt": "../raw/annotation/GPL570-55999.txt",
 "GeneCards_NETs.csv": "../raw/annotation/GeneCards_NETs.csv",
 "data_mrna_illumina_microarray.txt": "../raw/METABRIC/data_mrna_illumina_microarray.txt",
 "data_clinical_patient.txt": "../raw/METABRIC/data_clinical_patient.txt",
 "GSE96058-GPL11154_series_matrix.txt.gz": "../raw/GSE96058/GSE96058-GPL11154_series_matrix.txt.gz",
 "GSE96058-GPL18573_series_matrix.txt.gz": "../raw/GSE96058/GSE96058-GPL18573_series_matrix.txt.gz",
 "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz": "../raw/GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",
 "GDSC2_Expr (RMA Normalized and Log Transformed).rds": "../raw/drug/GDSC2_Expr (RMA Normalized and Log Transformed).rds",
 "GDSC2_Res.rds": "../raw/drug/GDSC2_Res.rds",
 "Breast_estimate_score.gct": "../raw/drug/Breast_estimate_score.gct",
}
for p in sorted(pathlib.Path("work/scripts").glob("*.R")):
    t = p.read_text(encoding="utf-8")
    orig = t
    # 注释掉 setwd 行
    t = re.sub(r"^\s*setwd\(.*\)\s*$", "# [patched] setwd removed; cwd is work/", t, flags=re.M)
    hits = []
    for old, new in MAP.items():
        # 只替换作为完整字符串字面量出现的原文件名
        pat = '"' + re.escape(old) + '"'
        if re.search(pat, t):
            t = re.sub(pat, '"' + new + '"', t)
            hits.append(old)
    if t != orig:
        p.write_text(t, encoding="utf-8")
    print(f"{p.name}: 路径替换 {len(hits)} 处 -> {hits}")
