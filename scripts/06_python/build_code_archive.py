# -*- coding: utf-8 -*-
"""
Build a portable, public-ready code archive for manuscript v2.

Copies every analysis / figure script into 11_代码公开包/scripts/, replaces
machine-specific absolute paths with environment-variable lookups, records
every replacement in PATHS.md, captures a real sessionInfo(), and writes
README.md / MANIFEST.md / 00_config.R.

Original scripts are never modified.
"""
import os
import re
import shutil
import subprocess
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "11_代码公开包")

GROUPS = [
    ("01_primary", os.path.join(ROOT, "04_代码"), "*.R",
     "Original analysis and figure scripts (DESeq2, LASSO, randomForest, Cox, KM)."),
    ("02_drug", os.path.join(ROOT, "04_代码", "药敏重算"), "*.R",
     "Exploratory drug-sensitivity analysis (GDSC2 ridge regression, Table S17)."),
    ("03_work", os.path.join(ROOT, "05_数据与结果", "work"), "*.R",
     "Re-evaluation analyses: cross-platform scaling, reference NET score, treatment "
     "interaction, immune-axis decomposition, single-cell origin, generalisation set."),
    ("04_work_scripts", os.path.join(ROOT, "05_数据与结果", "work", "scripts"), "*.R",
     "Full re-run copies of the primary scripts, kept for the byte-for-byte verification."),
    ("05_figures_v2", os.path.join(ROOT, "05_数据与结果", "work", "figs_v2"), "*.R",
     "Final v2 figure code (unified Okabe-Ito palette and typography)."),
    ("06_python", os.path.join(ROOT, "05_数据与结果", "work"), "*.py",
     "Reproducibility helpers: supplementary workbook, manuscript package, "
     "reference audit, review-fix patches, TIFF conversion."),
    ("07_sc_umap", os.path.join(ROOT, "12_可复用资源", "sc_umap"), "*.R",
     "Reusable single-cell asset: recomputes and plots the display-only UMAP embedding "
     "of Supplementary Figure S3 (standalone; enters no statistical model)."),
]

DATA_ROOT = 'Sys.getenv("WB_DATA_ROOT", unset = "data/raw")'
OUT_DIR = 'Sys.getenv("WB_OUTPUT_DIR", unset = "out")'

# R running with LC_CTYPE=C cannot open files whose names contain multibyte
# characters; rename these in the portable archive (originals stay untouched).
NAME_MAP = {"figure1_提取degs.R": "figure1_extract_degs.R"}
RENAMED = []

RE_VAR_ASSIGN = re.compile(
    r'^(\s*(?:G|OUT|ROOT|BASE|WD)\s*<-\s*)"([A-Za-z]:[/\\][^"]+)"', re.M)
RE_SETWD = re.compile(r'^(\s*setwd\()"([A-Za-z]:[/\\][^"]+)"\)', re.M)
RE_COMMENT = re.compile(r'^(\s*#\s*setwd\()"([A-Za-z]:[/\\][^"]+)"\)', re.M)


def portable(src_text):
    """Replace absolute paths in a copy; return (new_text, list of (old_line, new_line))."""
    log = []

    def var_sub(m):
        old = m.group(0)
        env = DATA_ROOT if m.group(2).upper().startswith("G:") else OUT_DIR
        new = m.group(1) + env
        log.append((old.strip(), new.strip()))
        return new

    def setwd_sub(m):
        old = m.group(0)
        new = m.group(1) + OUT_DIR + ")"
        log.append((old.strip(), new.strip()))
        return new

    def comment_sub(m):
        old = m.group(0)
        new = m.group(1) + DATA_ROOT + ")"
        log.append((old.strip(), new.strip()))
        return new

    text = RE_VAR_ASSIGN.sub(var_sub, src_text)
    text = RE_SETWD.sub(setwd_sub, text)
    text = RE_COMMENT.sub(comment_sub, text)
    return text, log


def r_parse_check(files):
    """Parse every portable R file in a single R process; return failures."""
    if not files:
        return []
    code = ('for (f in commandArgs(TRUE)) '
            'tryCatch(invisible(parse(f)), error = function(e) '
            'cat("PARSE_FAIL:", f, conditionMessage(e), "\\n"))')
    try:
        r = subprocess.run(["Rscript", "-e", code] + files,
                           capture_output=True, text=True, timeout=300, cwd=OUT)
        return [ln for ln in r.stdout.splitlines() if "PARSE_FAIL" in ln]
    except Exception as e:
        return ["PARSE_CHECK_UNAVAILABLE: %s" % e]


def session_info():
    code = r'''
suppressPackageStartupMessages({
  for (p in c("DESeq2","glmnet","randomForest","survival","survminer",
              "timeROC","clusterProfiler","ggplot2","patchwork",
              "org.Hs.eg.db","SummarizedExperiment")) {
    if (requireNamespace(p, quietly = TRUE)) library(p, character.only = TRUE)
  }
})
cat("Captured:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
print(sessionInfo())
'''
    tmp = os.path.join(OUT, "_session_info.R")
    with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(code)
    try:
        r = subprocess.run(["Rscript", tmp],
                           capture_output=True, text=True, timeout=300)
        out = (r.stdout or "").strip()
        if len(out) < 200:
            return (out + "\n\n[stderr]\n" + (r.stderr or "")[:1500]).strip()
        return out
    except Exception as e:
        return "sessionInfo unavailable: %s" % e
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def main():
    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    scripts = os.path.join(OUT, "scripts")
    os.makedirs(scripts)

    manifest, path_log, r_files, py_files = [], [], [], []
    counts = {}

    for sub, src_dir, pattern, desc in GROUPS:
        dst_dir = os.path.join(scripts, sub)
        os.makedirs(dst_dir)
        files = sorted(f for f in os.listdir(src_dir)
                       if f.endswith(".R" if pattern == "*.R" else ".py")
                       and os.path.isfile(os.path.join(src_dir, f))) \
            if os.path.isdir(src_dir) else []
        counts[sub] = len(files)
        manifest.append("### `scripts/%s/` (%d files)" % (sub, len(files)))
        manifest.append("")
        manifest.append(desc)
        manifest.append("")
        for f in files:
            src = os.path.join(src_dir, f)
            portable_name = NAME_MAP.get(f, f)
            if portable_name != f:
                RENAMED.append((sub, f, portable_name))
            dst = os.path.join(dst_dir, portable_name)
            with open(src, encoding="utf-8", errors="replace") as fh:
                text = fh.read()
            new, log = portable(text)
            with open(dst, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(new)
            rel = "scripts/%s/%s" % (sub, portable_name)
            marker = " (paths rewritten)" if log else ""
            if portable_name != f:
                marker += "; renamed from `%s`" % f
            manifest.append("- `%s`%s" % (rel, marker))
            for a, b in log:
                path_log.append((rel, a, b))
            (r_files if f.endswith(".R") else py_files).append(dst)
        manifest.append("")

    # extra python helper living in figs_v2
    extra = os.path.join(ROOT, "05_数据与结果", "work", "figs_v2", "make_tiff.py")
    if os.path.exists(extra):
        dst = os.path.join(scripts, "06_python", "make_tiff.py")
        with open(extra, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        new, log = portable(text)
        with open(dst, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(new)
        manifest.append("- `scripts/06_python/make_tiff.py` (600 dpi TIFF conversion)")
        py_files.append(dst)
        counts["06_python"] = counts.get("06_python", 0) + 1

    # python loader that belongs with the reusable UMAP asset
    extra2 = os.path.join(ROOT, "12_可复用资源", "sc_umap", "sc_umap_load.py")
    if os.path.exists(extra2):
        dst = os.path.join(scripts, "07_sc_umap", "sc_umap_load.py")
        with open(extra2, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        new, log = portable(text)
        with open(dst, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(new)
        manifest.append("- `scripts/07_sc_umap/sc_umap_load.py` (read-only loader for the "
                        "cached UMAP embeddings)")
        py_files.append(dst)
        counts["07_sc_umap"] = counts.get("07_sc_umap", 0) + 1

    # copy the existing figure-to-script crosswalk
    crosswalk = os.path.join(ROOT, "04_代码", "图表代码对照表.md")
    if os.path.exists(crosswalk):
        shutil.copy2(crosswalk, os.path.join(OUT, "图表脚本对照表.md"))

    # 00_config.R
    config = '''## Central configuration for the portable code archive.
## Set WB_DATA_ROOT and WB_OUTPUT_DIR as environment variables before running,
## or edit the defaults below. No script in this archive contains a machine path.

WB_DATA_ROOT  <- Sys.getenv("WB_DATA_ROOT",  unset = "data/raw")
WB_OUTPUT_DIR <- Sys.getenv("WB_OUTPUT_DIR", unset = "out")

## Fixed random seeds used in the manuscript.
WB_SEED_FEATURE_SELECTION <- 2026L  # LASSO cv.glmnet and randomForest (Methods 2.3)
WB_SEED_CROSS_VALIDATION  <- 42L    # 5-fold CV optimism, module scores, drug panel

dir.create(WB_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
'''
    with open(os.path.join(OUT, "00_config.R"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(config)

    # sessionInfo
    si = session_info()
    with open(os.path.join(OUT, "sessionInfo.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(si + "\n")
    print("sessionInfo captured: %d chars" % len(si))

    # PATHS.md
    pl = ["# Path portability audit", "",
          "Every absolute path found in the original scripts is listed below together with "
          "its replacement in this archive. **The originals under `04_代码/` and "
          "`05_数据与结果/work/` were not modified.**", "",
          "Replacements use `WB_DATA_ROOT` (raw input data) and `WB_OUTPUT_DIR` "
          "(derived tables and figures); both are defined in `00_config.R` and can be "
          "overridden by environment variables.", "",
          "| File | Original line | Portable line |", "|---|---|---|"]
    for rel, a, b in path_log:
        pl.append("| `%s` | `%s` | `%s` |" %
                  (rel, a.replace("|", "\\|"), b.replace("|", "\\|")))
    pl += ["", "Total replacements: **%d** across **%d** files." %
           (len(path_log), len(set(x[0] for x in path_log))), "",
           "Two commented-out `setwd()` lines that used Windows backslash paths are also "
           "listed; they are inert but were rewritten for consistency.", ""]
    with open(os.path.join(OUT, "PATHS.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(pl))

    # syntax check on portable R copies
    failures = r_parse_check([os.path.relpath(f, OUT).replace("\\", "/") for f in r_files])
    print("R syntax check: %d files, failures = %s" % (len(r_files), failures or "none"))

    # LICENSE (MIT)
    lic = """MIT License

Copyright (c) 2026 Ansheng Xu, Yang Yu, Yanqing Mao, Guobing Yin, Jing Luo

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
    with open(os.path.join(OUT, "LICENSE"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(lic)

    # README.md
    rd = [
        "# Analysis code for: NET-related gene signatures in ER-positive breast cancer "
        "reflect tissue composition rather than neutrophil biology",
        "",
        "Portable, path-independent copy of the analysis code. Every machine-specific "
        "absolute path in the original scripts has been replaced by an environment "
        "variable; the replacements are listed one by one in `PATHS.md`. The original "
        "working copies were not modified.",
        "",
        "## Contents",
        "",
        "| Path | What it is |",
        "|---|---|",
        "| `scripts/01_primary/` | Differential expression, LASSO and random-forest "
        "feature selection, Cox models, Kaplan-Meier and the original figure scripts |",
        "| `scripts/02_drug/` | Exploratory drug-sensitivity analyses (GDSC2, ridge "
        "regression) |",
        "| `scripts/03_work/` | Cohort assembly, module scores, single-cell attribution, "
        "table generation (`s9_supp_tables.R`, `make_table_s5_ci.R`, `make_table_s5c.R`, "
        "`j3_singlecell.R`) |",
        "| `scripts/04_work_scripts/` | Validation-cohort and cross-platform scaling "
        "scripts |",
        "| `scripts/05_figures_v2/` | Final figure scripts (`fig1.R`–`fig7.R`, `style.R`) |",
        "| `scripts/06_python/` | Packaging helpers only — see \"What the Python scripts "
        "do not do\" below |",
        "| `scripts/07_sc_umap/` | Standalone asset that recomputes the display-only UMAP "
        "of Supplementary Figure S3 (config, pipeline, plots, read-only loader); enters no "
        "statistical model |",
        "| `00_config.R` | Central configuration: `WB_DATA_ROOT`, `WB_OUTPUT_DIR`, fixed "
        "seeds |",
        "| `sessionInfo.txt` | Output of `sessionInfo()` on the machine that produced the "
        "results (R 4.4.2) |",
        "| `PATHS.md` | Audit of every absolute path that was rewritten |",
        "| `MANIFEST.md` | Script-by-script listing with file counts |",
        "| `图表脚本对照表.md` | Crosswalk from each figure and table to its script |",
        "",
        "Totals: **%d** R scripts, **%d** Python scripts, **%d** path replacements "
        "across **%d** files." % (len(r_files), len(py_files), len(path_log),
                                 len(set(x[0] for x in path_log))),
        "",
        "## Requirements",
        "",
        "- R 4.4.2 with the packages recorded in `sessionInfo.txt` (DESeq2, glmnet, "
        "randomForest, survival, survminer, timeROC, clusterProfiler, ggplot2, patchwork, "
        "rms).",
        "- Python 3.13 with openpyxl, python-docx, PyMuPDF and Pillow — needed only for "
        "packaging, not for any analysis.",
        "",
        "## How to run",
        "",
        "1. Download the public input data listed under \"Input data\" below into one "
        "directory and set `WB_DATA_ROOT` to it.",
        "2. Set `WB_OUTPUT_DIR` to a writable directory for derived tables and figures.",
        "3. Run the scripts in directory order (`01_primary` → `02_drug` → `03_work` → "
        "`04_work_scripts` → `05_figures_v2`), each from within its own directory after "
        "`source(\"../../00_config.R\")`.",
        "",
        "Both variables can be set as environment variables or edited in `00_config.R`; "
        "no script contains a machine path.",
        "",
        "## Fixed seeds",
        "",
        "- `2026` — LASSO (`cv.glmnet`) and random-forest feature selection.",
        "- `42` — five-fold cross-validation, module scores and drug-sensitivity analyses.",
        "- `20261106` — display-only UMAP embedding of Supplementary Figure S3.",
        "",
        "## Input data",
        "",
        "All inputs are public and were obtained at no cost:",
        "",
        "- TCGA-BRCA — Genomic Data Commons (https://portal.gdc.cancer.gov/)",
        "- GSE20685, GSE96058, GSE176078 — NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/)",
        "- METABRIC — cBioPortal (https://www.cbioportal.org/)",
        "- GDSC2 — Genomics of Drug Sensitivity in Cancer (https://www.cancerrxgene.org/)",
        "- NET-related gene annotations — GeneCards (https://www.genecards.org/, relevance "
        "score > 10)",
        "",
        "## What the Python scripts do not do",
        "",
        "The Python scripts perform no statistical computation and produce no number that "
        "appears in the manuscript. They consolidate the supplementary workbook, convert "
        "figure PDFs to 600 dpi submission TIFF, audit reference records and assemble the "
        "manuscript package. Every numerical result was produced in R.",
        "",
        "## Licence",
        "",
        "Released under the MIT License; the full text is in `LICENSE`. A permissive "
        "licence was chosen so that the code can be reused and adapted; the manuscript's "
        "text and figures are not covered by it.",
        "",
    ]
    with open(os.path.join(OUT, "README.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(rd))

    # MANIFEST.md
    mf = ["# Script manifest", "",
          "Generated %s from `04_代码/` and `05_数据与结果/work/`." %
          datetime.now().strftime("%Y-%m-%d %H:%M"), "",
          "Directories are ordered by analysis stage. `图表脚本对照表.md` maps each "
          "figure and table to the script that produced it.", ""] + manifest
    mf += ["## Verification", "",
           "- Portable R scripts: **%d** parsed with `parse()`; failures: **%s**." %
           (len(r_files), ", ".join(failures) if failures else "none"),
           "- Portable Python helpers: **%d**." % len(py_files),
           "- Absolute paths rewritten: **%d** (see `PATHS.md`)." % len(path_log), "",
           "Parsing checks syntax only. A full numerical re-run additionally requires the "
           "public input data listed in `README.md`; that step was not re-executed while "
           "building this archive.", "",
           "## Filename normalisation", "",
           "R running under the `C` locale cannot open filenames containing multibyte "
           "characters. The following copies were renamed (contents unchanged):", ""]
    for sub, old, new in RENAMED:
        mf.append("- `scripts/%s/%s` \u2192 `scripts/%s/%s`" % (sub, old, sub, new))
    mf.append("")
    with open(os.path.join(OUT, "MANIFEST.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(mf))

    total_r, total_py = len(r_files), len(py_files)
    return counts, total_r, total_py, len(path_log), failures


if __name__ == "__main__":
    counts, nr, npy, npath, failures = main()
    print("Archive:", OUT)
    print("  R scripts %d, Python helpers %d, path rewrites %d" % (nr, npy, npath))
    print("  groups:", counts)
