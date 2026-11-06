# -*- coding: utf-8 -*-
"""
Build a preliminary, self-contained manuscript package for manuscript v2.

Inputs (read-only):
  01_手稿/_md/manuscript_v2_tissue_composition.md
  07_投稿材料/Title_Page_and_Declarations.md
  02_主图/v2_组织成分版/Figure*.pdf
  07_投稿材料/figures_tiff/FigureN.tif
  03_补充材料/Supplementary_*.md / .xlsx

Outputs:
  10_初稿合订_v2/
      Manuscript_v2_draft.docx        editable master (text + legends + 7 figures)
      Manuscript_v2_draft.pdf         reading copy, same content
      Figure_Plates.pdf               7 main figures in one file
      Figures/                       7 submission TIFF + 7 vector PDF
      Supplementary/                 consolidated workbook + legends + notes
      README_初稿清单.md              package manifest and outstanding items

No journal-specific formatting is applied. Re-run after any source edit.
"""
import html
import os
import re
import shutil
import tempfile
from datetime import datetime
from pathlib import Path

import pymupdf
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

MS = ROOT / "01_手稿/_md/manuscript_v2_tissue_composition.md"
TP = ROOT / "07_投稿材料/Title_Page_and_Declarations.md"
FIG_PDF = ROOT / "02_主图/v2_组织成分版"
FIG_TIF = ROOT / "07_投稿材料/figures_tiff"
SUPP_DIR = ROOT / "03_补充材料"
OUT = ROOT / "10_初稿合订_v2"

FIGMAP = {
    1: "Figure1_ABCDEF.pdf",
    2: "Figure2_ABCD.pdf",
    3: "Figure3_ABCDEFG.pdf",
    4: "Figure4_ABCDEF.pdf",
    5: "Figure5_ABCD.pdf",
    6: "Figure6_ABCDEF.pdf",
    7: "Figure7_ABCD.pdf",
}
SUPP_FILES = [
    "Supplementary_Tables_S1-S17.xlsx",
    "Supplementary_Materials_v2.md",
    "Supplementary_Table_Legends.md",
    "Supplementary_Figure_Legends.md",
    "Supplementary Figure S3.pdf",
]

CITE = re.compile(r"\[\d+(?:[–-]\d+)?\]")
TODO = re.compile(r"\[[A-Z][^\]]{1,240}\]")
TOKEN = re.compile(r"(\*\*.+?\*\*|\[\d+(?:[–-]\d+)?\]|\[[A-Z][^\]]{1,240}\])")
LEGEND = re.compile(r"^\*\*Figure (\d+)\.")


# --------------------------------------------------------------------------
# Content model
# --------------------------------------------------------------------------
def read(path):
    return path.read_text(encoding="utf-8")


def parse_md(text):
    """Coarse markdown block parser; returns (kind, text) tuples."""
    out = []
    for line in text.splitlines():
        s = line.strip()
        if not s or s == "---":
            continue
        if s.startswith("### "):
            out.append(("h3", s[4:].strip()))
        elif s.startswith("## "):
            out.append(("h2", s[3:].strip()))
        elif s.startswith("# "):
            out.append(("h1", s[2:].strip()))
        elif s.startswith("> "):
            out.append(("note", s[2:].strip()))
        elif re.match(r"^\[\d+\] ", s):
            out.append(("ref", s))
        else:
            out.append(("p", s))
    return out


def title_page_blocks():
    sec = read(TP).split("## 1. Title page", 1)[1].split("\n---", 1)[0]
    out = []
    for line in sec.splitlines():
        s = line.strip()
        if not s or s == "---" or s.startswith("|") or s.startswith(">"):
            continue
        if s.startswith("**Counts at a glance**"):
            continue
        out.append(("p", s))
    return out


def declaration_blocks():
    sec = read(TP).split("## 2. Declarations", 1)[1].split("\n---", 1)[0]
    out = [("h2", "Declarations")]
    skip = False
    for kind, text in parse_md(sec):
        if kind == "h3" and text == "Data availability":
            skip = True
            out.append(("p", "**Data availability.** See the Data availability section above; "
                             "the two statements must be merged into one at journal formatting."))
            continue
        if kind == "h3":
            skip = False
        if skip or kind == "note":
            continue
        out.append((kind, text))
    return out


def build_blocks():
    blocks = [
        ("todo_notice",
         "PRELIMINARY DRAFT — %d %s. This consolidated file is for author review only. "
         "Items shown in [SQUARE BRACKETS] still require author input. Not for submission."
         % (datetime.now().day, datetime.now().strftime("%B %Y"))),
    ]
    blocks += title_page_blocks()
    blocks.append(("pagebreak", ""))

    ms = parse_md(read(MS))
    for kind, text in ms:
        if kind == "h1":
            continue  # title is on the title page
        if kind == "h2" and text == "References":
            blocks += declaration_blocks()
        if kind == "p" and LEGEND.match(text):
            kind = "legend"
        blocks.append((kind, text))

    blocks.append(("pagebreak", ""))
    blocks.append(("h2", "Figures"))
    blocks.append(("p", "Each figure is reproduced on its own page below. Standalone submission files "
                        "(600 dpi RGB/LZW TIFF) and vector PDFs are provided in the Figures folder."))
    for n in range(1, 8):
        if n > 1:
            blocks.append(("pagebreak", ""))
        blocks.append(("figlabel", "Figure %d" % n))
        blocks.append(("image", "Figure%d.jpg" % n))
    return blocks


# --------------------------------------------------------------------------
# Inline formatting
# --------------------------------------------------------------------------
def segments(text, cites=True, todos=True):
    """Yield (text, style) where style is plain/bold/sup/todo."""
    for part in TOKEN.split(text):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**") and len(part) > 4:
            yield part[2:-2], "bold"
        elif CITE.fullmatch(part) and cites:
            yield part, "sup"
        elif TODO.fullmatch(part) and todos:
            yield part, "todo"
        else:
            yield part, "plain"


def inline_html(text, cites=True, todos=True):
    out = []
    for seg, style in segments(text, cites, todos):
        e = html.escape(seg, quote=False)
        if style == "bold":
            out.append("<strong>%s</strong>" % e)
        elif style == "sup":
            out.append("<sup>%s</sup>" % e)
        elif style == "todo":
            out.append('<span class="todo">%s</span>' % e)
        else:
            out.append(e)
    return "".join(out)


CSS = """
body { font-family: "Times New Roman", serif; font-size: 10.5pt; line-height: 1.35; }
h1 { font-size: 16pt; text-align: center; line-height: 1.25; margin: 0 0 14pt 0; }
h2 { font-size: 12.5pt; margin: 16pt 0 8pt 0; }
h3 { font-size: 11pt; margin: 12pt 0 4pt 0; }
p  { margin: 0 0 7pt 0; text-align: justify; }
.ref { padding-left: 18pt; text-indent: -18pt; font-size: 9.5pt; margin-bottom: 4pt; }
.legend { font-size: 10pt; }
.note { color: #555555; font-size: 9.5pt; }
.todo { color: #b35000; font-weight: bold; }
.todo_notice { color: #b35000; font-size: 9.5pt; text-align: center; margin-bottom: 14pt; }
.figlabel { text-align: center; font-weight: bold; font-size: 11pt; margin: 0 0 6pt 0; }
img { width: 100%; height: auto; }
.pagebreak { page-break-before: always; }
"""


def blocks_to_html(blocks):
    out = []
    for kind, text in blocks:
        if kind == "pagebreak":
            out.append('<div class="pagebreak"></div>')
        elif kind == "image":
            out.append('<img src="%s">' % text)
        elif kind == "todo_notice":
            out.append('<p class="todo_notice">%s</p>' % inline_html(text, False, True))
        elif kind in ("h1", "h2", "h3"):
            out.append("<%s>%s</%s>" % (kind, inline_html(text, False, False), kind))
        elif kind == "figlabel":
            out.append('<p class="figlabel">%s</p>' % html.escape(text))
        elif kind == "note":
            out.append('<p class="note">%s</p>' % inline_html(text, False, True))
        elif kind == "ref":
            out.append('<p class="ref">%s</p>' % inline_html(text, False, False))
        elif kind == "legend":
            out.append('<p class="legend">%s</p>' % inline_html(text, True, False))
        else:
            out.append('<p>%s</p>' % inline_html(text, True, True))
    return "\n".join(out)


# --------------------------------------------------------------------------
# DOCX
# --------------------------------------------------------------------------
def set_cell_font(style, size, bold=False):
    style.font.name = "Times New Roman"
    style.font.size = Pt(size)
    style.font.bold = bold
    style._element.get_or_add_rPr().get_or_add_rFonts().set(qn("w:eastAsia"), "SimSun")


def add_runs(par, text, cites=True, todos=True):
    for seg, style in segments(text, cites, todos):
        r = par.add_run(seg)
        if style == "bold":
            r.bold = True
        elif style == "sup":
            r.font.superscript = True
        elif style == "todo":
            r.bold = True
            r.font.color.rgb = RGBColor(0xB3, 0x50, 0x00)
    return par


def build_docx(blocks, png_dir, path):
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(2.3)
    sec.top_margin = sec.bottom_margin = Cm(2.2)

    normal = doc.styles["Normal"]
    set_cell_font(normal, 10.5)
    normal.paragraph_format.line_spacing = 1.3
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    for name, size in (("Heading 1", 14), ("Heading 2", 12.5), ("Heading 3", 11)):
        st = doc.styles[name]
        set_cell_font(st, size, True)
        st.font.color.rgb = RGBColor(0, 0, 0)
        st.paragraph_format.space_before = Pt(12)
        st.paragraph_format.space_after = Pt(5)
        st.paragraph_format.keep_with_next = True

    footer = sec.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    fr = footer.add_run("Preliminary draft — 5 November 2026 — not for submission  |  page ")
    fr.font.size = Pt(8)
    fr.font.color.rgb = RGBColor(0x77, 0x77, 0x77)
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    footer._p.append(fld)

    for kind, text in blocks:
        if kind == "pagebreak":
            doc.add_page_break()
        elif kind == "image":
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_after = Pt(2)
            p.add_run().add_picture(str(png_dir / text), width=Cm(15.4))
        elif kind == "todo_notice":
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            for r in p.runs:
                pass
            add_runs(p, text, False, True)
            for r in p.runs:
                r.font.size = Pt(9)
                r.font.color.rgb = RGBColor(0xB3, 0x50, 0x00)
        elif kind in ("h1", "h2", "h3"):
            p = doc.add_paragraph(style="Heading %s" % kind[1])
            add_runs(p, text, False, False)
        elif kind == "figlabel":
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            r = p.add_run(text)
            r.bold = True
        elif kind == "note":
            p = doc.add_paragraph()
            add_runs(p, text, False, True)
            for r in p.runs:
                r.font.size = Pt(9)
                r.italic = True
        elif kind == "ref":
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.8)
            p.paragraph_format.first_line_indent = Cm(-0.8)
            p.paragraph_format.space_after = Pt(3)
            add_runs(p, text, False, False)
            for r in p.runs:
                r.font.size = Pt(9.5)
        elif kind == "legend":
            p = doc.add_paragraph()
            p.paragraph_format.space_after = Pt(9)
            add_runs(p, text, True, False)
        else:
            p = doc.add_paragraph()
            add_runs(p, text, True, True)
    doc.save(str(path))


# --------------------------------------------------------------------------
# PDF (reading copy) and figure plates
# --------------------------------------------------------------------------
def rectfn(rect_num, filled):
    page = pymupdf.paper_rect("a4")
    return page, page + (50, 52, -50, -52), pymupdf.Matrix(1, 1)


def add_page_numbers(path):
    doc = pymupdf.open(str(path))
    for i, page in enumerate(doc, 1):
        page.insert_text((page.rect.width / 2 - 8, page.rect.height - 26), str(i),
                         fontname="helv", fontsize=9, color=(0.5, 0.5, 0.5))
    doc.saveIncr()
    n = len(doc)
    doc.close()
    return n


def write_story(html_text, css, path, archive=None):
    story = pymupdf.Story(html=html_text, user_css=css, archive=archive)
    writer = pymupdf.DocumentWriter(str(path))
    story.write(writer, rectfn)
    writer.close()
    return add_page_numbers(path)


def build_pdf(blocks, png_dir, path):
    return write_story(blocks_to_html(blocks), CSS, path,
                       pymupdf.Archive(str(png_dir)))


SUPP_CSS = CSS + """
.bullet { padding-left: 14pt; text-indent: -10pt; }
h1 { font-size: 14pt; text-align: left; margin: 6pt 0 12pt 0; }
"""


def supp_html():
    out = []
    sources = ["Supplementary_Materials_v2.md", "Supplementary_Table_Legends.md",
               "Supplementary_Figure_Legends.md"]
    for idx, name in enumerate(sources):
        if idx:
            out.append('<div class="pagebreak"></div>')
        for line in read(SUPP_DIR / name).splitlines():
            s = line.strip().replace("\\*", "*")
            if not s or s == "---":
                continue
            if s.startswith("## "):
                out.append("<h2>%s</h2>" % html.escape(s[3:], quote=False))
            elif s.startswith("# "):
                out.append("<h1>%s</h1>" % html.escape(s[2:], quote=False))
            elif s.startswith("> "):
                out.append('<p class="note">%s</p>' % inline_html(s[2:], False, False))
            elif s.startswith("- ") or re.match(r"^\d+\. ", s):
                body = s[2:] if s.startswith("- ") else s
                out.append('<p class="bullet">%s</p>' % inline_html(body, False, False))
            else:
                out.append("<p>%s</p>" % inline_html(s, False, False))
    return "\n".join(out)


def build_supp_pdf(path):
    n = write_story(supp_html(), SUPP_CSS, path)
    s3 = SUPP_DIR / "Supplementary Figure S3.pdf"
    if s3.exists():
        doc = pymupdf.open(str(path))
        src = pymupdf.open(str(s3))
        doc.insert_pdf(src)
        # The S3 plates are full-bleed figure pages. A centred page number
        # would land on the panel E axis labels, so leave these two inserted
        # plates unnumbered; the text pages keep theirs.
        tmp = str(path) + ".tmp"
        doc.save(tmp, garbage=4, deflate=True)
        doc.close(); src.close()
        Path(tmp).replace(path)
        chk = pymupdf.open(str(path)); n = len(chk); chk.close()
    return n


def build_plates(path):
    out = pymupdf.open()
    toc = []
    for n, name in FIGMAP.items():
        src = pymupdf.open(str(FIG_PDF / name))
        toc.append([1, "Figure %d" % n, len(out) + 1])
        out.insert_pdf(src)
        src.close()
    out.set_toc(toc)
    out.save(str(path), garbage=4, deflate=True)
    pages = len(out)
    out.close()
    return pages


def render_pngs(png_dir, zoom=2.2):
    """Rasterise figure plates to quality-88 JPEG (reading copy / docx embedding)."""
    from PIL import Image
    for n, name in FIGMAP.items():
        doc = pymupdf.open(str(FIG_PDF / name))
        pix = doc[0].get_pixmap(matrix=pymupdf.Matrix(zoom, zoom), alpha=False)
        im = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
        im.save(png_dir / ("Figure%d.jpg" % n), format="JPEG", quality=88,
                optimize=True, dpi=(200, 200))
        doc.close()


# --------------------------------------------------------------------------
# Package assembly / manifest
# --------------------------------------------------------------------------
def collect_files():
    (OUT / "Figures").mkdir(parents=True, exist_ok=True)
    (OUT / "Supplementary").mkdir(parents=True, exist_ok=True)
    for n, name in FIGMAP.items():
        shutil.copy2(FIG_TIF / ("Figure%d.tif" % n), OUT / "Figures" / ("Figure%d.tif" % n))
        shutil.copy2(FIG_PDF / name, OUT / "Figures" / ("Figure%d.pdf" % n))
    for f in SUPP_FILES:
        shutil.copy2(SUPP_DIR / f, OUT / "Supplementary" / f)


MANIFEST = """# 初稿合订包 · 使用说明

**生成日期**：2026-11-05　|　**状态**：初步合订稿，供作者通读与传阅，**不可直接投稿**
**复现脚本**：`05_数据与结果/work/build_manuscript_package.py`（源文件改动后重跑即可）

## 一、这个文件夹是什么

把散在 `01_手稿/`、`02_主图/`、`03_补充材料/`、`07_投稿材料/` 的 v2 稿件材料，
按「正文 + 图 + 补充」三块收进一个目录。主文档已把 Title page、摘要、正文、声明、
全部编号参考文献、7 条图注和 7 张主图合在一起，**打开一个文件就能看完整篇稿件**。

```
10_初稿合订_v2/
├── Manuscript_v2_draft.docx      可编辑母版（Word / WPS）
├── Manuscript_v2_draft.pdf       阅读版（内容与 docx 相同，便于传阅/打印）
├── Figure_Plates.pdf             7 张主图单册（带书签）
├── Figures/                     Figure1–7.tif（600 dpi 投稿用）+ .pdf（矢量）
├── Supplementary/
│   ├── Supplementary_Tables_S1-S17.xlsx   数据表合订（46 工作表）
│   ├── Supplementary_Materials_v2.pdf     补充材料阅读版（表目录 + 说明 + S1–S4 图注 + S4 图页）
│   ├── Supplementary Figure S3.pdf       重算 UMAP 补充图（2 页，A4）
│   └── Supplementary_*.md                 可编辑的说明文本
└── README_初稿清单.md            本文件
```

## 二、三块内容与源文件的对应

| 部分 | 本包文件 | 来源 | 状态 |
|---|---|---|---|
| 正文 | `Manuscript_v2_draft.docx` / `.pdf` | `01_手稿/_md/manuscript_v2_tissue_composition.md` + `07_投稿材料/Title_Page_and_Declarations.md` | 正文事实层已核查；Title page 作者信息待补 |
| 主图 | `Figures/Figure1–7.tif`（+ `.pdf`）、`Figure_Plates.pdf` | `02_主图/v2_组织成分版/` 与 `07_投稿材料/figures_tiff/` | 7 张，600 dpi RGB/LZW TIFF，Okabe-Ito 统一配色 |
| 补充 | `Supplementary/Supplementary_Tables_S1-S17.xlsx` | `03_补充材料/` | 46 个工作表 = README + 45 张表（S1–S17 及子面板，含 S5c 比例风险检验与 S17a–f 药敏全表） |
| 补充阅读版 | `Supplementary/Supplementary_Materials_v2.pdf` | 由三份补充 md 排版 + Figure S4 图页装配 | 表目录 + 数据口径说明 + S1–S4 图注 + S4 图页，不用打开 Excel 也能通读 |

补充目录另附 `Supplementary_Materials_v2.md`（合订文本）、`Supplementary_Table_Legends.md`（纯表说明）、
`Supplementary_Figure_Legends.md`（S1–S3 图注）与 `Supplementary Figure S3.pdf`（重算 UMAP 与构成/marker/基因表达图）。

## 三、合订时顺手修掉的三处自相矛盾

上一轮把 Results 3.6 的「28 类细胞」改为原始数据真实的 9 个 major / 29 个 minor 注释后，
摘要与图注还留着旧口径。本次合订前一并改正，避免初稿内部打架：

1. 摘要「PLG was undetected in all 28 cell types」→「essentially undetected across all nine major cell types (≤0.2% of cells in any class)」。
2. Figure 6 图注 (A)「PLG was not detected in any compartment」→「essentially undetected (≤0.2% of cells in any compartment)」。
3. Figure 6 图注 (D) 曾改为「≤0.3% of cells」；第三轮核实后已改为 ER+ 与全细胞两套精确检出率，见下文第七节。

正文 `01_手稿/_md/manuscript_v2_tissue_composition.md` 已同步修改。

另外把补充合订文本里一句过期的「摘要 277 词」改成「约 303 词、按空格分词」，
并同步修了生成脚本 `build_supplementary_tables.py`，下次重跑不会再带回旧值。

## 四、需要你补齐 / 决定的事项

**Title page 仍留着 `[SQUARE BRACKETS]` 占位符**（docx 与 pdf 里以橙色标出），共 8 项：
作者名单与顺序、单位全称与编号、通讯作者信息（含 ORCID）、基金、利益冲突、preprint/会议摘要、
作者贡献（CRediT）、致谢。另需决定：代码是否公开托管取 DOI、三张旧补充图是否保留。

**目标期刊未定**，因此本包刻意没有做：期刊模板排版、摘要压缩（现约 370 词）、参考文献样式转换、
Cover letter、TRIPOD 版本确认。这些按你之前说的「等定具体期刊再说」。

## 五、投稿前已知待办（不影响通读）

1. **Data availability 有两套措辞**：正文里是「向通讯作者合理索取」，`Title_Page_and_Declarations.md`
   的声明草稿里有一版更详细、列出各数据库 URL 与补充表的写法。本包为避免重复，声明段只保留一句交叉引用；
   定刊格式时二选一合并。若改走 Zenodo/GitHub 公开托管，**两处都要同步改**。
2. **三张旧补充图 S1–S3 未被 v2 正文引用**（上一版遗留），本包未收入。定刊时决定保留还是删除。
3. **两份旧文件不可投稿**：`05_数据与结果/Supplementary Materials.docx`（旧 S1–S4 目录）、
   `05_数据与结果/work/Supplementary_Tables.xlsx`（旧 4 sheet）。本包未收录，原文件保留未动。
4. 补充图注文件里的 Figure S1–S3 章节在合订文本中已标为 `archived, not cited in v2`。
5. 词数按空格分词（2026-11-06 第七轮 C 档修订后重算）：摘要约 385 词、正文约 5,389 词；不同计数器结果有别，以目标期刊计数器为准。

## 六、校验记录

生成脚本自动检查：题录编号连续且全部被正文引用、7 条图注与 7 张图一一对应、
7 份 TIFF / 7 份 PDF 齐全、补充工作簿 46 个工作表、旧硬伤表述零残留。

## 七、2026-11-05 第三轮逐行审稿意见核实后的修订

18 条审稿意见（去重后 17 条）逐条回溯原始代码与结果表后完成以下修订：

1. 统一为「6 个已发表 + 本文 = 7 个签名 / 41 loci」；19.5%（8/41）是全体口径，21.6%（8/37）是仅已发表口径。摘要、Methods 2.9、Results 3.7、Discussion、Conclusion、Figure 7 图注与 Table S13a 标题同步。
2. Figure 1D 补齐 broadly expressed (3)，正文列出 12+20+29+15+14+14+3+25=132。
3. 增量价值拆成三个模型重述（Table S10a 的 0.725→0.727 / 0.657 / 0.765 与 Table S10b 的 0.734 / 0.657 / 0.765），删去不实的 "identical"；Results 3.4 与 Figure 4E 图注写明 n=327 / 1,445 / 2,840，并说明 Table S9 的 1,498 / 2,935 属于另一套完整病例口径。
4. Methods 2.1 改为「使用数据源已归一化矩阵」，不虚构 limma/affy 预处理步骤；补 GeneCards 引用 [40]，现共 40 条题录。
5. Methods 2.5 补全 11 个模块的枚举并说明两套模块体系；Methods 2.9 列出 9 个泛化模块；Figure 7D 图注给出映射，明确 NETosis 是 6 基因效应模块而非 22 基因参考评分。
6. Methods 2.7 在正文列出 22 个基因并改引 Table S14b；2.3 改写随机森林种子表述（单次 `set.seed(2026)`、跨包版本排序差异）；2.11 与 Code availability 补 10 个 Python 3.13 脚本的用途与版本。
7. 摘要 / Results 3.6 / Figure 6D 图注统一 all-cell 与 ER+ 检出率口径；Figure 6C 的 PLG 条置灰并标 N/A（Figure 6、7 已重画，7 张 TIFF 已重生成）；Figure 4D 图注解释 S0–S3 与 0.690 / 0.654 / 0.666 的关系；Table S13a 新增 Reference 列对应 [14]–[19]；Table S10b 补回 METABRIC CombinedScore 行并注明来源。
8. Discussion 改为「四个基因均属组织组成指标」，并说明 PLG 反映血液含量与血管灌注而非肿瘤驻留细胞。

**仍未解决的阻塞项**：P1-1 的 `[REPOSITORY URL]` / `[VERSION]` / `[DOI]` 需上传代码包取 DOI 后回填；8 项作者信息占位符；目标期刊与摘要压缩。

## 八、2026-11-05 第四轮复审（P3-19 / P3-20）

1. P3-19：Methods 2.1 原先只给完整协变量子集的样本量，与 Results 3.4 / Table S9 的表达+生存队列并列时产生歧义。现于 2.1 同时列出两套完整病例定义及各自事件数（GSE20685 327/83；METABRIC 1,498/879 与 1,445/848；GSE96058 2,935/264 与 2,840/247），并注明摘要等处默认指协变量完整子集。Results 3.4、Figure 4E 图注补上较大队列的事件数；Table S9 表注与补充材料 Notes 同步。**未采用审稿人建议的写法**，因为它把 848 事件挂在 1,498 例队列上，而 848 属于 1,445 例子集（1,498 例为 879 事件），GSE96058 同理。
2. P3-20：Results 3.6 补一句纯算术说明——ELANE 的 ER+ 检出率（0.13%）略高于全细胞率（0.11%），是因为 ER− 细胞检出率更低（约 0.10%），0.02 个百分点处于检出下限内。**未采用「ER+ 中性粒相关检出更高」的生物学措辞**，避免对十余个阳性细胞的差异做过度解读。

## 九、2026-11-05 综合评价核实后的四处文字补充

对第三方综合评价逐条回溯后，只采纳不涉及数据与图表的四处补充：

1. **跨亚型局限**：泛化分析在 1,111 例未筛选的 TCGA-BRCA 全亚型中评分，而签名在 ER+ 中构建，正文此前未说明。Limitations 补一句，注明相关性模式是类别层面比较而非 ER+ 特异效应估计。评价误称「作者已说明这一点」，实为真实缺口。
2. **偏相关标注与一处既有错配**：Results 3.6 的 −0.377 / −0.268 / −0.269 现标明是 Bcell_noCD79A 模块在 TCGA / GSE20685 / METABRIC 三队列的 `rho_adjPurity`（Table S12d）。同时修正正文原写「控制上皮纯度**和增殖**」与数值不符的问题——这三个值来自只控制纯度的列，已删去「和增殖」，与 Methods 2.5 一致；Table S12d 表注补列名说明。
3. **模板遵循声明**：Methods 2.3 补一句，说明流程原样复现已发表模板、未做刻意弱化选择，回应「作者自建较弱签名」的可能质疑。
4. **3.4 拆段**：在 "Two further observations point to measurement rather than modelling" 处断开，把三个拟合的区分与跨平台测量证据分成两段。

**未采纳**：删 Methods 2.5 的 11 模块枚举（Table S12a 已有 96 行完整列表，且与第三轮 P2-12 相冲突）、压缩药敏/治疗交互段（会切断 S17 / S15 的正文锚点）、为「正向对照」补造实验数据、补充材料新增核查清单（会改变 45 sheet 结构）。本轮未改任何数值、图或数据表。
## 十、单细胞措辞与注释体系说明（A 档，2026-11-06）

核查发现：本文单细胞部分复用 GSE176078 作者发布的 curated 标签做基因来源归因，全项目无 QC / 归一化 / HVG / 聚类 / UMAP 代码；GEO 处理文件的 `metadata.csv` 也不含 UMAP 坐标（原图谱 UMAP 见 Wu et al. 2021 Fig.1）。更关键的是，该注释体系的 9 个 major 与 29 个 minor 标签中**没有任何中性粒细胞 / 粒细胞类别**，因此「九大类中找不到 NET 效应基因」部分源于注释粒度，正文此前只写了「10x 低估中性粒细胞」。

本轮（A 档，仅文字）已执行：

1. 标题 `with single-cell profiling` → `with single-cell attribution`，Title page 与代码包 README 同步；章节 2.6 改名为 `Single-cell attribution using published annotations`。
2. Methods 2.6 明写复用作者标签、未做 QC / 重归一化 / HVG / 聚类 / 降维、未生成 UMAP 或 tSNE，并指向文献 [32] 的原图谱；同时列出髓系区室的实际构成，说明不存在中性粒细胞类（除 10x 低估外的第二个原因）。
3. Results 3.6、Figure 6 图注 (A)、Discussion Limitations 同步标注该注释体系边界，Limitations 明写未做 de novo 聚类或 UMAP。
4. 补充材料 Notes 新增 `Single-cell annotation audit`：9 大类 / 29 minor / 亚型细胞数、谱系 marker 验证值；Table S11a 表注补口径说明。
5. `j3_singlecell.R` 追加可复现的审计段，输出 `j3_annotation_audit.csv` 与 `j3_marker_verification.csv`（实测 9 / 29 / 3 类、中性粒细胞标签 0 个），公开代码包副本同步。

**A 档第 4 项（新增细胞类型构成 + marker 验证图）本轮未做**：该图与下一步「是否重跑 UMAP」共用同一张补充图，先画会返工；且 `03_补充材料/` 的 S1–S3 编号已被三张旧补图占用、去留未定，新增补充图需连带改合订包的图片装配逻辑。待决定是否做 C 档 UMAP 后一并处理。本轮**未改动任何数值、图或 7 图结构**，补充工作簿仍 45 sheet。

## 十一、C 档 UMAP 重算、亚型口径修正与补充图 S4（2026-11-06）

第十节记录的遗留项在本轮一次性闭环。核查发现 GSE176078 并非纯 ER+ 数据集（GEO 官方 11 ER+ / 5 HER2+ / 10 TNBC，26 例 100,064 细胞），且 GEO 样本级特征另有 2 例 `HER2+/ER+`（CID3586、CID4066）与 1 例标签冲突（CID3963：GEO 标 ER+、处理文件标 TNBC）。四种亚型口径的敏感性重算显示四个 hub 基因的主导细胞类型不变，结论稳健，但正文必须写清分母。

1. **正文口径**：Methods 2.6 明写两套细胞集合（hub 归因 / marker / Figure 6D 用 11 例 ER+ 38,241 细胞；Figure 1D 候选池区室分类用全部 26 例），并说明 ER+ 子集不含 HER2 共阳性、存在单样本标签冲突。摘要、Results 3.6、Discussion Limitations、Figure 1 / Figure 6 图注同步。
2. **UMAP 措辞反转**：A 档写的「未生成 UMAP」已改为「为展示重算了 UMAP，不重新聚类、不重新注释、不进入任何统计模型，坐标非原作者坐标」，指向新增的 Supplementary Figure S3（原编号 S4）。
3. **补充图 S3**：`03_补充材料/Supplementary Figure S3.pdf`（A4 纵向 2 页，UMAP 按 9 大类/亚型/ER+ 高亮 + 细胞构成 + marker 验证 + 12 个基因表达，4 x 3 排列，覆盖特征表全部基因）；原 S3（nomogram + DCA）因正文无对应方法与结果已移除，本图由原 S4 顺延编号；`Supplementary_Figure_Legends.md` 同步图注；合订脚本的 `supp_html()` 纳入第三份 md，`build_supp_pdf()` 会把该图页追加进补充阅读版。
4. **可复用产物**：`12_可复用资源/sc_umap/` 提供坐标表（csv.gz）、12 基因表达表、轻量 R 对象、JSON manifest（含 SHA-256 cache_key）、sessionInfo、config/pipeline/plots/load 脚本与 REUSE.md；公开代码包 `scripts/07_umap/` 同步脚本副本，`00_config.R` 增加 `WB_SEED_UMAP`。本机无 pyarrow/pandas，故主格式定为 csv.gz + Python 标准库加载器，未生成 Parquet（已在 REUSE.md 如实说明）。
5. **计数同步**：公开代码包由 78 R + 10 Python 变为 82 R + 11 Python；正文 2.11 与 Title page、代码包 README/MANIFEST 同步。词数重算为摘要 385 / 正文 5,389。

本轮**未改动任何主图、数值或补充工作簿的 45 个 sheet**；新增内容仅为补充图 S4、正文口径文字与可复用产物。（此句描述第十一节那一轮；紧接的第十二节记录了后续改动。）

## 十二、单细胞归一化口径统一（2026-11-06 追加）

复核 C 档 UMAP 时发现 `j3_singlecell.R` 用 20 个目标基因的 counts 之和作 CP10K 分母，与 UMAP 流程的全转录组 `LogNormalize` 不同口径；已改用 `metadata$nCount_RNA` 从 26 份 mtx 完整重跑，新表与 UMAP 表达表逐细胞校验一致（最大差 3.6×10⁻¹⁵）。

1. **主图 Figure 6 重画**：A 面板色标 0–5.4→0–1.7、B 面板对数轴 4e-4–25→6e-5–3、倍数标签更新为 19×/510×/53×/462×；`Figure6_ABCDEF.pdf` 与 600 dpi `Figure6.tif` 已同步到 `02_主图/`、`07_投稿材料/`、`10_初稿合订_v2/Figures/`。
2. **补充图 S4 重出**：E 面板改为 6 个谱系 marker（删去 MKI67），副标题、图注与补充材料 Notes 同步；PDF 三处副本哈希一致，压缩后 19.60 MB。
3. **数值同步**：Results 3.6、Methods 2.6、Figure S4 图注、Table S11a 说明与 45 sheet 工作簿的 `Table S11a` 均值列全部更新；检出率、Figure 6D、Table S11b、S4 第 2 页百分比及全部 bulk / 生存分析数值未变。
4. **排版修正**：补充阅读版中 S4 图页的页码原压在 E 面板 x 轴刻度上，现改为图页不编号（文字页 1–12 仍保留页码）。
5. 修正前的 4 张 j3 csv 与 `j3_allcells.rds` 备份于 `99_归档/归一化修正前_20261106/`；完整前后对照见 `改稿说明_v2.md` 第十七节。

"""


def validate(blocks):
    ms = read(MS)
    refs = re.findall(r"^\[(\d+)\]", ms, re.M)
    body = ms.split("## References")[0]
    cited = set()
    for m in re.finditer(r"\[(\d+)(?:[–-](\d+))?\]", body):
        a = int(m.group(1))
        b = int(m.group(2) or a)
        cited.update(range(a, b + 1))
    problems = []
    n_refs = len(refs)
    expected = set(range(1, n_refs + 1))
    if refs != [str(i) for i in range(1, n_refs + 1)]:
        problems.append("参考文献编号不连续: %s" % refs[:5])
    if cited != expected:
        problems.append("悬空/未引用编号: %s" % sorted(expected - cited))
    legends = [t for k, t in blocks if k == "legend"]
    if len(legends) != 7:
        problems.append("图注数 %d != 7" % len(legends))
    for n in range(1, 8):
        for p in (FIG_TIF / ("Figure%d.tif" % n), OUT / "Figures" / ("Figure%d.tif" % n),
                  OUT / "Figures" / ("Figure%d.pdf" % n)):
            if not p.exists() or p.stat().st_size < 1000:
                problems.append("缺文件 %s" % p.name)
    for f in SUPP_FILES:
        if not (OUT / "Supplementary" / f).exists():
            problems.append("缺补充文件 %s" % f)
    if not (OUT / "Supplementary" / "Supplementary_Materials_v2.pdf").exists():
        problems.append("缺补充材料阅读版 PDF")
    for bad in ("674188", "up to 30%", "28 cell types", "28 annotated",
                "PLG was not detected in any compartment",
                "R 4.4.3", "ridge regression (pRRophetic[33])",
                "reported for completeness only",
                "available from the corresponding author on reasonable request",
                "Only ER+ tumours were analysed", "no UMAP/tSNE embedding was generated",
                "no de novo clustering or UMAP embedding was generated"):
        if bad in ms:
            problems.append("旧硬伤残留: %s" % bad)
    import openpyxl
    wb = openpyxl.load_workbook(OUT / "Supplementary" / SUPP_FILES[0], read_only=True)
    n_sheets = len(wb.sheetnames)
    wb.close()
    EXPECTED_SHEETS = 46  # README + 45 data sheets (S1-S17, incl. S5c and S17a-S17f)
    if n_sheets != EXPECTED_SHEETS:
        problems.append("补充工作簿 sheet 数 %d != %d" % (n_sheets, EXPECTED_SHEETS))
    return problems, len(refs), len(legends), n_sheets


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    blocks = build_blocks()
    collect_files()

    tmp = Path(tempfile.mkdtemp(prefix="wb_manuscript_"))
    try:
        render_pngs(tmp)
        build_docx(blocks, tmp, OUT / "Manuscript_v2_draft.docx")
        n_pages = build_pdf(blocks, tmp, OUT / "Manuscript_v2_draft.pdf")
        n_plates = build_plates(OUT / "Figure_Plates.pdf")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    n_supp_pages = build_supp_pdf(OUT / "Supplementary" / "Supplementary_Materials_v2.pdf")
    (OUT / "README_初稿清单.md").write_text(MANIFEST, encoding="utf-8")
    problems, n_refs, n_legends, n_sheets = validate(blocks)

    print("输出目录:", OUT)
    for p in sorted(OUT.iterdir()):
        if p.is_file():
            print("  %-42s %8.1f KB" % (p.name, p.stat().st_size / 1024))
    for sub in ("Figures", "Supplementary"):
        files = sorted((OUT / sub).iterdir())
        print("  %-14s %d 个文件, %.2f MB" %
              (sub + "/", len(files), sum(f.stat().st_size for f in files) / 1048576))
    print("主文档 PDF %d 页, 图版 %d 页, 补充阅读版 %d 页 | 题录 %d 条 / 图注 %d 条 / 补充表 %d 个 sheet" %
          (n_pages, n_plates, n_supp_pages, n_refs, n_legends, n_sheets))
    print("校验:", "全部通过" if not problems else "发现问题 -> " + "; ".join(problems))
    if problems:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
