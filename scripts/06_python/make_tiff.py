"""
生成投稿用 TIFF（600 dpi, RGB, LZW）
输入：02_主图/v2_组织成分版/*.pdf  +  02_主图/figure1.tif, figure2.tif
输出：07_投稿材料/figures_tiff/FigureN.tif
"""
import os, shutil
import pymupdf
from PIL import Image

SRC = os.path.join("02_主图", "v2_组织成分版")
OUT = os.path.join("07_投稿材料", "figures_tiff")
os.makedirs(OUT, exist_ok=True)

DPI = 600
MAP = {
    "Figure1_ABCDEF.pdf": "Figure1.tif",
    "Figure2_ABCD.pdf": "Figure2.tif",
    "Figure3_ABCDEFG.pdf": "Figure3.tif",
    "Figure4_ABCDEF.pdf": "Figure4.tif",
    "Figure5_ABCD.pdf": "Figure5.tif",
    "Figure6_ABCDEF.pdf": "Figure6.tif",
    "Figure7_ABCD.pdf": "Figure7.tif",
}

for pdf, tif in MAP.items():
    p = os.path.join(SRC, pdf)
    doc = pymupdf.open(p)
    page = doc[0]
    zoom = DPI / 72.0
    pix = page.get_pixmap(matrix=pymupdf.Matrix(zoom, zoom), alpha=False)
    im = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    dst = os.path.join(OUT, tif)
    im.save(dst, format="TIFF", compression="tiff_lzw", dpi=(DPI, DPI))
    mb = os.path.getsize(dst) / 1024 / 1024
    print(f"{tif}: {pix.width}x{pix.height} px  "
          f"({pix.width/DPI*25.4:.0f}x{pix.height/DPI*25.4:.0f} mm)  {mb:.1f} MB")
    doc.close()


