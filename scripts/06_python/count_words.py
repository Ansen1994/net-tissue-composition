# -*- coding: utf-8 -*-
"""
稿件字数统计（唯一定义，避免不同文件各报一个数）

口径
----
摘要     : "## Abstract" 到 "**Keywords:**" 之间，按空白切分（space-delimited）。
           同时给出「不含四个小标题」与「含四个小标题」两个值。
正文     : "## 1. Introduction" 到 "## Data availability" 之前，
           去掉所有 markdown 标题行与 ** 加粗标记，按空白切分。
           含 Methods 内的加粗小节标题（2.1、2.2 …），不含图注、参考文献、表格。
图注     : "## Figure legends" 之后的全部段落。
参考文献 : "## References" 到 "## Figure legends" 之间，以 [n] 开头的行数。

用法: python count_words.py [稿件.md]
默认读取 ../01_手稿/_md/manuscript_v2_tissue_composition.md
"""
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "..", "..", "01_手稿", "_md",
                       "manuscript_v2_tissue_composition.md")


def plain(t):
    t = re.sub(r"^#+ .*$", "", t, flags=re.M)
    t = re.sub(r"\*\*", "", t)
    return t


def wc(t):
    return len(plain(t).split())


def section(lines, start, end):
    return "\n".join(lines[start:end])


def main(path):
    lines = io.open(path, encoding="utf-8").read().split("\n")

    def idx(prefix):
        for i, l in enumerate(lines):
            if l.startswith(prefix):
                return i
        return None

    a = idx("## Abstract")
    b = idx("## 1. Introduction")
    c = idx("## Data availability")
    d = idx("## References")
    e = idx("## Figure legends")

    # ---- 摘要 ----
    abs_txt = section(lines, a + 1, b)
    abs_txt = abs_txt.split("**Keywords:**")[0]
    heads = ["Background", "Methods", "Results", "Conclusions"]
    per = {}
    for h in heads:
        m = re.search(r"\*\*%s:\*\*(.*?)(?=\n\n|\Z)" % h, abs_txt, re.S)
        per[h] = wc(m.group(1)) if m else 0
    n_head_words = sum(per.values())
    n_with_heads = n_head_words + len(heads)   # 四个小标题各 1 词

    # ---- 正文 ----
    body = section(lines, b, c)

    # ---- 参考文献 ----
    ref_lines = [l for l in lines[d:e] if re.match(r"^\[\d+\] ", l.strip())]

    # ---- 图注 ----
    legends = section(lines, e, len(lines))

    print("文件:", os.path.relpath(path, HERE))
    print()
    print("摘要")
    for h in heads:
        print("  %-12s %4d" % (h, per[h]))
    print("  %-12s %4d   <- 不含小标题（投稿填报用）" % ("合计", n_head_words))
    print("  %-12s %4d   <- 含四个小标题" % ("合计", n_with_heads))
    print()
    print("正文（Introduction–Conclusion，不含标题行/图注/参考文献）  %d 词" % wc(body))
    print("参考文献  %d 条" % len(ref_lines))
    print("图注  %d 词" % wc(legends))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT)
