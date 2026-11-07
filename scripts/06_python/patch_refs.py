# -*- coding: utf-8 -*-
"""
参考文献按正文首次出现顺序重编（旧号 -> 新号）

映射由 _ref_map.json 提供（脚本生成，见同目录 _ref_map.json）。
变动项：
  [23]->[25] [24]->[26] [25]->[27] [26]->[28] [27]->[29] [28]->[30]
  [29]->[31] [30]->[32] [31]->[33] [32]->[24] [33]->[37]
  [37]->[39] [38]->[40] [39]->[38] [40]->[23]
  其余 25 条编号不变（含 [14]-[19] 六个签名文献与 [34]-[36] GDSC）。

同时把 References 段落按新编号重新排序。

用法: python patch_refs.py
"""
import io
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

MAP = {int(k): v for k, v in json.load(
    io.open(os.path.join(HERE, "_ref_map.json"), encoding="utf-8")).items()}

TARGETS = [
    "01_手稿/_md/manuscript_v2_tissue_composition.md",
    "07_投稿材料/Title_Page_and_Declarations.md",
]

RANGE = re.compile(r"\[(\d+)(\s*[–-]\s*)(\d+)\]")
SINGLE = re.compile(r"\[(\d+)\]")


def conv_range(m):
    a, b, sep = int(m.group(1)), int(m.group(3)), m.group(2)
    new = [MAP[n] for n in range(a, b + 1)]
    if all(new[i] + 1 == new[i + 1] for i in range(len(new) - 1)):
        return "[%d%s%d]" % (new[0], sep, new[-1])
    return "".join("[%d]" % n for n in new)


def convert(text):
    text = RANGE.sub(conv_range, text)
    return SINGLE.sub(lambda m: "[%d]" % MAP[int(m.group(1))], text)


def reorder_refs(text):
    """References 段按新编号排序。"""
    marker = "## References"
    if marker not in text:
        return text, 0
    i = text.index(marker) + len(marker)
    j = text.find("\n## ", i)
    if j == -1:
        j = len(text)
    body = text[i:j]

    refs = []
    for l in body.split("\n"):
        m = re.match(r"^\[(\d+)\] ", l.strip())
        if m:
            refs.append((int(m.group(1)), l.strip()))
    refs.sort(key=lambda x: x[0])
    nums = [n for n, _ in refs]
    assert nums == list(range(1, len(refs) + 1)), "编号不连续: %s" % nums

    rebuilt = "\n\n" + "\n".join(t for _, t in refs) + "\n\n"
    return text[:i] + rebuilt + text[j:].lstrip("\n"), len(refs)


def main():
    changed = {}
    for rel in TARGETS:
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            print("[skip]", rel)
            continue
        s = io.open(p, encoding="utf-8").read()
        before = s
        s = convert(s)
        nref = 0
        if rel.endswith(".md") and "## References" in s:
            s, nref = reorder_refs(s)
        if s != before:
            io.open(p, "w", encoding="utf-8", newline="\n").write(s)
            changed[rel] = nref
        print("[ok] %-60s refs=%s" % (rel, nref or "-"))
    print("\n共修改 %d 个文件" % len(changed))


if __name__ == "__main__":
    main()
