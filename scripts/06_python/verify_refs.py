# -*- coding: utf-8 -*-
"""
Reviewer-style reference audit for manuscript v2.

Step 1: parse the 37 numbered references and every in-text citation.
Step 2: query Crossref (by title) and Europe PMC (core result: metadata + abstract).
Output: ref_audit_data.json, ref_audit_citations.txt, ref_audit_summary.txt

Raw evidence only; the verdicts are written by us in the audit report.
"""
import difflib
import html
import json
import os
import re
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MS = os.path.join(ROOT, "01_手稿", "_md", "manuscript_v2_tissue_composition.md")
OUT_JSON = os.path.join(HERE, "ref_audit_data.json")
OUT_CITE = os.path.join(HERE, "ref_audit_citations.txt")
OUT_SUM = os.path.join(HERE, "ref_audit_summary.txt")

UA = "WorkBuddyReferenceAudit/1.0 (mailto:ref.audit@example.org)"


def get(url, tries=2, timeout=25):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            if i == tries - 1:
                return "ERR:" + str(e)
            time.sleep(1.0)
    return "ERR"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def sim(a, b):
    return difflib.SequenceMatcher(None, norm(a), norm(b)).ratio()


def extract_title(line):
    """'Authors. Title. Journal. YEAR;vol(issue):pages.' -> 'Title'."""
    ys = list(re.finditer(r"\.\s*\d{4};", line))
    if not ys:
        return line
    prefix = line[:ys[-1].start()].strip().rstrip(".")
    seg = re.split(r"\.\s+", prefix)
    if len(seg) >= 3:
        return seg[-2].strip()
    if len(seg) == 2:
        return seg[-1].strip()
    return prefix


def parse_manuscript():
    s = open(MS, encoding="utf-8").read()
    refs = {}
    for line in s.splitlines():
        m = re.match(r"^\[(\d+)\]\s+(.+)$", line.strip())
        if m:
            refs[int(m.group(1))] = m.group(2).strip()
    body = s.split("## References")[0]
    ctx = {}
    for m in re.finditer(r"\[(\d+(?:\s*[,–-]\s*\d+)*)\]", body):
        inner = m.group(1)
        nums = []
        for part in re.split(r"[,]", inner):
            pm = re.match(r"\s*(\d+)\s*[–-]\s*(\d+)\s*$", part)
            if pm:
                nums.extend(range(int(pm.group(1)), int(pm.group(2)) + 1))
            elif part.strip().isdigit():
                nums.append(int(part.strip()))
        frag = re.sub(r"\s+", " ", body[max(0, m.start() - 200):
                                       min(len(body), m.end() + 150)])
        for n in nums:
            ctx.setdefault(n, []).append(frag)
    return refs, ctx


def crossref(title):
    q = urllib.parse.urlencode({
        "query.title": title, "rows": 4,
        "select": "DOI,title,author,container-title,published,volume,issue,page,type,score,abstract",
        "mailto": "ref.audit@example.org",
    })
    raw = get("https://api.crossref.org/works?" + q)
    if raw.startswith("ERR:"):
        return {"error": raw}
    try:
        items = json.loads(raw)["message"]["items"]
    except Exception as e:
        return {"error": "parse " + str(e)}
    cands = []
    for it in items:
        cands.append({
            "doi": it.get("DOI"),
            "title": (it.get("title") or [""])[0],
            "authors": ["%s %s" % (a.get("given", ""), a.get("family", ""))
                        for a in (it.get("author") or [])],
            "journal": (it.get("container-title") or [""])[0],
            "year": (it.get("published", {}).get("date-parts") or [[None]])[0][0],
            "volume": it.get("volume"), "issue": it.get("issue"),
            "page": it.get("page"), "type": it.get("type"),
            "abstract": re.sub(r"<[^>]+>", " ", it.get("abstract", ""))[:1500],
            "score": it.get("score"),
        })
    cands.sort(key=lambda c: sim(title, c["title"]), reverse=True)
    return {"candidates": cands}


def epmc(title):
    q = urllib.parse.urlencode({
        "query": 'TITLE:"%s" AND SRC:MED' % title.replace('"', ""),
        "format": "json", "resultType": "core", "pageSize": 2,
    })
    raw = get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + q)
    if raw.startswith("ERR:"):
        return {"error": raw}
    try:
        res = json.loads(raw)["resultList"]["result"]
    except Exception as e:
        return {"error": "parse " + str(e)}
    if not res:
        return {"pmid": None}
    d = max(res, key=lambda x: sim(title, x.get("title", "")))
    return {
        "pmid": d.get("pmid"), "pmcid": d.get("pmcid"),
        "title": d.get("title"), "authors": d.get("authorString"),
        "journal": d.get("journalTitle"), "journal_abbr": d.get("journalTitle"),
        "year": d.get("pubYear"), "volume": d.get("journalVolume"),
        "issue": d.get("issue"), "pages": d.get("pageInfo"),
        "doi": d.get("doi"), "pubtype": d.get("pubType"),
        "abstract": re.sub(r"<[^>]+>", " ", html.unescape(d.get("abstractText") or "")),
        "first_date": d.get("firstPublicationDate"),
        "cited_by": d.get("citedByCount"),
    }


def write_outputs(data):
    json.dump(data, open(OUT_JSON, "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    with open(OUT_CITE, "w", encoding="utf-8") as fh:
        for n in sorted(data):
            fh.write("=== [%s] 正文引用 %d 处 ===\n" % (n, len(data[n]["contexts"])))
            for c in data[n]["contexts"]:
                fh.write("  " + c + "\n")
    write_summary(data)
    print("saved", OUT_JSON, OUT_CITE, OUT_SUM)


def write_summary(data):
    with open(OUT_SUM, "w", encoding="utf-8") as fh:
        for n in sorted(data, key=int):
            d = data[n]
            b, e = d["best"], d["epmc"]
            fh.write("=" * 100 + "\n[%d] 稿件: %s\n" % (int(n), d["manuscript"]))
            fh.write("     提取标题: %s\n" % d["parsed_title"])
            fh.write("     Crossref: %s | %s | %s %s;%s(%s):%s | doi %s | 相似度 %.2f\n" %
                     (b.get("title", ""), ", ".join(b.get("authors", [])[:3]),
                      b.get("journal", ""), b.get("year", ""), b.get("volume", ""),
                      b.get("issue", ""), b.get("page", ""), b.get("doi", ""),
                      sim(d["parsed_title"], b.get("title", ""))))
            fh.write("     EuropePMC: %s | %s | %s %s;%s(%s):%s | PMID %s | doi %s\n" %
                     (e.get("title", ""), (e.get("authors") or "")[:110],
                      e.get("journal", ""), e.get("year", ""), e.get("volume", ""),
                      e.get("issue", ""), e.get("pages", ""), e.get("pmid"),
                      e.get("doi", "")))
            fh.write("     摘要: %s\n" %
                     (e.get("abstract") or b.get("abstract") or "（无）")[:1400])
            fh.write("     引用处 %d, 相似度 %.2f\n" %
                     (len(d["contexts"]), sim(d["parsed_title"], b.get("title", ""))))


def main():
    if "--summary-only" in __import__("sys").argv:
        data = json.load(open(OUT_JSON, encoding="utf-8"))
        write_outputs(data)
        return
    refs, ctx = parse_manuscript()
    print("题录 %d 条，正文引用编号 %d 个" % (len(refs), len(ctx)), flush=True)
    data = {}
    for n in sorted(refs):
        line = refs[n]
        title = extract_title(line)
        cr = crossref(title)
        cands = cr.get("candidates", [])
        best = cands[0] if cands else {}
        ep = epmc(title)
        data[n] = {"manuscript": line, "parsed_title": title,
                   "crossref": cr, "best": best, "epmc": ep,
                   "contexts": ctx.get(n, [])}
        print("[%2d] title=%-48s CR=%-34s EPMC=%-9s abst=%4d" %
              (n, title[:48], (best.get("journal") or "NONE")[:34],
               ep.get("pmid") or "NONE", len(ep.get("abstract") or "")), flush=True)
        time.sleep(0.15)
    write_outputs(data)


if __name__ == "__main__":
    main()
