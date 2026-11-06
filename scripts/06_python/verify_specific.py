# -*- coding: utf-8 -*-
"""Spot checks by exact DOI / PMID for the references that fuzzy matching missed."""
import json
import urllib.parse
import urllib.request

UA = "WorkBuddyReferenceAudit/1.0 (mailto:ref.audit@example.org)"


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    return urllib.request.urlopen(req, timeout=25).read().decode("utf-8", "replace")


def cr(doi):
    try:
        d = json.loads(get("https://api.crossref.org/works/" + doi))["message"]
        return {"title": (d.get("title") or [""])[0],
                "journal": (d.get("container-title") or [""])[0],
                "year": (d.get("published", {}).get("date-parts") or [[None]])[0][0],
                "volume": d.get("volume"), "issue": d.get("issue"),
                "page": d.get("page"), "article-number": d.get("article-number"),
                "doi": d.get("DOI"),
                "authors": ["%s %s" % (a.get("given", ""), a.get("family", ""))
                            for a in (d.get("author") or [])][:5]}
    except Exception as e:
        return {"error": str(e)}


def epmc(query):
    u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
        {"query": query, "format": "json", "resultType": "core", "pageSize": 2})
    try:
        r = json.loads(get(u))["resultList"]["result"]
        return [{"title": x.get("title"), "authors": x.get("authorString"),
                 "journal": x.get("journalInfo", {}).get("journal", {}).get("title"),
                 "year": x.get("pubYear"), "volume": x.get("journalInfo", {}).get("volume"),
                 "issue": x.get("journalInfo", {}).get("issue"),
                 "pages": x.get("pageInfo"), "doi": x.get("doi"),
                 "pmid": x.get("pmid")} for x in r]
    except Exception as e:
        return [{"error": str(e)}]


CASES = [
    ("[5] 稿件写的编号 674188", "10.3389/fimmu.2021.674188", 'DOI:"10.3389/fimmu.2021.674188"'),
    ("[5] 权威记录编号 674192", "10.3389/fimmu.2021.674192", 'DOI:"10.3389/fimmu.2021.674192"'),
    ("[37] Ann Intern Med 版 TRIPOD", "10.7326/M14-0697", 'EXT_ID:25560714 AND SRC:MED'),
    ("[37] BJ C 同步版", "10.1038/bjc.2014.639", 'DOI:"10.1038/bjc.2014.639"'),
    ("[2] Harbeck Nat Rev Dis Primers", "10.1038/s41572-019-0011-6", 'DOI:"10.1038/s41572-019-0011-6"'),
    ("[23] DESeq2 Genome Biol 版", "10.1186/s13059-014-0550-8", 'DOI:"10.1186/s13059-014-0550-8"'),
    ("[31] Cell Rep 版", "10.1016/j.celrep.2016.12.019", 'DOI:"10.1016/j.celrep.2016.12.019"'),
]

for label, doi, q in CASES:
    print("=" * 95)
    print(label, "| DOI", doi)
    print("  Crossref :", json.dumps(cr(doi), ensure_ascii=False)[:600])
    print("  EuropePMC:", json.dumps(epmc(q), ensure_ascii=False)[:700])
