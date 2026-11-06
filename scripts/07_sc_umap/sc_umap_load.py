"""sc_umap_load.py - load-only helper, standard library only (no pandas required).

    import sys; sys.path.insert(0, "<dir>")
    from sc_umap_load import load_sc_umap
    cells = load_sc_umap(subtype="ER+")
    print(len(cells), cells[0]["UMAP_1"])

Set the WB_SC_UMAP_DIR environment variable, or pass dir=, to point at this
folder from another project. Returns a list of dicts. If pandas is available,
pass as_dataframe=True to get a DataFrame instead.
"""

import csv
import gzip
import json
import os


def load_sc_umap(dir=None, subtype=None, as_dataframe=False, validate=True):
    dir = dir or os.environ.get("WB_SC_UMAP_DIR", "out/sc_umap")
    emb_path = os.path.join(dir, "sc_umap_embeddings.csv.gz")
    man_path = os.path.join(dir, "sc_umap_manifest.json")
    if not os.path.exists(emb_path):
        raise FileNotFoundError(
            "Embeddings not found: %s\nSet WB_SC_UMAP_DIR or pass dir= to the sc_umap folder." % emb_path)

    with gzip.open(emb_path, "rt", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in ("UMAP_1", "UMAP_2", "nCount_RNA", "nFeature_RNA", "percent_mito"):
            if k in r and r[k] != "":
                r[k] = float(r[k])

    manifest = None
    if os.path.exists(man_path):
        with open(man_path, "r", encoding="utf-8") as fh:
            manifest = json.load(fh)

    if validate:
        need = {"cell_id", "sample", "subtype", "major", "minor", "UMAP_1", "UMAP_2"}
        missing = need - set(rows[0].keys() if rows else [])
        if missing:
            raise ValueError("Malformed embeddings table, missing: %s" % ", ".join(sorted(missing)))
        key = (manifest or {}).get("cache_key", "n/a")[:12]
        print("[load_sc_umap] %d cells, %d samples, %d major classes; cache_key=%s"
              % (len(rows), len({r["sample"] for r in rows}),
                 len({r["major"] for r in rows}), key))

    if subtype is not None:
        wanted = {subtype} if isinstance(subtype, str) else set(subtype)
        rows = [r for r in rows if r["subtype"] in wanted]

    if as_dataframe:
        import pandas as pd
        return pd.DataFrame(rows), manifest
    return rows, manifest
