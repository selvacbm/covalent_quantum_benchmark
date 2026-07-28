#!/usr/bin/env python3
"""Fetch assay descriptions for retained benchmark records.

This is especially required for the KRAS_G12C label because its ChEMBL target
entry is generic KRas.  Assay text and the cited paper must explicitly support
the mutation before the record is described as KRAS(G12C).
"""
import argparse
import csv
import time
from pathlib import Path

import requests

BASE = "https://www.ebi.ac.uk/chembl/api/data"


def get_json(session, url, attempts=5):
    error = None
    for attempt in range(attempts):
        try:
            response = session.get(url, timeout=120)
            response.raise_for_status()
            return response.json()
        except Exception as exc:
            error = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"ChEMBL request failed after {attempts} attempts: {error}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(args.benchmark.open()))
    session = requests.Session()
    session.headers.update({"User-Agent": "CovalentScope-QM/0.1 (assay-review)"})
    cache, output = {}, []
    for index, row in enumerate(rows, 1):
        assay_id = row.get("assay_chembl_id", "")
        if assay_id not in cache:
            try:
                cache[assay_id] = get_json(session, f"{BASE}/assay/{assay_id}.json")
            except RuntimeError as exc:
                print(f"WARNING: {assay_id}: {exc}", flush=True)
                cache[assay_id] = {}
        assay = cache[assay_id]
        description = assay.get("description", "")
        output.append({
            "record_id": row.get("record_id", ""),
            "target_label": row.get("target_label", ""),
            "molecule_chembl_id": row.get("molecule_chembl_id", ""),
            "assay_chembl_id": assay_id,
            "document_chembl_id": row.get("document_chembl_id", ""),
            "assay_description": description,
            "assay_type": assay.get("assay_type", ""),
            "assay_organism": assay.get("assay_organism", ""),
            "assay_target_chembl_id": assay.get("target_chembl_id", ""),
            "variant_or_mechanism_confirmed": "",
            "evidence_location": "",
            "review_status": "pending",
        })
        if index % 50 == 0 or index == len(rows):
            print(f"assays={index}/{len(rows)}", flush=True)
    fields = list(output[0]) if output else ["record_id"]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)
    print(f"Wrote {len(output)} assay-review rows to {args.out}")


if __name__ == "__main__":
    main()
