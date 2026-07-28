#!/usr/bin/env python3
"""Fetch bibliographic metadata for the manual electrophile-evidence review.

This script does not infer covalent mechanism.  It simply turns ChEMBL
document identifiers into a source-review sheet containing title, DOI and
PubMed ID, so that the experimental paper can be checked efficiently.
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
    raise RuntimeError(f"Request failed after {attempts} attempts: {error}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("queue", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    queue = list(csv.DictReader(args.queue.open()))
    session = requests.Session()
    session.headers.update({"User-Agent": "CovalentScope-QM/0.1 (source-review)"})
    cache, output = {}, []
    for index, row in enumerate(queue, 1):
        document_id = row.get("document_chembl_id", "")
        if document_id not in cache:
            try:
                source = get_json(session, f"{BASE}/document/{document_id}.json")
            except RuntimeError as exc:
                print(f"WARNING: {document_id}: {exc}", flush=True)
                source = {}
            cache[document_id] = source
        source = cache[document_id]
        output.append({
            **row,
            "title": source.get("title", ""),
            "doi": source.get("doi", ""),
            "pubmed_id": source.get("pubmed_id", ""),
            "journal": source.get("journal", ""),
            "year": source.get("year", ""),
            "covalent_mechanism_evidence": "",
            "evidence_location": "",
            "reviewer": "",
            "review_date": "",
            "review_status": "pending",
        })
        if index % 20 == 0 or index == len(queue):
            print(f"review_sources={index}/{len(queue)}", flush=True)
    fields = [
        "target_label", "document_chembl_id", "warhead_class", "n_records", "note",
        "title", "doi", "pubmed_id", "journal", "year",
        "covalent_mechanism_evidence", "evidence_location", "reviewer", "review_date", "review_status",
    ]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{field: row.get(field, "") for field in fields} for row in output])
    print(f"Wrote {len(output)} source-review rows to {args.out}")


if __name__ == "__main__":
    main()
