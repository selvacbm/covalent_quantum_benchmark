#!/usr/bin/env python3
"""Retrieve provenance-preserving activity records for configured ChEMBL targets."""
import argparse
import csv
import time
from pathlib import Path

import requests

BASE = "https://www.ebi.ac.uk/chembl/api/data"


def get_json(session, url, retries=5):
    last_error = None
    for attempt in range(retries):
        try:
            response = session.get(url, timeout=120)
            response.raise_for_status()
            return response.json()
        except Exception as exc:  # ChEMBL intermittently returns HTTP 500
            last_error = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"ChEMBL request failed after {retries} attempts: {url}; {last_error}")


def read_targets(path):
    with path.open(newline="") as handle:
        targets = list(csv.DictReader(handle))
    required = {"target_label", "target_chembl_id"}
    if not targets or not required.issubset(targets[0]):
        raise ValueError("Target CSV must contain target_label and target_chembl_id columns")
    return targets


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--max-pages", type=int, default=30)
    parser.add_argument("--page-size", type=int, default=100)
    args = parser.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    targets = read_targets(args.targets)
    session = requests.Session()
    session.headers.update({"User-Agent": "CovalentScope-QM/0.1 (reproducible-benchmark)"})
    activities = []

    for target in targets:
        target_id = target["target_chembl_id"].strip()
        label = target["target_label"].strip()
        print(f"Fetching {label} ({target_id})", flush=True)
        for page in range(args.max_pages):
            offset = page * args.page_size
            url = f"{BASE}/activity.json?target_chembl_id__exact={target_id}&limit={args.page_size}&offset={offset}"
            try:
                payload = get_json(session, url)
            except RuntimeError as exc:
                print(f"WARNING: stopping {label} at page {page + 1}: {exc}", flush=True)
                break
            batch = payload.get("activities", [])
            if not batch:
                break
            for row in batch:
                row["target_label"] = label
                activities.append(row)
            print(f"  page={page + 1} target_rows={len(batch)} total={len(activities)}", flush=True)
            if len(batch) < args.page_size:
                break

    structure_cache = {}
    molecule_ids = sorted({row.get("molecule_chembl_id", "") for row in activities if row.get("molecule_chembl_id")})
    for index, molecule_id in enumerate(molecule_ids, 1):
        try:
            molecule = get_json(session, f"{BASE}/molecule/{molecule_id}.json")
            structure_cache[molecule_id] = (molecule.get("molecule_structures") or {}).get("canonical_smiles", "")
        except RuntimeError as exc:
            print(f"WARNING: no structure for {molecule_id}: {exc}", flush=True)
        if index % 25 == 0 or index == len(molecule_ids):
            print(f"  structures={index}/{len(molecule_ids)}", flush=True)

    fields = [
        "target_label", "molecule_chembl_id", "target_chembl_id", "assay_chembl_id",
        "document_chembl_id", "standard_type", "standard_relation", "standard_value",
        "standard_units", "pchembl_value", "activity_comment", "canonical_smiles",
    ]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in activities:
            item = {field: row.get(field, "") for field in fields}
            item["canonical_smiles"] = row.get("canonical_smiles", "") or structure_cache.get(row.get("molecule_chembl_id", ""), "")
            writer.writerow(item)
    print(f"Wrote {len(activities)} raw activity records to {args.out}")


if __name__ == "__main__":
    main()
