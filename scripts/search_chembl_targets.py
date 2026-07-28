#!/usr/bin/env python3
"""Create a reviewable list of ChEMBL target candidates for expansion.

Target IDs are never guessed: review the output and retain only the human,
single-protein targets that match the intended pharmacological target.
"""
import argparse
import csv
import time
from pathlib import Path
from urllib.parse import quote

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
    parser.add_argument("--queries", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--limit", type=int, default=20)
    args = parser.parse_args()
    if args.limit < 1:
        raise SystemExit("--limit must be positive")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update({"User-Agent": "CovalentScope-QM/0.1 (target-discovery)"})
    output = []
    for query_row in csv.DictReader(args.queries.open()):
        label = query_row["target_label"].strip()
        query = query_row["query"].strip()
        urls = [
            f"{BASE}/target/search.json?q={quote(query)}&limit={args.limit}",
            f"{BASE}/target.json?pref_name__icontains={quote(query)}&limit={args.limit}",
        ]
        payload = None
        for url in urls:
            try:
                payload = get_json(session, url)
                break
            except RuntimeError as exc:
                print(f"WARNING: {label}: {exc}", flush=True)
        if not payload:
            continue
        candidates = payload.get("targets") or payload.get("target") or []
        seen = set()
        for rank, target in enumerate(candidates, 1):
            target_id = target.get("target_chembl_id", "")
            if not target_id or target_id in seen:
                continue
            seen.add(target_id)
            output.append({
                "requested_label": label,
                "query": query,
                "rank": rank,
                "target_chembl_id": target_id,
                "pref_name": target.get("pref_name", ""),
                "target_type": target.get("target_type", ""),
                "organism": target.get("organism", ""),
                "review_decision": "",
                "review_note": "Retain only the intended human single-protein target.",
            })
        print(f"{label}: {len(seen)} target candidates", flush=True)
    fields = [
        "requested_label", "query", "rank", "target_chembl_id", "pref_name",
        "target_type", "organism", "review_decision", "review_note",
    ]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)
    print(f"Wrote {len(output)} target candidates to {args.out}")


if __name__ == "__main__":
    main()
