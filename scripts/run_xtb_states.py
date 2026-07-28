#!/usr/bin/env python3
"""Fault-tolerant parallel GFN2-xTB state calculations.

Each calculation runs in its own directory, avoiding xTB temporary-file
collisions. A timeout is recorded as a failed state, rather than aborting the
whole experiment. This command deliberately uses a GFN2 single-point
calculation after RDKit/MMFF geometry generation; calculation mode is retained
in the output and must be reported.
"""
import argparse
import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def state_succeeded(out_path):
    return out_path.exists() and "normal termination" in out_path.read_text(errors="ignore")


def run_one(row, outdir, timeout):
    state_id = row["state_id"]
    state_dir = outdir / state_id
    state_dir.mkdir(parents=True, exist_ok=True)
    out_path = state_dir / "xtb.out"
    if state_succeeded(out_path):
        return state_id, "cached_success"
    command = [
        "xtb", str(Path(row["xyz_path"]).resolve()), "--gfn2", "--chrg",
        str(row["formal_charge"]), "--uhf", "0", "--iterations", "300",
    ]
    try:
        with out_path.open("w") as handle:
            finished = subprocess.run(
                command,
                stdout=handle,
                stderr=subprocess.STDOUT,
                cwd=state_dir,
                timeout=timeout,
                check=False,
            )
        return state_id, "success" if finished.returncode == 0 and state_succeeded(out_path) else "failed"
    except subprocess.TimeoutExpired:
        with out_path.open("a") as handle:
            handle.write("\nTIMEOUT: calculation exceeded configured wall time\n")
        return state_id, "timeout"
    except FileNotFoundError:
        return state_id, "xtb_not_found"
    except Exception as exc:
        with out_path.open("a") as handle:
            handle.write(f"\nPYTHON_ERROR: {exc}\n")
        return state_id, "error"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("states", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=420)
    args = parser.parse_args()
    if args.workers < 1 or args.timeout < 30:
        raise SystemExit("workers must be >= 1 and timeout must be >= 30 seconds")
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(args.states.open()))
    if not rows:
        raise SystemExit("No states found")
    status_rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one, row, args.outdir, args.timeout): row for row in rows}
        for index, future in enumerate(as_completed(futures), 1):
            row = futures[future]
            try:
                state_id, status = future.result()
            except Exception as exc:
                state_id, status = row["state_id"], f"executor_error:{exc}"
            status_rows.append({"state_id": state_id, "status": status})
            print(f"{index}/{len(rows)} {state_id}: {status}", flush=True)
    with (args.outdir / "run_status.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["state_id", "status"])
        writer.writeheader()
        writer.writerows(sorted(status_rows, key=lambda row: row["state_id"]))


if __name__ == "__main__":
    main()
