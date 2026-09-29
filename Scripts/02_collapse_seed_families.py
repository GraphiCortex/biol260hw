#!/usr/bin/env python3
"""
02_collapse_seed_families.py

Collapse mature miRNAs with the same nt2-8 seed into seed families, rank those
families by the best sequence-level support found in the circRNA scan, and
separately report any canonical sites that cross the back-splice junction.

Inputs
------
1) Summary CSV produced by 01_scan_cdc14b_circrna_mirna.py
2) Site-level CSV produced by 01_scan_cdc14b_circrna_mirna.py

Outputs
-------
seed_family_summary.csv
sequence_shortlist.csv
junction_spanning_sites.csv

This script uses only the Python standard library.
"""

from __future__ import annotations
import argparse
import csv
from collections import defaultdict
from pathlib import Path

def as_int(x):
    return 0 if x in (None, "") else int(float(x))

def as_float(x):
    return 0.0 if x in (None, "") else float(x)

def read_csv(path):
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))

def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True, type=Path)
    ap.add_argument("--sites", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--min-score", type=int, default=4)
    args = ap.parse_args()

    if not args.summary.exists():
        raise SystemExit(f"ERROR: summary file not found: {args.summary}")
    if not args.sites.exists():
        raise SystemExit(f"ERROR: site-level file not found: {args.sites}")

    summary = read_csv(args.summary)
    sites = read_csv(args.sites)

    grouped = defaultdict(list)
    for row in summary:
        grouped[(row["circRNA"], row["seed_2_8"])].append(row)

    family_rows = []
    for (circrna, seed), members in grouped.items():
        members_sorted = sorted(
            members,
            key=lambda r: (
                -as_int(r["weighted_seed_score"]),
                -as_int(r["n_8mer"]),
                -as_int(r["n_canonical_sites"]),
                r["miRNA"],
            ),
        )
        rep = members_sorted[0]
        family_rows.append({
            "circRNA": circrna,
            "seed_2_8": seed,
            "representative_miRNA": rep["miRNA"],
            "family_members": ";".join(sorted(r["miRNA"] for r in members)),
            "n_mature_miRNAs_in_family": len(members),
            "n_8mer": as_int(rep["n_8mer"]),
            "n_7mer_m8": as_int(rep["n_7mer_m8"]),
            "n_7mer_A1": as_int(rep["n_7mer_A1"]),
            "n_canonical_sites": as_int(rep["n_canonical_sites"]),
            "n_junction_spanning_sites": as_int(rep["n_junction_spanning_sites"]),
            "site_density_per_kb": as_float(rep["site_density_per_kb"]),
            "weighted_seed_score": as_int(rep["weighted_seed_score"]),
        })

    family_rows.sort(
        key=lambda r: (
            r["circRNA"],
            -r["weighted_seed_score"],
            -r["n_8mer"],
            -r["n_canonical_sites"],
            r["representative_miRNA"],
        )
    )

    current = None
    rank = 0
    for row in family_rows:
        if row["circRNA"] != current:
            current = row["circRNA"]
            rank = 1
        else:
            rank += 1
        row["seed_family_rank"] = rank

    fields = [
        "seed_family_rank","circRNA","seed_2_8","representative_miRNA",
        "family_members","n_mature_miRNAs_in_family","n_8mer","n_7mer_m8",
        "n_7mer_A1","n_canonical_sites","n_junction_spanning_sites",
        "site_density_per_kb","weighted_seed_score",
    ]

    shortlist = [r.copy() for r in family_rows if r["weighted_seed_score"] >= args.min_score]

    junction_rows = []
    for row in sites:
        val = str(row.get("crosses_backsplice_junction", "")).strip().lower()
        if val in {"true","1","yes","y"}:
            junction_rows.append(row)

    junction_rows.sort(
        key=lambda r: (r["circRNA"], r["miRNA"], as_int(r["site_start_1based"]))
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    family_path = args.outdir / "seed_family_summary.csv"
    shortlist_path = args.outdir / "sequence_shortlist.csv"
    junction_path = args.outdir / "junction_spanning_sites.csv"

    write_csv(family_path, family_rows, fields)
    write_csv(shortlist_path, shortlist, fields)

    site_fields = list(sites[0].keys()) if sites else [
        "circRNA","circRNA_length_nt","miRNA","miRNA_sequence","seed_2_8",
        "site_type","site_start_1based","site_end_1based_circular",
        "crosses_backsplice_junction","target_site_5to3",
    ]
    write_csv(junction_path, junction_rows, site_fields)

    print(f"[ok] mature miRNAs in input: {len(summary)}")
    print(f"[ok] unique nt2-8 seed families: {len(family_rows)}")
    print(f"[ok] sequence shortlist (score >= {args.min_score}): {len(shortlist)} families")
    print(f"[ok] junction-spanning canonical sites: {len(junction_rows)}")
    print(f"[ok] wrote: {family_path}")
    print(f"[ok] wrote: {shortlist_path}")
    print(f"[ok] wrote: {junction_path}")
    print()
    print("NOTE: seed-family ranking is a transparent prioritization heuristic,")
    print("not experimental proof of miRNA sponging.")

if __name__ == "__main__":
    main()
