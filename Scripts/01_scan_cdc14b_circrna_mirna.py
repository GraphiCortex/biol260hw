#!/usr/bin/env python3
"""
01_scan_cdc14b_circrna_mirna.py

Purpose
-------
Reconstruct the miRNA-site screening step of the BIOL 260 circRNA exercise in a
transparent, reproducible way.

The original CircInteractome uses TargetScan logic to identify candidate miRNA
sites on mature circRNAs. This script implements the three strongest canonical
TargetScan seed classes:

    8mer      : match to miRNA nts 2-8 + downstream A
    7mer-m8   : match to miRNA nts 2-8
    7mer-A1   : match to miRNA nts 2-7 + downstream A

Because a circRNA is circular, the script ALSO searches sites that cross the
back-splice junction. This is important and is often missed by naive linear
sequence scans.

This is a sequence-level prediction, NOT evidence that a circRNA actually
sponges a miRNA in vivo.

Usage
-----
Option A: let the script download a pinned miRBase v22.1 human mature-miRNA FASTA

    python 01_scan_cdc14b_circrna_mirna.py \
        --circ-fasta cdc14b_circ_0087640.fasta \
        --out cdc14b_circ_0087640_mirna_scan.csv

Option B: supply your own miRNA FASTA (full miRBase mature.fa is okay; only hsa-
records are kept)

    python 01_scan_cdc14b_circrna_mirna.py \
        --circ-fasta cdc14b_circ_0087640.fasta \
        --mirna-fasta mature.fa \
        --out cdc14b_circ_0087640_mirna_scan.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

DEFAULT_MIRBASE_URL = (
    "https://raw.githubusercontent.com/KechrisLab/miR-MaGiC/master/"
    "resources/reference_sequences/miRBase_v22.1_mature_sequences_Homo_sapiens.fasta"
)

SITE_WEIGHT = {
    "8mer": 3,
    "7mer-m8": 2,
    "7mer-A1": 1,
}

DNA_COMPLEMENT = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def reverse_complement(seq: str) -> str:
    return seq.translate(DNA_COMPLEMENT)[::-1].upper().replace("U", "T")


def normalize_dna(seq: str) -> str:
    seq = "".join(seq.split()).upper().replace("U", "T")
    bad = set(seq) - set("ACGTN")
    if bad:
        raise ValueError(f"Unexpected sequence characters: {sorted(bad)}")
    return seq


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, List[str]] = {}
    current = None
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(f"Duplicate FASTA ID: {current}")
                records[current] = []
            else:
                if current is None:
                    raise ValueError("FASTA sequence encountered before a header.")
                records[current].append(line)
    return {k: normalize_dna("".join(v)) for k, v in records.items()}


def download_default_mirbase(dest: Path) -> None:
    print(f"[download] miRBase v22.1 mature human FASTA -> {dest}")
    try:
        urllib.request.urlretrieve(DEFAULT_MIRBASE_URL, dest)
    except Exception as exc:
        raise RuntimeError(
            "Could not download the pinned miRBase FASTA. "
            "Download a mature-miRNA FASTA manually and pass it with "
            "--mirna-fasta."
        ) from exc


def circular_subseq(seq: str, start: int, length: int) -> str:
    L = len(seq)
    return "".join(seq[(start + j) % L] for j in range(length))


def site_crosses_junction(start: int, length: int, L: int) -> bool:
    return start + length > L


def scan_one_mirna(circ: str, mirna: str) -> List[Tuple[str, int, int, bool, str]]:
    """
    Return canonical sites as:
      (site_type, start_1based, end_1based_circular, crosses_junction, target_site)

    Nested 7mer components of an 8mer are not double-counted.
    """
    L = len(circ)
    mirna = normalize_dna(mirna)
    if len(mirna) < 8:
        return []

    motif_2_8 = reverse_complement(mirna[1:8])
    motif_2_7 = reverse_complement(mirna[1:7])

    eightmer_starts = set()
    sites: List[Tuple[str, int, int, bool, str]] = []

    for s in range(L):
        w8 = circular_subseq(circ, s, 8)
        if w8 == motif_2_8 + "A":
            eightmer_starts.add(s)
            sites.append(
                ("8mer", s + 1, ((s + 7) % L) + 1,
                 site_crosses_junction(s, 8, L), w8)
            )

    for s in range(L):
        if s in eightmer_starts:
            continue
        w7 = circular_subseq(circ, s, 7)
        if w7 == motif_2_8:
            sites.append(
                ("7mer-m8", s + 1, ((s + 6) % L) + 1,
                 site_crosses_junction(s, 7, L), w7)
            )

    for s in range(L):
        prev = (s - 1) % L
        if prev in eightmer_starts:
            continue
        w7 = circular_subseq(circ, s, 7)
        if w7 == motif_2_7 + "A":
            sites.append(
                ("7mer-A1", s + 1, ((s + 6) % L) + 1,
                 site_crosses_junction(s, 7, L), w7)
            )

    class_order = {"8mer": 0, "7mer-m8": 1, "7mer-A1": 2}
    sites.sort(key=lambda x: (class_order[x[0]], x[1]))
    return sites


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--circ-fasta", required=True, type=Path,
                    help="FASTA containing one or more mature circRNA sequences.")
    ap.add_argument("--mirna-fasta", type=Path, default=None,
                    help="Mature miRNA FASTA. If omitted, pinned miRBase v22.1 "
                         "human FASTA is downloaded.")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output CSV with one row per circRNA-miRNA pair.")
    ap.add_argument("--sites-out", type=Path, default=None,
                    help="Optional site-level CSV. Default: <out stem>_sites.csv")
    ap.add_argument("--top", type=int, default=100,
                    help="Number of ranked miRNAs per circRNA to keep in summary "
                         "(default: 100; use 0 for all).")
    args = ap.parse_args()

    if not args.circ_fasta.exists():
        sys.exit(f"ERROR: circRNA FASTA not found: {args.circ_fasta}")

    mirna_path = args.mirna_fasta
    if mirna_path is None:
        mirna_path = args.out.parent / "miRBase_v22.1_mature_Homo_sapiens.fasta"
        if not mirna_path.exists():
            download_default_mirbase(mirna_path)

    if not mirna_path.exists():
        sys.exit(f"ERROR: miRNA FASTA not found: {mirna_path}")

    circ_records = read_fasta(args.circ_fasta)
    mirna_records = {
        k: v for k, v in read_fasta(mirna_path).items()
        if k.lower().startswith("hsa-")
    }

    if not mirna_records:
        sys.exit("ERROR: no human miRNA records (IDs starting with 'hsa-') found.")

    summary_rows = []
    site_rows = []

    for circ_id, circ_seq in circ_records.items():
        L = len(circ_seq)
        rows_this_circ = []

        for mirna_id, mirna_seq in mirna_records.items():
            sites = scan_one_mirna(circ_seq, mirna_seq)
            if not sites:
                continue

            counts = {k: 0 for k in SITE_WEIGHT}
            for stype, start, end, crosses, target_site in sites:
                counts[stype] += 1
                site_rows.append({
                    "circRNA": circ_id,
                    "circRNA_length_nt": L,
                    "miRNA": mirna_id,
                    "miRNA_sequence": mirna_seq.replace("T", "U"),
                    "seed_2_8": mirna_seq[1:8].replace("T", "U"),
                    "site_type": stype,
                    "site_start_1based": start,
                    "site_end_1based_circular": end,
                    "crosses_backsplice_junction": crosses,
                    "target_site_5to3": target_site,
                })

            weighted = sum(counts[k] * SITE_WEIGHT[k] for k in SITE_WEIGHT)
            site_count = sum(counts.values())
            junction_count = sum(1 for x in sites if x[3])
            rows_this_circ.append({
                "circRNA": circ_id,
                "circRNA_length_nt": L,
                "miRNA": mirna_id,
                "miRNA_sequence": mirna_seq.replace("T", "U"),
                "seed_2_8": mirna_seq[1:8].replace("T", "U"),
                "n_8mer": counts["8mer"],
                "n_7mer_m8": counts["7mer-m8"],
                "n_7mer_A1": counts["7mer-A1"],
                "n_canonical_sites": site_count,
                "n_junction_spanning_sites": junction_count,
                "site_density_per_kb": round(site_count / L * 1000.0, 4),
                "weighted_seed_score": weighted,
            })

        rows_this_circ.sort(
            key=lambda r: (
                -r["weighted_seed_score"],
                -r["n_8mer"],
                -r["n_canonical_sites"],
                r["miRNA"],
            )
        )

        if args.top > 0:
            rows_this_circ = rows_this_circ[:args.top]

        for rank, row in enumerate(rows_this_circ, start=1):
            row["rank_within_circRNA"] = rank
            summary_rows.append(row)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    sites_out = args.sites_out or args.out.with_name(args.out.stem + "_sites.csv")

    summary_fields = [
        "rank_within_circRNA", "circRNA", "circRNA_length_nt", "miRNA",
        "miRNA_sequence", "seed_2_8", "n_8mer", "n_7mer_m8", "n_7mer_A1",
        "n_canonical_sites", "n_junction_spanning_sites", "site_density_per_kb",
        "weighted_seed_score",
    ]
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=summary_fields)
        w.writeheader()
        w.writerows(summary_rows)

    site_fields = [
        "circRNA", "circRNA_length_nt", "miRNA", "miRNA_sequence", "seed_2_8",
        "site_type", "site_start_1based", "site_end_1based_circular",
        "crosses_backsplice_junction", "target_site_5to3",
    ]
    with sites_out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=site_fields)
        w.writeheader()
        w.writerows(site_rows)

    print(f"[ok] circRNAs scanned: {len(circ_records)}")
    print(f"[ok] human miRNAs scanned: {len(mirna_records)}")
    print(f"[ok] ranked summary: {args.out}")
    print(f"[ok] site-level output: {sites_out}")
    print()
    print("IMPORTANT: these are sequence predictions, not experimental proof of sponging.")


if __name__ == "__main__":
    main()
