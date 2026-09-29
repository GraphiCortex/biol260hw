#!/usr/bin/env python3
"""
03_extract_targetscan_targets_v2.py

Download TargetScanHuman 8.0 data and extract predicted mRNA targets for the
miRNA seed families shortlisted in Stage 2.

Why this script?
----------------
The BIOL 260 homework explicitly asks for TargetScan downstream targets.
Rather than manually copy three genes from the web interface, this script
downloads the TargetScanHuman tables and keeps the complete target space for
our shortlisted miRNA families. We can then choose biologically coherent
targets transparently.

Inputs
------
sequence_shortlist.csv   (output from 02_collapse_seed_families.py)

Outputs
-------
targetscan_family_map.csv
targetscan_targets_all.csv
targetscan_top_targets.csv
unmapped_seed_families.csv

Notes
-----
- Uses TargetScanHuman 8.0 "all predictions" summary table so that poorly
  conserved/non-conserved miRNA families are not silently discarded.
- More-negative cumulative weighted context++ scores predict stronger
  repression in the TargetScan7-style ranking.
- TargetScan 8.0 reports essentially the same target predictions/UTR profiles
  as release 7.x, while also offering newer ranking models on the website.
- Predictions are computational, not proof of direct regulation.

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path


TARGETSCAN_BASES = [
    "https://www.targetscan.org/vert_80/vert_80_data_download",
    "http://www.targetscan.org/vert_80/vert_80_data_download",
]

FAMILY_ZIP = "miR_Family_Info.txt.zip"
SUMMARY_ZIP = "Summary_Counts.all_predictions.txt.zip"


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (s or "").strip().lower()).strip()


def clean_mirna(s: str) -> str:
    return (s or "").strip().lower().replace("_", "-")


def clean_seed(s: str) -> str:
    return (s or "").strip().upper().replace("T", "U")


def find_col(fieldnames, candidates, required=True):
    """
    Flexible header lookup. Each candidate is a substring after normalization.
    """
    normalized = {f: norm(f) for f in fieldnames}
    cand_norm = [norm(c) for c in candidates]

    # exact normalized match first
    for c in cand_norm:
        for original, n in normalized.items():
            if n == c:
                return original

    # substring match second
    for c in cand_norm:
        for original, n in normalized.items():
            if c in n:
                return original

    if required:
        raise KeyError(
            f"Could not find any of {candidates} in headers:\n{fieldnames}"
        )
    return None


def read_csv(path: Path):
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _validate_zip(path: Path) -> None:
    if not path.exists() or path.stat().st_size < 100:
        raise RuntimeError(f"Downloaded file is empty or too small: {path}")
    if not zipfile.is_zipfile(path):
        # This catches HTML error pages saved with a .zip extension.
        try:
            preview = path.read_text(encoding="utf-8", errors="replace")[:300]
        except Exception:
            preview = "<binary/non-text content>"
        raise RuntimeError(
            f"Downloaded file is not a valid ZIP: {path}\n"
            f"First bytes/text: {preview!r}"
        )


def download_with_curl(url: str, dest: Path):
    """
    TargetScan's web server can terminate Python/OpenSSL TLS handshakes early
    on some Windows systems. Windows 10/11 ships curl.exe, which is much more
    tolerant here. -k is intentional for this public data download only.
    """
    curl = shutil.which("curl.exe") or shutil.which("curl")
    if not curl:
        raise RuntimeError("curl/curl.exe not found")

    tmp = dest.with_suffix(dest.suffix + ".part")
    if tmp.exists():
        tmp.unlink()

    print(f"[curl] {url}")
    cmd = [
        curl,
        "-L",                    # follow redirects
        "--fail",               # fail on HTTP >= 400
        "--retry", "5",
        "--retry-delay", "2",
        "--retry-all-errors",
        "--connect-timeout", "30",
        "--max-time", "600",
        "--http1.1",
        "-k",                    # tolerate TargetScan TLS/certificate quirks
        "-A", "Mozilla/5.0 BIOL260-reproducible-analysis",
        "-o", str(tmp),
        url,
    ]
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if proc.returncode != 0:
        if tmp.exists():
            tmp.unlink()
        raise RuntimeError(
            f"curl failed with exit code {proc.returncode}: "
            f"{proc.stderr.strip()}"
        )

    tmp.replace(dest)
    _validate_zip(dest)


def download_with_urllib(url: str, dest: Path):
    print(f"[urllib] {url}")
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0 BIOL260-reproducible-analysis"},
    )
    tmp = dest.with_suffix(dest.suffix + ".part")
    if tmp.exists():
        tmp.unlink()
    try:
        with urllib.request.urlopen(req, timeout=120) as r, tmp.open("wb") as out:
            shutil.copyfileobj(r, out)
        tmp.replace(dest)
        _validate_zip(dest)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise


def download_with_fallback(filename: str, dest: Path):
    # Reuse a previously completed manual or scripted download.
    if dest.exists():
        try:
            _validate_zip(dest)
            print(f"[cached] {dest}")
            return
        except Exception:
            print(f"[warning] removing invalid cached file: {dest}")
            dest.unlink()

    errors = []
    for base in TARGETSCAN_BASES:
        url = f"{base}/{filename}"

        # 1) Prefer curl on Windows/Linux because TargetScan has known TLS
        #    compatibility problems with Python urllib on some machines.
        try:
            download_with_curl(url, dest)
            return
        except Exception as exc:
            errors.append(f"CURL {url}: {exc}")

        # 2) Keep urllib as a secondary route.
        try:
            download_with_urllib(url, dest)
            return
        except Exception as exc:
            errors.append(f"URLLIB {url}: {exc}")

    raise RuntimeError(
        "Could not download TargetScan data automatically.\n\n"
        + "\n".join(errors)
        + "\n\nFallback: download these two files manually from the "
          "TargetScanHuman 8.0 data-download page:\n"
          "  miR_Family_Info.txt.zip\n"
          "  Summary_Counts.all_predictions.txt.zip\n"
          "and place them in:\n"
          f"  {dest.parent}\n"
          "Then rerun the SAME command. The script will detect and reuse them."
    )


def open_first_text_member(zip_path: Path, preferred_fragment: str):
    zf = zipfile.ZipFile(zip_path)
    members = [
        x for x in zf.namelist()
        if not x.endswith("/") and x.lower().endswith((".txt", ".tsv"))
    ]
    if not members:
        raise RuntimeError(f"No text table found inside {zip_path}")

    preferred = [
        x for x in members if preferred_fragment.lower() in x.lower()
    ]
    member = preferred[0] if preferred else members[0]
    raw = zf.open(member, "r")
    text = io.TextIOWrapper(raw, encoding="utf-8-sig", newline="")
    return zf, text, member


def to_float(x, default=float("inf")):
    try:
        if x is None or str(x).strip() in {"", "NA", "N/A"}:
            return default
        return float(x)
    except Exception:
        return default


def to_int(x, default=0):
    try:
        if x is None or str(x).strip() in {"", "NA", "N/A"}:
            return default
        return int(float(x))
    except Exception:
        return default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shortlist", required=True, type=Path,
                    help="sequence_shortlist.csv from Stage 2.")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--top-per-family", type=int, default=100,
                    help="Top unique gene symbols retained per TargetScan family.")
    args = ap.parse_args()

    if not args.shortlist.exists():
        sys.exit(f"ERROR: shortlist file not found: {args.shortlist}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    dldir = args.outdir / "downloads"
    dldir.mkdir(parents=True, exist_ok=True)

    fam_zip = dldir / FAMILY_ZIP
    sum_zip = dldir / SUMMARY_ZIP
    download_with_fallback(FAMILY_ZIP, fam_zip)
    download_with_fallback(SUMMARY_ZIP, sum_zip)

    shortlist = read_csv(args.shortlist)

    # ---------------------------------------------------------------
    # Read TargetScan miRNA-family table
    # ---------------------------------------------------------------
    zf, fh, member = open_first_text_member(fam_zip, "miR_Family_Info")
    fam_reader = csv.DictReader(fh, delimiter="\t")
    fam_fields = fam_reader.fieldnames or []

    col_ts_family = find_col(fam_fields, ["miR family", "miRNA family"])
    col_seed = find_col(fam_fields, ["Seed+m8", "seed m8", "seed"], required=False)
    col_species = find_col(fam_fields, ["Species ID"], required=False)
    col_mirbase = find_col(fam_fields, ["MiRBase ID", "miRBase"], required=False)

    ts_rows = []
    for row in fam_reader:
        if col_species and str(row.get(col_species, "")).strip() not in {"9606", "9606.0"}:
            continue
        ts_rows.append(row)
    fh.close()
    zf.close()

    by_mirna = defaultdict(set)
    by_seed = defaultdict(set)

    for row in ts_rows:
        fam = (row.get(col_ts_family) or "").strip()
        if not fam:
            continue
        if col_mirbase:
            mid = clean_mirna(row.get(col_mirbase, ""))
            if mid:
                by_mirna[mid].add(fam)
        if col_seed:
            seed = clean_seed(row.get(col_seed, ""))
            if seed:
                by_seed[seed].add(fam)

    family_map_rows = []
    wanted_ts_families = set()
    seed_to_ts = {}

    for srow in shortlist:
        seed = clean_seed(srow["seed_2_8"])
        names = []
        for x in (srow.get("family_members") or "").split(";"):
            if x.strip():
                names.append(clean_mirna(x))
        rep = clean_mirna(srow.get("representative_miRNA", ""))
        if rep and rep not in names:
            names.append(rep)

        exact = set()
        for n in names:
            exact |= by_mirna.get(n, set())

        seed_matches = by_seed.get(seed, set())
        mapped = exact if exact else seed_matches

        seed_to_ts[seed] = sorted(mapped)
        wanted_ts_families |= mapped

        family_map_rows.append({
            "seed_family_rank": srow.get("seed_family_rank", ""),
            "seed_2_8": seed,
            "representative_miRNA": srow.get("representative_miRNA", ""),
            "family_members": srow.get("family_members", ""),
            "weighted_seed_score": srow.get("weighted_seed_score", ""),
            "mapping_method": "MiRBase_ID" if exact else ("seed_2_8" if seed_matches else "UNMAPPED"),
            "TargetScan_families": ";".join(sorted(mapped)),
            "n_TargetScan_families": len(mapped),
        })

    map_fields = [
        "seed_family_rank", "seed_2_8", "representative_miRNA",
        "family_members", "weighted_seed_score", "mapping_method",
        "TargetScan_families", "n_TargetScan_families"
    ]
    write_csv(args.outdir / "targetscan_family_map.csv",
              family_map_rows, map_fields)

    unmapped = [r for r in family_map_rows if r["n_TargetScan_families"] == 0]
    write_csv(args.outdir / "unmapped_seed_families.csv",
              unmapped, map_fields)

    print(f"[ok] shortlisted seed families: {len(shortlist)}")
    print(f"[ok] mapped TargetScan families: {len(wanted_ts_families)}")
    print(f"[ok] unmapped seed families: {len(unmapped)}")

    # ---------------------------------------------------------------
    # Stream TargetScan all-predictions summary and retain only the
    # TargetScan families we need.
    # ---------------------------------------------------------------
    zf, fh, member = open_first_text_member(sum_zip, "Summary_Counts")
    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames or []

    col_family = find_col(fields, ["miRNA family", "miR family"])
    col_species2 = find_col(fields, ["Species ID"], required=False)
    col_gene = find_col(fields, ["Gene Symbol"])
    col_gene_id = find_col(fields, ["Gene ID"], required=False)
    col_transcript = find_col(fields, ["Transcript ID"], required=False)
    col_cwcs = find_col(
        fields,
        ["Cumulative weighted context++ score",
         "cumulative weighted context score"],
        required=False,
    )
    col_cons = find_col(fields, ["Total num conserved sites"], required=False)
    col_noncons = find_col(fields, ["Total num nonconserved sites"], required=False)
    col_pct = find_col(fields, ["Aggregate PCT"], required=False)

    # reverse map TargetScan family -> our shortlisted seed-family rows
    ts_to_short = defaultdict(list)
    for row in family_map_rows:
        for fam in (row["TargetScan_families"] or "").split(";"):
            if fam:
                ts_to_short[fam].append(row)

    kept = []
    for row in reader:
        fam = (row.get(col_family) or "").strip()
        if fam not in wanted_ts_families:
            continue
        if col_species2 and str(row.get(col_species2, "")).strip() not in {"9606", "9606.0"}:
            continue

        for smap in ts_to_short[fam]:
            kept.append({
                "seed_family_rank": smap["seed_family_rank"],
                "seed_2_8": smap["seed_2_8"],
                "representative_miRNA": smap["representative_miRNA"],
                "TargetScan_family": fam,
                "Gene_Symbol": row.get(col_gene, ""),
                "Gene_ID": row.get(col_gene_id, "") if col_gene_id else "",
                "Transcript_ID": row.get(col_transcript, "") if col_transcript else "",
                "Cumulative_weighted_context_score":
                    row.get(col_cwcs, "") if col_cwcs else "",
                "Total_num_conserved_sites":
                    row.get(col_cons, "") if col_cons else "",
                "Total_num_nonconserved_sites":
                    row.get(col_noncons, "") if col_noncons else "",
                "Aggregate_PCT":
                    row.get(col_pct, "") if col_pct else "",
            })

    fh.close()
    zf.close()

    target_fields = [
        "seed_family_rank", "seed_2_8", "representative_miRNA",
        "TargetScan_family", "Gene_Symbol", "Gene_ID", "Transcript_ID",
        "Cumulative_weighted_context_score",
        "Total_num_conserved_sites",
        "Total_num_nonconserved_sites",
        "Aggregate_PCT",
    ]

    kept.sort(
        key=lambda r: (
            to_int(r["seed_family_rank"], 10**9),
            r["TargetScan_family"],
            to_float(r["Cumulative_weighted_context_score"]),
            r["Gene_Symbol"],
        )
    )
    write_csv(args.outdir / "targetscan_targets_all.csv",
              kept, target_fields)

    # ---------------------------------------------------------------
    # Top unique gene symbols per shortlisted family.
    # Keep the transcript with the most favorable (most negative) CWCS.
    # ---------------------------------------------------------------
    grouped = defaultdict(list)
    for row in kept:
        key = (row["seed_2_8"], row["TargetScan_family"])
        grouped[key].append(row)

    top_rows = []
    for key, rows in grouped.items():
        rows.sort(key=lambda r: (
            to_float(r["Cumulative_weighted_context_score"]),
            -to_int(r["Total_num_conserved_sites"]),
            -to_float(r["Aggregate_PCT"], default=-1.0),
            r["Gene_Symbol"],
        ))

        seen = set()
        rank = 0
        for row in rows:
            gene = (row["Gene_Symbol"] or "").strip()
            if not gene or gene in seen:
                continue
            seen.add(gene)
            rank += 1
            out = dict(row)
            out["TargetScan_gene_rank"] = rank
            top_rows.append(out)
            if rank >= args.top_per_family:
                break

    top_fields = ["TargetScan_gene_rank"] + target_fields
    top_rows.sort(key=lambda r: (
        to_int(r["seed_family_rank"], 10**9),
        r["TargetScan_family"],
        to_int(r["TargetScan_gene_rank"], 10**9),
    ))
    write_csv(args.outdir / "targetscan_top_targets.csv",
              top_rows, top_fields)

    print(f"[ok] matching TargetScan target rows: {len(kept)}")
    print(f"[ok] wrote top {args.top_per_family} unique genes per mapped family")
    print(f"[ok] output directory: {args.outdir}")
    print()
    print("IMPORTANT:")
    print("TargetScan predictions are computational hypotheses, not proof of")
    print("direct miRNA-mediated repression. We will integrate these results")
    print("with breast-cancer literature rather than treating rank alone as truth.")


if __name__ == "__main__":
    main()
