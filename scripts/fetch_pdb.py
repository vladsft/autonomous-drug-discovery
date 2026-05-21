#!/usr/bin/env python3
"""Fetch PDB files from RCSB into the local data/processed/ tree.

Used by `batch_cloud_run.py` to populate the input set before pods are
provisioned. Idempotent: a target already present locally is reported and
skipped (compare-by-size is intentional — we never want to re-download a
file we already have).

Validation: a downloaded file must (a) be > 1 KB, (b) start with one of the
PDB header records (HEADER / TITLE / ATOM / REMARK). RCSB serves a tiny HTML
error page for unknown codes; this check catches that and reports an error
without writing an obviously-broken PDB into the tree.

Usage:
    python scripts/fetch_pdb.py --targets 1M17 2HYY 6P3D
    python scripts/fetch_pdb.py --targets KRAS_G12C --out custom/dir
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUT = REPO_ROOT / "autonomous_drug_discovery" / "data" / "processed"

RCSB_URL = "https://files.rcsb.org/download/{target}.pdb"
# PDB codes are 4 chars: digit + 3 alphanumerics. We're more permissive here
# (any 4-char alphanumeric) so the script also accepts research aliases like
# "8P1L" or future RCSB extensions to 5-char codes.
TARGET_PATTERN = re.compile(r"^[A-Za-z0-9]{4,8}$")

# RCSB returns plain text PDB; first non-blank line should start with one of:
PDB_HEADER_RECORDS = ("HEADER", "TITLE ", "ATOM  ", "REMARK", "OBSLTE",
                      "EXPDTA", "MODEL ", "CRYST1")


def validate_target(target: str) -> bool:
    return bool(TARGET_PATTERN.match(target))


def is_valid_pdb_payload(data: bytes) -> bool:
    """Trip-wire against RCSB's HTML error page being saved as a PDB."""
    if len(data) < 1024:
        return False
    try:
        text = data[:4096].decode("ascii", errors="replace")
    except UnicodeDecodeError:
        return False
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        return any(stripped.startswith(rec) for rec in PDB_HEADER_RECORDS)
    return False


def fetch_one(target: str, out_dir: Path, force: bool) -> int:
    if not validate_target(target):
        print(f"[fetch] {target}: invalid format (expected 4-8 alphanumeric chars)", file=sys.stderr)
        return 1
    dest = out_dir / f"{target}.pdb"
    if dest.exists() and not force:
        print(f"[fetch] {target}: present already ({dest.stat().st_size // 1024} KB) — skipping")
        return 0

    url = RCSB_URL.format(target=target)
    print(f"[fetch] {target}: GET {url}")
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = resp.read()
    except urllib.error.HTTPError as e:
        print(f"[fetch] {target}: HTTP {e.code} — {e.reason}", file=sys.stderr)
        return 1
    except urllib.error.URLError as e:
        print(f"[fetch] {target}: network error — {e.reason}", file=sys.stderr)
        return 1

    if not is_valid_pdb_payload(data):
        print(f"[fetch] {target}: payload doesn't look like a PDB "
              f"({len(data)} bytes) — refusing to write", file=sys.stderr)
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    print(f"[fetch] {target}: wrote {dest} ({len(data) // 1024} KB)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--targets", nargs="+", required=True,
                    help="PDB codes (4-8 alphanumeric chars each)")
    ap.add_argument("--out", default=str(DEFAULT_OUT), type=Path,
                    help="Output directory (default: data/processed/)")
    ap.add_argument("--force", action="store_true",
                    help="Re-download even if the file already exists")
    args = ap.parse_args()

    failed: list[str] = []
    for t in args.targets:
        if fetch_one(t, args.out, args.force) != 0:
            failed.append(t)

    if failed:
        print(f"\n[fetch] FAILED: {' '.join(failed)}", file=sys.stderr)
        return 1
    print(f"\n[fetch] all {len(args.targets)} target(s) ready in {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
