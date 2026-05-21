#!/usr/bin/env python3
"""Regenerate the static chemist dashboard from local telemetry.

This is what `make dashboard` calls. Auto-discovers, per target and per
generator backend, the latest successful campaign in telemetry.db and assembles
a multi-backend `professor_demo.js` so the dashboard can toggle between
generators. The heavy lifting (RDKit drawing, ADMET flags, composite score)
lives in `build_demo_dataset.py`; this script is the glue that finds the
inputs and calls into it.

Usage:
    python scripts/regenerate_dashboard.py                       # target=1M17, all backends found
    python scripts/regenerate_dashboard.py --target 2HYY         # different target
    python scripts/regenerate_dashboard.py --data-dir data       # explicit data dir
    python scripts/regenerate_dashboard.py --out dashboard       # explicit dashboard dir

For each backend, the script looks for an AiZynthFinder output at
`<data_dir>/outputs/aizynth_<backend>.json` (the convention this repo already
uses for the EGFR demo). If absent, that backend's molecules are emitted
without synthesis annotations — the dashboard still loads.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PKG_ROOT = REPO_ROOT / "autonomous_drug_discovery"

# build_demo_dataset.py lives next to us; reuse its dataset assembly.
sys.path.insert(0, str(REPO_ROOT / "scripts"))
import build_demo_dataset as bdd  # noqa: E402

# A target's "known drug" caption on the dashboard. Add new targets here as we
# extend Phase 1; missing entries fall back to a generic caption.
TARGET_METADATA: dict[str, dict] = {
    "1M17": {
        "name": "EGFR (Epidermal Growth Factor Receptor)",
        "disease": "Non-small cell lung cancer",
        "known_drug": "Erlotinib",
    },
    "2HYY": {
        "name": "BCR-ABL kinase",
        "disease": "Chronic myeloid leukemia",
        "known_drug": "Imatinib",
    },
    "6P3D": {
        "name": "BRAF V600E",
        "disease": "Melanoma",
        "known_drug": "Ponatinib",
    },
    "8P1L": {
        "name": "8P1L (research target)",
        "disease": "Internal validation",
        "known_drug": None,
    },
}

# Order matters: the first backend with a complete campaign becomes the
# dashboard's default tab. RDKit ships in every image, so put it first.
BACKEND_ORDER = ["rdkit", "targetdiff", "pocket2mol"]


def latest_successful_campaign(
    db_path: Path, target: str, backend: str
) -> str | None:
    """Return the most recent campaign_id whose generation stage succeeded
    for (target, backend), or None if no such campaign exists."""
    if not db_path.exists():
        return None
    conn = sqlite3.connect(str(db_path))
    try:
        # Match the manifest by stem-prefix so both `<TARGET>_manifest.json`
        # and per-pocket variants like `<TARGET>_pocket2_manifest.json` count.
        rows = conn.execute(
            """SELECT campaign_id, parameters
               FROM runs
               WHERE module_name = '02_generation'
                 AND status = 'success'
                 AND input_path LIKE ?
               ORDER BY started_at DESC""",
            (f"%/{target}%_manifest.json",),
        ).fetchall()
    finally:
        conn.close()

    for cid, params_json in rows:
        try:
            params = json.loads(params_json) if params_json else {}
        except json.JSONDecodeError:
            continue
        if params.get("mode") == backend:
            return cid
    return None


# Map backend → list of acceptable AiZynth filenames inside data/outputs/. The
# bare `aizynth_<backend>.json` form is the canonical name going forward; the
# legacy aliases below are the names used by the original EGFR demo (2026-05-10).
# When you add a new backend, prefer the canonical form and skip the alias.
LEGACY_AIZYNTH_ALIASES: dict[str, list[str]] = {
    "rdkit": ["aizynth_rdkit.json", "aizynth_top10.json"],
    "pocket2mol": ["aizynth_pocket2mol.json", "aizynth_p2m.json"],
    "targetdiff": ["aizynth_targetdiff.json", "aizynth_td.json"],
}


def discover_aizynth(data_dir: Path, backend: str, campaign_id: str) -> Path | None:
    """Find an AiZynthFinder result JSON for this campaign, if any.

    Preference order:
      1. Per-campaign:  data/<campaign>/aizynth_top10.json (canonical, future)
      2. Global:        data/outputs/aizynth_<backend>.json (+ legacy aliases)
    """
    per_campaign = data_dir / campaign_id / "aizynth_top10.json"
    if per_campaign.exists():
        return per_campaign
    candidates = LEGACY_AIZYNTH_ALIASES.get(backend, [f"aizynth_{backend}.json"])
    for name in candidates:
        path = data_dir / "outputs" / name
        if path.exists():
            return path
    # build_demo_dataset.py expects a JSON file; an empty list keeps it happy
    # and the dashboard handles missing synthesis annotations gracefully.
    return None


def _build_empty_aizynth(tmp_path: Path) -> Path:
    """Write a sentinel empty AiZynth file so build_demo_dataset can be reused."""
    tmp_path.write_text("[]")
    return tmp_path


def assemble(
    target: str,
    data_dir: Path,
    out_dir: Path,
    max_per_backend: int,
    backends_requested: list[str],
) -> int:
    db_path = data_dir / "telemetry.db"
    if not db_path.exists():
        print(f"[regen] telemetry.db not found at {db_path}", file=sys.stderr)
        return 1

    # The campaign data lives under data/, but bdd is hard-coded to look
    # under autonomous_drug_discovery/data — repoint it at the actual dir.
    bdd.CAMPAIGN_BASE = data_dir

    backends_data: dict[str, dict] = {}
    insertion_order: list[str] = []
    empty_path = data_dir / ".empty_aizynth.json"
    empty_built = False

    for backend in backends_requested:
        cid = latest_successful_campaign(db_path, target, backend)
        if cid is None:
            print(f"[regen] no successful {backend} campaign for {target}; skipping")
            continue
        az_path = discover_aizynth(data_dir, backend, cid)
        if az_path is None:
            if not empty_built:
                _build_empty_aizynth(empty_path)
                empty_built = True
            az_path = empty_path
        print(f"[regen] {backend}: campaign={cid} aizynth={az_path.name}")
        try:
            backends_data[backend] = bdd.build_backend_dataset(
                backend, cid, az_path, max_molecules=max_per_backend
            )
            insertion_order.append(backend)
        except FileNotFoundError as e:
            print(f"[regen] {backend}: skipped — {e}", file=sys.stderr)

    if empty_built:
        empty_path.unlink(missing_ok=True)

    if not backends_data:
        print(f"[regen] no usable campaigns for target {target}", file=sys.stderr)
        return 1

    meta = TARGET_METADATA.get(target, {
        "name": target, "disease": "—", "known_drug": None,
    })
    doc = {
        "target": {"pdb": target, **meta},
        "generated_at": __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc
        ).isoformat(),
        "pipeline": {
            "screening_backend": "MolScore + ADMET-AI (104 properties)",
            "docking_backend": "AutoDock Vina (production)",
            "synthesis_backend": "AiZynthFinder (USPTO + ZINC)",
        },
        "backends": backends_data,
        "backend_order": insertion_order,
        "default_backend": insertion_order[0],
        "generator_descriptions": {
            name: bdd.GENERATOR_DESCRIPTIONS.get(name, "") for name in insertion_order
        },
        "admet_pass_rules": {
            k: {"kind": v[0], "threshold": v[1]}
            for k, v in bdd.ADMET_PASS_RULES.items()
        },
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "professor_demo.json"
    js_path = out_dir / "professor_demo.js"
    json_path.write_text(json.dumps(doc, indent=2))
    with js_path.open("w") as f:
        f.write("window.PROFESSOR_DEMO_DATA = ")
        json.dump(doc, f)
        f.write(";\n")

    print(f"[regen] wrote {json_path}")
    print(f"[regen] wrote {js_path}")
    print(f"[regen] backends: {', '.join(insertion_order)}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", default="1M17",
                    help="PDB stem (default: 1M17 — the EGFR demo target)")
    ap.add_argument("--data-dir", default=str(PKG_ROOT / "data"), type=Path,
                    help="Root data directory (must contain telemetry.db)")
    ap.add_argument("--out", default=str(REPO_ROOT / "dashboard"), type=Path,
                    help="Dashboard output directory")
    ap.add_argument("--max-per-backend", type=int, default=30,
                    help="Cap molecules per backend (default 30)")
    ap.add_argument("--backends", nargs="+", default=BACKEND_ORDER,
                    help=f"Backends to include (default: {' '.join(BACKEND_ORDER)})")
    args = ap.parse_args()
    return assemble(
        target=args.target,
        data_dir=Path(args.data_dir),
        out_dir=Path(args.out),
        max_per_backend=args.max_per_backend,
        backends_requested=list(args.backends),
    )


if __name__ == "__main__":
    sys.exit(main())
