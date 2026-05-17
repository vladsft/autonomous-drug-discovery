"""
Module 01: Ingestion — Pocket Detection (P2Rank / fpocket).

Two backends:
  - p2rank (default): ML-based pocket prediction, 10-20% better recall than fpocket.
  - fpocket: Geometry-based detection, fallback if Java/P2Rank unavailable.

Produces a manifest.json (pocket data) and run_metadata.json (telemetry).

Input contract:  .pdb file path
Output contract: {stem}_manifest.json + run_metadata.json
"""

import os
import sys
import argparse
import csv
import subprocess
import shutil
import json
import hashlib
import re
import traceback
from datetime import datetime, timezone
from pathlib import Path

# Paths
MODULE_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = MODULE_DIR.parent.parent
FPOCKET_BIN = Path(os.environ.get("FPOCKET_BIN", os.path.expanduser("~/fpocket/bin/fpocket")))
P2RANK_BIN = Path(os.environ.get("P2RANK_BIN", os.path.expanduser("~/p2rank_2.5.1/prank")))

sys.path.insert(0, str(PROJECT_ROOT))
from telemetry import TelemetryDB


def _file_hash(filepath):
    """Compute SHA256 of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def _clean_pdb(pdb_in: Path, pdb_out: Path) -> Path:
    """Clean a PDB for pocket detection.

    Drops HETATM records (waters, ligands, cofactors), keeps only the first
    MODEL for NMR ensembles, and strips alternate-location atoms other than A.
    Returns the output path.
    """
    in_first_model = True
    saw_model = False
    with open(pdb_in) as fin, open(pdb_out, "w") as fout:
        for line in fin:
            rec = line[:6]
            if rec == "MODEL ":
                if saw_model:
                    in_first_model = False
                saw_model = True
                continue
            if rec == "ENDMDL":
                in_first_model = False
                continue
            if not in_first_model:
                continue
            if rec == "HETATM":
                continue
            if rec == "ATOM  ":
                alt_loc = line[16]
                if alt_loc not in (" ", "A"):
                    continue
            fout.write(line)
    return pdb_out


def _compute_pocket_center(pocket_pdb_path: Path) -> list[float] | None:
    """Return [cx, cy, cz] centroid of atoms in a pocket PDB file."""
    coords_x, coords_y, coords_z, n = 0.0, 0.0, 0.0, 0
    with open(pocket_pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    coords_x += float(line[30:38])
                    coords_y += float(line[38:46])
                    coords_z += float(line[46:54])
                    n += 1
                except ValueError:
                    continue
    if n == 0:
        return None
    return [round(coords_x / n, 3), round(coords_y / n, 3), round(coords_z / n, 3)]


def _parse_pocket_scores(info_txt_path):
    """Parse fpocket _info.txt file and return dict of {pocket_number: druggability_score}."""
    scores = {}
    current_pocket = None
    with open(info_txt_path, "r") as f:
        for line in f:
            m = re.match(r"^Pocket\s+(\d+)\s*:", line)
            if m:
                current_pocket = int(m.group(1))
                continue
            if current_pocket is not None:
                m = re.match(r"\s+Druggability Score\s*:\s+([\d.]+)", line)
                if m:
                    scores[current_pocket] = float(m.group(1))
                    current_pocket = None
    return scores


_P2RANK_REQUIRED_COLS = (
    "name", "rank", "score", "probability",
    "center_x", "center_y", "center_z", "residue_ids",
)


def _parse_p2rank_predictions(predictions_csv):
    """Parse P2Rank predictions CSV. Returns list of dicts sorted by rank.

    Robust to column reordering and optional columns; raises a clear error
    if a required column is missing so we fail loud rather than silently.
    """
    pockets = []
    with open(predictions_csv, newline="") as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        if reader.fieldnames:
            reader.fieldnames = [h.strip() for h in reader.fieldnames]
        for raw in reader:
            row = {(k or "").strip(): (v.strip() if isinstance(v, str) else v)
                   for k, v in raw.items()}
            missing = [c for c in _P2RANK_REQUIRED_COLS if c not in row]
            if missing:
                raise RuntimeError(
                    f"P2Rank predictions CSV missing required columns: {missing}. "
                    f"Found columns: {list(row)}"
                )
            pockets.append({
                "name": row["name"],
                "rank": int(row["rank"]),
                "score": float(row["score"]),
                "probability": float(row["probability"]),
                "center_x": float(row["center_x"]),
                "center_y": float(row["center_y"]),
                "center_z": float(row["center_z"]),
                "residue_ids": row["residue_ids"],
                "surf_atoms": int(row["surf_atoms"]) if row.get("surf_atoms") else None,
            })
    pockets.sort(key=lambda p: p["rank"])
    return pockets


def run_p2rank(pdb_file, output_dir, db_path=None, campaign_id=None, clean=False):
    """Run P2Rank on a PDB file with full telemetry capture."""
    pdb_path = Path(pdb_file).resolve()
    out_path = Path(output_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).isoformat()

    db = None
    run_id = None
    parameters = {"backend": "p2rank", "p2rank_binary": str(P2RANK_BIN), "clean": bool(clean)}

    if db_path and campaign_id:
        db = TelemetryDB(db_path)
        run_id = db.start_run(
            campaign_id=campaign_id,
            module_name="01_ingestion",
            input_path=str(pdb_path),
            parameters=parameters,
        )

    try:
        if not pdb_path.exists():
            raise FileNotFoundError(f"Input PDB file {pdb_path} does not exist.")
        if not P2RANK_BIN.exists():
            raise FileNotFoundError(f"P2Rank binary not found at {P2RANK_BIN}")

        input_hash = _file_hash(pdb_path)

        # Optionally clean: drop HETATM/water/alt locs; write next to output dir.
        run_pdb = pdb_path
        if clean:
            cleaned = out_path / f"{pdb_path.stem}_clean.pdb"
            _clean_pdb(pdb_path, cleaned)
            run_pdb = cleaned
            print(f"[Ingestion] Cleaned PDB written to {cleaned}")

        p2rank_out = out_path / f"{pdb_path.stem}_p2rank"

        print(f"[Ingestion] Running P2Rank on {run_pdb.name}...")

        cmd = [
            str(P2RANK_BIN), "predict",
            "-f", str(run_pdb),
            "-o", str(p2rank_out),
            "-visualizations", "0",
        ]
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Parse predictions (filename uses the PDB we actually ran against)
        predictions_csv = p2rank_out / f"{run_pdb.name}_predictions.csv"
        if not predictions_csv.exists():
            # Fall back to the original PDB filename if P2Rank used that
            predictions_csv = p2rank_out / f"{pdb_path.name}_predictions.csv"
        if not predictions_csv.exists():
            raise RuntimeError(f"P2Rank predictions not found in {p2rank_out}")

        pockets = _parse_p2rank_predictions(predictions_csv)
        print(f"[Ingestion] P2Rank found {len(pockets)} pockets")

        # Fail loud on no pockets — downstream stages can't do anything useful.
        if not pockets:
            raise RuntimeError(
                f"P2Rank found 0 pockets in {pdb_path.name}. "
                "Downstream generation/docking cannot proceed without a pocket."
            )

        # Build manifest — write best pocket residues as a PDB for downstream compatibility
        best_pocket = pockets[0]
        best_pocket_pdb = None

        print(f"[Ingestion] Best pocket: {best_pocket['name']} "
              f"(score: {best_pocket['score']:.2f}, probability: {best_pocket['probability']:.3f})")

        # Extract pocket atoms from PDB for docking box centroid
        residue_ids = best_pocket["residue_ids"].split()
        pocket_residues = set()
        for rid in residue_ids:
            parts = rid.split("_")
            if len(parts) == 2:
                try:
                    pocket_residues.add((parts[0], int(parts[1])))
                except ValueError:
                    continue

        # Write pocket PDB (atoms belonging to pocket residues)
        best_pocket_pdb_path = p2rank_out / f"{pdb_path.stem}_pocket1_atm.pdb"
        with open(run_pdb) as fin, open(best_pocket_pdb_path, "w") as fout:
            for line in fin:
                if line.startswith(("ATOM", "HETATM")):
                    chain = line[21]
                    try:
                        resnum = int(line[22:26].strip())
                    except ValueError:
                        continue
                    if (chain, resnum) in pocket_residues:
                        fout.write(line)
        best_pocket_pdb = str(best_pocket_pdb_path)

        manifest = {
            "input_pdb": str(pdb_path),
            "run_pdb": str(run_pdb),
            # receptor_pdb is what downstream docking should use: the cleaned
            # structure when --clean was applied, else the original. Keeping it
            # distinct from input_pdb means cleaning actually reaches docking.
            "receptor_pdb": str(run_pdb),
            "pocket_backend": "p2rank",
            "p2rank_out_dir": str(p2rank_out),
            "pockets_found": len(pockets),
            "best_pocket": best_pocket_pdb,
            "best_pocket_score": best_pocket["score"],
            "best_pocket_probability": best_pocket["probability"],
            "best_pocket_center": [
                best_pocket["center_x"], best_pocket["center_y"], best_pocket["center_z"]
            ],
        }

        manifest_path = out_path / f"{pdb_path.stem}_manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=4)

        metadata = {
            "module": "01_ingestion",
            "timestamp": timestamp,
            "backend": "p2rank",
            "input_pdb": str(pdb_path),
            "input_hash_sha256": input_hash,
            "pockets_found": len(pockets),
            "best_pocket": best_pocket_pdb,
            "manifest_path": str(manifest_path),
            "status": "success",
        }
        with open(out_path / "run_metadata.json", "w") as f:
            json.dump(metadata, f, indent=2)

        if db and run_id:
            db.complete_run(run_id, "success", str(manifest_path))

        print(f"[Ingestion] Manifest written to {manifest_path}")

    except Exception as e:
        error_msg = traceback.format_exc()
        print(f"[Ingestion] FAILED: {e}")
        if db and run_id:
            db.complete_run(run_id, "failed", error_trace=error_msg)
        fail_metadata = {
            "module": "01_ingestion", "timestamp": timestamp,
            "input_pdb": str(pdb_path), "status": "failed",
            "error": str(e), "traceback": error_msg,
        }
        with open(out_path / "run_metadata.json", "w") as f:
            json.dump(fail_metadata, f, indent=2)
        sys.exit(1)
    finally:
        if db:
            db.close()


def _get_fpocket_version():
    """Attempt to get fpocket version string."""
    try:
        result = subprocess.run(
            [str(FPOCKET_BIN)], capture_output=True, text=True, timeout=5
        )
        # fpocket prints version info to stderr or stdout when called without args
        output = result.stdout + result.stderr
        for line in output.split("\n"):
            if "version" in line.lower() or "fpocket" in line.lower():
                return line.strip()
    except Exception:
        pass
    return None


def run_fpocket(pdb_file, output_dir, db_path=None, campaign_id=None, clean=False):
    """Run fpocket on a PDB file with full telemetry capture.

    Args:
        pdb_file: Path to input PDB file.
        output_dir: Directory for output files.
        db_path: Optional telemetry database path.
        campaign_id: Optional campaign identifier.
        clean: If True, drop HETATM/waters/alt locs before running fpocket.
    """
    pdb_path = Path(pdb_file).resolve()
    out_path = Path(output_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).isoformat()

    # Telemetry setup
    db = None
    run_id = None
    parameters = {
        "backend": "fpocket",
        "fpocket_binary": str(FPOCKET_BIN),
        "fpocket_version": _get_fpocket_version(),
        "clean": bool(clean),
    }

    if db_path and campaign_id:
        db = TelemetryDB(db_path)
        run_id = db.start_run(
            campaign_id=campaign_id,
            module_name="01_ingestion",
            input_path=str(pdb_path),
            parameters=parameters,
        )

    try:
        if not pdb_path.exists():
            raise FileNotFoundError(f"Input PDB file {pdb_path} does not exist.")

        if not FPOCKET_BIN.exists():
            raise FileNotFoundError(f"fpocket binary not found at {FPOCKET_BIN}")

        input_hash = _file_hash(pdb_path)

        # Optionally clean the PDB; regardless, always stage a copy in the
        # output dir because fpocket writes its results next to the input file.
        target_pdb_path = out_path / pdb_path.name
        if clean:
            _clean_pdb(pdb_path, target_pdb_path)
            print(f"[Ingestion] Cleaned PDB written to {target_pdb_path}")
        elif pdb_path.resolve() != target_pdb_path.resolve():
            shutil.copy2(pdb_path, target_pdb_path)

        print(f"[Ingestion] Running fpocket on {target_pdb_path.name}...")

        cmd = [str(FPOCKET_BIN), "-f", str(target_pdb_path)]

        subprocess.check_call(cmd, cwd=out_path, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Analyze output
        out_folder_name = f"{target_pdb_path.stem}_out"
        out_folder_path = out_path / out_folder_name

        if not out_folder_path.exists():
            raise RuntimeError(f"Expected output folder {out_folder_path} not found.")

        print(f"[Ingestion] fpocket finished. Results in {out_folder_path}")

        pockets_dir = out_folder_path / "pockets"
        if not pockets_dir.exists():
            raise RuntimeError(
                f"fpocket produced no pockets/ subdir under {out_folder_path}."
            )

        pockets = list(pockets_dir.glob("pocket*_atm.pdb"))
        if not pockets:
            raise RuntimeError(
                f"fpocket found 0 pockets in {pdb_path.name}. "
                "Downstream generation/docking cannot proceed without a pocket."
            )

        # Parse druggability scores from fpocket info file; fall back to pocket index
        info_txt = out_folder_path / f"{target_pdb_path.stem}_info.txt"
        best_score = None
        if info_txt.exists():
            pocket_scores = _parse_pocket_scores(info_txt)
            pockets.sort(
                key=lambda p: pocket_scores.get(
                    int(p.stem.replace("pocket", "").replace("_atm", "")), -1,
                ),
                reverse=True,
            )
            best_num = int(pockets[0].stem.replace("pocket", "").replace("_atm", ""))
            best_score = pocket_scores.get(best_num)
            print(f"[Ingestion] Best pocket: pocket{best_num} (Druggability Score: {best_score})")
        else:
            pockets.sort(
                key=lambda p: int(p.stem.replace("pocket", "").replace("_atm", ""))
            )
            print("[Ingestion] Warning: _info.txt not found, using fpocket default pocket ordering")

        best_pocket_pdb = str(pockets[0])
        best_pocket_center = _compute_pocket_center(pockets[0])

        manifest = {
            "input_pdb": str(pdb_path),
            "run_pdb": str(target_pdb_path),
            # receptor_pdb: the structure downstream docking should use (cleaned
            # when --clean was applied, else the staged copy).
            "receptor_pdb": str(target_pdb_path),
            "pocket_backend": "fpocket",
            "fpocket_out_dir": str(out_folder_path),
            "pockets_found": len(pockets),
            "best_pocket": best_pocket_pdb,
            "best_pocket_score": best_score,
            "best_pocket_center": best_pocket_center,
        }

        manifest_path = out_path / f"{target_pdb_path.stem}_manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=4)

        # Write run_metadata.json
        metadata = {
            "module": "01_ingestion",
            "timestamp": timestamp,
            "input_pdb": str(pdb_path),
            "input_hash_sha256": input_hash,
            "fpocket_version": parameters["fpocket_version"],
            "pockets_found": manifest["pockets_found"],
            "best_pocket": manifest["best_pocket"],
            "manifest_path": str(manifest_path),
            "status": "success",
        }
        metadata_path = out_path / "run_metadata.json"
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2)

        # Complete telemetry
        if db and run_id:
            db.complete_run(run_id, "success", str(manifest_path))

        print(f"[Ingestion] Manifest written to {manifest_path}")
        print(f"[Ingestion] Metadata written to {metadata_path}")

    except Exception as e:
        error_msg = traceback.format_exc()
        print(f"[Ingestion] FAILED: {e}")

        if db and run_id:
            db.complete_run(run_id, "failed", error_trace=error_msg)

        # Write failure metadata
        fail_metadata = {
            "module": "01_ingestion",
            "timestamp": timestamp,
            "input_pdb": str(pdb_path),
            "status": "failed",
            "error": str(e),
            "traceback": error_msg,
        }
        fail_path = out_path / "run_metadata.json"
        with open(fail_path, "w") as f:
            json.dump(fail_metadata, f, indent=2)

        sys.exit(1)

    finally:
        if db:
            db.close()


def main():
    parser = argparse.ArgumentParser(description="Module 01: Ingestion (P2Rank / fpocket)")
    parser.add_argument("--pdb", required=True, help="Input PDB file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--db_path", default=None, help="Path to telemetry database")
    parser.add_argument("--campaign_id", default=None, help="Campaign ID for telemetry")
    parser.add_argument("--backend", choices=["p2rank", "fpocket"], default="p2rank",
                        help="Pocket detection backend (default: p2rank)")
    parser.add_argument("--clean", action="store_true", help="Clean PDB before processing")
    args = parser.parse_args()

    if args.backend == "p2rank":
        run_p2rank(args.pdb, args.output_dir, args.db_path, args.campaign_id, clean=args.clean)
    else:
        run_fpocket(args.pdb, args.output_dir, args.db_path, args.campaign_id, clean=args.clean)


if __name__ == "__main__":
    main()
