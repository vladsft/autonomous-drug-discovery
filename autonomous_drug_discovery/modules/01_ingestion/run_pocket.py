"""
Module 01: Ingestion — fpocket Wrapper.

Wraps the fpocket binary to identify binding pockets on a protein PDB.
Produces a manifest.json (pocket data) and run_metadata.json (telemetry).

Input contract:  .pdb file path
Output contract: {stem}_manifest.json + run_metadata.json
"""

import os
import sys
import argparse
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
FPOCKET_BIN = Path("/home/vladsft/fpocket/bin/fpocket")

sys.path.insert(0, str(PROJECT_ROOT))
from telemetry import TelemetryDB


def _file_hash(filepath):
    """Compute SHA256 of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


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
        return "unknown"
    except Exception:
        return "unknown"


def run_fpocket(pdb_file, output_dir, db_path=None, campaign_id=None):
    """Run fpocket on a PDB file with full telemetry capture.

    Args:
        pdb_file: Path to input PDB file.
        output_dir: Directory for output files.
        db_path: Optional telemetry database path.
        campaign_id: Optional campaign identifier.
    """
    pdb_path = Path(pdb_file).resolve()
    out_path = Path(output_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).isoformat()

    # Telemetry setup
    db = None
    run_id = None
    parameters = {
        "fpocket_binary": str(FPOCKET_BIN),
        "fpocket_version": _get_fpocket_version(),
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

        print(f"[Ingestion] Running fpocket on {pdb_path.name}...")

        # Copy PDB to output directory for fpocket (skip if already there)
        target_pdb_name = pdb_path.name
        target_pdb_path = out_path / target_pdb_name
        if pdb_path.resolve() != target_pdb_path.resolve():
            shutil.copy2(pdb_path, target_pdb_path)

        cmd = [str(FPOCKET_BIN), "-f", str(target_pdb_path)]

        subprocess.check_call(cmd, cwd=out_path, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Analyze output
        out_folder_name = f"{target_pdb_path.stem}_out"
        out_folder_path = out_path / out_folder_name

        if not out_folder_path.exists():
            raise RuntimeError(f"Expected output folder {out_folder_path} not found.")

        print(f"[Ingestion] fpocket finished. Results in {out_folder_path}")

        # Create manifest
        pockets_dir = out_folder_path / "pockets"
        manifest = {
            "input_pdb": str(target_pdb_path),
            "fpocket_out_dir": str(out_folder_path),
            "pockets_found": 0,
            "best_pocket": None,
        }

        if pockets_dir.exists():
            pockets = list(pockets_dir.glob("pocket*_atm.pdb"))
            manifest["pockets_found"] = len(pockets)
            if pockets:
                # Parse druggability scores from fpocket info file and select best pocket
                info_txt = out_folder_path / f"{target_pdb_path.stem}_info.txt"
                if info_txt.exists():
                    pocket_scores = _parse_pocket_scores(info_txt)
                    # Sort by druggability score (descending), fall back to pocket number
                    pockets.sort(
                        key=lambda p: pocket_scores.get(
                            int(p.stem.replace("pocket", "").replace("_atm", "")),
                            -1,
                        ),
                        reverse=True,
                    )
                    best_num = int(pockets[0].stem.replace("pocket", "").replace("_atm", ""))
                    best_score = pocket_scores.get(best_num, "N/A")
                    print(f"[Ingestion] Best pocket: pocket{best_num} (Druggability Score: {best_score})")
                else:
                    # Fallback: sort numerically by pocket number (fpocket default ranking)
                    pockets.sort(
                        key=lambda p: int(p.stem.replace("pocket", "").replace("_atm", ""))
                    )
                    print("[Ingestion] Warning: _info.txt not found, using fpocket default pocket ordering")
                manifest["best_pocket"] = str(pockets[0])

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
    parser = argparse.ArgumentParser(description="Module 01: Ingestion (fpocket)")
    parser.add_argument("--pdb", required=True, help="Input PDB file")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--db_path", default=None, help="Path to telemetry database")
    parser.add_argument("--campaign_id", default=None, help="Campaign ID for telemetry")
    parser.add_argument("--clean", action="store_true", help="Clean PDB before processing")
    args = parser.parse_args()

    run_fpocket(args.pdb, args.output_dir, args.db_path, args.campaign_id)


if __name__ == "__main__":
    main()
