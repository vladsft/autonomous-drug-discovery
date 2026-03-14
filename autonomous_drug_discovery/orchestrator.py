"""
Orchestrator: The central CLI driver for the Adaptive Discovery Pipeline.

Manages data flow between isolated modules, generates campaign IDs,
and ensures all runs are captured in the telemetry database.

Usage:
    python orchestrator.py ingest target.pdb
    python orchestrator.py generate manifest.json
    python orchestrator.py screen candidates.sdf
    python orchestrator.py dock manifest.json
    python orchestrator.py run target.pdb [--mode simulation|production]
"""

import os
import sys
import uuid
import argparse
import subprocess
from pathlib import Path

from telemetry import TelemetryDB

# Configuration
BASE_DIR = Path(__file__).parent.resolve()
DATA_DIR = BASE_DIR / "data"
MODULES_DIR = BASE_DIR / "modules"
DEFAULT_DB_PATH = str(DATA_DIR / "telemetry.db")


def ensure_dirs():
    """Ensure data directories exist."""
    for subdir in ["raw", "processed", "candidates", "screened", "results"]:
        (DATA_DIR / subdir).mkdir(parents=True, exist_ok=True)


def get_python_cmd(env_name):
    """Return the command to run python in a specific conda environment."""
    return ["conda", "run", "-n", env_name, "python"]


def run_ingestion(pdb_path, db_path, campaign_id, clean_pdb=False):
    """Run the ingestion module (fpocket wrapper)."""
    print(f"\n{'='*60}")
    print(f"[Orchestrator] Stage 1: INGESTION — {pdb_path}")
    print(f"{'='*60}")

    script_path = MODULES_DIR / "01_ingestion" / "run_pocket.py"
    if not script_path.exists():
        print(f"Error: Ingestion module not found at {script_path}")
        return False

    cmd = get_python_cmd("base") + [
        str(script_path),
        "--pdb", str(pdb_path),
        "--output_dir", str(DATA_DIR / "processed"),
        "--db_path", db_path,
        "--campaign_id", campaign_id,
    ]
    if clean_pdb:
        cmd.append("--clean")

    try:
        subprocess.check_call(cmd)
        print("[Orchestrator] Ingestion complete.\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[Orchestrator] ERROR: Ingestion failed with code {e.returncode}")
        return False


def run_generation(manifest_path, db_path, campaign_id, mode="simulation"):
    """Run the generation module (TargetDiff wrapper)."""
    print(f"\n{'='*60}")
    print(f"[Orchestrator] Stage 2: GENERATION — mode={mode}")
    print(f"{'='*60}")

    script_path = MODULES_DIR / "02_generation" / "run_generation.py"
    if not script_path.exists():
        print(f"Error: Generation module not found at {script_path}")
        return False

    # Use targetdiff_env only for targetdiff mode; base env for everything else
    if mode == "targetdiff":
        cmd = get_python_cmd("targetdiff_env")
    else:
        cmd = get_python_cmd("base")

    cmd += [
        str(script_path),
        "--manifest", str(manifest_path),
        "--output_dir", str(DATA_DIR / "candidates"),
        "--db_path", db_path,
        "--campaign_id", campaign_id,
        "--mode", mode,
    ]

    try:
        subprocess.check_call(cmd)
        print("[Orchestrator] Generation complete.\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[Orchestrator] ERROR: Generation failed with code {e.returncode}")
        return False


def run_screening(sdf_path, db_path, campaign_id):
    """Run the screening module (RDKit fast triage)."""
    print(f"\n{'='*60}")
    print(f"[Orchestrator] Stage 3: SCREENING (Fast Triage)")
    print(f"{'='*60}")

    script_path = MODULES_DIR / "03_screening" / "run_screening.py"
    if not script_path.exists():
        print(f"Error: Screening module not found at {script_path}")
        return False

    cmd = get_python_cmd("base") + [
        str(script_path),
        "--input_sdf", str(sdf_path),
        "--output_dir", str(DATA_DIR / "screened"),
        "--db_path", db_path,
        "--campaign_id", campaign_id,
    ]

    try:
        subprocess.check_call(cmd)
        print("[Orchestrator] Screening complete.\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[Orchestrator] ERROR: Screening failed with code {e.returncode}")
        return False


def run_docking(manifest_path, db_path, campaign_id, mode="simulation"):
    """Run the docking module (Vina wrapper)."""
    print(f"\n{'='*60}")
    print(f"[Orchestrator] Stage 4: DOCKING — mode={mode}")
    print(f"{'='*60}")

    script_path = MODULES_DIR / "04_docking" / "run_docking.py"
    if not script_path.exists():
        print(f"Error: Docking module not found at {script_path}")
        return False

    cmd = get_python_cmd("base")
    cmd += [
        str(script_path),
        "--manifest", str(manifest_path),
        "--candidates_dir", str(DATA_DIR / "screened"),
        "--output_dir", str(DATA_DIR / "results"),
        "--db_path", db_path,
        "--campaign_id", campaign_id,
        "--mode", mode,
    ]

    try:
        subprocess.check_call(cmd)
        print("[Orchestrator] Docking complete.\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[Orchestrator] ERROR: Docking failed with code {e.returncode}")
        return False


def print_campaign_summary(db_path, campaign_id):
    """Print a summary of the campaign from telemetry."""
    try:
        db = TelemetryDB(db_path)
        summary = db.get_campaign_summary(campaign_id)
        db.close()

        print(f"\n{'='*60}")
        print(f"[Orchestrator] CAMPAIGN SUMMARY")
        print(f"{'='*60}")
        print(f"  Campaign ID:     {summary['campaign_id']}")
        print(f"  Runs by status:  {summary['runs_by_status']}")
        print(f"  Total molecules: {summary['total_molecules']}")
        if summary['triage_stats']:
            passed = summary['triage_stats'].get('1', 0)
            failed = summary['triage_stats'].get('0', 0)
            print(f"  Triage passed:   {passed}")
            print(f"  Triage rejected: {failed}")
        print(f"  Telemetry DB:    {db_path}")
        print(f"{'='*60}\n")
    except Exception as e:
        print(f"[Orchestrator] Warning: Could not print campaign summary: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Adaptive Discovery Orchestrator — Autonomous Drug Discovery Pipeline"
    )
    subparsers = parser.add_subparsers(dest="command", help="Pipeline command")

    # Global options
    parser.add_argument("--db_path", default=DEFAULT_DB_PATH,
                        help="Path to telemetry database")

    # Ingest command
    ingest_parser = subparsers.add_parser("ingest", help="Ingest a PDB file and identify pockets")
    ingest_parser.add_argument("pdb_file", help="Path to input PDB file")
    ingest_parser.add_argument("--clean", action="store_true", help="Clean PDB before processing")

    # Generate command
    gen_parser = subparsers.add_parser("generate", help="Generate molecules from a manifest")
    gen_parser.add_argument("manifest", help="Path to ingestion manifest.json")
    gen_parser.add_argument("--mode", choices=["simulation", "production"],
                            default="simulation", help="Execution mode")

    # Screen command
    screen_parser = subparsers.add_parser("screen", help="Screen molecules through fast triage")
    screen_parser.add_argument("input_sdf", help="Path to SDF file to screen")

    # Dock command
    dock_parser = subparsers.add_parser("dock", help="Dock generated molecules")
    dock_parser.add_argument("manifest", help="Path to ingestion manifest.json")
    dock_parser.add_argument("--mode", choices=["simulation", "production"],
                            default="simulation", help="Execution mode")

    # Full Pipeline command
    pipeline_parser = subparsers.add_parser("run", help="Run full pipeline")
    pipeline_parser.add_argument("pdb_file", help="Path to input PDB file")
    pipeline_parser.add_argument("--mode", choices=["simulation", "production"],
                                 default="simulation", help="Execution mode")

    args = parser.parse_args()

    ensure_dirs()

    # Generate a new campaign ID for tracking
    campaign_id = f"campaign_{uuid.uuid4().hex[:8]}"
    db_path = args.db_path

    print(f"[Orchestrator] Campaign ID: {campaign_id}")
    print(f"[Orchestrator] Telemetry DB: {db_path}")

    if args.command == "ingest":
        run_ingestion(args.pdb_file, db_path, campaign_id, args.clean)

    elif args.command == "generate":
        mode = getattr(args, "mode", "simulation")
        run_generation(args.manifest, db_path, campaign_id, mode)

    elif args.command == "screen":
        run_screening(args.input_sdf, db_path, campaign_id)

    elif args.command == "dock":
        mode = getattr(args, "mode", "simulation")
        run_docking(args.manifest, db_path, campaign_id, mode)

    elif args.command == "run":
        mode = getattr(args, "mode", "simulation")
        pdb_path = Path(args.pdb_file)

        # Map pipeline mode to per-stage modes
        if mode == "simulation":
            gen_mode, dock_mode = "simulation", "simulation"
        else:
            gen_mode, dock_mode = "rdkit", "production"

        print(f"[Orchestrator] Running FULL PIPELINE — gen={gen_mode}, dock={dock_mode}")

        if run_ingestion(pdb_path, db_path, campaign_id):
            manifest_path = DATA_DIR / "processed" / f"{pdb_path.stem}_manifest.json"
            if run_generation(manifest_path, db_path, campaign_id, gen_mode):
                sdf_path = DATA_DIR / "candidates" / "generated_molecules.sdf"
                if run_screening(sdf_path, db_path, campaign_id):
                    screened_sdf = DATA_DIR / "screened" / "screened_molecules.sdf"
                    run_docking(manifest_path, db_path, campaign_id, dock_mode)

        print_campaign_summary(db_path, campaign_id)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
