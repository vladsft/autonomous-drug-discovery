"""
Orchestrator: The central CLI driver for the Adaptive Discovery Pipeline.

Manages data flow between isolated modules, generates campaign IDs,
and ensures all runs are captured in the telemetry database.

Usage:
    python orchestrator.py ingest target.pdb
    python orchestrator.py generate manifest.json
    python orchestrator.py screen candidates.sdf
    python orchestrator.py dock manifest.json
    python orchestrator.py rank docking_results.csv [--screening_json screening_report.json]
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
    """Ensure data directories exist.

    `processed/` holds ingestion outputs (manifests + pocket PDBs) that are
    keyed by PDB stem and shared across campaigns. `candidates/`, `screened/`,
    and `results/` are fallback shared directories for single-stage invocations;
    full pipeline runs use per-campaign subdirectories instead.
    """
    for subdir in ["processed", "candidates", "screened", "results"]:
        (DATA_DIR / subdir).mkdir(parents=True, exist_ok=True)


def get_python_cmd(env_name):
    """Return the command to run python in a specific conda environment."""
    conda_bin = os.environ.get("CONDA_EXE", "conda")
    return [conda_bin, "run", "-n", env_name, "python"]


def run_ingestion(pdb_path, db_path, campaign_id, clean_pdb=False, pocket_backend="p2rank"):
    """Run the ingestion module (P2Rank or fpocket)."""
    print(f"\n{'='*60}")
    print(f"[Orchestrator] Stage 1: INGESTION — {pdb_path} (backend={pocket_backend})")
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
        "--backend", pocket_backend,
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


def run_generation(manifest_path, db_path, campaign_id, mode="simulation",
                   output_dir=None, num_samples=None):
    """Run the generation module.

    Dispatches to the appropriate conda env based on `mode`:
    RDKit/simulation stay in `base`, `targetdiff` → `targetdiff_env`,
    `pocket2mol` → `pocket2mol_env`.

    `num_samples`, when set, overrides the per-mode default campaign size.
    """
    print(f"\n{'='*60}")
    print(f"[Orchestrator] Stage 2: GENERATION — mode={mode}")
    print(f"{'='*60}")

    script_path = MODULES_DIR / "02_generation" / "run_generation.py"
    if not script_path.exists():
        print(f"Error: Generation module not found at {script_path}")
        return False

    out_dir = output_dir or str(DATA_DIR / "candidates")

    env_map = {"targetdiff": "targetdiff_env", "pocket2mol": "pocket2mol_env"}
    cmd = get_python_cmd(env_map.get(mode, "base"))

    cmd += [
        str(script_path),
        "--manifest", str(manifest_path),
        "--output_dir", out_dir,
        "--db_path", db_path,
        "--campaign_id", campaign_id,
        "--mode", mode,
    ]
    if num_samples is not None:
        cmd += ["--num_samples", str(num_samples)]

    try:
        subprocess.check_call(cmd)
        print("[Orchestrator] Generation complete.\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[Orchestrator] ERROR: Generation failed with code {e.returncode}")
        return False


def run_screening(sdf_path, db_path, campaign_id, output_dir=None):
    """Run the screening module (MolScore descriptors + ADMET-AI enrichment)."""
    print(f"\n{'='*60}")
    print(f"[Orchestrator] Stage 3: SCREENING (Fast Triage)")
    print(f"{'='*60}")

    script_path = MODULES_DIR / "03_screening" / "run_screening.py"
    if not script_path.exists():
        print(f"Error: Screening module not found at {script_path}")
        return False

    out_dir = output_dir or str(DATA_DIR / "screened")

    cmd = get_python_cmd("base") + [
        str(script_path),
        "--input_sdf", str(sdf_path),
        "--output_dir", out_dir,
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


def run_docking(manifest_path, db_path, campaign_id, mode="simulation",
                candidates_dir=None, output_dir=None):
    """Run the docking module (Vina wrapper)."""
    print(f"\n{'='*60}")
    print(f"[Orchestrator] Stage 4: DOCKING — mode={mode}")
    print(f"{'='*60}")

    script_path = MODULES_DIR / "04_docking" / "run_docking.py"
    if not script_path.exists():
        print(f"Error: Docking module not found at {script_path}")
        return False

    cand_dir = candidates_dir or str(DATA_DIR / "screened")
    out_dir = output_dir or str(DATA_DIR / "results")

    cmd = get_python_cmd("base")
    cmd += [
        str(script_path),
        "--manifest", str(manifest_path),
        "--candidates_dir", cand_dir,
        "--output_dir", out_dir,
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


def run_ranking(docking_csv, screening_json, db_path, campaign_id, output_dir=None):
    """Run the ranking module (multi-criteria final ranker)."""
    print(f"\n{'='*60}")
    print(f"[Orchestrator] Stage 5: RANKING")
    print(f"{'='*60}")

    script_path = MODULES_DIR / "05_ranking" / "run_ranking.py"
    if not script_path.exists():
        print(f"Error: Ranking module not found at {script_path}")
        return False

    out_dir = output_dir or str(DATA_DIR / "results")

    cmd = get_python_cmd("base") + [
        str(script_path),
        "--docking_csv", str(docking_csv),
        "--output_dir", out_dir,
        "--db_path", db_path,
        "--campaign_id", campaign_id,
    ]
    if screening_json:
        cmd += ["--screening_json", str(screening_json)]

    try:
        subprocess.check_call(cmd)
        print("[Orchestrator] Ranking complete.\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[Orchestrator] ERROR: Ranking failed with code {e.returncode}")
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
    ingest_parser.add_argument("--backend", choices=["p2rank", "fpocket"], default="p2rank",
                               help="Pocket detection backend (default: p2rank)")
    ingest_parser.add_argument("--clean", action="store_true",
                               help="Drop HETATM/waters/alt locs before detection")

    # Generate command
    gen_parser = subparsers.add_parser("generate", help="Generate molecules from a manifest")
    gen_parser.add_argument("manifest", help="Path to ingestion manifest.json")
    gen_parser.add_argument("--mode", choices=["simulation", "rdkit", "targetdiff", "pocket2mol", "production"],
                            default="simulation", help="Execution mode")
    gen_parser.add_argument("--num_samples", type=int, default=None,
                            help="Number of molecules to generate (overrides the per-mode default)")

    # Screen command
    screen_parser = subparsers.add_parser("screen", help="Screen molecules through fast triage")
    screen_parser.add_argument("input_sdf", help="Path to SDF file to screen")

    # Dock command
    dock_parser = subparsers.add_parser("dock", help="Dock generated molecules")
    dock_parser.add_argument("manifest", help="Path to ingestion manifest.json")
    dock_parser.add_argument("--mode", choices=["simulation", "triage", "production"],
                            default="simulation", help="Execution mode")

    # Rank command (Stage 5: multi-criteria ranker)
    rank_parser = subparsers.add_parser("rank", help="Multi-criteria ranking of docked candidates")
    rank_parser.add_argument("docking_csv", help="Path to docking_results.csv (Stage 4 output)")
    rank_parser.add_argument("--screening_json", default=None,
                             help="Optional screening_report.json for ADMET enrichment")

    # Full Pipeline command
    pipeline_parser = subparsers.add_parser("run", help="Run full pipeline")
    pipeline_parser.add_argument("pdb_file", help="Path to input PDB file")
    pipeline_parser.add_argument("--mode", choices=["simulation", "production", "targetdiff", "pocket2mol"],
                                 default="simulation", help="Execution mode")
    pipeline_parser.add_argument("--backend", choices=["p2rank", "fpocket"], default="p2rank",
                                 help="Pocket detection backend (default: p2rank)")
    pipeline_parser.add_argument("--num_samples", type=int, default=None,
                                 help="Number of molecules to generate (overrides the per-mode default)")
    pipeline_parser.add_argument("--clean", action="store_true",
                                 help="Drop HETATM/waters/alt locs before detection")

    args = parser.parse_args()

    ensure_dirs()

    # Generate a new campaign ID for tracking
    campaign_id = f"campaign_{uuid.uuid4().hex[:8]}"
    db_path = args.db_path

    print(f"[Orchestrator] Campaign ID: {campaign_id}")
    print(f"[Orchestrator] Telemetry DB: {db_path}")

    ok = True
    if args.command == "ingest":
        ok = run_ingestion(args.pdb_file, db_path, campaign_id,
                           clean_pdb=args.clean, pocket_backend=args.backend)

    elif args.command == "generate":
        mode = getattr(args, "mode", "simulation")
        ok = run_generation(args.manifest, db_path, campaign_id, mode,
                            num_samples=args.num_samples)

    elif args.command == "screen":
        ok = run_screening(args.input_sdf, db_path, campaign_id)

    elif args.command == "dock":
        mode = getattr(args, "mode", "simulation")
        ok = run_docking(args.manifest, db_path, campaign_id, mode)

    elif args.command == "rank":
        ok = run_ranking(args.docking_csv, args.screening_json,
                         db_path, campaign_id)

    elif args.command == "run":
        mode = getattr(args, "mode", "simulation")
        pdb_path = Path(args.pdb_file)

        # Map pipeline mode to per-stage modes
        if mode == "simulation":
            gen_mode, dock_mode = "simulation", "simulation"
        elif mode == "targetdiff":
            gen_mode, dock_mode = "targetdiff", "production"
        elif mode == "pocket2mol":
            gen_mode, dock_mode = "pocket2mol", "production"
        else:
            gen_mode, dock_mode = "rdkit", "production"

        # Per-campaign output directories prevent file collisions between runs
        campaign_dir = DATA_DIR / campaign_id
        gen_dir = str(campaign_dir / "candidates")
        screen_dir = str(campaign_dir / "screened")
        dock_dir = str(campaign_dir / "results")
        rank_dir = str(campaign_dir / "ranked")

        print(f"[Orchestrator] Running FULL PIPELINE — gen={gen_mode}, dock={dock_mode}")
        print(f"[Orchestrator] Campaign directory: {campaign_dir}")

        ok = run_ingestion(pdb_path, db_path, campaign_id,
                           clean_pdb=args.clean, pocket_backend=args.backend)
        if ok:
            manifest_path = DATA_DIR / "processed" / f"{pdb_path.stem}_manifest.json"
            ok = run_generation(manifest_path, db_path, campaign_id, gen_mode,
                                output_dir=gen_dir, num_samples=args.num_samples)
        if ok:
            sdf_path = Path(gen_dir) / "generated_molecules.sdf"
            ok = run_screening(sdf_path, db_path, campaign_id,
                               output_dir=screen_dir)
        if ok:
            ok = run_docking(manifest_path, db_path, campaign_id, dock_mode,
                             candidates_dir=screen_dir, output_dir=dock_dir)
        if ok:
            docking_csv = Path(dock_dir) / "docking_results.csv"
            screening_json = Path(screen_dir) / "screening_report.json"
            ok = run_ranking(docking_csv,
                             screening_json if screening_json.exists() else None,
                             db_path, campaign_id, output_dir=rank_dir)

        print_campaign_summary(db_path, campaign_id)
    else:
        parser.print_help()
        return 2

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main() or 0)
