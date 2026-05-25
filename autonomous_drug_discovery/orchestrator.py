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
    """Return the command to run python in a specific conda environment.

    If we are ALREADY running inside `env_name`, reuse this very interpreter
    (`sys.executable`) instead of spawning a nested `conda run -n env_name`.

    Why this matters: the container entrypoint launches the orchestrator with
    `conda run -n base python orchestrator.py`. A stage that also runs in `base`
    (docking, screening, ranking) would otherwise be spawned as a SECOND,
    nested `conda run -n base` — and that nested activation does not reliably
    restore the dynamic-loader paths that conda packages with compiled
    extensions need. `vina` links Boost `.so`s out of `$CONDA_PREFIX/lib`, so
    `import vina` then fails at runtime even though a *single* `conda run -n
    base` imports it fine (which the Docker build check proves). `sys.executable`
    already carries the parent's fully-activated environment, so the import
    works. Cross-env hops (e.g. base -> targetdiff_env) still go through
    `conda run`, which is correct.
    """
    if os.environ.get("CONDA_DEFAULT_ENV") == env_name:
        return [sys.executable]
    conda_bin = os.environ.get("CONDA_EXE", "conda")
    return [conda_bin, "run", "-n", env_name, "python"]


def stage_env(env_name):
    """Environment for a stage subprocess.

    When a stage runs in the env we're ALREADY in (base — see get_python_cmd),
    prepend that env's lib dir to LD_LIBRARY_PATH so the dynamic loader resolves
    `libstdc++.so.6` (and other compiled deps) from conda, not the older system
    copy under /lib/x86_64-linux-gnu.

    Why: a stage like docking imports `tdc`/`scipy` *before* `vina`. Those map
    the system libstdc++ into the process first; vina's compiled extension then
    binds to that already-loaded lib, which lacks the `CXXABI_1.3.15` symbol its
    conda build needs, and `import vina` aborts. LD_LIBRARY_PATH is read by the
    loader at exec time, so it must be set on the child's environment here —
    setting it from inside the child is too late once the system lib is mapped.

    This is deliberately NOT applied to cross-env launches (e.g. base ->
    targetdiff_env via `conda run`): prepending base's lib there could shadow
    targetdiff_env's own CUDA/torch libraries. The CONDA_DEFAULT_ENV guard
    confines the override to same-env (base) stages, which is exactly where the
    sys.executable launch path runs.
    """
    env = os.environ.copy()
    if os.environ.get("CONDA_DEFAULT_ENV") == env_name:
        lib = os.path.join(sys.prefix, "lib")
        prev = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = lib + (os.pathsep + prev if prev else "")
    return env


def check_stage_output(path, description, min_bytes=1):
    """Tripwire: verify a stage actually produced a non-trivial output file.

    A module can exit 0 yet leave an empty or missing artefact (a broken
    receptor once produced a docking CSV of all-zero scores). Returns True if
    `path` exists and is at least `min_bytes`; otherwise prints why and returns
    False so the pipeline stops instead of feeding garbage to the next stage.
    """
    p = Path(path)
    if not p.exists():
        print(f"[Orchestrator] TRIPWIRE: {description} missing — expected {p}")
        return False
    if p.stat().st_size < min_bytes:
        print(f"[Orchestrator] TRIPWIRE: {description} is empty — {p}")
        return False
    return True


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
        subprocess.check_call(cmd, env=stage_env("base"))
        print("[Orchestrator] Ingestion complete.\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[Orchestrator] ERROR: Ingestion failed with code {e.returncode}")
        return False


def run_generation(manifest_path, db_path, campaign_id, mode="simulation",
                   output_dir=None, num_samples=None, device="auto"):
    """Run the generation module.

    Dispatches to the appropriate conda env based on `mode`:
    RDKit/simulation stay in `base`, `targetdiff` → `targetdiff_env`,
    `pocket2mol` → `pocket2mol_env`.

    `num_samples`, when set, overrides the per-mode default campaign size.
    `device` ("auto"/"cuda"/"cpu") selects the compute device for the
    GPU-capable backends; "auto" detects a GPU.
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
    gen_env = env_map.get(mode, "base")
    cmd = get_python_cmd(gen_env)

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
    if device:
        cmd += ["--device", str(device)]

    try:
        subprocess.check_call(cmd, env=stage_env(gen_env))
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
        subprocess.check_call(cmd, env=stage_env("base"))
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
        subprocess.check_call(cmd, env=stage_env("base"))
        print("[Orchestrator] Docking complete.\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[Orchestrator] ERROR: Docking failed with code {e.returncode}")
        return False


def run_ranking(docking_csv, screening_json, db_path, campaign_id, output_dir=None,
                aizynth_config=None):
    """Run the ranking module (multi-criteria final ranker).

    `aizynth_config`, when set, points the ranker at an AiZynthFinder config so
    the top-N candidates get retrosynthetic feasibility scores; otherwise the
    synthesis term stays neutral.
    """
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
    if aizynth_config:
        cmd += ["--aizynth_config", str(aizynth_config)]

    try:
        subprocess.check_call(cmd, env=stage_env("base"))
        print("[Orchestrator] Ranking complete.\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[Orchestrator] ERROR: Ranking failed with code {e.returncode}")
        return False


def run_validation(holo_pdb, db_path, campaign_id, ligand_resname,
                   ligand_smiles=None, exhaustiveness=16, decoys=None):
    """Run retrospective docking validation (re-dock a known ligand)."""
    print(f"\n{'='*60}")
    print(f"[Orchestrator] DOCKING VALIDATION — {holo_pdb} (ligand={ligand_resname})")
    print(f"{'='*60}")

    script_path = MODULES_DIR / "04_docking" / "validate_docking.py"
    if not script_path.exists():
        print(f"Error: Validation module not found at {script_path}")
        return False

    cmd = get_python_cmd("base") + [
        str(script_path),
        "--holo_pdb", str(holo_pdb),
        "--ligand_resname", ligand_resname,
        "--output_dir", str(DATA_DIR / "validation"),
        "--db_path", db_path,
        "--campaign_id", campaign_id,
        "--exhaustiveness", str(exhaustiveness),
    ]
    if ligand_smiles:
        cmd += ["--ligand_smiles", ligand_smiles]
    if decoys:
        cmd += ["--decoys", str(decoys)]

    try:
        subprocess.check_call(cmd, env=stage_env("base"))
        print("[Orchestrator] Validation PASSED (pose RMSD within threshold).\n")
        return True
    except subprocess.CalledProcessError as e:
        # validate_docking.py exits 2 specifically when re-docking misses the
        # RMSD bar — that's a meaningful "docking setup is not trustworthy"
        # result, not a crash.
        if e.returncode == 2:
            print("[Orchestrator] Validation FAILED: re-docking did not "
                  "reproduce the crystal pose within the RMSD threshold.")
        else:
            print(f"[Orchestrator] ERROR: Validation failed with code {e.returncode}")
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
    gen_parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto",
                            help="Compute device for GPU-capable backends "
                                 "(targetdiff/pocket2mol); 'auto' detects a GPU")

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
    rank_parser.add_argument("--aizynth_config", default=os.environ.get("AIZYNTH_CONFIG"),
                             help="AiZynthFinder config.yml — enables retrosynthetic "
                                  "scoring of the top-N candidates (default: $AIZYNTH_CONFIG)")

    # Validate command (retrospective docking validation)
    validate_parser = subparsers.add_parser(
        "validate", help="Retrospective docking validation — re-dock a known ligand")
    validate_parser.add_argument("holo_pdb", help="Holo PDB (protein + co-crystal ligand)")
    validate_parser.add_argument("--ligand_resname", required=True,
                                 help="3-letter HETATM residue name of the ligand")
    validate_parser.add_argument("--ligand_smiles", default=None,
                                 help="Ligand SMILES — assigns correct bond orders "
                                      "to the crystal pose (recommended)")
    validate_parser.add_argument("--exhaustiveness", type=int, default=16,
                                 help="Vina exhaustiveness for the re-docking")
    validate_parser.add_argument("--decoys", default=None,
                                 help="Optional decoy SMILES file for an "
                                      "enrichment check")

    # Full Pipeline command
    pipeline_parser = subparsers.add_parser("run", help="Run full pipeline")
    pipeline_parser.add_argument("pdb_file", help="Path to input PDB file")
    pipeline_parser.add_argument("--mode",
                                 choices=["simulation", "production", "rdkit", "targetdiff", "pocket2mol"],
                                 default="simulation",
                                 help="Execution mode (`rdkit` and `production` are aliases)")
    pipeline_parser.add_argument("--backend", choices=["p2rank", "fpocket"], default="p2rank",
                                 help="Pocket detection backend (default: p2rank)")
    pipeline_parser.add_argument("--num_samples", type=int, default=None,
                                 help="Number of molecules to generate (overrides the per-mode default)")
    pipeline_parser.add_argument("--clean", action="store_true",
                                 help="Drop HETATM/waters/alt locs before detection")
    pipeline_parser.add_argument("--aizynth_config", default=os.environ.get("AIZYNTH_CONFIG"),
                                 help="AiZynthFinder config.yml — enables retrosynthetic "
                                      "scoring in Stage 5 (default: $AIZYNTH_CONFIG)")
    pipeline_parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto",
                                 help="Compute device for GPU-capable generation "
                                      "(targetdiff/pocket2mol); 'auto' detects a GPU")

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
                            num_samples=args.num_samples, device=args.device)

    elif args.command == "screen":
        ok = run_screening(args.input_sdf, db_path, campaign_id)

    elif args.command == "dock":
        mode = getattr(args, "mode", "simulation")
        ok = run_docking(args.manifest, db_path, campaign_id, mode)

    elif args.command == "rank":
        ok = run_ranking(args.docking_csv, args.screening_json,
                         db_path, campaign_id,
                         aizynth_config=args.aizynth_config)

    elif args.command == "validate":
        ok = run_validation(args.holo_pdb, db_path, campaign_id,
                            ligand_resname=args.ligand_resname,
                            ligand_smiles=args.ligand_smiles,
                            exhaustiveness=args.exhaustiveness,
                            decoys=args.decoys)

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

        manifest_path = DATA_DIR / "processed" / f"{pdb_path.stem}_manifest.json"
        sdf_path = Path(gen_dir) / "generated_molecules.sdf"
        screened_sdf = Path(screen_dir) / "screened_molecules.sdf"
        docking_csv = Path(dock_dir) / "docking_results.csv"

        ok = run_ingestion(pdb_path, db_path, campaign_id,
                           clean_pdb=args.clean, pocket_backend=args.backend)
        if ok:
            ok = check_stage_output(manifest_path, "ingestion manifest")
        if ok:
            ok = run_generation(manifest_path, db_path, campaign_id, gen_mode,
                                output_dir=gen_dir, num_samples=args.num_samples,
                                device=args.device)
        if ok:
            ok = check_stage_output(sdf_path, "generated molecules SDF", min_bytes=64)
        if ok:
            ok = run_screening(sdf_path, db_path, campaign_id,
                               output_dir=screen_dir)
        if ok:
            ok = check_stage_output(screened_sdf, "screened molecules SDF")
        if ok:
            ok = run_docking(manifest_path, db_path, campaign_id, dock_mode,
                             candidates_dir=screen_dir, output_dir=dock_dir)
        if ok:
            ok = check_stage_output(docking_csv, "docking results CSV", min_bytes=64)
        if ok:
            screening_json = Path(screen_dir) / "screening_report.json"
            ok = run_ranking(docking_csv,
                             screening_json if screening_json.exists() else None,
                             db_path, campaign_id, output_dir=rank_dir,
                             aizynth_config=args.aizynth_config)

        print_campaign_summary(db_path, campaign_id)
    else:
        parser.print_help()
        return 2

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main() or 0)
