"""
Module 02: Generation — Structure-Aware Molecule Generator.

Three generation backends:
  - simulation: Stub SDF for pipeline testing.
  - rdkit: Fragment-based combinatorial generation using RDKit. Produces
    diverse drug-like molecules with 3D conformers, optionally filtered by
    pocket volume. No GPU required.
  - targetdiff: TargetDiff diffusion model (requires separate conda env + checkpoint).

Input contract:  manifest.json (from ingestion module)
Output contract: generated_molecules.sdf + run_metadata.json
"""

import os
import sys
import argparse
import subprocess
import json
import hashlib
import random
import traceback
from datetime import datetime, timezone
from pathlib import Path

# Paths
MODULE_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = MODULE_DIR.parent.parent
TARGETDIFF_REPO = MODULE_DIR / "targetdiff"

sys.path.insert(0, str(PROJECT_ROOT))
from telemetry import TelemetryDB

# Check for RDKit
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdFMCS
    from rdkit.Chem import BRICS
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


def _get_git_commit(repo_path):
    """Get the current git commit hash of a repository."""
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_path), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        return result.stdout.strip() if result.returncode == 0 else None
    except Exception:
        return None


# Default generation parameters
DEFAULT_PARAMS = {
    "num_samples": 100,
    "sampling_steps": 1000,
    "noise_schedule": "polynomial_2",
    "batch_size": 16,
    "device": "cpu",
}


# ---------------------------------------------------------------------------
# Fragment library for RDKit-based generation
# ---------------------------------------------------------------------------

# Drug-like scaffolds — common cores found in approved drugs
SCAFFOLDS = [
    "c1ccc2[nH]ccc2c1",          # indole
    "c1ccc2nc[nH]c2c1",          # benzimidazole
    "c1cnc2ccccc2n1",             # quinazoline
    "c1ccc2ncccc2c1",             # quinoline
    "c1ccncc1",                   # pyridine
    "c1cnc[nH]1",                 # imidazole
    "c1ccoc1",                    # furan
    "c1ccsc1",                    # thiophene
    "c1cc[nH]c1",                 # pyrrole
    "c1ccc(cc1)c1ccccc1",        # biphenyl
    "c1ccc2c(c1)cccc2",          # naphthalene
    "C1CCCCC1",                   # cyclohexane
    "C1CCNCC1",                   # piperidine
    "C1CNCCN1",                   # piperazine
    "C1CCOCC1",                   # morpholine (tetrahydro-1,4-oxazine)
    "c1ccc(cc1)O",               # phenol
    "c1ccc(cc1)N",               # aniline
    "C1CC1",                      # cyclopropane
]

# Substituents — small functional groups to attach at open positions
SUBSTITUENTS = [
    "O",           # hydroxyl
    "N",           # amine
    "C(=O)O",     # carboxylic acid
    "C(=O)N",     # amide
    "F",           # fluorine
    "Cl",          # chlorine
    "C",           # methyl
    "CC",          # ethyl
    "OC",          # methoxy
    "C#N",         # nitrile
    "C(F)(F)F",   # trifluoromethyl
    "S(=O)(=O)N", # sulfonamide
    "NC(=O)",     # reverse amide
    "OCC",         # ethoxy
    "C(=O)",      # carbonyl
]

# Linkers — connect scaffolds to make larger molecules
LINKERS = [
    "",            # direct bond
    "C",           # methylene
    "CC",          # ethylene
    "NC(=O)",     # amide linker
    "C(=O)N",     # reverse amide
    "O",           # ether
    "NC",          # aminomethyl
    "OC",          # oxymethyl
]


def run_generation_simulation(manifest, out_path, parameters):
    """Simulation mode: generate a stub SDF for pipeline testing."""
    print("[Generation] (SIMULATION MODE) Generating stub molecules...")

    simulated_mol = out_path / "generated_molecules.sdf"
    out_path.mkdir(parents=True, exist_ok=True)

    sdf_content = """
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124   -0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2124   -0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
>  <molecule_id>
mol_0000

$$$$
"""
    with open(simulated_mol, "w") as f:
        f.write(sdf_content)

    print(f"[Generation] Stub SDF written to {simulated_mol}")
    return str(simulated_mol)


# ---------------------------------------------------------------------------
# RDKit Fragment-Based Generation
# ---------------------------------------------------------------------------

def _pocket_volume_radius(pocket_pdb_path):
    """Estimate the pocket radius from atom coordinates."""
    try:
        coords = []
        with open(pocket_pdb_path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    coords.append((x, y, z))
        if not coords:
            return None
        n = len(coords)
        cx = sum(c[0] for c in coords) / n
        cy = sum(c[1] for c in coords) / n
        cz = sum(c[2] for c in coords) / n
        max_dist = max(
            ((x - cx)**2 + (y - cy)**2 + (z - cz)**2)**0.5
            for x, y, z in coords
        )
        return max_dist
    except Exception:
        return None


def _estimate_heavy_atom_range(pocket_pdb_path):
    """Estimate a reasonable heavy-atom count range based on pocket size.

    Rule of thumb: ~10 Angstrom radius pocket fits ~15-25 heavy atoms.
    """
    radius = _pocket_volume_radius(pocket_pdb_path)
    if radius is None:
        return 10, 35  # default range

    # Rough empirical: 1.5 heavy atoms per Angstrom of radius
    center = int(radius * 1.5)
    lo = max(8, center - 8)
    hi = min(50, center + 8)
    return lo, hi


def _assemble_molecule(scaffolds, substituents, linkers, rng):
    """Assemble a random molecule from fragments.

    Strategy: pick 1-2 scaffolds, connect with a linker, decorate with
    1-3 substituents. This mimics medicinal chemistry fragment assembly.
    """
    # Pick primary scaffold
    scaffold_smi = rng.choice(scaffolds)

    # Optionally attach a second scaffold via a linker (30% chance)
    if rng.random() < 0.3:
        scaffold2 = rng.choice(scaffolds)
        linker = rng.choice(linkers)
        # Use SMILES concatenation — RDKit will try to parse it
        combined = f"{scaffold_smi}{linker}{scaffold2}"
    else:
        combined = scaffold_smi

    mol = Chem.MolFromSmiles(combined)
    if mol is None:
        return None

    # Try BRICS decomposition and reassembly for more diversity
    try:
        frags = list(BRICS.BRICSDecompose(mol))
        if len(frags) >= 2:
            rebuilt = list(BRICS.BRICSBuild(
                [Chem.MolFromSmiles(f) for f in frags if Chem.MolFromSmiles(f) is not None]
            ))
            if rebuilt:
                candidate = rng.choice(rebuilt[:10])  # pick from first few
                candidate = Chem.MolFromSmiles(Chem.MolToSmiles(candidate))
                if candidate is not None:
                    mol = candidate
    except Exception:
        pass  # BRICS can fail on some scaffolds, that's fine

    # Add 1-3 substituents by replacing Hs on the scaffold
    n_subs = rng.randint(1, 3)
    current_smi = Chem.MolToSmiles(mol)
    for _ in range(n_subs):
        sub = rng.choice(substituents)
        # Simple attachment: wrap scaffold + substituent
        trial_smi = f"{current_smi}.{sub}"
        trial = Chem.MolFromSmiles(trial_smi)
        if trial is not None:
            # Try to create a bond between the fragments
            try:
                combo = AllChem.CombineMols(mol, Chem.MolFromSmiles(sub))
                editable = Chem.RWMol(combo)
                # Find atoms that can accept a bond
                n_atoms_main = mol.GetNumAtoms()
                # Pick a random atom from main mol and sub mol
                main_idx = rng.randint(0, n_atoms_main - 1)
                sub_idx = rng.randint(n_atoms_main, editable.GetNumAtoms() - 1)
                editable.AddBond(main_idx, sub_idx, Chem.BondType.SINGLE)
                new_mol = editable.GetMol()
                try:
                    Chem.SanitizeMol(new_mol)
                    mol = new_mol
                    current_smi = Chem.MolToSmiles(mol)
                except Exception:
                    pass  # invalid bond, skip this substituent
            except Exception:
                pass

    return mol


def _is_drug_like_quick(mol, min_ha, max_ha):
    """Fast pre-filter before full screening. Rejects obvious junk."""
    ha = mol.GetNumHeavyAtoms()
    if ha < min_ha or ha > max_ha:
        return False
    mw = Descriptors.MolWt(mol)
    if mw < 100 or mw > 600:
        return False
    logp = Descriptors.MolLogP(mol)
    if logp < -3 or logp > 6:
        return False
    # Reject molecules with too many rings (likely polymeric junk)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 6:
        return False
    return True


def run_generation_rdkit(manifest, out_path, parameters):
    """RDKit fragment-based generation: produces diverse drug-like molecules."""
    if not HAS_RDKIT:
        raise ImportError("RDKit is required for rdkit generation mode.")

    RDLogger.DisableLog("rdApp.*")

    pocket_pdb = manifest.get("best_pocket")
    num_samples = parameters.get("num_samples", 100)
    seed = parameters.get("seed", 42)
    rng = random.Random(seed)

    # Estimate pocket-appropriate molecule size
    if pocket_pdb and Path(pocket_pdb).exists():
        min_ha, max_ha = _estimate_heavy_atom_range(pocket_pdb)
        radius = _pocket_volume_radius(pocket_pdb)
        print(f"[Generation] Pocket radius: {radius:.1f} A → targeting {min_ha}-{max_ha} heavy atoms")
    else:
        min_ha, max_ha = 10, 35
        print(f"[Generation] No pocket info, using default range: {min_ha}-{max_ha} heavy atoms")

    out_path.mkdir(parents=True, exist_ok=True)
    output_sdf = out_path / "generated_molecules.sdf"

    # Generate molecules with deduplication
    seen_smiles = set()
    molecules = []
    attempts = 0
    max_attempts = num_samples * 20  # give up after this many tries

    print(f"[Generation] Generating {num_samples} drug-like molecules...")

    while len(molecules) < num_samples and attempts < max_attempts:
        attempts += 1
        mol = _assemble_molecule(SCAFFOLDS, SUBSTITUENTS, LINKERS, rng)
        if mol is None:
            continue

        # Canonicalize and deduplicate
        try:
            smi = Chem.MolToSmiles(mol)
        except Exception:
            continue
        if smi in seen_smiles:
            continue

        # Quick drug-likeness filter
        if not _is_drug_like_quick(mol, min_ha, max_ha):
            continue

        # Generate 3D conformer
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, randomSeed=seed + attempts)
        if result != 0:
            # Retry with different params
            result = AllChem.EmbedMolecule(
                mol, AllChem.ETKDGv3(), randomSeed=seed + attempts
            )
            if result != 0:
                continue

        # Optimize geometry
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass

        mol = Chem.RemoveHs(mol)
        mol.SetProp("molecule_id", f"mol_{len(molecules):04d}")
        mol.SetProp("_Name", f"mol_{len(molecules):04d}")
        mol.SetProp("smiles", smi)

        seen_smiles.add(smi)
        molecules.append(mol)

        if len(molecules) % 25 == 0:
            print(f"[Generation]   {len(molecules)}/{num_samples} molecules generated ({attempts} attempts)")

    # Write SDF
    writer = Chem.SDWriter(str(output_sdf))
    for mol in molecules:
        writer.write(mol)
    writer.close()

    print(f"[Generation] {len(molecules)} molecules generated in {attempts} attempts")
    print(f"[Generation] Output: {output_sdf}")
    return str(output_sdf)


def run_generation_targetdiff(manifest, out_path, parameters):
    """TargetDiff mode: execute TargetDiff diffusion model inference."""
    if not TARGETDIFF_REPO.exists():
        raise FileNotFoundError(
            f"TargetDiff repository not found at {TARGETDIFF_REPO}. "
            f"Clone it: git clone https://github.com/guanjq/targetdiff.git {TARGETDIFF_REPO}"
        )

    pocket_pdb = Path(manifest.get("best_pocket"))
    if not pocket_pdb.exists():
        raise FileNotFoundError(f"Best pocket file {pocket_pdb} not found.")

    print(f"[Generation] (TARGETDIFF MODE) Running TargetDiff on {pocket_pdb}...")

    sample_script = TARGETDIFF_REPO / "scripts" / "sample_for_pocket.py"
    if not sample_script.exists():
        raise FileNotFoundError(f"TargetDiff sample script not found at {sample_script}")

    out_path.mkdir(parents=True, exist_ok=True)

    cmd = [
        "conda", "run", "-n", "targetdiff_env", "python",
        str(sample_script),
        "--pdb_path", str(pocket_pdb),
        "--result_path", str(out_path),
        "--num_samples", str(parameters.get("num_samples", 100)),
    ]

    subprocess.check_call(cmd)

    output_sdf = out_path / "generated_molecules.sdf"
    if not output_sdf.exists():
        raise RuntimeError(f"Expected output {output_sdf} was not produced by TargetDiff.")

    print(f"[Generation] Molecules generated at {output_sdf}")
    return str(output_sdf)


def run_generation(manifest_path, output_dir, mode="simulation",
                   db_path=None, campaign_id=None):
    """Run molecule generation with full telemetry.

    Args:
        manifest_path: Path to ingestion manifest.json.
        output_dir: Directory for output files.
        mode: "simulation", "rdkit", or "targetdiff".
        db_path: Optional telemetry database path.
        campaign_id: Optional campaign identifier.
    """
    manifest_path = Path(manifest_path).resolve()
    out_path = Path(output_dir).resolve()
    timestamp = datetime.now(timezone.utc).isoformat()

    if not manifest_path.exists():
        print(f"Error: Manifest file {manifest_path} not found.")
        sys.exit(1)

    with open(manifest_path, "r") as f:
        manifest = json.load(f)

    parameters = {**DEFAULT_PARAMS, "mode": mode}
    git_commit = _get_git_commit(TARGETDIFF_REPO) if TARGETDIFF_REPO.exists() else None

    db = None
    run_id = None
    if db_path and campaign_id:
        db = TelemetryDB(db_path)
        run_id = db.start_run(
            campaign_id=campaign_id,
            module_name="02_generation",
            input_path=str(manifest_path),
            parameters=parameters,
            git_commit=git_commit,
        )

    try:
        if mode == "simulation":
            output_sdf = run_generation_simulation(manifest, out_path, parameters)
        elif mode == "rdkit":
            output_sdf = run_generation_rdkit(manifest, out_path, parameters)
        elif mode == "targetdiff":
            output_sdf = run_generation_targetdiff(manifest, out_path, parameters)
        # Backwards compat: "production" maps to "rdkit" (was targetdiff stub)
        elif mode == "production":
            output_sdf = run_generation_rdkit(manifest, out_path, parameters)
        else:
            raise ValueError(f"Unknown mode: {mode}")

        metadata = {
            "module": "02_generation",
            "timestamp": timestamp,
            "mode": mode,
            "manifest_path": str(manifest_path),
            "pocket_pdb": manifest.get("best_pocket"),
            "output_sdf": output_sdf,
            "parameters": parameters,
            "targetdiff_git_commit": git_commit,
            "status": "success",
        }
        metadata_path = out_path / "run_metadata.json"
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2)

        if db and run_id:
            db.complete_run(run_id, "success", output_sdf)

        print(f"[Generation] Metadata written to {metadata_path}")

    except Exception as e:
        error_msg = traceback.format_exc()
        print(f"[Generation] FAILED: {e}")

        if db and run_id:
            db.complete_run(run_id, "failed", error_trace=error_msg)

        fail_metadata = {
            "module": "02_generation",
            "timestamp": timestamp,
            "mode": mode,
            "manifest_path": str(manifest_path),
            "status": "failed",
            "error": str(e),
            "traceback": error_msg,
        }
        fail_path = out_path / "run_metadata.json"
        out_path.mkdir(parents=True, exist_ok=True)
        with open(fail_path, "w") as f:
            json.dump(fail_metadata, f, indent=2)

        sys.exit(1)

    finally:
        if db:
            db.close()


def main():
    parser = argparse.ArgumentParser(description="Module 02: Generation")
    parser.add_argument("--manifest", required=True, help="Input manifest from ingestion")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--mode", choices=["simulation", "rdkit", "targetdiff", "production"],
                        default="simulation", help="Execution mode")
    parser.add_argument("--db_path", default=None, help="Path to telemetry database")
    parser.add_argument("--campaign_id", default=None, help="Campaign ID for telemetry")
    args = parser.parse_args()

    run_generation(args.manifest, args.output_dir, args.mode, args.db_path, args.campaign_id)


if __name__ == "__main__":
    main()
