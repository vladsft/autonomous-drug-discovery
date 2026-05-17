"""
Module 02: Generation — Structure-Aware Molecule Generator.

Four generation backends (plus a `production` legacy alias for `rdkit`):
  - simulation: Stub SDF for pipeline testing.
  - rdkit: Fragment-based combinatorial generation using RDKit. Produces
    diverse drug-like molecules with 3D conformers, optionally filtered by
    pocket volume. No GPU required.
  - pocket2mol: Pocket2Mol autoregressive generator (requires `pocket2mol_env`
    + checkpoint). Pocket-conditioned, ~11x faster than TargetDiff.
  - targetdiff: TargetDiff diffusion model (requires `targetdiff_env`
    + checkpoint). Highest-fidelity 3D, slowest.

Input contract:  manifest.json (from ingestion module)
Output contract: generated_molecules.sdf + run_metadata.json
"""

import os
import sys
import shutil
import argparse
import subprocess
import json
import random
import traceback
from datetime import datetime, timezone
from pathlib import Path

# Paths
MODULE_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = MODULE_DIR.parent.parent
TARGETDIFF_REPO = MODULE_DIR / "targetdiff"
POCKET2MOL_REPO = MODULE_DIR / "pocket2mol"

sys.path.insert(0, str(PROJECT_ROOT))
from telemetry import TelemetryDB

# Check for RDKit
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
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


# Per-mode generation defaults. Only the keys relevant to a given backend
# are logged to telemetry, so parameter audits aren't polluted with
# TargetDiff-specific knobs for RDKit runs (or vice versa).
_DEFAULT_PARAMS_BY_MODE = {
    "simulation": {"num_samples": 1},
    "rdkit": {"num_samples": 100, "seed": 42},
    "targetdiff": {
        "num_samples": 3,
        "sampling_steps": 1000,
        "noise_schedule": "polynomial_2",
        "batch_size": 16,
        "device": "cpu",
        "seed": 2021,
    },
    "pocket2mol": {
        "num_samples": 1,
        # device=cpu because pocket2mol_env's pytorch 1.10/cu113 doesn't support
        # Blackwell sm_120 (RTX 5060). Until the env is rebuilt on torch 2.4+/cu121
        # (plan.md Step 9), GPU mode JIT-thrashes and never produces molecules.
        "device": "cpu",
        "seed": 2020,
        "beam_size": 20,  # repo default 300; small for fast init on CPU
    },
}


def _default_params(mode: str) -> dict:
    """Return the defaults for a given generation mode."""
    if mode == "production":  # legacy alias → rdkit
        mode = "rdkit"
    return dict(_DEFAULT_PARAMS_BY_MODE.get(mode, _DEFAULT_PARAMS_BY_MODE["rdkit"]))


def _resolve_device(requested):
    """Resolve a device request to a concrete 'cuda' or 'cpu'.

    'auto' (or None) picks 'cuda' when a usable GPU is visible, else 'cpu'.
    An explicit 'cuda'/'cpu' is returned unchanged.
    """
    if requested not in (None, "auto"):
        return requested
    try:
        import torch
        return "cuda" if torch.cuda.is_available() else "cpu"
    except Exception:
        return "cpu"


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


def _find_attachment_points(mol):
    """Find atoms that can accept a new single bond (have implicit Hs).

    Returns list of atom indices that have at least one implicit hydrogen,
    meaning they can form a new bond without violating valence rules.
    """
    attachment = []
    for atom in mol.GetAtoms():
        if atom.GetNumImplicitHs() > 0:
            attachment.append(atom.GetIdx())
    return attachment


def _attach_fragment(mol, frag_smi, rng):
    """Attach a fragment to the molecule at a valid attachment point.

    Uses RDKit's RWMol to combine molecules and form a proper bond between
    atoms that have available valence (implicit Hs). Returns None if no
    valid attachment is possible.
    """
    frag = Chem.MolFromSmiles(frag_smi)
    if frag is None:
        return None

    main_points = _find_attachment_points(mol)
    frag_points = _find_attachment_points(frag)

    if not main_points or not frag_points:
        return None

    # Try a few random attachment point pairs
    attempts = min(5, len(main_points) * len(frag_points))
    for _ in range(attempts):
        main_idx = rng.choice(main_points)
        frag_idx = rng.choice(frag_points)

        combo = AllChem.CombineMols(mol, frag)
        editable = Chem.RWMol(combo)
        # Fragment atom indices are offset by the number of atoms in mol
        editable.AddBond(main_idx, frag_idx + mol.GetNumAtoms(), Chem.BondType.SINGLE)
        try:
            new_mol = editable.GetMol()
            Chem.SanitizeMol(new_mol)
            return new_mol
        except Exception:
            continue

    return None


def _assemble_molecule(scaffolds, substituents, linkers, rng):
    """Assemble a random molecule from fragments.

    Strategy: pick 1-2 scaffolds, connect with a linker, decorate with
    1-3 substituents. Uses proper RDKit bond formation to ensure valid
    valence at all attachment points.
    """
    # Pick primary scaffold
    scaffold_smi = rng.choice(scaffolds)
    mol = Chem.MolFromSmiles(scaffold_smi)
    if mol is None:
        return None

    # Optionally attach a second scaffold via a linker (30% chance)
    if rng.random() < 0.3:
        scaffold2_smi = rng.choice(scaffolds)
        linker_smi = rng.choice(linkers)

        if linker_smi:
            # Attach linker to scaffold, then second scaffold to linker
            linked = _attach_fragment(mol, linker_smi, rng)
            if linked is not None:
                mol = linked
                linked2 = _attach_fragment(mol, scaffold2_smi, rng)
                if linked2 is not None:
                    mol = linked2
        else:
            # Direct bond between scaffolds
            linked = _attach_fragment(mol, scaffold2_smi, rng)
            if linked is not None:
                mol = linked

    # Try BRICS decomposition and reassembly for more diversity
    try:
        frags = list(BRICS.BRICSDecompose(Chem.MolToSmiles(mol)))
        if len(frags) >= 2:
            frag_mols = [Chem.MolFromSmiles(f) for f in frags]
            frag_mols = [f for f in frag_mols if f is not None]
            if len(frag_mols) >= 2:
                rebuilt = list(BRICS.BRICSBuild(frag_mols))
                if rebuilt:
                    candidate = rng.choice(rebuilt[:10])
                    candidate = Chem.MolFromSmiles(Chem.MolToSmiles(candidate))
                    if candidate is not None:
                        mol = candidate
    except Exception:
        pass  # BRICS can fail on some scaffolds, that's fine

    # Add 1-3 substituents at valid attachment points
    n_subs = rng.randint(1, 3)
    for _ in range(n_subs):
        sub_smi = rng.choice(substituents)
        decorated = _attach_fragment(mol, sub_smi, rng)
        if decorated is not None:
            mol = decorated

    return mol


# Loose pre-screening bounds. Intentionally wider than default_scoring_config.json
# so the generator doesn't silently shrink the candidate pool before the real
# screening step gets to apply its own (tunable) thresholds.
_PRESCREEN_MW_RANGE = (80, 650)
_PRESCREEN_LOGP_RANGE = (-3, 7)
_PRESCREEN_MAX_RINGS = 6


def _is_drug_like_quick(mol, min_ha, max_ha):
    """Fast pre-filter. Rejects molecules clearly outside drug-like envelope.

    Tuned to be strictly looser than the Stage-3 screening config so we never
    eliminate a molecule the user-configurable screener would have kept.
    """
    ha = mol.GetNumHeavyAtoms()
    if ha < min_ha or ha > max_ha:
        return False
    mw = Descriptors.MolWt(mol)
    if mw < _PRESCREEN_MW_RANGE[0] or mw > _PRESCREEN_MW_RANGE[1]:
        return False
    logp = Descriptors.MolLogP(mol)
    if logp < _PRESCREEN_LOGP_RANGE[0] or logp > _PRESCREEN_LOGP_RANGE[1]:
        return False
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > _PRESCREEN_MAX_RINGS:
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
    """TargetDiff mode: execute TargetDiff diffusion model inference.

    TargetDiff's sample_for_pocket.py requires a config YAML (with checkpoint path)
    as a positional arg, and writes individual SDF files into a sdf/ subdirectory.
    This wrapper handles config generation and consolidates the output into a single
    generated_molecules.sdf for downstream pipeline consumption.
    """
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

    # Locate checkpoint — use the pretrained model shipped with the repo
    checkpoint = TARGETDIFF_REPO / "pretrained_models" / "pretrained_diffusion.pt"
    if not checkpoint.exists():
        raise FileNotFoundError(
            f"TargetDiff checkpoint not found at {checkpoint}. "
            f"Download it to {checkpoint}."
        )

    out_path.mkdir(parents=True, exist_ok=True)

    # Generate a sampling config YAML for this run
    num_samples = parameters.get("num_samples", 100)
    sampling_config = out_path / "sampling_config.yml"
    import yaml
    config_data = {
        "model": {
            "checkpoint": str(checkpoint),
        },
        "sample": {
            "seed": parameters.get("seed", 2021),
            "num_samples": num_samples,
            "num_steps": parameters.get("sampling_steps", 1000),
            "pos_only": False,
            "center_pos_mode": "protein",
            "sample_num_atoms": "prior",
        },
    }
    with open(sampling_config, "w") as f:
        yaml.dump(config_data, f)

    # Wipe stale outputs so a rerun can't silently include molecules from a
    # previous campaign when we glob the SDF directory.
    td_result_path = out_path / "targetdiff_raw"
    if td_result_path.exists():
        shutil.rmtree(td_result_path)

    # num_samples is set in the YAML above; do not pass --num_samples here,
    # otherwise the CLI flag silently overrides the config (order-dependent).
    cmd = [
        os.environ.get("CONDA_EXE", "conda"), "run", "-n", "targetdiff_env", "python",
        str(sample_script),
        str(sampling_config),
        "--pdb_path", str(pocket_pdb),
        "--result_path", str(td_result_path),
        "--device", parameters["device"],
        "--batch_size", str(parameters.get("batch_size", 16)),
    ]

    # The script does `import utils.misc` relative to the repo root, but
    # Python sets sys.path[0] to the script's directory (scripts/), not cwd.
    # Inject PYTHONPATH so the in-repo `utils/`, `models/`, `datasets/` packages
    # are importable.
    env = os.environ.copy()
    existing_pp = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = (
        f"{TARGETDIFF_REPO}{os.pathsep}{existing_pp}" if existing_pp else str(TARGETDIFF_REPO)
    )

    subprocess.check_call(cmd, cwd=str(TARGETDIFF_REPO), env=env)

    # Consolidate individual SDF files into a single generated_molecules.sdf
    sdf_dir = td_result_path / "sdf"
    output_sdf = out_path / "generated_molecules.sdf"

    if not HAS_RDKIT:
        raise ImportError("RDKit is required to consolidate TargetDiff output.")

    RDLogger.DisableLog("rdApp.*")
    individual_sdfs = sorted(sdf_dir.glob("*.sdf")) if sdf_dir.exists() else []
    if not individual_sdfs:
        raise RuntimeError(
            f"TargetDiff produced no SDF files in {sdf_dir}. "
            f"Check {td_result_path} for errors."
        )

    writer = Chem.SDWriter(str(output_sdf))
    mol_count = 0
    for sdf_file in individual_sdfs:
        supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
        for mol in supplier:
            if mol is not None:
                mol.SetProp("molecule_id", f"mol_{mol_count:04d}")
                mol.SetProp("_Name", f"mol_{mol_count:04d}")
                try:
                    mol.SetProp("smiles", Chem.MolToSmiles(mol))
                except Exception:
                    pass
                writer.write(mol)
                mol_count += 1
    writer.close()

    print(f"[Generation] {mol_count} molecules consolidated from {len(individual_sdfs)} TargetDiff outputs")
    print(f"[Generation] Output: {output_sdf}")
    return str(output_sdf)


def _compute_pocket_centroid(pocket_pdb_path):
    """Return (cx, cy, cz) centroid of atoms in a pocket PDB file."""
    coords = []
    with open(pocket_pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not coords:
        raise ValueError(f"No atoms found in pocket PDB: {pocket_pdb_path}")
    n = len(coords)
    return sum(c[0] for c in coords) / n, sum(c[1] for c in coords) / n, sum(c[2] for c in coords) / n


def run_generation_pocket2mol(manifest, out_path, parameters):
    """Pocket2Mol mode: autoregressive pocket-conditioned molecule generation.

    Pocket2Mol generates molecules atom-by-atom conditioned on the binding pocket
    geometry. ~11x faster than TargetDiff (7s/molecule vs 78s on GPU).

    Requires:
      - pocket2mol_env conda environment (envs/env_pocket2mol.yml)
      - Pocket2Mol repo cloned at modules/02_generation/pocket2mol/
      - Pretrained checkpoint at pocket2mol/ckpt/pretrained.pt
        (download: https://drive.google.com/drive/folders/1KfdOczjUPITPhIvCuBmnj4xFTV-iI2xB)
    """
    if not POCKET2MOL_REPO.exists():
        raise FileNotFoundError(
            f"Pocket2Mol repo not found at {POCKET2MOL_REPO}. "
            f"Clone it: git clone https://github.com/pengxingang/Pocket2Mol.git {POCKET2MOL_REPO}"
        )

    checkpoint = POCKET2MOL_REPO / "ckpt" / "pretrained_Pocket2Mol.pt"
    if not checkpoint.exists():
        raise FileNotFoundError(
            f"Pocket2Mol checkpoint not found at {checkpoint}. "
            f"Download from https://drive.google.com/drive/folders/1KfdOczjUPITPhIvCuBmnj4xFTV-iI2xB"
        )

    sample_script = POCKET2MOL_REPO / "sample_for_pdb.py"
    if not sample_script.exists():
        raise FileNotFoundError(f"Pocket2Mol sample script not found at {sample_script}")

    # Get the full protein PDB (not the pocket PDB — Pocket2Mol needs the full structure)
    input_pdb = Path(manifest.get("input_pdb", ""))
    if not input_pdb.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb}")

    # Get pocket center: P2Rank provides it directly; fpocket requires centroid computation
    center = manifest.get("best_pocket_center")
    if center is not None:
        cx, cy, cz = center
    else:
        pocket_pdb = manifest.get("best_pocket")
        if not pocket_pdb or not Path(pocket_pdb).exists():
            raise ValueError("No pocket center or pocket PDB available in manifest.")
        cx, cy, cz = _compute_pocket_centroid(pocket_pdb)

    print(f"[Generation] (POCKET2MOL MODE) Pocket center: ({cx:.2f}, {cy:.2f}, {cz:.2f})")
    print(f"[Generation] Running Pocket2Mol on {input_pdb.name}...")

    out_path.mkdir(parents=True, exist_ok=True)
    p2m_raw = out_path / "pocket2mol_raw"
    # Wipe stale SDFs so we don't mistakenly sweep in molecules from a prior run.
    if p2m_raw.exists():
        shutil.rmtree(p2m_raw)
    p2m_raw.mkdir(parents=True, exist_ok=True)

    num_samples = parameters.get("num_samples", 100)

    # Write a run-specific config so we can set the absolute checkpoint path
    # and num_samples without mutating the repo's config file.
    import yaml
    base_config_path = POCKET2MOL_REPO / "configs" / "sample_for_pdb.yml"
    with open(base_config_path) as f:
        run_config = yaml.safe_load(f)

    run_config["model"]["checkpoint"] = str(checkpoint)
    run_config["sample"]["num_samples"] = num_samples
    if "beam_size" in parameters:
        run_config["sample"]["beam_size"] = parameters["beam_size"]

    run_config_path = out_path / "pocket2mol_run_config.yml"
    with open(run_config_path, "w") as f:
        yaml.dump(run_config, f)

    # Center passed as a plain comma-separated string — no leading space needed
    # when using subprocess list args (the leading-space trick is shell-only).
    center_str = f"{cx:.3f},{cy:.3f},{cz:.3f}"

    device = parameters.get("device", "cpu")

    conda_bin = os.environ.get("CONDA_EXE", "conda")
    cmd = [
        conda_bin, "run", "-n", "pocket2mol_env", "python",
        str(sample_script),
        "--pdb_path", str(input_pdb),
        "--center", center_str,
        "--outdir", str(p2m_raw),
        "--config", str(run_config_path),
        "--device", device,
    ]

    subprocess.check_call(cmd, cwd=str(POCKET2MOL_REPO))

    # Consolidate all SDF outputs into a single generated_molecules.sdf
    if not HAS_RDKIT:
        raise ImportError("RDKit is required to consolidate Pocket2Mol output.")

    RDLogger.DisableLog("rdApp.*")
    output_sdf = out_path / "generated_molecules.sdf"
    # Pocket2Mol writes to {outdir}/{timestamp}/SDF/*.sdf
    sdf_files = sorted(p2m_raw.glob("**/SDF/*.sdf"))

    if not sdf_files:
        raise RuntimeError(
            f"Pocket2Mol produced no SDF files in {p2m_raw}. "
            "Check for initialization failures (known issue on some targets)."
        )

    writer = Chem.SDWriter(str(output_sdf))
    mol_count = 0
    fail_count = 0
    for sdf_file in sdf_files:
        supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
        for mol in supplier:
            if mol is None:
                fail_count += 1
                continue
            try:
                smi = Chem.MolToSmiles(mol)
            except Exception:
                fail_count += 1
                continue
            mol.SetProp("molecule_id", f"mol_{mol_count:04d}")
            mol.SetProp("_Name", f"mol_{mol_count:04d}")
            mol.SetProp("smiles", smi)
            mol.SetProp("generator", "pocket2mol")
            writer.write(mol)
            mol_count += 1
    writer.close()

    if mol_count == 0:
        raise RuntimeError(
            f"No valid molecules produced by Pocket2Mol. "
            f"{fail_count} reconstruction failures. "
            "Try re-running — diffusion models have stochastic failure modes."
        )

    print(f"[Generation] Pocket2Mol: {mol_count} valid molecules ({fail_count} reconstruction failures)")
    print(f"[Generation] Output: {output_sdf}")
    return str(output_sdf)


def run_generation(manifest_path, output_dir, mode="simulation",
                   db_path=None, campaign_id=None, num_samples=None, device="auto"):
    """Run molecule generation with full telemetry.

    Args:
        manifest_path: Path to ingestion manifest.json.
        output_dir: Directory for output files.
        mode: "simulation", "rdkit", or "targetdiff".
        db_path: Optional telemetry database path.
        campaign_id: Optional campaign identifier.
        num_samples: Optional override for the per-mode default campaign size.
        device: Compute device for GPU-capable backends — "auto" (detect a
            GPU), "cuda", or "cpu". Ignored by the CPU-only backends.
    """
    manifest_path = Path(manifest_path).resolve()
    out_path = Path(output_dir).resolve()
    timestamp = datetime.now(timezone.utc).isoformat()

    if not manifest_path.exists():
        print(f"Error: Manifest file {manifest_path} not found.")
        sys.exit(1)

    with open(manifest_path, "r") as f:
        manifest = json.load(f)

    parameters = {**_default_params(mode), "mode": mode}
    if num_samples is not None:
        parameters["num_samples"] = num_samples
    # Device selection. Only the pocket-conditioned backends carry a "device"
    # default. An explicit --device always wins; a bare "auto" upgrades
    # TargetDiff's legacy cpu default to the detected GPU, while Pocket2Mol
    # keeps its intentional cpu default (its cu113 env predates Blackwell).
    if device not in (None, "auto"):
        parameters["device"] = device
    elif mode == "targetdiff":
        parameters["device"] = _resolve_device("auto")
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
        elif mode == "pocket2mol":
            output_sdf = run_generation_pocket2mol(manifest, out_path, parameters)
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
    parser.add_argument("--mode", choices=["simulation", "rdkit", "targetdiff", "pocket2mol", "production"],
                        default="simulation", help="Execution mode")
    parser.add_argument("--db_path", default=None, help="Path to telemetry database")
    parser.add_argument("--campaign_id", default=None, help="Campaign ID for telemetry")
    parser.add_argument("--num_samples", type=int, default=None,
                        help="Number of molecules to generate (overrides the per-mode default)")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto",
                        help="Compute device for GPU-capable backends "
                             "(targetdiff/pocket2mol); 'auto' detects a GPU")
    args = parser.parse_args()

    run_generation(args.manifest, args.output_dir, args.mode, args.db_path,
                   args.campaign_id, num_samples=args.num_samples, device=args.device)


if __name__ == "__main__":
    main()
