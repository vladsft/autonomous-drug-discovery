"""
Module 04: Docking — AutoDock Vina / TDC Oracle Wrapper.

Three execution modes:
  - simulation: Dummy scores for pipeline testing.
  - triage:     TDC Docking Oracle (Vina under the hood, SMILES-based, fast).
  - production: Full manual Vina pipeline (PDB→PDBQT, box definition, per-ligand docking).

Input contract:  manifest.json + candidates SDF directory
Output contract: docking_results.csv + run_metadata.json
"""

import sys
import os
import csv
import argparse
import subprocess
import json
import traceback
from datetime import datetime, timezone
from pathlib import Path

# Paths
MODULE_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = MODULE_DIR.parent.parent

sys.path.insert(0, str(PROJECT_ROOT))
from telemetry import TelemetryDB

# Check for TDC availability
try:
    from tdc import Oracle as TDCOracle
    HAS_TDC = True
except ImportError:
    HAS_TDC = False

# Check for RDKit (needed for SDF reading)
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

# Check for Vina Python API + Meeko (production docking)
try:
    from vina import Vina
    HAS_VINA = True
except ImportError:
    HAS_VINA = False

try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    HAS_MEEKO = True
except ImportError:
    HAS_MEEKO = False


# Default docking parameters
DEFAULT_PARAMS = {
    "exhaustiveness": 8,
    "num_modes": 9,
    "energy_range": 3,
    "box_size": [20, 20, 20],
}

# Metal ions worth keeping in the receptor: for metalloenzymes (e.g. the
# cap-snatching endonuclease) the catalytic metals ARE the pharmacophore, and
# stripping them as "ligands" guts the binding site. Keyed by element symbol;
# value is the AutoDock atom type.
_RECEPTOR_METALS = {
    "MG": "Mg", "MN": "Mn", "ZN": "Zn", "CA": "Ca", "FE": "Fe",
    "NA": "Na", "K": "K", "CO": "Co", "NI": "Ni", "CU": "Cu",
}


def _receptor_pdb_from_manifest(manifest: dict) -> Path:
    """Resolve the receptor PDB to dock against.

    Prefers the dedicated `receptor_pdb` key (set by ingestion — this is the
    cleaned structure when --clean was used), falling back to `input_pdb` for
    manifests written before that field existed.
    """
    receptor = manifest.get("receptor_pdb") or manifest.get("input_pdb")
    if not receptor:
        raise ValueError("Manifest has neither 'receptor_pdb' nor 'input_pdb'.")
    path = Path(receptor)
    if not path.exists():
        raise FileNotFoundError(f"Receptor PDB {path} not found.")
    return path


# ---------------------------------------------------------------------------
# Mode: Simulation
# ---------------------------------------------------------------------------

def run_docking_simulation(manifest, candidates_dir, out_path, parameters, db=None, run_id=None):
    """Simulation mode: generate dummy docking results."""
    print("[Docking] (SIMULATION MODE) Generating dummy docking scores...")

    out_path.mkdir(parents=True, exist_ok=True)

    simulated = [
        {"ligand_id": "mol_0000", "smiles": "c1ccccc1", "affinity": -9.5},
        {"ligand_id": "mol_0001", "smiles": "CCO", "affinity": -8.2},
        {"ligand_id": "mol_0002", "smiles": "CC(=O)O", "affinity": -7.1},
    ]

    results_file = out_path / "docking_results.csv"
    with open(results_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ligand_id", "smiles", "affinity"])
        for r in simulated:
            writer.writerow([r["ligand_id"], r["smiles"], r["affinity"]])

    if db and run_id:
        db.log_molecules_batch(run_id, [
            {
                "molecule_id": r["ligand_id"],
                "smiles": r["smiles"],
                "docking_score": r["affinity"],
            }
            for r in simulated
        ])

    print(f"[Docking] Simulated results: {results_file}")
    return str(results_file)


# ---------------------------------------------------------------------------
# Mode: Triage (TDC Oracle)
# ---------------------------------------------------------------------------

def run_docking_triage(manifest, candidates_dir, out_path, parameters, db=None, run_id=None):
    """Triage mode: fast SMILES-based docking via TDC Oracle."""
    if not HAS_TDC or not HAS_RDKIT:
        missing = []
        if not HAS_TDC:
            missing.append(
                "PyTDC (note: not installable on Python 3.13+; "
                "requires base env on Python ≤3.11)"
            )
        if not HAS_RDKIT:
            missing.append("RDKit (conda install -c conda-forge rdkit)")
        print(f"[Docking] Triage mode requires: {', '.join(missing)}.")
        print(f"[Docking] Falling back to simulation mode — "
              "use --mode production for real docking without PyTDC.")
        return run_docking_simulation(manifest, candidates_dir, out_path, parameters, db, run_id)

    receptor_pdb = _receptor_pdb_from_manifest(manifest)

    # Find candidate SDF
    sdf_file = _find_candidates_sdf(candidates_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    print(f"[Docking] (TRIAGE MODE) TDC Oracle scoring against {receptor_pdb.name}...")
    print(f"[Docking] Candidates: {sdf_file}")

    # Initialize TDC docking oracle
    oracle = TDCOracle(name='Docking_Score', target_pdb=str(receptor_pdb))

    # Read molecules from SDF
    supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
    results = []
    mol_batch = []

    for idx, mol in enumerate(supplier):
        mol_id = _mol_id_from_mol(mol, idx)
        if mol is None:
            results.append({"ligand_id": mol_id, "smiles": None, "affinity": None, "error": "invalid_mol"})
            continue

        try:
            smiles = Chem.MolToSmiles(mol)
            score = oracle(smiles)
            results.append({"ligand_id": mol_id, "smiles": smiles, "affinity": float(score)})
            mol_batch.append({
                "molecule_id": mol_id,
                "smiles": smiles,
                "docking_score": float(score),
            })

            print(f"  {mol_id}: {smiles[:40]}... → {score:.2f} kcal/mol")

        except Exception as e:
            results.append({"ligand_id": mol_id, "smiles": None, "affinity": None, "error": str(e)})
            print(f"  {mol_id}: ERROR — {e}")

    if db and run_id and mol_batch:
        db.log_molecules_batch(run_id, mol_batch)

    # Sort by affinity (best = most negative)
    scored = [r for r in results if r["affinity"] is not None]
    scored.sort(key=lambda r: r["affinity"])

    # Write results
    results_file = out_path / "docking_results.csv"
    with open(results_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ligand_id", "smiles", "affinity"])
        for r in scored:
            writer.writerow([r["ligand_id"], r["smiles"], r["affinity"]])

    print(f"[Docking] {len(scored)} molecules scored. Best: {scored[0]['affinity']:.2f} kcal/mol" if scored else "[Docking] No molecules scored.")
    print(f"[Docking] Results: {results_file}")
    return str(results_file)


# ---------------------------------------------------------------------------
# Mode: Production (Manual Vina)
# ---------------------------------------------------------------------------

def run_docking_production(manifest, candidates_dir, out_path, parameters, db=None, run_id=None):
    """Production mode: full Vina docking pipeline using Python API + Meeko.

    Creates one Vina instance, loads the receptor and computes grid maps ONCE,
    then iterates ligands. This is ~50x faster than per-ligand setup when
    docking 100 molecules against the same target.
    """
    if not HAS_VINA or not HAS_MEEKO or not HAS_RDKIT:
        missing = []
        if not HAS_VINA:
            missing.append("vina (conda install -c conda-forge vina)")
        if not HAS_MEEKO:
            missing.append("meeko (pip install meeko)")
        if not HAS_RDKIT:
            missing.append("rdkit")
        raise ImportError(f"Production docking requires: {', '.join(missing)}")

    receptor_pdb = _receptor_pdb_from_manifest(manifest)

    sdf_file = _find_candidates_sdf(candidates_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    print(f"[Docking] (PRODUCTION MODE) Full Vina docking against {receptor_pdb.name}...")

    # 1. Prepare receptor ONCE.
    receptor_pdbqt = out_path / f"{receptor_pdb.stem}_receptor.pdbqt"
    print(f"[Docking] Preparing receptor PDBQT for {receptor_pdb.name}...")
    n_receptor_atoms = _prepare_receptor_pdbqt(receptor_pdb, receptor_pdbqt)
    print(f"[Docking] Receptor PDBQT: {n_receptor_atoms} atoms.")

    # 2. Box center: prefer ingestion's pre-computed centroid.
    if manifest.get("best_pocket_center"):
        box_center = manifest["best_pocket_center"]
    elif manifest.get("best_pocket"):
        box_center = _compute_pocket_centroid(manifest["best_pocket"])
    else:
        raise ValueError("Manifest has neither best_pocket_center nor best_pocket.")
    # Box size: scale to the actual pocket extent rather than a fixed cube, so a
    # small pocket isn't given room for ligands to drift and a large one isn't
    # truncated. Falls back to the configured default if the pocket PDB is absent.
    box_size = _pocket_box_size(manifest.get("best_pocket"),
                                parameters.get("box_size", [20, 20, 20]))
    parameters["box_size"] = box_size
    print(f"[Docking] Box center: {box_center}, size: {box_size}")

    # 3. Create Vina instance and load receptor + maps ONCE.
    v = Vina(sf_name="vina")
    v.set_receptor(str(receptor_pdbqt))
    v.compute_vina_maps(center=box_center, box_size=box_size)

    exhaustiveness = parameters.get("exhaustiveness", 8)
    n_poses = parameters.get("num_modes", 9)

    # 4. Iterate ligands; reuse the same Vina/receptor/maps across all of them.
    supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
    results = []
    mol_batch = []

    for idx, mol in enumerate(supplier):
        mol_id = _mol_id_from_mol(mol, idx)
        if mol is None:
            results.append({"ligand_id": mol_id, "smiles": None, "affinity": None, "error": "invalid_mol"})
            continue

        smiles = Chem.MolToSmiles(mol)
        try:
            affinity = _dock_one_ligand(v, mol, exhaustiveness, n_poses, out_path, mol_id)
            results.append({"ligand_id": mol_id, "smiles": smiles, "affinity": affinity})
            mol_batch.append({
                "molecule_id": mol_id,
                "smiles": smiles,
                "docking_score": affinity,
            })
            print(f"  {mol_id}: {smiles[:50]} → {affinity:.2f} kcal/mol")

        except Exception as e:
            results.append({"ligand_id": mol_id, "smiles": smiles, "affinity": None, "error": str(e)})
            mol_batch.append({
                "molecule_id": mol_id,
                "smiles": smiles,
                "docking_score": None,
                "stage_eliminated": f"docking_error: {e}",
            })
            print(f"  {mol_id}: ERROR — {e}")

    # Batch-write telemetry (one transaction instead of per-ligand commits).
    if db and run_id and mol_batch:
        db.log_molecules_batch(run_id, mol_batch)

    # Sort by affinity (most negative = best)
    scored = [r for r in results if r["affinity"] is not None]
    scored.sort(key=lambda r: r["affinity"])

    # Guard: real docking always produces a spread of affinities. If every
    # ligand scored identically (classically all 0.0 — an empty/garbage
    # receptor), the run is invalid; fail loudly instead of writing a CSV of
    # noise that downstream ranking would treat as real.
    if len(scored) >= 2:
        spread = max(r["affinity"] for r in scored) - min(r["affinity"] for r in scored)
        if spread < 1e-6:
            raise RuntimeError(
                f"All {len(scored)} ligands docked to an identical affinity "
                f"({scored[0]['affinity']}). This indicates a broken receptor "
                f"or grid setup, not a real result — refusing to report it."
            )

    results_file = out_path / "docking_results.csv"
    with open(results_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ligand_id", "smiles", "affinity"])
        for r in scored:
            writer.writerow([r["ligand_id"], r["smiles"], r["affinity"]])

    if scored:
        print(f"[Docking] {len(scored)} molecules scored. Best: {scored[0]['affinity']:.2f} kcal/mol")
    else:
        print("[Docking] No molecules scored successfully.")
    print(f"[Docking] Results: {results_file}")
    return str(results_file)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_candidates_sdf(candidates_dir):
    """Find the candidates SDF file in the given directory."""
    candidates_path = Path(candidates_dir).resolve()
    for name in ["screened_molecules.sdf", "generated_molecules.sdf"]:
        sdf = candidates_path / name
        if sdf.exists():
            return sdf
    raise FileNotFoundError(f"No candidate SDF found in {candidates_path}")


def _compute_pocket_centroid(pocket_pdb_path):
    """Compute pocket centroid for Vina box placement."""
    try:
        coords = []
        with open(pocket_pdb_path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    coords.append((x, y, z))
        if coords:
            n = len(coords)
            return [round(sum(c[0] for c in coords) / n, 2),
                    round(sum(c[1] for c in coords) / n, 2),
                    round(sum(c[2] for c in coords) / n, 2)]
    except Exception:
        pass
    return [0.0, 0.0, 0.0]


# --- Atom typing and charge assignment for the receptor PDBQT ---------------
#
# AutoDock Vina atom types: C, A (aromatic C), N, NA (H-bond acceptor N),
# OA (H-bond acceptor O), S, SA (H-bond acceptor S), HD (polar H), F, Cl, Br,
# I, P. We also apply simple integer charges on known charged side chains
# (ASP/GLU -1, LYS/ARG +1, HIS +1 on NE2/ND1 when protonated) so that the
# electrostatics term isn't entirely zero for the receptor.

_AROMATIC_CARBONS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "CD2", "CE1"},
}

# Simple residue-level charges distributed over the key side-chain atoms.
# We keep it integer-valued because we don't have Gasteiger calculations.
_CHARGED_ATOMS = {
    ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5,
    ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5,
    ("LYS", "NZ"): +1.0,
    ("ARG", "NH1"): +0.5, ("ARG", "NH2"): +0.5,
    # HIS is often neutral at physiological pH; leave it at 0.
}


def _autodock_type(atom_name: str, residue_name: str, element: str) -> str:
    el = element.strip().upper()
    name = atom_name.strip().upper()

    if el == "C":
        if name in _AROMATIC_CARBONS.get(residue_name, ()):
            return "A"
        return "C"
    if el == "N":
        # Aromatic sp2 nitrogen acceptors in HIS
        if residue_name == "HIS" and name in ("ND1", "NE2"):
            return "NA"
        # TRP NE1 is a donor (H-bond donor); leave as N so HD gets the donor mark.
        return "N"
    if el == "O":
        return "OA"
    if el == "S":
        if residue_name == "CYS" and name == "SG":
            return "SA"
        if residue_name == "MET" and name == "SD":
            return "S"
        return "S"
    if el == "H":
        if name.startswith(("HN", "HE", "HD", "HH", "HG", "HZ")):
            return "HD"
        return "H"
    if el == "F":
        return "F"
    if el == "CL":
        return "Cl"
    if el == "BR":
        return "Br"
    if el == "I":
        return "I"
    if el == "P":
        return "P"
    return el


def _pdbqt_atom_line(serial, name, resname, chain, resnum, x, y, z, charge, ad_type):
    """Format one PDBQT ATOM record."""
    return (
        f"ATOM  {serial:5d} {name:<4s} {resname:>3s} "
        f"{chain:1s}{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}    "
        f"{charge:+6.3f} {ad_type:<2s}"
    )


def _prepare_receptor_pdbqt(receptor_pdb, output_pdbqt) -> int:
    """Convert receptor PDB to PDBQT with AutoDock atom typing + basic charges.

    Waters and organic ligands are stripped, but catalytic metal ions are
    deliberately KEPT — for a metalloenzyme the metals are part of the binding
    site, and discarding them silently changes the docking problem. Alt-loc
    atoms other than A/blank are skipped. Charges are zero for uncharged atoms
    and integer-valued on classic charged side chains (ASP/GLU/LYS/ARG).

    Returns the number of ATOM records written, and raises RuntimeError if that
    is zero — an empty receptor produces meaningless (typically all-0.0)
    docking scores, so it must fail loudly here rather than downstream.

    Note: this does not add hydrogens or compute Gasteiger charges. For
    publication-grade work, pre-protonate the structure (pdb2pqr / reduce);
    this routine is a dependency-light approximation, not a substitute.
    """
    import gemmi

    structure = gemmi.read_structure(str(receptor_pdb))

    # Collect catalytic metal ions BEFORE remove_ligands_and_waters() drops them.
    metals = []
    if len(structure) > 0:
        for chain in structure[0]:
            for residue in chain:
                resname = residue.name.strip().upper()
                if resname not in _RECEPTOR_METALS:
                    continue
                for atom in residue:
                    metals.append((chain.name, residue.seqid.num, resname, atom))
    if metals:
        print(f"[Docking] Keeping {len(metals)} metal ion(s) in the receptor: "
              f"{sorted({m[2] for m in metals})}")

    structure.remove_ligands_and_waters()

    lines = [
        f"REMARK  Name = {receptor_pdb}",
        "REMARK                            x       y       z     vdW  Elec       q    Type",
        "REMARK                         _______ _______ _______ _____ _____    ______ ____",
    ]
    n_atoms = 0

    if len(structure) > 0:
        for chain in structure[0]:  # first model only (NMR ensembles)
            for residue in chain:
                for atom in residue:
                    # Skip alt-loc atoms that aren't the primary one.
                    # gemmi >=0.7 reports "no altloc" as a NUL byte ("\x00"),
                    # which str.strip() does not remove — drop it explicitly so
                    # ordinary atoms aren't mistaken for secondary alt-locs.
                    alt = (atom.altloc or "").replace("\x00", "").strip()
                    if alt and alt.upper() != "A":
                        continue
                    element = atom.element.name.strip()
                    ad_type = _autodock_type(atom.name, residue.name, element)
                    charge = _CHARGED_ATOMS.get(
                        (residue.name, atom.name.strip().upper()), 0.0)
                    lines.append(_pdbqt_atom_line(
                        atom.serial, atom.name, residue.name, chain.name,
                        residue.seqid.num, atom.pos.x, atom.pos.y, atom.pos.z,
                        charge, ad_type))
                    n_atoms += 1

    # Re-append the metal ions with their formal charge and element AD type.
    for chain_name, resnum, resname, atom in metals:
        ad_type = _RECEPTOR_METALS[resname]
        lines.append(_pdbqt_atom_line(
            atom.serial, atom.name, resname, chain_name, resnum,
            atom.pos.x, atom.pos.y, atom.pos.z, 2.0, ad_type))
        n_atoms += 1

    lines.append("END")

    if n_atoms == 0:
        raise RuntimeError(
            f"Receptor preparation produced 0 atoms from {receptor_pdb}. "
            "The PDBQT would be empty and every ligand would score 0.0 — "
            "check the input structure and gemmi parsing."
        )

    with open(output_pdbqt, "w") as f:
        f.write("\n".join(lines) + "\n")
    return n_atoms


def _pocket_box_size(pocket_pdb, default, padding=8.0, lo=16.0, hi=30.0):
    """Vina box dimensions scaled to the pocket's atomic extent.

    Returns [x, y, z] = per-axis span of the pocket atoms + `padding`, each
    clamped to [lo, hi]. Falls back to `default` when the pocket PDB is missing
    or unreadable so docking never breaks on box sizing.
    """
    if not pocket_pdb or not Path(pocket_pdb).exists():
        return list(default)
    xs, ys, zs = [], [], []
    try:
        with open(pocket_pdb) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
    except (ValueError, OSError):
        return list(default)
    if not xs:
        return list(default)
    box = []
    for axis in (xs, ys, zs):
        span = (max(axis) - min(axis)) + padding
        box.append(round(min(hi, max(lo, span)), 1))
    return box


def _prepare_ligand_pdbqt(mol):
    """Convert an RDKit mol to PDBQT string using Meeko."""
    preparator = MoleculePreparation()
    # Ensure 3D coords exist
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, randomSeed=42)
    mol_with_h = Chem.AddHs(mol, addCoords=True)
    if mol_with_h.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol_with_h, randomSeed=42)
    preparator.prepare(mol_with_h)
    pdbqt_string = PDBQTWriterLegacy.write_string(preparator.setup)[0]
    return pdbqt_string


def _mol_id_from_mol(mol, idx: int) -> str:
    """Preserve upstream molecule_id property if present; else derive from index."""
    if mol is not None and mol.HasProp("molecule_id"):
        return mol.GetProp("molecule_id")
    return f"mol_{idx:04d}"


def _dock_one_ligand(v, mol, exhaustiveness, n_poses, work_dir, mol_id):
    """Dock a single ligand against a pre-configured Vina instance.

    The receptor and grid maps must already be loaded into `v` — we only
    swap the ligand. This is dramatically faster than rebuilding the grid
    maps for every molecule.
    """
    ligand_pdbqt_str = _prepare_ligand_pdbqt(mol)
    v.set_ligand_from_string(ligand_pdbqt_str)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

    energies = v.energies()
    best_affinity = float(energies[0][0])  # first pose, total energy

    out_pdbqt = work_dir / f"docked_{mol_id}.pdbqt"
    v.write_poses(str(out_pdbqt), n_poses=1, overwrite=True)

    return best_affinity


# ---------------------------------------------------------------------------
# Main Entry Point
# ---------------------------------------------------------------------------

def run_docking(manifest_path, candidates_dir, output_dir, mode="simulation",
                db_path=None, campaign_id=None):
    """Run docking with full telemetry."""
    manifest_path = Path(manifest_path).resolve()
    candidates_path = Path(candidates_dir).resolve()
    out_path = Path(output_dir).resolve()
    timestamp = datetime.now(timezone.utc).isoformat()

    if not manifest_path.exists():
        print(f"Error: Manifest {manifest_path} not found.")
        sys.exit(1)

    with open(manifest_path) as f:
        manifest = json.load(f)

    parameters = {**DEFAULT_PARAMS, "mode": mode}

    db = None
    run_id = None
    if db_path and campaign_id:
        db = TelemetryDB(db_path)
        run_id = db.start_run(
            campaign_id=campaign_id,
            module_name="04_docking",
            input_path=str(manifest_path),
            parameters=parameters,
        )

    try:
        if mode == "simulation":
            results_path = run_docking_simulation(manifest, candidates_path, out_path, parameters, db, run_id)
        elif mode == "triage":
            results_path = run_docking_triage(manifest, candidates_path, out_path, parameters, db, run_id)
        elif mode == "production":
            results_path = run_docking_production(manifest, candidates_path, out_path, parameters, db, run_id)
        else:
            raise ValueError(f"Unknown mode: {mode}")

        metadata = {
            "module": "04_docking", "timestamp": timestamp, "mode": mode,
            "manifest_path": str(manifest_path), "candidates_dir": str(candidates_path),
            "results_file": results_path, "parameters": parameters,
            "tdc_available": HAS_TDC, "status": "success",
        }
        with open(out_path / "run_metadata.json", "w") as f:
            json.dump(metadata, f, indent=2)

        if db and run_id:
            db.complete_run(run_id, "success", results_path)

        print(f"[Docking] Metadata written to {out_path / 'run_metadata.json'}")

    except Exception as e:
        error_msg = traceback.format_exc()
        print(f"[Docking] FAILED: {e}")

        if db and run_id:
            db.complete_run(run_id, "failed", error_trace=error_msg)

        out_path.mkdir(parents=True, exist_ok=True)
        with open(out_path / "run_metadata.json", "w") as f:
            json.dump({
                "module": "04_docking", "timestamp": timestamp, "mode": mode,
                "status": "failed", "error": str(e), "traceback": error_msg,
            }, f, indent=2)
        sys.exit(1)

    finally:
        if db:
            db.close()


def main():
    parser = argparse.ArgumentParser(description="Module 04: Docking (Vina / TDC Oracle)")
    parser.add_argument("--manifest", required=True, help="Input manifest")
    parser.add_argument("--candidates_dir", required=True, help="Candidates directory")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--mode", choices=["simulation", "triage", "production"],
                        default="simulation", help="Execution mode")
    parser.add_argument("--db_path", default=None, help="Telemetry DB path")
    parser.add_argument("--campaign_id", default=None, help="Campaign ID")
    args = parser.parse_args()

    run_docking(args.manifest, args.candidates_dir, args.output_dir,
                args.mode, args.db_path, args.campaign_id)


if __name__ == "__main__":
    main()
