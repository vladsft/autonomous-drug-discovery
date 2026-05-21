"""
Retrospective docking validation.

Before any de-novo docking score can be trusted, the docking setup must be
shown to work on a known answer. This module does that two ways:

  1. Re-docking (always): take a holo structure with a co-crystallised ligand,
     dock that ligand back into its own pocket, and measure the RMSD between
     the top docked pose and the crystal pose. < 2.0 A is the accepted bar for
     "the protocol can reproduce a known binding mode".

  2. Decoy enrichment (optional, --decoys): dock the true active alongside a
     set of decoy molecules and report where the active ranks. If the active
     does not beat most decoys, docking scores are not separating binders from
     non-binders and downstream rankings are noise.

Input contract:  a holo .pdb (protein + HETATM ligand) + the ligand resname.
Output contract: validation_report.json (+ run_metadata.json).
"""

import sys
import json
import argparse
import traceback
from datetime import datetime, timezone
from pathlib import Path

MODULE_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = MODULE_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(MODULE_DIR))

from telemetry import TelemetryDB  # noqa: E402
import run_docking as rd  # noqa: E402

RMSD_PASS_THRESHOLD = 2.0  # Angstrom — standard re-docking success criterion


def _split_holo_pdb(holo_pdb: Path, ligand_resname: str, out_dir: Path):
    """Split a holo PDB into a receptor PDB and the ligand's HETATM block.

    Returns (receptor_pdb_path, ligand_pdb_block). The receptor keeps protein
    ATOM records (and any other HETATM, e.g. metals); only the named ligand is
    pulled out.
    """
    receptor_lines, ligand_lines = [], []
    resn = ligand_resname.strip().upper()
    with open(holo_pdb) as f:
        for line in f:
            rec = line[:6]
            if rec in ("ATOM  ", "HETATM"):
                this_resn = line[17:20].strip().upper()
                if rec == "HETATM" and this_resn == resn:
                    ligand_lines.append(line)
                else:
                    receptor_lines.append(line)
            elif rec in ("TER   ", "END   ", "ENDMDL"):
                receptor_lines.append(line)
    if not ligand_lines:
        raise ValueError(
            f"No HETATM records for ligand '{resn}' found in {holo_pdb.name}."
        )
    receptor_pdb = out_dir / f"{holo_pdb.stem}_receptor.pdb"
    receptor_pdb.write_text("".join(receptor_lines) + "END\n")
    return receptor_pdb, "".join(ligand_lines)


def _crystal_ligand_mol(ligand_pdb_block: str, ligand_smiles: str | None):
    """Build an RDKit mol of the crystal ligand pose.

    When `ligand_smiles` is given, bond orders are assigned from it (PDB bond
    perception is unreliable); otherwise RDKit's distance-based perception is
    used, which is good enough for heavy-atom RMSD.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromPDBBlock(ligand_pdb_block, removeHs=True, sanitize=True)
    if mol is None:
        raise ValueError("RDKit could not parse the crystal ligand block.")
    if ligand_smiles:
        template = Chem.MolFromSmiles(ligand_smiles)
        if template is not None:
            try:
                mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
            except Exception as e:
                print(f"[Validate] Bond-order assignment from SMILES failed "
                      f"({e}); using perceived bonds.")
    return mol


def _centroid_and_extent(mol):
    """Return (centroid[xyz], box_size[xyz]) for the molecule's conformer."""
    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    zs = [conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())]
    centroid = [round(sum(xs) / len(xs), 3),
                round(sum(ys) / len(ys), 3),
                round(sum(zs) / len(zs), 3)]
    box = [round(min(30.0, max(16.0, (max(a) - min(a)) + 10.0)), 1)
           for a in (xs, ys, zs)]
    return centroid, box


def _docked_pose_mol(docked_pdbqt: Path):
    """Read the top pose from a Vina output PDBQT into an RDKit mol (no Hs)."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem

    pdbqt_mol = PDBQTMolecule.from_file(str(docked_pdbqt), skip_typing=True)
    rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    mol = next((m for m in rdkit_mols if m is not None), None)
    if mol is None:
        raise ValueError(f"Could not reconstruct an RDKit mol from {docked_pdbqt}.")
    return Chem.RemoveHs(mol)


def _pose_rmsd(docked_mol, crystal_mol) -> float:
    """Symmetry-corrected heavy-atom RMSD between docked and crystal poses.

    Both poses are already in the receptor coordinate frame, so no alignment is
    applied — this is the true positional error, not a shape comparison.
    """
    from rdkit.Chem import rdMolAlign, Chem

    probe = Chem.RemoveHs(docked_mol)
    ref = Chem.RemoveHs(crystal_mol)
    return rdMolAlign.CalcRMS(probe, ref)


def validate_redocking(holo_pdb, ligand_resname, output_dir,
                       ligand_smiles=None, exhaustiveness=16, decoys=None,
                       db_path=None, campaign_id=None):
    """Re-dock a co-crystallised ligand and report pose RMSD (+ decoy enrichment).

    Returns the validation report dict.
    """
    holo_pdb = Path(holo_pdb).resolve()
    out_path = Path(output_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).isoformat()

    if not (rd.HAS_VINA and rd.HAS_MEEKO and rd.HAS_RDKIT):
        raise ImportError("Validation needs vina + meeko + rdkit in this env.")

    db = run_id = None
    if db_path and campaign_id:
        db = TelemetryDB(db_path)
        run_id = db.start_run(campaign_id, "04_docking_validate",
                              str(holo_pdb),
                              {"ligand_resname": ligand_resname,
                               "exhaustiveness": exhaustiveness})

    try:
        from vina import Vina

        # 1. Split holo into receptor + crystal ligand.
        receptor_pdb, ligand_block = _split_holo_pdb(holo_pdb, ligand_resname, out_path)
        crystal_mol = _crystal_ligand_mol(ligand_block, ligand_smiles)
        centroid, box = _centroid_and_extent(crystal_mol)
        print(f"[Validate] Crystal ligand '{ligand_resname}': "
              f"{crystal_mol.GetNumHeavyAtoms()} heavy atoms, "
              f"box centre {centroid}, size {box}")

        # 2. Prepare receptor + Vina maps once.
        receptor_pdbqt = out_path / f"{holo_pdb.stem}_receptor.pdbqt"
        n_atoms = rd._prepare_receptor_pdbqt(receptor_pdb, receptor_pdbqt)
        print(f"[Validate] Receptor PDBQT: {n_atoms} atoms.")
        v = Vina(sf_name="vina")
        v.set_receptor(str(receptor_pdbqt))
        v.compute_vina_maps(center=centroid, box_size=box)

        # 3. Re-dock the crystal ligand and score pose RMSD.
        affinity = rd._dock_one_ligand(v, crystal_mol, exhaustiveness, 9,
                                       out_path, "redock_active")
        docked_mol = _docked_pose_mol(out_path / "docked_redock_active.pdbqt")
        rmsd = _pose_rmsd(docked_mol, crystal_mol)
        passed = rmsd < RMSD_PASS_THRESHOLD
        print(f"[Validate] Re-docking: affinity {affinity:.2f} kcal/mol, "
              f"pose RMSD {rmsd:.2f} A → {'PASS' if passed else 'FAIL'} "
              f"(threshold {RMSD_PASS_THRESHOLD} A)")

        report = {
            "timestamp": timestamp,
            "holo_pdb": str(holo_pdb),
            "ligand_resname": ligand_resname,
            "redocking": {
                "affinity_kcal_mol": round(affinity, 3),
                "pose_rmsd_angstrom": round(rmsd, 3),
                "rmsd_threshold": RMSD_PASS_THRESHOLD,
                "passed": passed,
            },
        }

        # 4. Optional decoy enrichment.
        if decoys:
            report["decoy_enrichment"] = _decoy_enrichment(
                v, crystal_mol, affinity, decoys, exhaustiveness, out_path)

        report_path = out_path / "validation_report.json"
        report_path.write_text(json.dumps(report, indent=2))
        with open(out_path / "run_metadata.json", "w") as f:
            json.dump({"module": "04_docking_validate", "timestamp": timestamp,
                       "status": "success", "report": str(report_path)}, f, indent=2)
        if db and run_id:
            db.complete_run(run_id, "success", str(report_path))
        print(f"[Validate] Report: {report_path}")
        return report

    except Exception as e:
        if db and run_id:
            db.complete_run(run_id, "failed", error_trace=traceback.format_exc())
        raise
    finally:
        if db:
            db.close()


def _decoy_enrichment(v, active_mol, active_affinity, decoy_file,
                      exhaustiveness, out_path):
    """Dock decoy SMILES against the same maps; report the active's rank.

    A trustworthy docking setup ranks the true active above most decoys. The
    active's percentile (fraction of all molecules it beats) is the headline
    number — > 0.9 is healthy, near 0.5 means docking is not discriminating.
    """
    from rdkit import Chem

    smiles = [s.strip() for s in Path(decoy_file).read_text().splitlines()
              if s.strip() and not s.startswith("#")]
    decoy_scores = []
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        try:
            aff = rd._dock_one_ligand(v, mol, exhaustiveness, 9,
                                      out_path, f"decoy_{i:03d}")
            decoy_scores.append(aff)
        except Exception as e:
            print(f"[Validate] decoy {i} ({smi}) failed: {e}")
    if not decoy_scores:
        return {"error": "no decoys docked successfully"}

    # Lower (more negative) affinity = better. The active beats a decoy when
    # its affinity is more negative.
    beaten = sum(1 for d in decoy_scores if active_affinity < d)
    percentile = round(beaten / len(decoy_scores), 3)
    return {
        "n_decoys": len(decoy_scores),
        "active_affinity": round(active_affinity, 3),
        "decoy_affinity_mean": round(sum(decoy_scores) / len(decoy_scores), 3),
        "active_beats_fraction": percentile,
        "healthy": percentile >= 0.9,
    }


def main():
    parser = argparse.ArgumentParser(description="Retrospective docking validation")
    parser.add_argument("--holo_pdb", required=True,
                        help="Holo PDB: protein + co-crystallised ligand")
    parser.add_argument("--ligand_resname", required=True,
                        help="3-letter HETATM residue name of the ligand")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--ligand_smiles", default=None,
                        help="SMILES of the ligand — used to assign correct "
                             "bond orders to the crystal pose (recommended)")
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--decoys", default=None,
                        help="Optional text file of decoy SMILES, one per line")
    parser.add_argument("--db_path", default=None)
    parser.add_argument("--campaign_id", default=None)
    args = parser.parse_args()

    try:
        report = validate_redocking(
            args.holo_pdb, args.ligand_resname, args.output_dir,
            ligand_smiles=args.ligand_smiles, exhaustiveness=args.exhaustiveness,
            decoys=args.decoys, db_path=args.db_path, campaign_id=args.campaign_id)
    except Exception as e:
        print(f"[Validate] FAILED: {e}")
        sys.exit(1)
    # Non-zero exit when re-docking fails its RMSD bar, so CI / scripts can gate.
    sys.exit(0 if report["redocking"]["passed"] else 2)


if __name__ == "__main__":
    main()
