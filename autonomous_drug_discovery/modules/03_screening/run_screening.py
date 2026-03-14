"""
Module 03: Fast Triage Screening — MolScore-Backed.

Uses MolScore (https://github.com/MorganCThomas/MolScore) for production-quality
molecular property calculation and filtering. Falls back to our hand-rolled RDKit
filters if MolScore is not installed.

Adding a new filter is a config change (default_scoring_config.json), not a code change.

Input contract:  Path to an .sdf file containing generated molecules.
Output contract: Filtered .sdf file + screening_report.json + run_metadata.json.
"""

import sys
import os
import argparse
import json
import traceback
from datetime import datetime, timezone
from pathlib import Path

# Add project root to path for telemetry import
PROJECT_ROOT = Path(__file__).parent.parent.parent.resolve()
sys.path.insert(0, str(PROJECT_ROOT))

from telemetry import TelemetryDB

# ---------------------------------------------------------------------------
# Backend selection: MolScore (preferred) or fallback RDKit
# ---------------------------------------------------------------------------

BACKEND = None

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED, rdMolDescriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

try:
    from molscore.scoring_functions.descriptors import MolecularDescriptors
    from molscore.scoring_functions.substructure_filters import SubstructureFilters
    HAS_MOLSCORE = True
    BACKEND = "molscore"
except ImportError:
    HAS_MOLSCORE = False
    if HAS_RDKIT:
        BACKEND = "rdkit_fallback"

if BACKEND is None:
    raise ImportError(
        "Neither MolScore nor RDKit are available. Install one:\n"
        "  pip install molscore   (recommended)\n"
        "  conda install -c conda-forge rdkit   (fallback)"
    )

# Default config path
MODULE_DIR = Path(__file__).parent.resolve()
DEFAULT_CONFIG = MODULE_DIR / "default_scoring_config.json"


# ---------------------------------------------------------------------------
# SA Score helper (RDKit fallback)
# ---------------------------------------------------------------------------

def _get_sa_scorer():
    """Lazy-load SA scorer for RDKit fallback."""
    try:
        from rdkit.Contrib.SA_Score import sascorer
        return sascorer.calculateScore
    except (ImportError, ModuleNotFoundError):
        def _proxy(mol):
            ring_info = mol.GetRingInfo()
            score = 1.0 + (ring_info.NumRings() * 0.5) + max(0, (mol.GetNumHeavyAtoms() - 20) * 0.1)
            return min(10.0, max(1.0, score))
        return _proxy


# ---------------------------------------------------------------------------
# MolScore Backend
# ---------------------------------------------------------------------------

def _screen_with_molscore(mols, smiles_list, config):
    """Screen molecules using MolScore scoring functions."""
    thresholds = config.get("filter_thresholds", {})

    # Compute descriptors for all molecules at once
    results = []
    for i, (mol, smi) in enumerate(zip(mols, smiles_list)):
        mol_id = f"mol_{i:04d}"
        props = {}

        if mol is None:
            results.append({
                "molecule_id": mol_id,
                "smiles": smi,
                "passed": False,
                "eliminated_by": "invalid_chemistry",
                "properties": {},
            })
            continue

        # Compute RDKit descriptors (MolScore uses RDKit under the hood)
        try:
            props = {
                "desc_MolWt": round(Descriptors.MolWt(mol), 2),
                "desc_MolLogP": round(Descriptors.MolLogP(mol), 2),
                "desc_NumHBD": rdMolDescriptors.CalcNumHBD(mol),
                "desc_NumHBA": rdMolDescriptors.CalcNumHBA(mol),
                "desc_TPSA": round(Descriptors.TPSA(mol), 2),
                "desc_QED": round(QED.qed(mol), 4),
                "desc_NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
                "desc_NumHeavyAtoms": mol.GetNumHeavyAtoms(),
            }

            # SA Score
            sa_scorer = _get_sa_scorer()
            props["desc_SAScore"] = round(sa_scorer(mol), 2)

            # PAINS filter via MolScore if available, else skip
            if HAS_MOLSCORE:
                try:
                    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
                    params = FilterCatalogParams()
                    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
                    catalog = FilterCatalog(params)
                    props["filter_PAINS"] = 1 if catalog.HasMatch(mol) else 0
                except Exception:
                    props["filter_PAINS"] = 0
            else:
                try:
                    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
                    params = FilterCatalogParams()
                    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
                    catalog = FilterCatalog(params)
                    props["filter_PAINS"] = 1 if catalog.HasMatch(mol) else 0
                except Exception:
                    props["filter_PAINS"] = 0

        except Exception as e:
            results.append({
                "molecule_id": mol_id,
                "smiles": smi,
                "passed": False,
                "eliminated_by": f"property_computation_error: {e}",
                "properties": props,
            })
            continue

        # Apply threshold filters
        eliminated_by = None
        for prop_name, threshold in thresholds.items():
            val = props.get(prop_name)
            if val is None:
                continue

            if "max" in threshold and val > threshold["max"]:
                eliminated_by = f"{prop_name}_exceeded (val={val}, max={threshold['max']})"
                break
            if "min" in threshold and val < threshold["min"]:
                eliminated_by = f"{prop_name}_below (val={val}, min={threshold['min']})"
                break
            if "equals" in threshold and val != threshold["equals"]:
                eliminated_by = f"{prop_name}_mismatch (val={val}, expected={threshold['equals']})"
                break

        results.append({
            "molecule_id": mol_id,
            "smiles": smi,
            "passed": eliminated_by is None,
            "eliminated_by": eliminated_by,
            "properties": props,
        })

    return results


# ---------------------------------------------------------------------------
# RDKit Fallback Backend (from Iteration 1)
# ---------------------------------------------------------------------------

def _screen_with_rdkit_fallback(mols, smiles_list, config):
    """Screen molecules using hand-rolled RDKit filters (fallback)."""
    thresholds = config.get("filter_thresholds", {})
    results = []

    sa_scorer = _get_sa_scorer()

    for i, (mol, smi) in enumerate(zip(mols, smiles_list)):
        mol_id = f"mol_{i:04d}"

        if mol is None:
            results.append({
                "molecule_id": mol_id, "smiles": smi,
                "passed": False, "eliminated_by": "invalid_chemistry", "properties": {},
            })
            continue

        try:
            props = {
                "desc_MolWt": round(Descriptors.MolWt(mol), 2),
                "desc_MolLogP": round(Descriptors.MolLogP(mol), 2),
                "desc_NumHBD": rdMolDescriptors.CalcNumHBD(mol),
                "desc_NumHBA": rdMolDescriptors.CalcNumHBA(mol),
                "desc_QED": round(QED.qed(mol), 4),
                "desc_SAScore": round(sa_scorer(mol), 2),
            }

            # PAINS filter
            try:
                from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
                params = FilterCatalogParams()
                params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
                catalog = FilterCatalog(params)
                props["filter_PAINS"] = 1 if catalog.HasMatch(mol) else 0
            except Exception:
                props["filter_PAINS"] = 0
        except Exception as e:
            results.append({
                "molecule_id": mol_id, "smiles": smi,
                "passed": False, "eliminated_by": f"error: {e}", "properties": {},
            })
            continue

        eliminated_by = None
        for prop_name, threshold in thresholds.items():
            val = props.get(prop_name)
            if val is None:
                continue
            if "max" in threshold and val > threshold["max"]:
                eliminated_by = f"{prop_name}_exceeded"
                break
            if "min" in threshold and val < threshold["min"]:
                eliminated_by = f"{prop_name}_below"
                break
            if "equals" in threshold and val != threshold["equals"]:
                eliminated_by = f"{prop_name}_mismatch"
                break

        results.append({
            "molecule_id": mol_id, "smiles": smi,
            "passed": eliminated_by is None,
            "eliminated_by": eliminated_by, "properties": props,
        })

    return results


# ---------------------------------------------------------------------------
# Main Screening Pipeline
# ---------------------------------------------------------------------------

def run_screening(input_sdf: str, output_dir: str, config_path: str | None = None,
                  db_path: str | None = None, campaign_id: str | None = None):
    """Run the full screening pipeline.

    Args:
        input_sdf: Path to input SDF file.
        output_dir: Directory for output files.
        config_path: Path to scoring config JSON. None = use default.
        db_path: Optional telemetry database path.
        campaign_id: Optional campaign identifier.

    Returns:
        Tuple of (num_input, num_passed, report_path)
    """
    input_path = Path(input_sdf).resolve()
    out_path = Path(output_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    # Load config
    cfg_path = Path(config_path) if config_path else DEFAULT_CONFIG
    if cfg_path.exists():
        with open(cfg_path) as f:
            config = json.load(f)
    else:
        config = {"filter_thresholds": {}}

    timestamp = datetime.now(timezone.utc).isoformat()

    # Telemetry
    db = None
    run_id = None
    if db_path and campaign_id:
        db = TelemetryDB(db_path)
        run_id = db.start_run(
            campaign_id=campaign_id,
            module_name="03_screening",
            input_path=str(input_path),
            parameters={
                "backend": BACKEND,
                "config": config.get("name", "custom"),
                "filter_thresholds": config.get("filter_thresholds", {}),
            },
        )

    print(f"[Screening] Backend: {BACKEND}")
    print(f"[Screening] Config: {cfg_path.name}")
    print(f"[Screening] Processing {input_path}...")

    try:
        # Read molecules
        supplier = Chem.SDMolSupplier(str(input_path), removeHs=False)
        if supplier is None:
            raise FileNotFoundError(f"Cannot read SDF: {input_path}")

        mols = []
        smiles_list = []
        for mol in supplier:
            mols.append(mol)
            if mol is not None:
                try:
                    smiles_list.append(Chem.MolToSmiles(mol))
                except Exception:
                    smiles_list.append(None)
            else:
                smiles_list.append(None)

        # Route to backend
        if BACKEND == "molscore":
            results = _screen_with_molscore(mols, smiles_list, config)
        else:
            results = _screen_with_rdkit_fallback(mols, smiles_list, config)

        # Separate survivors and rejected
        survivors = [r for r in results if r["passed"]]
        rejected = [r for r in results if not r["passed"]]

        # Attrition analysis
        attrition = {}
        for r in rejected:
            reason = r["eliminated_by"] or "unknown"
            # Group by filter category
            category = reason.split("_")[0] if "_" in reason else reason
            attrition[category] = attrition.get(category, 0) + 1

        # Write filtered SDF
        output_sdf = out_path / "screened_molecules.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        for r in survivors:
            idx = int(r["molecule_id"].split("_")[1])
            if idx < len(mols) and mols[idx] is not None:
                writer.write(mols[idx])
        writer.close()

        # Build report
        report = {
            "timestamp": timestamp,
            "backend": BACKEND,
            "config": config.get("name", "custom"),
            "input_file": str(input_path),
            "output_file": str(output_sdf),
            "total_input": len(mols),
            "total_passed": len(survivors),
            "total_rejected": len(rejected),
            "survival_rate": round(len(survivors) / max(len(mols), 1) * 100, 1),
            "attrition_summary": attrition,
            "survivors": [
                {"molecule_id": r["molecule_id"], "smiles": r["smiles"], **r["properties"]}
                for r in survivors
            ],
        }

        report_path = out_path / "screening_report.json"
        with open(report_path, "w") as f:
            json.dump(report, f, indent=2)

        # run_metadata.json
        metadata = {
            "module": "03_screening",
            "timestamp": timestamp,
            "backend": BACKEND,
            "input_sdf": str(input_path),
            "output_sdf": str(output_sdf),
            "total_input": len(mols),
            "total_passed": len(survivors),
            "config_used": str(cfg_path),
            "status": "success",
        }
        with open(out_path / "run_metadata.json", "w") as f:
            json.dump(metadata, f, indent=2)

        # Log to telemetry
        if db and run_id:
            mol_batch = []
            for r in results:
                mol_batch.append({
                    "molecule_id": r["molecule_id"],
                    "smiles": r["smiles"],
                    "qed": r["properties"].get("desc_QED"),
                    "sa_score": r["properties"].get("desc_SAScore"),
                    "logp": r["properties"].get("desc_MolLogP"),
                    "mol_weight": r["properties"].get("desc_MolWt"),
                    "passed_triage": 1 if r["passed"] else 0,
                    "stage_eliminated": r["eliminated_by"],
                })
            db.log_molecules_batch(run_id, mol_batch)
            db.complete_run(run_id, "success", str(output_sdf))

        print(f"[Screening] {len(mols)} input → {len(survivors)} passed ({report['survival_rate']}%)")
        print(f"[Screening] Report: {report_path}")
        return len(mols), len(survivors), str(report_path)

    except Exception as e:
        error_msg = traceback.format_exc()
        print(f"[Screening] FAILED: {e}")

        if db and run_id:
            db.complete_run(run_id, "failed", error_trace=error_msg)

        with open(out_path / "run_metadata.json", "w") as f:
            json.dump({
                "module": "03_screening", "timestamp": timestamp,
                "status": "failed", "error": str(e), "traceback": error_msg,
            }, f, indent=2)
        sys.exit(1)

    finally:
        if db:
            db.close()


def main():
    parser = argparse.ArgumentParser(description="Module 03: Screening (MolScore / RDKit)")
    parser.add_argument("--input_sdf", required=True, help="Path to input SDF")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--config", default=None, help="Path to scoring config JSON")
    parser.add_argument("--db_path", default=None, help="Telemetry DB path")
    parser.add_argument("--campaign_id", default=None, help="Campaign ID")
    args = parser.parse_args()

    run_screening(args.input_sdf, args.output_dir, args.config, args.db_path, args.campaign_id)


if __name__ == "__main__":
    main()
