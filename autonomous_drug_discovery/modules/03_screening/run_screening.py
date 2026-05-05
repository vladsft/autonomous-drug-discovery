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

# Try multiple import paths for MolScore descriptors — the package has
# reorganised modules across versions. We only promote BACKEND to "molscore"
# if we can actually instantiate a MolScore descriptor callable.
HAS_MOLSCORE = False
_molscore_descriptor_fn = None
try:
    # Newer MolScore (>=0.9) exposes RDKitDescriptors / MolecularDescriptors
    # under scoring_functions.descriptors.
    from molscore.scoring_functions import descriptors as _ms_desc  # type: ignore
    if hasattr(_ms_desc, "RDKitDescriptors"):
        _molscore_descriptor_fn = _ms_desc.RDKitDescriptors
    elif hasattr(_ms_desc, "MolecularDescriptors"):
        _molscore_descriptor_fn = _ms_desc.MolecularDescriptors
    if _molscore_descriptor_fn is not None:
        HAS_MOLSCORE = True
        BACKEND = "molscore"
except ImportError:
    pass

if BACKEND is None and HAS_RDKIT:
    BACKEND = "rdkit"

if BACKEND is None:
    raise ImportError(
        "Neither MolScore nor RDKit are available. Install at least one:\n"
        "  conda install -c conda-forge rdkit   (required baseline)\n"
        "  pip install molscore                  (optional, preferred)"
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
# Property computation + filter application (shared across backends)
# ---------------------------------------------------------------------------

def _pains_catalog():
    """Return a cached PAINS filter catalog, or None if unavailable."""
    try:
        from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        return FilterCatalog(params)
    except Exception:
        return None


def _compute_properties(mol, sa_scorer, pains_catalog, molscore_fn=None) -> dict:
    """Compute all descriptors + PAINS flag for a single molecule.

    If `molscore_fn` is provided (a MolScore descriptor callable), it is used
    as the source of truth for the RDKit descriptors it exposes. Missing values
    fall back to direct RDKit calls so the property dict is always complete.
    """
    props: dict = {}

    if molscore_fn is not None:
        try:
            ms_out = molscore_fn(mols=[mol], prefix="desc_")
            if isinstance(ms_out, list) and ms_out:
                # MolScore returns one dict per input molecule; values may be
                # prefixed (e.g. desc_MolWt). Normalise keys to our schema.
                for k, v in ms_out[0].items():
                    props[k if k.startswith("desc_") else f"desc_{k}"] = v
        except Exception:
            # Any MolScore failure → fall through to RDKit.
            props = {}

    # RDKit descriptors (fills in anything MolScore didn't provide).
    props.setdefault("desc_MolWt", round(Descriptors.MolWt(mol), 2))
    props.setdefault("desc_MolLogP", round(Descriptors.MolLogP(mol), 2))
    props.setdefault("desc_NumHBD", rdMolDescriptors.CalcNumHBD(mol))
    props.setdefault("desc_NumHBA", rdMolDescriptors.CalcNumHBA(mol))
    props.setdefault("desc_TPSA", round(Descriptors.TPSA(mol), 2))
    props.setdefault("desc_QED", round(QED.qed(mol), 4))
    props.setdefault("desc_NumRotatableBonds", Descriptors.NumRotatableBonds(mol))
    props.setdefault("desc_NumHeavyAtoms", mol.GetNumHeavyAtoms())
    props.setdefault("desc_SAScore", round(sa_scorer(mol), 2))

    if pains_catalog is not None:
        try:
            props["filter_PAINS"] = 1 if pains_catalog.HasMatch(mol) else 0
        except Exception:
            props["filter_PAINS"] = 0
    else:
        props["filter_PAINS"] = 0

    return props


def _violations(props: dict, thresholds: dict) -> list[str]:
    """Return a list of threshold violations. Empty list = passes all filters."""
    reasons = []
    for prop_name, spec in thresholds.items():
        val = props.get(prop_name)
        if val is None:
            continue
        if "max" in spec and val > spec["max"]:
            reasons.append(f"{prop_name}_exceeded (val={val}, max={spec['max']})")
        if "min" in spec and val < spec["min"]:
            reasons.append(f"{prop_name}_below (val={val}, min={spec['min']})")
        if "equals" in spec and val != spec["equals"]:
            reasons.append(f"{prop_name}_mismatch (val={val}, expected={spec['equals']})")
    return reasons


def _mol_id_from_supplier(mol, idx: int) -> str:
    """Preserve upstream mol_id if present in the SDF; otherwise assign by index.

    Pocket2Mol and TargetDiff tag molecules with `molecule_id` during
    consolidation; preserving it keeps telemetry joins across stages correct.
    """
    if mol is not None and mol.HasProp("molecule_id"):
        return mol.GetProp("molecule_id")
    return f"mol_{idx:04d}"


def _screen_molecules(mols, smiles_list, config, backend):
    """Apply property computation + threshold filters to all molecules.

    Returns a list of per-molecule result dicts. Each dict has keys:
    molecule_id, smiles, passed, eliminated_by, violations (full list), properties.
    """
    thresholds = config.get("filter_thresholds", {})
    sa_scorer = _get_sa_scorer()
    pains_catalog = _pains_catalog()

    # Instantiate MolScore descriptor object once if available.
    molscore_fn = None
    if backend == "molscore" and _molscore_descriptor_fn is not None:
        try:
            molscore_fn = _molscore_descriptor_fn(prefix="desc_")
        except Exception:
            molscore_fn = None

    results = []
    for i, (mol, smi) in enumerate(zip(mols, smiles_list)):
        mol_id = _mol_id_from_supplier(mol, i)

        if mol is None:
            results.append({
                "molecule_id": mol_id,
                "smiles": smi,
                "passed": False,
                "eliminated_by": "invalid_chemistry",
                "violations": ["invalid_chemistry"],
                "properties": {},
            })
            continue

        try:
            props = _compute_properties(mol, sa_scorer, pains_catalog, molscore_fn)
        except Exception as e:
            results.append({
                "molecule_id": mol_id,
                "smiles": smi,
                "passed": False,
                "eliminated_by": f"property_computation_error: {e}",
                "violations": [f"property_computation_error: {e}"],
                "properties": {},
            })
            continue

        violations = _violations(props, thresholds)
        results.append({
            "molecule_id": mol_id,
            "smiles": smi,
            "passed": not violations,
            "eliminated_by": violations[0] if violations else None,
            "violations": violations,
            "properties": props,
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

        results = _screen_molecules(mols, smiles_list, config, BACKEND)

        # Separate survivors and rejected
        survivors = [r for r in results if r["passed"]]
        rejected = [r for r in results if not r["passed"]]

        # Attrition analysis — count every violation, so a molecule that fails
        # three filters contributes to three attrition buckets.
        attrition: dict[str, int] = {}
        for r in rejected:
            for v in (r["violations"] or ["unknown"]):
                category = v.split("_")[0] if "_" in v else v
                attrition[category] = attrition.get(category, 0) + 1

        # Write filtered SDF — zip results with the input mol list by position.
        output_sdf = out_path / "screened_molecules.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        for r, mol in zip(results, mols):
            if r["passed"] and mol is not None:
                # Ensure the preserved molecule_id is stamped into the SDF so
                # downstream stages (docking, telemetry) see the same key.
                mol.SetProp("molecule_id", r["molecule_id"])
                writer.write(mol)
        writer.close()

        # ADMET-AI enrichment for survivors.
        # Run in small chunks so one poisoned SMILES can't sink the whole batch.
        admet_results: dict[str, dict] = {}
        admet_errors: dict[str, str] = {}
        ADMET_KEY_PROPS = [
            "hERG", "AMES", "DILI", "CYP2D6_Veith", "CYP3A4_Veith",
            "Caco2_Wang", "HIA_Hou", "BBB_Martins", "Clearance_Hepatocyte_AZ",
            "Bioavailability_Ma", "LD50_Zhu",
        ]
        try:
            from admet_ai import ADMETModel
            admet_model = ADMETModel()
            survivor_smiles = [r["smiles"] for r in survivors if r["smiles"]]
            if survivor_smiles:
                print(f"[Screening] Running ADMET-AI on {len(survivor_smiles)} survivors...")
                chunk_size = 32
                for start in range(0, len(survivor_smiles), chunk_size):
                    chunk = survivor_smiles[start:start + chunk_size]
                    try:
                        preds = admet_model.predict(smiles=chunk)
                        for i, smi in enumerate(chunk):
                            row = preds.iloc[i]
                            admet_results[smi] = {
                                k: round(float(v), 4) if isinstance(v, (float, int)) else v
                                for k, v in row.items()
                            }
                    except Exception as chunk_err:
                        # Retry one-by-one so we can pinpoint and skip the bad SMILES.
                        for smi in chunk:
                            try:
                                single = admet_model.predict(smiles=[smi])
                                row = single.iloc[0]
                                admet_results[smi] = {
                                    k: round(float(v), 4) if isinstance(v, (float, int)) else v
                                    for k, v in row.items()
                                }
                            except Exception as single_err:
                                admet_errors[smi] = str(single_err)
                        print(f"[Screening] ADMET-AI: chunk starting at {start} failed ({chunk_err}); "
                              f"fell back to per-molecule, {len(admet_errors)} failures so far")

                for r in survivors:
                    smi = r["smiles"]
                    if smi in admet_results:
                        for prop in ADMET_KEY_PROPS:
                            if prop in admet_results[smi]:
                                r["properties"][f"admet_{prop}"] = admet_results[smi][prop]
                    elif smi in admet_errors:
                        r["properties"]["admet_error"] = admet_errors[smi]
                print(f"[Screening] ADMET-AI: {len(admet_results)} annotated, "
                      f"{len(admet_errors)} failed")
        except ImportError:
            print("[Screening] ADMET-AI not installed, skipping ADMET predictions")

        # Build report
        report = {
            "timestamp": timestamp,
            "backend": BACKEND,
            "admet_ai": len(admet_results) > 0,
            "admet_failures": len(admet_errors),
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
            "rejections": [
                {"molecule_id": r["molecule_id"], "smiles": r["smiles"],
                 "violations": r["violations"]}
                for r in rejected
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
                # Join the full violations list so expert tuning sees every
                # filter that fired, not only the first.
                reasons = r["violations"] or ([] if r["passed"] else ["unknown"])
                mol_batch.append({
                    "molecule_id": r["molecule_id"],
                    "smiles": r["smiles"],
                    "qed": r["properties"].get("desc_QED"),
                    "sa_score": r["properties"].get("desc_SAScore"),
                    "logp": r["properties"].get("desc_MolLogP"),
                    "mol_weight": r["properties"].get("desc_MolWt"),
                    "passed_triage": 1 if r["passed"] else 0,
                    "stage_eliminated": "; ".join(reasons) if reasons else None,
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
