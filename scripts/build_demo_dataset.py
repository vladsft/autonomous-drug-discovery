"""Build a multi-backend demo dataset for the chemist dashboard.

Reads one campaign per backend (RDKit, Pocket2Mol, …), runs AiZynthFinder on
the top-N by docking score, pre-renders 2D structures, and writes a single
JSON + JS bridge keyed by backend so the dashboard can toggle between them.

Usage:
  python scripts/build_demo_dataset.py \
      --backend rdkit:campaign_23355c05:.../aizynth_top10.json \
      --backend pocket2mol:campaign_xxxxxxxx:.../aizynth_p2m_top10.json \
      --output_json data/outputs/professor_demo.json \
      --output_js   dashboard/professor_demo.js
"""

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN_BASE = PROJECT_ROOT / "autonomous_drug_discovery" / "data"


def render_svg(smiles: str, w: int = 280, h: int = 200) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    AllChem.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 1.4
    opts.padding = 0.06
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().replace(
        "<?xml version='1.0' encoding='iso-8859-1'?>\n", ""
    )


def annotate_tree_with_svgs(node, w=100, h=70):
    if not node:
        return
    if node.get("type") == "mol" and node.get("smiles"):
        node["svg"] = render_svg(node["smiles"], w=w, h=h)
    for child in node.get("children") or []:
        annotate_tree_with_svgs(child, w=w, h=h)


ADMET_KEYS = {
    "hERG": "admet_hERG", "AMES": "admet_AMES", "DILI": "admet_DILI",
    "CYP2D6": "admet_CYP2D6_Veith", "CYP3A4": "admet_CYP3A4_Veith",
    "Caco2": "admet_Caco2_Wang", "HIA": "admet_HIA_Hou", "BBB": "admet_BBB_Martins",
    "Clearance": "admet_Clearance_Hepatocyte_AZ",
    "Bioavailability": "admet_Bioavailability_Ma", "LD50": "admet_LD50_Zhu",
}

ADMET_PASS_RULES = {
    "hERG": ("max", 0.5),
    "AMES": ("max", 0.5),
    "DILI": ("max", 0.5),
}


def build_backend_dataset(backend_name: str, campaign_id: str, aizynth_json: Path,
                          max_molecules: int | None = None) -> dict:
    """Assemble per-backend dataset from a campaign + AiZynth output."""
    campaign_dir = CAMPAIGN_BASE / campaign_id
    screening_path = campaign_dir / "screened" / "screening_report.json"
    docking_path = campaign_dir / "results" / "docking_results.csv"

    if not screening_path.exists():
        raise FileNotFoundError(f"Screening report missing: {screening_path}")
    if not docking_path.exists():
        raise FileNotFoundError(f"Docking results missing: {docking_path}")

    with open(screening_path) as f:
        screening = json.load(f)
    with open(docking_path) as f:
        docking = list(csv.DictReader(f))
    with open(aizynth_json) as f:
        aizynth = json.load(f)

    survivors_by_id = {m["molecule_id"]: m for m in screening["survivors"]}
    rejections_by_id = {m["molecule_id"]: m for m in screening.get("rejections", [])}
    docking_by_id = {r["ligand_id"]: r for r in docking}
    aizynth_by_id = {r["ligand_id"]: r for r in aizynth}

    all_ids = sorted(set(survivors_by_id) | set(rejections_by_id) | set(docking_by_id))

    molecules = []
    for mid in all_ids:
        surv = survivors_by_id.get(mid)
        rej = rejections_by_id.get(mid)
        src = surv or rej or {}
        smiles = src.get("smiles") or docking_by_id.get(mid, {}).get("smiles", "")

        record = {
            "molecule_id": mid,
            "smiles": smiles,
            "svg": render_svg(smiles) if smiles else "",
            "screening_passed": bool(surv),
            "rejected_reason": rej.get("reason") if rej else None,
            "properties": {
                "mol_weight": src.get("desc_MolWt"),
                "logp": src.get("desc_MolLogP"),
                "qed": src.get("desc_QED"),
                "sa_score": src.get("desc_SAScore"),
                "tpsa": src.get("desc_TPSA"),
                "hbd": src.get("desc_NumHBD"),
                "hba": src.get("desc_NumHBA"),
                "rotatable_bonds": src.get("desc_NumRotatableBonds"),
                "heavy_atoms": src.get("desc_NumHeavyAtoms"),
                "pains_alerts": src.get("filter_PAINS"),
            },
            "admet": {label: src.get(key) for label, key in ADMET_KEYS.items()},
            "admet_flags": {},
            "docking_score": (
                float(docking_by_id[mid]["affinity"]) if mid in docking_by_id else None
            ),
        }

        for label, (kind, thresh) in ADMET_PASS_RULES.items():
            v = record["admet"].get(label)
            if v is None:
                record["admet_flags"][label] = None
            elif kind == "max":
                record["admet_flags"][label] = bool(v <= thresh)
            else:
                record["admet_flags"][label] = bool(v >= thresh)

        az = aizynth_by_id.get(mid)
        if az:
            tree = az.get("route_tree")
            annotate_tree_with_svgs(tree)
            record["synthesis"] = {
                "evaluated": True,
                "route_found": az.get("route_found", False),
                "n_steps": az.get("n_steps", 0),
                "n_routes": az.get("n_routes", 0),
                "precursors": az.get("precursors", []),
                "route_tree": tree,
                "error": az.get("error"),
            }
        else:
            record["synthesis"] = {
                "evaluated": False, "route_found": None, "n_steps": None,
                "n_routes": None, "precursors": [], "route_tree": None,
            }

        molecules.append(record)

    docking_scores = [m["docking_score"] for m in molecules if m["docking_score"] is not None]
    best = min(docking_scores) if docking_scores else 0
    worst = max(docking_scores) if docking_scores else 0

    for m in molecules:
        ds = m["docking_score"]
        dock_norm = 0.0
        if ds is not None and worst != best:
            dock_norm = (worst - ds) / (worst - best)
        flags = [v for v in m["admet_flags"].values() if v is not None]
        admet_norm = sum(1 for v in flags if v) / len(flags) if flags else 0.5
        synth = m["synthesis"]
        synth_norm = (1.0 if synth["route_found"] else 0.0) if synth["evaluated"] else 0.5
        m["composite_score"] = round(
            0.5 * dock_norm + 0.3 * admet_norm + 0.2 * synth_norm, 4
        )

    molecules.sort(key=lambda m: m["composite_score"], reverse=True)

    if max_molecules is not None and len(molecules) > max_molecules:
        molecules = molecules[:max_molecules]

    for i, m in enumerate(molecules, start=1):
        m["final_rank"] = i

    summary = {
        "total_generated": len(molecules),
        "passed_screening": sum(1 for m in molecules if m["screening_passed"]),
        "synthesis_evaluated": sum(1 for m in molecules if m["synthesis"]["evaluated"]),
        "synthesis_routes_found": sum(
            1 for m in molecules if m["synthesis"].get("route_found")
        ),
        "best_docking_score": min(docking_scores) if docking_scores else None,
        "median_docking_score": (
            sorted(docking_scores)[len(docking_scores) // 2] if docking_scores else None
        ),
        "mean_qed": round(
            sum(m["properties"]["qed"] for m in molecules
                if m["properties"]["qed"] is not None)
            / max(1, sum(1 for m in molecules if m["properties"]["qed"] is not None)),
            3,
        ),
        "mean_sa": round(
            sum(m["properties"]["sa_score"] for m in molecules
                if m["properties"]["sa_score"] is not None)
            / max(1, sum(1 for m in molecules if m["properties"]["sa_score"] is not None)),
            3,
        ),
    }

    return {
        "backend": backend_name,
        "campaign_id": campaign_id,
        "screening_backend": screening.get("backend", "molscore"),
        "summary": summary,
        "molecules": molecules,
    }


def parse_backend_arg(spec: str) -> tuple[str, str, Path]:
    """Parse a `name:campaign_id:aizynth_path` triple."""
    parts = spec.split(":", 2)
    if len(parts) != 3:
        raise ValueError(f"Bad --backend spec '{spec}'; expected name:campaign_id:aizynth_path")
    return parts[0], parts[1], Path(parts[2])


GENERATOR_DESCRIPTIONS = {
    "rdkit": (
        "RDKit fragment-based combinatorial generator. Assembles drug-like "
        "scaffolds + functional-group substituents + linkers into 3D-embedded, "
        "MMFF-optimised molecules. Pocket-aware only via molecule-size targeting "
        "(uses the pocket radius to pick a heavy-atom range). Runs on CPU in "
        "milliseconds per molecule. Strong baseline for diversity, but the "
        "geometry it produces ignores the 3D shape of the binding site."
    ),
    "pocket2mol": (
        "Pocket2Mol autoregressive pocket-conditioned generator (graph neural "
        "network, ICML 2022). Builds molecules atom-by-atom inside the binding "
        "pocket, conditioning each step on the protein context. Runs on GPU "
        "(~7 sec/molecule) and produces molecules whose 3D shape is explicitly "
        "complementary to the pocket geometry. Has a known beam-search "
        "'initialization failed' edge case — request more samples than needed."
    ),
    "targetdiff": (
        "TargetDiff E(3)-equivariant diffusion generator (ICLR 2023). Denoises "
        "from random noise across 1000 steps, conditioned on the pocket shape. "
        "Highest-fidelity 3D, slowest (~30-60 sec/molecule on GPU). Has a ~50% "
        "reconstruction failure rate; the survivors tend to score competitively "
        "but with inflated docking due to larger-molecule bias (ICLR 2025)."
    ),
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--backend", action="append", required=True,
                        help="Format: name:campaign_id:aizynth_json_path")
    parser.add_argument("--output_json", required=True, type=Path)
    parser.add_argument("--output_js", required=True, type=Path)
    parser.add_argument("--max_per_backend", type=int, default=30,
                        help="Cap molecules per backend (default 30)")
    args = parser.parse_args()

    backends_data: dict[str, dict] = {}
    insertion_order: list[str] = []
    for spec in args.backend:
        name, campaign_id, az_path = parse_backend_arg(spec)
        print(f"[Build] {name}: campaign={campaign_id}, aizynth={az_path}")
        backends_data[name] = build_backend_dataset(
            name, campaign_id, az_path, max_molecules=args.max_per_backend
        )
        insertion_order.append(name)

    # Preserve user-supplied order so the first --backend arg becomes the default.
    backend_names = insertion_order
    doc = {
        "target": {
            "pdb": "1M17",
            "name": "EGFR (Epidermal Growth Factor Receptor)",
            "disease": "Non-small cell lung cancer",
            "known_drug": "Erlotinib",
        },
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "pipeline": {
            "screening_backend": "MolScore + ADMET-AI (104 properties)",
            "docking_backend": "AutoDock Vina (production)",
            "synthesis_backend": "AiZynthFinder (USPTO + ZINC)",
        },
        "backends": backends_data,
        "backend_order": backend_names,
        "default_backend": backend_names[0] if backend_names else None,
        "generator_descriptions": {
            name: GENERATOR_DESCRIPTIONS.get(name, "") for name in backend_names
        },
        "admet_pass_rules": {
            k: {"kind": v[0], "threshold": v[1]} for k, v in ADMET_PASS_RULES.items()
        },
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_json, "w") as f:
        json.dump(doc, f, indent=2)

    args.output_js.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_js, "w") as f:
        f.write("window.PROFESSOR_DEMO_DATA = ")
        json.dump(doc, f)
        f.write(";\n")

    print(f"\n[Build] Wrote {args.output_json}")
    print(f"[Build] Wrote {args.output_js}")
    print(f"[Build] Backends: {', '.join(backend_names)}")
    for name, d in backends_data.items():
        s = d["summary"]
        print(f"  {name:12s} mols={s['total_generated']:3d} "
              f"screening_pass={s['passed_screening']:3d} "
              f"synth={s['synthesis_routes_found']}/{s['synthesis_evaluated']} "
              f"best_dock={s['best_docking_score']}")


if __name__ == "__main__":
    main()
