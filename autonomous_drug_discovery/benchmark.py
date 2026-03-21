"""
Benchmark Comparison Script — M2 Validation.

Queries the telemetry database for all completed production campaigns,
compares results across validation targets (EGFR, BCR-ABL, BRAF V600E),
and prints a human-readable report.

Usage:
    conda run -n base python benchmark.py
    conda run -n base python benchmark.py --db_path data/telemetry.db
"""

import sys
import json
import sqlite3
import argparse
from pathlib import Path
from datetime import datetime


DEFAULT_DB = Path(__file__).parent / "data" / "telemetry.db"

# Known validation targets and their expected properties
KNOWN_TARGETS = {
    "1M17": {
        "name": "EGFR kinase",
        "disease": "Non-small cell lung cancer",
        "known_drug": "Erlotinib",
        "known_drug_affinity_range": (-9.5, -7.0),  # literature Kd/Ki estimates
    },
    "2HYY": {
        "name": "ABL1 kinase (BCR-ABL)",
        "disease": "Chronic myeloid leukemia",
        "known_drug": "Imatinib",
        "known_drug_affinity_range": (-10.0, -7.5),
    },
    "6P3D": {
        "name": "BRAF V600E kinase",
        "disease": "Melanoma",
        "known_drug": "Ponatinib",
        "known_drug_affinity_range": (-10.5, -8.0),
    },
}


def get_campaigns(db_path):
    """Get all campaigns with their targets and stages."""
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    campaigns = {}
    rows = conn.execute("""
        SELECT campaign_id, module_name, status, input_path, started_at, completed_at
        FROM runs ORDER BY started_at
    """).fetchall()

    for row in rows:
        cid = row["campaign_id"]
        if cid not in campaigns:
            campaigns[cid] = {"stages": {}, "input_path": None, "target": None, "started": row["started_at"]}
        campaigns[cid]["stages"][row["module_name"]] = row["status"]
        if row["input_path"] and row["module_name"] == "01_ingestion":
            campaigns[cid]["input_path"] = row["input_path"]
            # Extract PDB ID from path
            pdb_name = Path(row["input_path"]).stem
            campaigns[cid]["target"] = pdb_name

    conn.close()
    return campaigns


def get_campaign_molecules(db_path, campaign_id):
    """Get molecule scores for a campaign."""
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row

    rows = conn.execute("""
        SELECT ms.smiles, ms.qed, ms.sa_score, ms.logp, ms.mol_weight,
               ms.docking_score, ms.passed_triage, ms.stage_eliminated,
               r.module_name
        FROM molecule_scores ms
        JOIN runs r ON ms.run_id = r.run_id
        WHERE r.campaign_id = ?
        ORDER BY ms.docking_score
    """, (campaign_id,)).fetchall()

    conn.close()
    return [dict(r) for r in rows]


def analyze_campaign(db_path, campaign_id, target_pdb):
    """Analyze a single campaign's results."""
    molecules = get_campaign_molecules(db_path, campaign_id)

    # Separate screening vs docking data
    screened = [m for m in molecules if m["module_name"] == "03_screening"]
    docked = [m for m in molecules if m["module_name"] == "04_docking" and m["docking_score"] is not None]

    # Screening stats
    total_screened = len(screened)
    passed_screening = sum(1 for m in screened if m["passed_triage"] == 1)

    # Docking stats
    docking_scores = sorted([m["docking_score"] for m in docked])
    best_score = docking_scores[0] if docking_scores else None
    worst_score = docking_scores[-1] if docking_scores else None
    median_score = docking_scores[len(docking_scores) // 2] if docking_scores else None

    # Drug-likeness of docked molecules (from screening stage)
    qed_values = [m["qed"] for m in screened if m["qed"] is not None and m["passed_triage"] == 1]
    sa_values = [m["sa_score"] for m in screened if m["sa_score"] is not None and m["passed_triage"] == 1]

    # Count strong binders (better than -7 kcal/mol is generally considered decent)
    strong_binders = sum(1 for s in docking_scores if s < -7.0)
    very_strong = sum(1 for s in docking_scores if s < -9.0)

    return {
        "campaign_id": campaign_id,
        "target": target_pdb,
        "total_screened": total_screened,
        "passed_screening": passed_screening,
        "survival_rate": round(passed_screening / max(total_screened, 1) * 100, 1),
        "total_docked": len(docked),
        "best_docking_score": best_score,
        "worst_docking_score": worst_score,
        "median_docking_score": median_score,
        "strong_binders_lt_neg7": strong_binders,
        "very_strong_lt_neg9": very_strong,
        "mean_qed": round(sum(qed_values) / max(len(qed_values), 1), 3) if qed_values else None,
        "mean_sa_score": round(sum(sa_values) / max(len(sa_values), 1), 2) if sa_values else None,
        "docking_scores": docking_scores,
    }


def print_report(analyses):
    """Print a formatted benchmark comparison report."""
    print("\n" + "=" * 80)
    print("  BENCHMARK REPORT — M2 Validation Against Known Targets")
    print("  Generated:", datetime.now().strftime("%Y-%m-%d %H:%M"))
    print("=" * 80)

    for a in analyses:
        target_info = KNOWN_TARGETS.get(a["target"], {})
        print(f"\n{'─' * 80}")
        print(f"  TARGET: {a['target']} — {target_info.get('name', 'Unknown')}")
        print(f"  Disease: {target_info.get('disease', 'Unknown')}")
        print(f"  Known drug: {target_info.get('known_drug', 'Unknown')}")
        print(f"  Campaign: {a['campaign_id']}")
        print(f"{'─' * 80}")

        print(f"\n  SCREENING")
        print(f"    Molecules screened:    {a['total_screened']}")
        print(f"    Passed screening:      {a['passed_screening']}")
        print(f"    Survival rate:         {a['survival_rate']}%")

        print(f"\n  DOCKING")
        print(f"    Molecules docked:      {a['total_docked']}")
        if a["best_docking_score"] is not None:
            print(f"    Best score:            {a['best_docking_score']:.2f} kcal/mol")
            print(f"    Median score:          {a['median_docking_score']:.2f} kcal/mol")
            print(f"    Worst score:           {a['worst_docking_score']:.2f} kcal/mol")
            print(f"    Strong binders (<-7):  {a['strong_binders_lt_neg7']}")
            print(f"    Very strong (<-9):     {a['very_strong_lt_neg9']}")

            known_range = target_info.get("known_drug_affinity_range")
            if known_range:
                in_range = sum(1 for s in a["docking_scores"] if known_range[0] <= s <= known_range[1])
                print(f"    In known drug range:   {in_range} molecules score between {known_range[0]} and {known_range[1]} kcal/mol")
        else:
            print(f"    (no docking results)")

        print(f"\n  DRUG-LIKENESS (screened survivors)")
        if a["mean_qed"] is not None:
            print(f"    Mean QED:              {a['mean_qed']}  (0-1, higher=more drug-like)")
            print(f"    Mean SA Score:         {a['mean_sa_score']}  (1-10, lower=easier to make)")
        else:
            print(f"    (no screening data)")

    # Cross-target comparison
    if len(analyses) > 1:
        print(f"\n{'=' * 80}")
        print(f"  CROSS-TARGET COMPARISON")
        print(f"{'=' * 80}")
        print(f"\n  {'Target':<10} {'Screened':>10} {'Survived':>10} {'Rate':>8} {'Docked':>8} {'Best':>10} {'Strong':>8}")
        print(f"  {'─' * 10} {'─' * 10} {'─' * 10} {'─' * 8} {'─' * 8} {'─' * 10} {'─' * 8}")
        for a in analyses:
            best = f"{a['best_docking_score']:.2f}" if a['best_docking_score'] else "N/A"
            print(f"  {a['target']:<10} {a['total_screened']:>10} {a['passed_screening']:>10} {a['survival_rate']:>7}% {a['total_docked']:>8} {best:>10} {a['strong_binders_lt_neg7']:>8}")

    print(f"\n{'=' * 80}")
    print(f"  INTERPRETATION GUIDE")
    print(f"{'=' * 80}")
    print("""
  Docking scores (kcal/mol) — more negative = stronger predicted binding:
    > -5:     Weak or no binding (probably not useful)
    -5 to -7: Moderate binding (worth investigating)
    -7 to -9: Strong binding (promising candidate)
    < -9:     Very strong binding (excellent candidate)

  QED (Quantitative Estimate of Drug-likeness):
    0.0-0.3:  Poor drug-likeness
    0.3-0.6:  Moderate
    0.6-1.0:  Good drug-likeness

  SA Score (Synthetic Accessibility):
    1-3:  Easy to synthesize
    3-6:  Moderate difficulty
    6-10: Very hard to synthesize

  Survival rate: What fraction of generated molecules pass basic drug safety
  filters. 40-60% is typical for fragment-based generation. Too high (>80%)
  may mean filters are too loose; too low (<20%) means generation is producing
  too much junk.

  A successful validation means: the pipeline finds molecules that score in
  the same ballpark as known drugs for that target, even though it has no
  knowledge of those drugs. This shows the pipeline can independently identify
  promising chemistry.
""")


def main():
    parser = argparse.ArgumentParser(description="M2 Benchmark Comparison")
    parser.add_argument("--db_path", default=str(DEFAULT_DB), help="Telemetry database path")
    args = parser.parse_args()

    db_path = Path(args.db_path)
    if not db_path.exists():
        print(f"Error: Database not found at {db_path}")
        sys.exit(1)

    campaigns = get_campaigns(db_path)

    # Find the most recent successful production campaign for each target
    target_campaigns = {}
    for cid, info in campaigns.items():
        target = info.get("target")
        if target is None:
            continue
        # Must have completed all 4 stages successfully
        stages = info["stages"]
        if stages.get("04_docking") == "success":
            # Keep the most recent one per target
            if target not in target_campaigns or info["started"] > target_campaigns[target]["started"]:
                target_campaigns[target] = {"campaign_id": cid, **info}

    if not target_campaigns:
        print("No completed production campaigns found in the database.")
        print("\nAvailable campaigns:")
        for cid, info in campaigns.items():
            print(f"  {cid}: target={info.get('target', '?')}, stages={info['stages']}")
        sys.exit(1)

    print(f"Found {len(target_campaigns)} completed target(s): {', '.join(target_campaigns.keys())}")

    analyses = []
    for target, info in sorted(target_campaigns.items()):
        print(f"Analyzing {target} (campaign {info['campaign_id']})...")
        analysis = analyze_campaign(db_path, info["campaign_id"], target)
        analyses.append(analysis)

    print_report(analyses)


if __name__ == "__main__":
    main()
