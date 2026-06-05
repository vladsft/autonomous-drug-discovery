#!/usr/bin/env python3
"""Pharmacophore bridge — harness TargetDiff as binding-mode hypotheses.

The cascade docks BOTH backends. TargetDiff's molecules are usually unmakeable
(~0% routable), so they never win the synthesis-aware ranking on their own — but
their docked poses encode WHERE and HOW the pocket likes to be engaged. This
bridge turns those poses into a target: it scores every (makeable) candidate by
how well it REPRODUCES TargetDiff's binding-mode pharmacophore, and folds that
into the final ranking. Net effect: makeable molecules that hit the same pocket
interactions TargetDiff discovered rise to the top — TargetDiff's GPU work shapes
the selection instead of being discarded.

Because both backends dock into the same receptor, their poses share one
coordinate frame, so feature overlap is a direct geometric comparison — no
alignment needed. When a docked pose isn't available for a molecule, the bridge
falls back to a 2D pharmacophore-fingerprint similarity (topology proxy).

Inputs:
  --ranked        ranked_candidates.json from the normal Stage-5 rank (BOTH
                  backends present, each with scores.{docking,admet,synthesis}).
  --docked-poses  docked_poses.sdf from Stage 4 (mols w/ molecule_id, backend,
                  3D coords). Optional — bridge falls back to 2D without it.
  --candidates-sdf  merged generated_molecules.sdf (molecule_id + backend) —
                  used to tag backend when poses are missing.
  --output        augmented ranked_candidates.json (pharmacophore_match added,
                  re-ranked).

Runs in `base` (RDKit). See plan.md "Pharmacophore-guided cascade".
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def _feature_factory():
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    fdef = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
    return ChemicalFeatures.BuildFeatureFactory(fdef)


# Pharmacophore families we compare on (the interaction-relevant ones).
_FAMILIES = {"Donor", "Acceptor", "Aromatic", "Hydrophobe",
             "LumpedHydrophobe", "PosIonizable", "NegIonizable"}


def _features_3d(mol, ff):
    """[(family, (x,y,z))] for a mol that has a 3D conformer."""
    if mol is None or mol.GetNumConformers() == 0:
        return []
    conf = mol.GetConformer()
    out = []
    for f in ff.GetFeaturesForMol(mol):
        fam = f.GetFamily()
        if fam not in _FAMILIES:
            continue
        # feature position is the centroid of its atoms
        p = f.GetPos(conf.GetId()) if hasattr(f, "GetPos") else conf.GetAtomPosition(f.GetAtomIds()[0])
        out.append((fam, (p.x, p.y, p.z)))
    return out


def _overlap_3d(cand_feats, seed_feats, tol=2.0):
    """Fraction of SEED features reproduced by a candidate feature of the same
    family within `tol` Å. In [0,1]; rewards reproducing the binding mode."""
    if not seed_feats:
        return 0.0
    matched = 0
    for sfam, (sx, sy, sz) in seed_feats:
        for cfam, (cx, cy, cz) in cand_feats:
            if cfam != sfam:
                continue
            if (sx - cx) ** 2 + (sy - cy) ** 2 + (sz - cz) ** 2 <= tol * tol:
                matched += 1
                break
    return matched / len(seed_feats)


def _sig_factory(ff):
    from rdkit.Chem.Pharm2D.SigFactory import SigFactory
    sf = SigFactory(ff, minPointCount=2, maxPointCount=3)
    sf.SetBins([(0, 2), (2, 5), (5, 8), (8, 100)])
    sf.Init()
    return sf


def _pharm2d_fp(smiles, ff, sig):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Pharm2D import Generate
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    AllChem.Compute2DCoords(m)
    return Generate.Gen2DFingerprint(m, sig)


def _tanimoto(fp1, fp2):
    from rdkit import DataStructs
    return float(DataStructs.TanimotoSimilarity(fp1, fp2))


def _backend_of(mol_id, sdf_backend, ranked_backend=None):
    """Best-effort backend tag: explicit prop, else id prefix heuristic."""
    b = sdf_backend or ranked_backend
    if b:
        return b.lower()
    mid = mol_id.lower()
    if mid.startswith("td") or "targetdiff" in mid:
        return "td"
    if mid.startswith("rdkit"):
        return "rdkit"
    return "unknown"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ranked", required=True, help="Stage-5 ranked_candidates.json")
    ap.add_argument("--output", required=True, help="augmented ranked json out")
    ap.add_argument("--docked-poses", default=None, help="docked_poses.sdf (Stage 4)")
    ap.add_argument("--candidates-sdf", default=None, help="merged generated SDF (backend tags)")
    ap.add_argument("--seeds-top-k", type=int, default=5,
                    help="how many top-docked TargetDiff poses define the binding-mode hypotheses")
    ap.add_argument("--feature-tol", type=float, default=2.0,
                    help="Å tolerance for 3D feature-overlap matching")
    ap.add_argument("--w-dock", type=float, default=0.35)
    ap.add_argument("--w-admet", type=float, default=0.20)
    ap.add_argument("--w-synth", type=float, default=0.15)
    ap.add_argument("--w-pharm", type=float, default=0.30)
    args = ap.parse_args()

    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    ranked = json.load(open(args.ranked))
    cands = ranked.get("candidates", [])
    if not cands:
        print("[bridge] no candidates in ranked file; nothing to do", file=sys.stderr)
        return 1

    # backend tags from the merged candidates SDF
    sdf_backend: dict[str, str] = {}
    if args.candidates_sdf and Path(args.candidates_sdf).exists():
        for m in Chem.SDMolSupplier(args.candidates_sdf, removeHs=False):
            if m is None or not m.HasProp("molecule_id"):
                continue
            if m.HasProp("backend"):
                sdf_backend[m.GetProp("molecule_id")] = m.GetProp("backend")

    # docked poses (3D) keyed by molecule_id
    poses: dict[str, object] = {}
    pose_backend: dict[str, str] = {}
    if args.docked_poses and Path(args.docked_poses).exists():
        for m in Chem.SDMolSupplier(args.docked_poses, removeHs=False):
            if m is None or not m.HasProp("molecule_id"):
                continue
            mid = m.GetProp("molecule_id")
            poses[mid] = m
            if m.HasProp("backend"):
                pose_backend[mid] = m.GetProp("backend")
    have_3d = len(poses) > 0
    print(f"[bridge] {len(cands)} ranked candidates; {len(poses)} docked poses "
          f"({'3D overlap' if have_3d else '2D fallback'} mode)")

    def backend(mid):
        return _backend_of(mid, sdf_backend.get(mid) or pose_backend.get(mid))

    # ── seeds: TargetDiff candidates by best docking affinity ──────────────
    td = [c for c in cands if backend(c["ligand_id"]) == "td"]
    td.sort(key=lambda c: c.get("docking_affinity", 0.0))  # most negative first
    seeds = td[: args.seeds_top_k]
    if not seeds:
        print("[bridge] WARNING: no TargetDiff seeds found — pharmacophore term "
              "will be 0 for all (check backend tags). Writing pass-through.")

    ff = _feature_factory()
    sig = _sig_factory(ff) if not have_3d or any(s["ligand_id"] not in poses for s in seeds) else None

    # Precompute seed representations
    seed_feats_3d = {s["ligand_id"]: _features_3d(poses.get(s["ligand_id"]), ff) for s in seeds}
    seed_fp_2d = {}
    if sig is not None:
        for s in seeds:
            seed_fp_2d[s["ligand_id"]] = _pharm2d_fp(s["smiles"], ff, sig)

    def pharm_match(cand):
        cid = cand["ligand_id"]
        best, best_seed = 0.0, None
        for s in seeds:
            sid = s["ligand_id"]
            score = 0.0
            if have_3d and cid in poses and seed_feats_3d.get(sid):
                score = _overlap_3d(_features_3d(poses[cid], ff), seed_feats_3d[sid],
                                    tol=args.feature_tol)
            elif sig is not None and seed_fp_2d.get(sid) is not None:
                fp = _pharm2d_fp(cand["smiles"], ff, sig)
                if fp is not None:
                    score = _tanimoto(fp, seed_fp_2d[sid])
            if score > best:
                best, best_seed = score, sid
        return best, best_seed

    # ── re-score ───────────────────────────────────────────────────────────
    wsum = args.w_dock + args.w_admet + args.w_synth + args.w_pharm
    for c in cands:
        sc = c.get("scores", {})
        pm, seed_id = pharm_match(c)
        c.setdefault("scores", {})["pharmacophore_match"] = round(pm, 4)
        c["pharmacophore_seed"] = seed_id
        c["backend"] = backend(c["ligand_id"])
        combined = (args.w_dock * sc.get("docking", 0.0)
                    + args.w_admet * sc.get("admet", 0.0)
                    + args.w_synth * sc.get("synthesis", 0.0)
                    + args.w_pharm * pm) / wsum
        c["combined_score"] = round(combined, 5)

    cands.sort(key=lambda c: c["combined_score"], reverse=True)
    for i, c in enumerate(cands, 1):
        c["final_rank"] = i

    ranked["weights"] = {"docking": args.w_dock, "admet": args.w_admet,
                         "synthesis": args.w_synth, "pharmacophore_match": args.w_pharm}
    ranked["pharmacophore_bridge"] = {
        "mode": "3d_overlap" if have_3d else "2d_fingerprint_fallback",
        "feature_tol_A": args.feature_tol,
        "seeds": [{"ligand_id": s["ligand_id"], "smiles": s["smiles"],
                   "docking_affinity": s.get("docking_affinity")} for s in seeds],
        "note": "TargetDiff top-docked poses used as binding-mode hypotheses; "
                "makeable candidates scored by pharmacophore fidelity to them.",
    }
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    json.dump(ranked, open(args.output, "w"), indent=2)
    print(f"[bridge] wrote pharmacophore-weighted ranking → {args.output} "
          f"({len(seeds)} seeds, top mol now {cands[0]['ligand_id']} "
          f"pharm={cands[0]['scores'].get('pharmacophore_match')})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
