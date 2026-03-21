# Autonomous Drug Discovery Pipeline

## What this is
End-to-end computational drug discovery. Input: protein PDB file. Output: ranked drug candidates with ADMET profiles. Stages: P2Rank pocket detection → RDKit molecule generation → RDKit+ADMET-AI screening → AutoDock Vina docking. All logged to SQLite telemetry.

Target users: academic labs and small biotechs who can't afford Schrodinger. See `reports/north_star.md` for vision, `compass_artifact_*.md` at root for strategic/market analysis.

## Python environments
- **base env** (`conda run -n base python`): RDKit, Vina, Meeko, gemmi, ADMET-AI, Java 17. Use this for everything except TargetDiff.
- **targetdiff_env** (`conda run -n targetdiff_env python`): PyTorch 1.13.1+cu117, PyG 2.5.2, Python 3.8. Working. Checkpoint downloaded. No GPU on this machine — CPU inference takes hours.
- **Never use** bare `python3` — system Python lacks all chemistry deps.

## Running the pipeline
```bash
cd /home/vladsft/agent-harness/autonomous_drug_discovery
conda run -n base python orchestrator.py run data/processed/1M17.pdb --mode production
conda run -n base python benchmark.py   # cross-target comparison report
```

## Key paths
```
autonomous_drug_discovery/
  orchestrator.py                    # CLI entrypoint
  telemetry.py                       # SQLite: tables runs, molecule_scores
  benchmark.py                       # M2 validation report
  modules/
    01_ingestion/run_pocket.py       # P2Rank (default) or fpocket via --backend
    02_generation/run_generation.py  # RDKit fragment-based or TargetDiff diffusion
    02_generation/targetdiff/        # TargetDiff repo + pretrained checkpoint
    03_screening/run_screening.py    # RDKit filters + ADMET-AI (104 properties)
    04_docking/run_docking.py        # AutoDock Vina + Meeko
  data/
    processed/*.pdb                  # Input protein structures
    telemetry.db                     # Source of truth for all results
reports/
  north_star.md                      # Vision doc, milestones M1-M5
  testing_guide.md                   # Plain-language guide for non-experts
```

## External binaries
- P2Rank: `/home/vladsft/p2rank_2.5.1/prank`
- fpocket: `/home/vladsft/fpocket/bin/fpocket`

## Project status (as of 2026-03-21)
- **M1 (Working Pipeline):** DONE
- **M2 (Validation):** DONE — 3 cancer targets validated against X-ray crystallography
- **M3 (Expert Review):** NOT STARTED — need medicinal chemistry advisor ("non-negotiable" per strategy doc)
- **M4-M5:** NOT STARTED

## Validation results (P2Rank + ADMET-AI)
| Target | PDB | Best Dock | Pocket Distance | Residue Overlap |
|--------|-----|-----------|----------------|-----------------|
| EGFR | 1M17 | -9.32 | 2.7 A | 82% |
| BCR-ABL | 2HYY | -12.56 | 2.7 A | 92% |
| BRAF V600E | 6P3D | -10.40 | 3.1 A | 89% |

## Caveats you must understand
- **Vina docking scores correlate weakly with real binding** (r=0.4-0.6). They're a ranking signal, not a measurement. Don't present them as proof a molecule works.
- **SBDD models can inflate Vina scores** by generating larger molecules. Always check molecular weight alongside docking score.
- **Concurrent pipeline runs collide** — shared output dirs with hardcoded filenames. Run targets sequentially.
- **Screening survival rates are high** (73-98%). Should be 40-60%. Thresholds in `default_scoring_config.json` need tightening with domain expert input.

## What to work on next
1. Run TargetDiff diffusion generation (env ready, slow on CPU)
2. Add GNINA re-scoring (needs GPU binary from github.com/gnina/gnina/releases)
3. Add AiZynthFinder retrosynthetic feasibility (`pip install aizynthfinder`, needs policy models)
4. Tighten screening thresholds
5. Per-campaign output directories in orchestrator (prevents file collisions)
6. Engage a medicinal chemistry advisor for M3
