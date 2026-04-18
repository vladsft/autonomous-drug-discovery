# Autonomous Drug Discovery Pipeline

## What this is
End-to-end computational drug discovery. Input: protein PDB file. Output: ranked drug candidates with ADMET profiles. Stages: P2Rank pocket detection → molecule generation (RDKit / TargetDiff / Pocket2Mol) → MolScore + ADMET-AI screening → AutoDock Vina docking. All logged to SQLite telemetry.

Target users: academic labs and small biotechs who can't afford Schrödinger. See `docs/north-star.md` for vision and strategy.

## Python environments
- **base** (`conda run -n base python`): RDKit, Vina, Meeko, gemmi, ADMET-AI, MolScore, Java 17, pyyaml. Default env for everything except deep-learning generation.
- **targetdiff_env** (`conda run -n targetdiff_env python`): PyTorch 1.13 + CUDA 11.7, PyG 2.5, Python 3.8. For `--mode targetdiff`.
- **pocket2mol_env** (`conda run -n pocket2mol_env python`): PyTorch 1.10.1 + CUDA 11.3, PyG, Python 3.8. For `--mode pocket2mol`.
- **Never use** bare `python3` — system Python lacks all chemistry deps.
- Miniconda installed at `~/miniconda3/`; conda binary at `~/miniconda3/bin/conda`. `CONDA_EXE` env var is respected by orchestrator.

## Running the pipeline
```bash
cd /home/vladsft/autonomous-drug-discovery/autonomous_drug_discovery
conda run -n base python orchestrator.py run data/processed/1M17.pdb --mode production
conda run -n base python orchestrator.py run data/processed/6P3D.pdb --mode pocket2mol
conda run -n base python benchmark.py   # cross-target comparison report
```

Modes: `simulation` (stub), `production` (RDKit + Vina), `targetdiff` (diffusion), `pocket2mol` (autoregressive).

## Key paths
```
autonomous_drug_discovery/
  orchestrator.py                    # CLI entrypoint
  telemetry.py                       # SQLite: tables runs, molecule_scores
  benchmark.py                       # M2 validation report
  modules/
    01_ingestion/run_pocket.py       # P2Rank (default) or fpocket via --backend
    02_generation/run_generation.py  # Modes: simulation, rdkit, targetdiff, pocket2mol
    02_generation/targetdiff/        # TargetDiff repo + checkpoint (pretrained_diffusion.pt)
    02_generation/pocket2mol/        # Pocket2Mol repo + checkpoint (pretrained_Pocket2Mol.pt)
    03_screening/run_screening.py    # MolScore (primary) or RDKit fallback + ADMET-AI (104 properties)
    04_docking/run_docking.py        # AutoDock Vina + Meeko; modes: simulation, triage, production
  envs/
    env_orchestrator.yml             # base env spec
    env_targetdiff.yml               # targetdiff_env spec
    env_pocket2mol.yml               # pocket2mol_env spec
  data/
    processed/*.pdb                  # Input protein structures
    telemetry.db                     # Source of truth for all results
docs/                                # north-star, pipeline-guide, installation, targetdiff-setup, etc.
```

## External binaries
- P2Rank: `/home/vladsft/p2rank_2.5.1/prank`
- fpocket: `/home/vladsft/fpocket/bin/fpocket`

## Project status (as of 2026-04-18)
- **M1 (Working Pipeline):** DONE
- **M2 (Validation):** DONE — 3 cancer targets validated against X-ray crystallography (fpocket and P2Rank)
- **M2.5 (TargetDiff POC):** DONE — standalone diffusion generation on BRAF V600E
- **M2.6 (Pocket2Mol integration):** DONE — wired into orchestrator as `--mode pocket2mol`, ~11× faster than TargetDiff
- **M3 (Expert Review):** NOT STARTED — need medicinal chemistry advisor ("non-negotiable" per strategy doc)
- **M4-M5:** NOT STARTED

## Validation results (P2Rank + ADMET-AI, RDKit generation)
| Target | PDB | Best Dock | Avg Dock | Pocket Distance | Residue Overlap |
|--------|-----|-----------|----------|----------------|-----------------|
| EGFR | 1M17 | -9.32 | -6.58 | 2.7 A | 82% |
| BCR-ABL | 2HYY | -12.59 | -9.25 | 2.7 A | 92% |
| BRAF V600E | 6P3D | -11.20 | -8.39 | 3.1 A | 89% |

## TargetDiff diffusion results (standalone POC, BRAF V600E pocket)
| Molecule | SMILES | MW | QED | Dock | Ligand Eff. |
|----------|--------|----|-----|------|-------------|
| Mol 1 | `C1=CN=CC=C(c2cccc(Nc3ccncc3)c2)C1` | 261 | 0.899 | -7.59 | 0.38 |
| Mol 2 | `COc1cnc(C(=O)NCc2cccc(C)n2)cn1` | 258 | 0.889 | -7.38 | 0.39 |

~12 min/molecule on CPU. 50% reconstruction failure rate (normal for diffusion). Now runnable through orchestrator via `--mode targetdiff`.

## Caveats you must understand
- **Vina docking scores correlate weakly with real binding** (r=0.4-0.6). Ranking signal, not measurement.
- **SBDD models inflate Vina scores** by generating larger molecules. Always check MW; use ligand efficiency for fair comparison.
- **Screening survival rates are high** (73-98%). Should be 40-60%. Thresholds in `default_scoring_config.json` need tightening with expert input.
- **Diffusion models (TargetDiff) have ~50% reconstruction failure** — generate more than needed.
- **Pocket2Mol is autoregressive, not diffusion** — uses pocket bounding box (default 23 Å), known "initialization failed" edge cases on some targets.
- **pytdc is not installed** on this system (incompatible with Python 3.13). Triage docking mode unavailable; production Vina is unaffected.

## What to work on next
1. Empirical validation: run Pocket2Mol on EGFR / BCR-ABL / BRAF and compare to RDKit+TargetDiff results
2. Tighten screening thresholds with domain expert input
3. Add GNINA CNN re-scoring (needs GPU binary from github.com/gnina/gnina/releases)
4. Add AiZynthFinder retrosynthetic feasibility (`pip install aizynthfinder`, needs policy models)
5. Docker containerization for reproducible deployment
6. Engage medicinal chemistry advisor for M3
