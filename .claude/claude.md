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
    02_generation/targetdiff/        # NOT checked in — user clones per docs/targetdiff-setup.md (needs pretrained_diffusion.pt)
    02_generation/pocket2mol/        # Pocket2Mol repo + checkpoint (pretrained_Pocket2Mol.pt) — present
    03_screening/run_screening.py    # MolScore (primary) or RDKit fallback + ADMET-AI (104 properties)
    04_docking/run_docking.py        # AutoDock Vina + Meeko; modes: simulation, triage, production
    05_ranking/run_ranking.py        # Multi-criteria final ranker — wired into orchestrator `run` + `rank` subcommand. Composite = 0.5·docking + 0.3·ADMET + 0.2·synthesis (synth = 0.5 placeholder until AiZynth slot is wired)
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

## Project status (as of 2026-05-10)
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
- **Pocket2Mol GPU mode is broken on this machine because of a PyTorch ↔ GPU architecture mismatch.** Root cause confirmed 2026-05-10:
  - RTX 5060 is **Blackwell, sm_120**, released late 2024.
  - `pocket2mol_env` pins `pytorch=1.10.1` + `cudatoolkit=11.3` (mirrors upstream Pocket2Mol's 2022 spec). CUDA 11.3 only knows architectures up to **sm_86 (Ampere)**.
  - The driver's PTX-JIT forward-compatibility path *technically* runs but is pathologically slow for sm_120 from sm_86-vintage PTX — minutes per kernel call. Pocket2Mol's autoregressive loop fires thousands of small kernel calls, each one re-triggering JIT.
  - Diagnostic: a trivial `torch.randn(2000,2000).cuda(); for _ in range(10): x @ x.T` test ran for 2:26 with no output before being killed. Same script in any modern env runs in <1 sec.
  - **CPU mode (`--device cpu`) works**: init completed in 4 sec on the official `4yhj.pdb` example (vs >10 min stalled on broken GPU), then proceeded to actual atom-by-atom sampling. Slower than a properly-configured GPU but functional.
  - **Same root cause likely affects the planned `targetdiff_env`** — its spec also pins pre-Blackwell PyTorch (1.13 + cu117). Once we recover the TargetDiff checkpoint, expect the same wall.
  - **Fix planned (see plan.md Step 9)**: rebuild both envs on PyTorch 2.4+ / CUDA 12.x.
- **TargetDiff is not reproducible on this machine.** The standalone POC referenced in `validation_results` (BRAF, 2 molecules, ~7 kcal/mol) was a one-off run on 2026-04-18; the `targetdiff_env`, the `targetdiff/` repo subdirectory, and `pretrained_diffusion.pt` are all gone, and the original Google Drive folder for the checkpoint now returns 404. The only surviving artefacts are the PNG visualisations in `reports/`. To re-run, you need a backup checkpoint or to contact the paper authors.
- **pytdc is not installed** on this system (incompatible with Python 3.13). Triage docking mode unavailable; production Vina is unaffected.
- **GPU available**: NVIDIA RTX 5060 (8 GB VRAM, CUDA 13.2) — TargetDiff/Pocket2Mol now feasible at full speed.

## What to work on next (4-phase roadmap — see `autonomous_drug_discovery/plan.md`)
1. **Phase 1 — Cloud GPU + more targets.** Recover TargetDiff checkpoint; rebuild `pocket2mol_env_v2` / `targetdiff_env_v2` on PyTorch 2.4 + cu121 (current envs don't run on RTX 5060); rent ~$5-10 of cloud GPU time; run cascade on 5-10 oncology kinase targets (KRAS G12C, JAK2, CDK4/6, ER-α, AR + the existing three).
2. **Phase 2 — Professor review (M3).** Walk advisor through dashboard, capture per-candidate annotations, calibrate filters, get input on the adaptive layer design.
3. **Phase 3 — Adaptive layer.** Ship in order: (a) Obsidian campaign emitter, (b) Bayesian (Thompson-sampling) strategy recommender, (c) Sonnet-in-the-loop agent for retrieval / mid-campaign sanity checks / report writing.
4. **Phase 4 — Bayesian evaluation.** Posteriors with credible intervals on lift vs deterministic baseline. Honest, even if null.

Already done: AiZynthFinder integration, Stage 5 ranker wiring, multi-backend dashboard, Pocket2Mol CPU-mode patch.
