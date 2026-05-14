# Autonomous Drug Discovery Pipeline

An end-to-end computational drug discovery system that takes a protein structure (PDB file) and produces a ranked list of candidate drug molecules. The pipeline automates pocket detection, molecule generation, drug-likeness screening, molecular docking, and multi-criteria ranking — with every step logged to a telemetry database for reproducibility and analysis.

## What it does

Given a protein target, the pipeline:

1. **Detects binding pockets** with P2Rank (ML-based, default) or fpocket (geometry fallback).
2. **Generates candidate molecules** with one of:
   - **RDKit** fragment-based combinatorial generation (default, fast, CPU)
   - **TargetDiff** E(3)-equivariant diffusion (pocket-conditioned, highest-fidelity 3D; needs GPU in practice)
   - **Pocket2Mol** autoregressive generator (pocket-conditioned, faster than TargetDiff; *currently deferred — pretrained checkpoint not recovered*)
3. **Screens candidates** against drug-likeness filters (Lipinski, QED, SA, PAINS) and ADMET-AI (104 properties).
4. **Docks survivors** into the pocket with AutoDock Vina.
5. **Ranks** via a composite score (docking + ADMET + optional AiZynthFinder synthesis feasibility) and writes a final scorecard.

Validated against three cancer targets with crystallographic ground truth:

| Target | Disease | Best Dock Score | Pocket Accuracy |
|---|---|---|---|
| EGFR (1M17) | Lung cancer | -9.32 kcal/mol | 2.7 Å from erlotinib, 82% residue overlap |
| BCR-ABL (2HYY) | Leukemia | -12.59 kcal/mol | 2.7 Å from imatinib, 92% residue overlap |
| BRAF V600E (6P3D) | Melanoma | -11.20 kcal/mol | 3.1 Å from ponatinib, 89% residue overlap |

## Quick start

The pipeline ships as a single Docker image with all conda environments, P2Rank, and (mirrored) model checkpoints baked in. Both contributors run the same image; the cloud GPU pod runs the same image.

> **Note on the Docker image:** The Dockerfile and CI workflow are part of the Phase 1 work currently in flight (see [`autonomous_drug_discovery/plan.md`](autonomous_drug_discovery/plan.md)). Until that lands, the image at `ghcr.io/<you>/agent-harness` does not yet exist; the legacy manual conda setup at the bottom of this README is the fallback. The commands below describe the target workflow.

### Prerequisites

- Docker 24+ (with the `compose` plugin)
- For local GPU runs: an NVIDIA GPU + the NVIDIA Container Toolkit
- For cloud GPU runs: a [RunPod](https://www.runpod.io/) account
- For state syncing: a [Cloudflare R2](https://www.cloudflare.com/products/r2/) bucket

### One-time setup

```bash
git clone <this-repo> && cd agent-harness
make bootstrap   # pulls Docker image, configures rclone for R2
```

### Run the full pipeline

```bash
# Local CPU run (RDKit, fast)
make run TARGET=1M17 MODE=production

# Local CPU sanity check
make run TARGET=1M17 MODE=simulation

# Cloud GPU run (TargetDiff, ~$0.30, ~30 min on an RTX 3090)
make cloud-run TARGET=2HYY MODE=targetdiff NUM=30
```

Output lands in `data/campaign_<id>/`. Sync to/from the shared R2 bucket with `make push` / `make pull`.

### View the dashboard

The static dashboard auto-deploys to GitHub Pages on every `main` commit. URL: `https://<you>.github.io/agent-harness/`. Locally:

```bash
make dashboard    # regenerate dashboard JSON from latest telemetry
open dashboard/index.html
```

## Documentation

| Doc | Question it answers |
|---|---|
| [`autonomous_drug_discovery/plan.md`](autonomous_drug_discovery/plan.md) | **Canonical architecture + roadmap.** Container architecture, operating principles, 4-phase roadmap, 7-day immediate action plan. Start here. |
| [`docs/north-star.md`](docs/north-star.md) | Vision and market positioning. Why this project exists, who it serves, what differentiates it. |
| [`docs/pipeline-guide.md`](docs/pipeline-guide.md) | Operational usage — `make` commands, per-stage commands, parameters, output structure. |
| [`docs/testing-guide.md`](docs/testing-guide.md) | Plain-language explanation of the science, validation experiments run to date, glossary. |
| [`docs/telemetry-guide.md`](docs/telemetry-guide.md) | SQLite schema, common queries, Python API. |

## Repository structure

```
.
├── autonomous_drug_discovery/       # Main application
│   ├── orchestrator.py              # CLI entrypoint — runs stages individually or as full pipeline
│   ├── agent_planner.py             # LLM-driven adaptive loop (Phase 3 scaffold)
│   ├── telemetry.py                 # SQLite telemetry database
│   ├── plan.md                      # CANONICAL architecture + roadmap
│   │
│   ├── modules/
│   │   ├── 01_ingestion/            # P2Rank / fpocket pocket detection
│   │   ├── 02_generation/           # RDKit / TargetDiff / Pocket2Mol molecule generators
│   │   ├── 03_screening/            # Lipinski + QED + SA + PAINS + ADMET-AI
│   │   ├── 04_docking/              # AutoDock Vina
│   │   └── 05_ranking/              # Multi-criteria composite scoring
│   │
│   ├── data/processed/              # Input PDBs and ingestion manifests
│   ├── envs/                        # Conda env specs (baked into Docker image at build time)
│   └── tests/
│
├── dashboard/                       # Static HTML + JSON dashboard
├── docs/                            # See table above
└── scripts/                         # bootstrap.sh, cloud_run.sh, regenerate_dashboard.py
```

## Architecture summary

The full architecture lives in [`plan.md`](autonomous_drug_discovery/plan.md). One-paragraph version:

A single Docker image (built by GitHub Actions, pushed to GHCR, weights mirrored to Hugging Face) runs identically on either contributor's laptop and on rented RunPod GPU pods. Campaign outputs and the telemetry SQLite DB live in a Cloudflare R2 bucket as the cross-machine source of truth. The chemist dashboard is a static site deployed to GitHub Pages and regenerated on every commit. A top-level `Makefile` is the user interface. The architecture is deliberately *not* Kubernetes, *not* a VM, *not* Modal — those are evaluated as post-demo work in `plan.md`.

## Status

- M1 (working pipeline) — done.
- M2 (validation against 3 crystal-structure targets) — done.
- **Phase 1 (containerise + cloud GPU + 20 quality results for the professor) — in flight; see `plan.md` Immediate Action Plan.**
- M3 (medicinal chemistry expert review) — Phase 2, post-demo.
- Adaptive Layer (Sonnet + Bayesian + Obsidian + Cascade) — Phase 3.
- Bayesian evaluation — Phase 4.

## Legacy / fallback: manual conda setup

While the Docker image is in flight, the manual conda workflow still works on machines that have the local checkpoint and environments already configured. It is **not** recommended for fresh setups — the Pocket2Mol Google Drive folder is permanently dead and the TargetDiff one is also gone (the dev box has a pre-takedown copy; see `plan.md`). For full historical install steps see the git history of `docs/installation.md` (deleted 2026-05-11) or the per-stage commands in [`docs/pipeline-guide.md`](docs/pipeline-guide.md).

The minimal legacy invocation looks like:

```bash
conda env create -f autonomous_drug_discovery/envs/env_orchestrator.yml
conda run -n base python autonomous_drug_discovery/orchestrator.py run \
    autonomous_drug_discovery/data/processed/1M17.pdb --mode simulation
```

Beyond simulation mode you also need P2Rank installed (`~/p2rank_2.5.1/prank`), and for TargetDiff/Pocket2Mol the relevant conda env *and* the pretrained checkpoint at the expected path.

## License

See `LICENSE`.
