# Autonomous Drug Discovery Pipeline

An end-to-end computational drug discovery system that takes a protein structure (PDB file) and produces a ranked list of candidate drug molecules. The pipeline automates pocket detection, molecule generation, drug-likeness screening, molecular docking, and multi-criteria ranking — with every step logged to a telemetry database for reproducibility and analysis.

## What it does

Given a protein target, the pipeline:

1. **Detects binding pockets** with P2Rank (ML-based, default) or fpocket (geometry fallback).
2. **Generates candidate molecules**. The recommended mode is **`cascade`** — RDKit + TargetDiff merged:
   - **RDKit** fragment-based combinatorial generation (CPU, fast; the *makeable* workhorse — ~30% of its hits have a real synthetic route)
   - **TargetDiff** E(3)-equivariant diffusion (pocket-conditioned, highest-fidelity 3D; needs a GPU). Reframed as a *binding-mode proposer* — its raw output is potent but rarely synthesizable on its own.
   - **Pocket2Mol** autoregressive generator — *dropped*: can't target NVIDIA Blackwell (CPU-only there) and adds little makeable matter the other two don't.
3. **Screens candidates** against drug-likeness filters (Lipinski, QED, SA, PAINS) and ADMET-AI (104 properties).
4. **Synthesizability gate (Stage 2.5)** — runs AiZynthFinder retrosynthesis (with an optional fast RAScore pre-filter). In single-backend modes it *filters* unmakeable molecules before docking; in `cascade` it runs as an *annotation* (so TargetDiff's poses are retained for step 6). Added after measuring that the 3D generators produce mostly unmakeable molecules (0/28 top kinase candidates had a route) — synthesizability is now a first-class signal, not a footnote.
5. **Docks survivors** (both backends, in `cascade`) into the pocket with AutoDock Vina; emits docked 3D poses.
6. **Ranks** via a composite score, and in `cascade` runs a **pharmacophore bridge (Stage 5.5)**: TargetDiff's top-docked poses become binding-mode hypotheses, and each makeable candidate is scored by how well its pose *reproduces* that pharmacophore — so makeable molecules that recreate TargetDiff's binding mode rise to the top. This is how TargetDiff's output is *harnessed* rather than discarded. Writes a final scorecard.

Validated against three cancer targets with crystallographic ground truth:

| Target | Disease | Best Dock Score | Pocket Accuracy |
|---|---|---|---|
| EGFR (1M17) | Lung cancer | -9.32 kcal/mol | 2.7 Å from erlotinib, 82% residue overlap |
| BCR-ABL (2HYY) | Leukemia | -12.59 kcal/mol | 2.7 Å from imatinib, 92% residue overlap |
| BRAF V600E (6P3D) | Melanoma | -11.20 kcal/mol | 3.1 Å from ponatinib, 89% residue overlap |

## Quick start

The pipeline ships as a single Docker image with all conda environments, P2Rank, and (mirrored) model checkpoints baked in. Both contributors run the same image; the cloud GPU pod runs the same image. CI builds the image on every `main` push and publishes it to `ghcr.io/vladsft/autonomous-drug-discovery:latest`.

> **Note on the Docker image:** The image is published by the `build` workflow on every `main` push (`:latest` and `:<sha>`). If `docker pull` fails immediately after a fresh clone, the most likely reason is that no commit has hit `main` yet on this fork — fall back to the [no-Docker path](#running-locally-without-docker) below until CI has run once.

### Prerequisites

- Docker 24+ (with the `compose` plugin)
- For local GPU runs: an NVIDIA GPU + the NVIDIA Container Toolkit
- For cloud GPU runs: a [RunPod](https://www.runpod.io/) account
- For state syncing: a [Cloudflare R2](https://www.cloudflare.com/products/r2/) bucket

### One-time setup

```bash
git clone https://github.com/vladsft/autonomous-drug-discovery.git
cd autonomous-drug-discovery
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

For 30-50 targets fire-and-forget: open the repo's **Actions** tab → **batch_cloud_run** → **Run workflow** → paste a list of PDB codes → submit. The action provisions concurrent RunPod pods, waits for all of them, regenerates the multi-target dashboard, and redeploys Pages — no laptop involvement after kickoff. See [`docs/pipeline-guide.md`](docs/pipeline-guide.md#batch-runs-cloud-hands-free--phase-15).

Output lands in `data/campaign_<id>/`. Sync to/from the shared R2 bucket with `make push` / `make pull`.

To run TargetDiff diffusion on a **local** GPU — the only option for NVIDIA Blackwell (RTX 50-series), which the cloud image's CUDA 11.7 stack cannot target — see [Running locally without Docker](#running-locally-without-docker).

### View the dashboard

The static dashboard auto-deploys to GitHub Pages on every `main` commit. URL: `https://vladsft.github.io/autonomous-drug-discovery/`. Locally:

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
└── scripts/                         # cloud_run.sh, pod_campaign.sh, telemetry + build helpers
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

## Running locally without Docker

Docker is the *portable* path, not the only one. The orchestrator runs each stage in a conda environment directly, so a machine with the environments installed can run the full pipeline — including TargetDiff diffusion on a local GPU — without building or pulling an image.

This is also the **only** way to use a GPU the cloud image's CUDA 11.7 stack cannot target — notably NVIDIA Blackwell (RTX 50-series, `sm_120`), whose architecture post-dates that image's PyTorch 1.13 / CUDA 11.7 environment.

### One-time setup

```bash
# Orchestrator + CPU stages (ingestion, RDKit, screening, docking, ranking)
conda env create -f autonomous_drug_discovery/envs/env_orchestrator.yml

# TargetDiff generation — pick the env spec by local GPU generation:
#   Ampere / Ada / Turing  ->  env_targetdiff.yml            (PyTorch 1.13 / cu117)
#   Blackwell (RTX 50xx)   ->  env_targetdiff_blackwell.yml  (PyTorch 2.8  / cu128)
conda env create -f autonomous_drug_discovery/envs/env_targetdiff_blackwell.yml

# Re-apply the pinned-submodule patches (NumPy + PyTorch 2.x compatibility)
scripts/apply_targetdiff_patches.sh
```

Both TargetDiff env specs install under the same conda env name (`targetdiff_env`); only the file differs, chosen per machine. You also need P2Rank installed (`~/p2rank_2.5.1/prank`) for pocket detection, and the TargetDiff checkpoint at `autonomous_drug_discovery/modules/02_generation/targetdiff/pretrained_models/pretrained_diffusion.pt` (mirrored at <https://huggingface.co/vladsft/agent-harness-weights>).

### Run — the `*-local` make targets (no Docker)

`make run` and `make dashboard` run *inside Docker*; for a no-Docker machine use their `*-local` counterparts, which call conda directly. Only `make pull` / `make push` (pure rclone) and these `*-local` targets work without Docker.

```bash
make run-local TARGET=1IEP MODE=cascade NUM=5      # full pipeline via conda
make dashboard-local                               # regenerate dashboard/ from telemetry
```

`MODE` accepts `cascade` (RDKit + TargetDiff), `rdkit`, `targetdiff`, `simulation`; `DEVICE` defaults to `auto` (detects a GPU, else CPU). Equivalent raw invocation:

```bash
conda run -n base python autonomous_drug_discovery/orchestrator.py run \
    autonomous_drug_discovery/data/processed/1IEP.pdb --mode cascade --device auto --num_samples 5
```

> **⚠️ TargetDiff on CPU is slow — ~12 min/molecule** (vs ~30 s on GPU). So on a CPU-only box: `rdkit` mode is ~2 min; a `cascade`/`targetdiff` run is only practical at a **small `NUM`** (e.g. `NUM=5` ≈ ~1 h end-to-end). Don't run `NUM=50` on CPU — that's ~10 h. The synthesizability gate also stays off unless you pass `--aizynth_config` (it needs AiZynthFinder + its ZINC stock); cascade + the pharmacophore bridge still run without it.

### No-Docker demo / parity checklist

To stand up this repo on another machine and both **run it** and **show the accumulated results**:

1. **Clone + submodules**
   ```bash
   git clone https://github.com/vladsft/autonomous-drug-discovery.git && cd autonomous-drug-discovery
   git submodule update --init --recursive && scripts/apply_targetdiff_patches.sh
   ```
2. **Copy `.env` manually** (scp / USB / password manager) — it holds the Cloudflare R2 + RunPod credentials and is gitignored. Without it, `make pull` can't reach R2. *This is the one step that can't be automated — secrets never travel through git or the data bucket.*
3. **Build the conda envs** (one-time, slow): `env_orchestrator.yml` (→ `base`) and an `env_targetdiff*.yml` (→ `targetdiff_env`) per the setup above. Install P2Rank + fetch the TargetDiff checkpoint.
4. **Pull data + telemetry from R2:** `make pull` — brings `telemetry.db`, every campaign's molecules/scores/rankings, and the `reports/` images into `data/`. (Heavy regenerable intermediates are excluded from the bucket; the demo-relevant outputs are all there.)
5. **Show the accumulated results** (no GPU needed): `make dashboard-local` then open `dashboard/index.html`; query history with `sqlite3 data/telemetry.db` (see [`docs/telemetry-guide.md`](docs/telemetry-guide.md)).
6. **One live run:** `make run-local TARGET=1IEP MODE=cascade NUM=5` (kick it off early — ~1 h on CPU — and walk through the dashboard while it runs).

For per-stage commands and parameters, see [`docs/pipeline-guide.md`](docs/pipeline-guide.md).

## License

See `LICENSE`.
