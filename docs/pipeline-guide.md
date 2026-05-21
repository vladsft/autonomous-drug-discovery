# Pipeline Usage Guide

## Overview

The pipeline takes a protein structure (PDB file) and produces a ranked list of candidate drug molecules through 5 stages:

```
PDB file  -->  [1. Pocket Detection]  -->  [2. Molecule Generation]  -->  [3. Screening]  -->  [4. Docking]  -->  [5. Ranking]  -->  Ranked candidates
```

Each stage reads the output of the previous stage and writes structured output (JSON manifests, SDF files, CSV results). Everything is logged to a SQLite telemetry database.

For setup, see the [README](../README.md). For architectural context (why Docker, why RunPod, why R2), see [`autonomous_drug_discovery/plan.md`](../autonomous_drug_discovery/plan.md).

## Running the Pipeline (Docker, recommended)

The repo ships a `Makefile` that wraps Docker. Once `make bootstrap` has run on a machine, everything below works.

### Full pipeline, local CPU

```bash
# Production (RDKit fragment-based generation, fast, no GPU needed)
make run TARGET=1M17 MODE=production

# Simulation (stub data, tests plumbing only)
make run TARGET=1M17 MODE=simulation
```

`TARGET` is the PDB stem; the file must exist at `data/processed/<TARGET>.pdb`. `MODE` is one of `simulation`, `production` (alias for `rdkit`), `rdkit`, `targetdiff`, or `pocket2mol`. Output lands in `data/campaign_<id>/`.

### Full pipeline, cloud GPU

For diffusion-based generation (TargetDiff), CPU is impractical (~12 min/mol). Use the cloud-run wrapper instead:

```bash
make cloud-run TARGET=2HYY MODE=targetdiff NUM=30
```

This provisions a RunPod GPU pod with the same Docker image, runs the full pipeline, syncs results to Cloudflare R2, and tears down the pod. Typical cost: ~$0.30 per run. See [`plan.md`](../autonomous_drug_discovery/plan.md#cloud-gpu-on-runpod-wrapped-in-a-script) for the architecture and `scripts/cloud_run.sh` for the implementation.

The cloud image's `targetdiff_env` is PyTorch 1.13 / CUDA 11.7, which runs on RunPod's Ampere GPUs (RTX 3090 / A40) unchanged.

### Batch runs (cloud, hands-free) — Phase 1.5

A `workflow_dispatch` from the GitHub Actions UI runs the full pipeline against a target list, no laptop involvement after kickoff. This section is the operator's manual.

#### One-time setup

Six repo secrets and one Pages toggle. Do these once, never again.

**Repo Settings → Secrets and variables → Actions → New repository secret:**

| Secret | Where to get it |
|---|---|
| `RUNPOD_API_KEY` | runpod.io → Settings → API Keys → "Create API Key" (Read + Write) |
| `RUNPOD_NETWORK_VOLUME_ID` | runpod.io → Storage → your Network Volume → ID field (the same one `make cloud-run` uses) |
| `R2_BUCKET` | The Cloudflare R2 bucket name, e.g. `agent-harness` |
| `R2_ACCESS_KEY_ID` | Cloudflare dashboard → R2 → Manage R2 API Tokens → Create API token |
| `R2_SECRET_ACCESS_KEY` | (shown once at token creation — save it) |
| `R2_ENDPOINT` | `https://<account-id>.r2.cloudflarestorage.com` (account ID is on the R2 overview page) |

**Repo Settings → Pages → Source = "GitHub Actions"** (one toggle, no further config).

Optional variables (Settings → Variables): `IMAGE` to override the GHCR tag, `RUNPOD_GPU_TYPE` to switch from the default 3090 (e.g. `NVIDIA RTX A6000`).

#### Firing a run

1. Open the repo on github.com → **Actions** tab → pick **batch_cloud_run** in the left sidebar.
2. Click the green **"Run workflow"** button on the right.
3. Fill the form (see field reference below), click the bottom green **"Run workflow"**.
4. The page refreshes; a new run appears at the top. Click into it to watch live logs.
5. Walk away. When the run finishes, you'll get the standard GitHub email notification ("workflow batch_cloud_run completed").

#### Field reference

| Field | Default | What to put |
|---|---|---|
| **targets** | (required) | Whitespace- or comma-separated PDB codes — e.g. `1M17 2HYY 6P3D 8P1L`. 4–8 alphanumeric characters each. `scripts/fetch_pdb.py` downloads any code not already cached. |
| **mode** | `targetdiff` | Which generator: `targetdiff` (diffusion, GPU-bound, ~30 min/target, best 3D fidelity), `rdkit` (fragment-based, CPU, ~2 min/target, broad chemical diversity), `pocket2mol` (pocket-aware GNN, GPU, ~7 s/mol — checkpoint not yet rehosted, so currently fails fast), `simulation` (stub, instant, plumbing-test only). |
| **num_samples** | `30` | Molecules generated per target. 30 is the validated default. Going to 100+ multiplies wallclock and cost roughly linearly. |
| **parallelism** | `5` | Concurrent RunPod pods. 5 is the sweet spot for the free tier. Bump to 10 if you have ample credit; >10 risks RunPod balance starvation mid-batch. |
| **force** | `false` | If `true`, ignores the 24-hour idempotency skip-list and re-runs every target. Use this when you change pipeline params; otherwise leave at `false` to avoid paying for work you already have. |

##### Recommended first-time run

Before running 30-50 targets, do a small live shake-out: `targets="1M17 2HYY"`, `mode=rdkit`, `num_samples=5`, `parallelism=2`. Costs cents, finishes in ~5 minutes, validates that all six secrets are wired correctly. The dashboard at `https://vladsft.github.io/autonomous-drug-discovery/` updates with the two targets when it lands.

#### Timeline of a real run (parallelism=5, 30 targets, targetdiff)

| T+      | What happens |
|---------|--------------|
| 0 s     | You click Run workflow. Action runner spins up (~30 s GitHub-side). |
| 30 s    | Python + rclone install; secrets exported to env. |
| 1 min   | `batch_cloud_run.py` starts. Pulls `telemetry.db` from R2; builds skip-list. Fetches any missing PDBs from RCSB. Queries RunPod balance. |
| 2 min   | First 5 RunPod pods provisioned. Each pulls its input PDB from R2, runs orchestrator end-to-end. |
| ~30 min | First wave of pods finishes (one wave at a time on `parallelism=5`). `.done` sentinels appear in R2. Dispatcher cycles in the next 5. |
| ~3 h    | Last sentinel lands. Dispatcher syncs all campaign output from R2, writes `data/batch_summary.json`. |
| ~3 h 2 m | `regenerate_dashboard.py --all-targets` runs, writes `dashboard/professor_demo.{js,json}`. |
| ~3 h 3 m | Commits `dashboard/` back to `main` with a "batch: refresh dashboard" message; pushes. |
| ~3 h 4 m | Inline Pages deploy. URL goes live. |
| ~3 h 5 m | GitHub emails you. |

If any pod fails its first attempt, the dispatcher retries it up to 3 times. Pods that fail all 3 attempts get a `.failed` sentinel in R2 and are reported in `batch_summary.json` as `outcome: failed`. The rest of the batch is unaffected.

#### Where the results live

Three places, each with a different audience:

1. **The chemist dashboard** — `https://vladsft.github.io/autonomous-drug-discovery/`
    The fastest path. Multi-target dropdown in the header; each target carries a per-backend tab (RDKit / TargetDiff / Pocket2Mol), a sortable molecule table, and a detail panel with SVG structure + drug-likeness + ADMET + synthesis route. This is what you'd show a medicinal chemist.

2. **Cloudflare R2 bucket** — `r2:<R2_BUCKET>/`
    The raw artefacts. For every campaign the pipeline produced, you can pull down:
    - `r2:<bucket>/campaign_<id>/candidates/generated_molecules.sdf` — raw 3D coordinates from the generator (open in PyMOL: `pymol generated_molecules.sdf`)
    - `r2:<bucket>/campaign_<id>/screened/screened_molecules.sdf` — survivors after drug-likeness filters
    - `r2:<bucket>/campaign_<id>/screened/screening_report.json` — per-molecule properties + ADMET + attrition counts
    - `r2:<bucket>/campaign_<id>/results/docking_results.csv` — sorted by Vina affinity
    - `r2:<bucket>/campaign_<id>/results/docked_mol_*.pdbqt` — the best binding pose per ligand (also openable in PyMOL)
    - `r2:<bucket>/campaign_<id>/ranked/ranked_candidates.json` — the final scorecard with composite score
    - `r2:<bucket>/telemetry.db` — SQLite DB of every run

    Sync locally with `make pull` (or directly: `rclone copy r2:<bucket>/campaign_<id> ./local-dir`).

3. **GitHub Action artifact** — `batch-summary` (attached to the workflow run)
    Machine-readable per-target outcome list. Download from the run summary page → "Artifacts" section. Schema:
    ```jsonc
    {
      "generated_at": "2026-05-21T20:00:00Z",
      "mode": "targetdiff", "num_samples": 30, "parallelism": 5,
      "results": [
        { "target": "1M17", "outcome": "done",   "sentinel_key": "1M17-targetdiff-abc1",
          "pod_id": "pod_xyz", "tail": "[pod] === done (rc=0) ===" },
        { "target": "2HYY", "outcome": "failed", "sentinel_key": "2HYY-targetdiff-def2",
          "pod_id": "pod_uvw", "tail": "[pod] FATAL: pipeline exceeded 85-minute fuse." },
        ...
      ]
    }
    ```

#### What the results actually look like

You get **per-molecule structured data plus aggregate per-backend statistics, with rendered 2D structures, for every target × generator combination that produced output.**

##### Per-molecule record (one entry per generated molecule)

```jsonc
{
  "molecule_id": "mol_0012",
  "smiles": "COc1nc2c(C(N)=O)cccc2nc1O",            // canonical SMILES
  "svg": "<svg …>",                                  // pre-rendered 2D structure
  "screening_passed": true,
  "rejected_reason": null,                           // populated when screening_passed=false

  // Drug-likeness, from RDKit. All filters in screening config feed this.
  "properties": {
    "mol_weight": 219.2,            // Da
    "logp": 0.44,                    // Crippen logP
    "qed": 0.7632,                   // 0-1, higher = more drug-like
    "sa_score": 2.0,                 // 1-10, lower = easier to synthesise
    "tpsa": 98.33,                   // topological polar surface area, Å²
    "hbd": 2,                        // H-bond donors
    "hba": 5,                        // H-bond acceptors
    "rotatable_bonds": 2,
    "heavy_atoms": 16,
    "pains_alerts": 0                // PAINS substructure count (0 = clean)
  },

  // ADMET-AI, 11 endpoints exposed (104 total available via the raw screening_report).
  "admet": {
    "hERG": 0.061,            // hERG channel block probability (lower = safer)
    "AMES": 0.2616,           // mutagenicity probability
    "DILI": 0.9419,           // drug-induced liver injury risk
    "CYP2D6": 0.0057,         // metabolism interference (lower better)
    "CYP3A4": 0.003,
    "Caco2": -4.9759,         // permeability (Caco-2)
    "HIA": 0.9969,            // human intestinal absorption
    "BBB": 0.8098,            // blood-brain barrier
    "Clearance": 43.6822,
    "Bioavailability": 0.9339,
    "LD50": 1.6498
  },
  "admet_flags": {              // simple pass/fail for the dashboard's red/green badges
    "hERG": true, "AMES": true, "DILI": false
  },

  // AutoDock Vina, kcal/mol, more negative = stronger binding.
  "docking_score": -1.815,

  // AiZynthFinder retrosynthesis, only run for top-N by docking + ADMET.
  "synthesis": {
    "evaluated": true,
    "route_found": true,
    "n_steps": 1,                 // synthesis steps
    "n_routes": 110,              // total alternative routes found
    "precursors": [               // SMILES of the building blocks
      "COC(=O)C(=O)O",
      "NC(=O)c1cccc(N)c1N"
    ],
    "route_tree": { … }           // full retrosynthetic tree, renderable in the dashboard
  },

  // Multi-criteria final score: 0.5·dock + 0.3·ADMET + 0.2·synthesis, normalised.
  "composite_score": 0.8421,
  "final_rank": 1
}
```

##### Per-backend summary (one entry per target × backend)

```jsonc
{
  "total_generated": 30,           // molecules out of generator
  "passed_screening": 28,          // survived all Lipinski/QED/SA/PAINS filters
  "synthesis_evaluated": 15,       // top-N that AiZynth saw
  "synthesis_routes_found": 8,     // of those, how many had at least one route
  "best_docking_score": -9.47,     // most negative Vina score
  "median_docking_score": -6.82,
  "mean_qed": 0.62,
  "mean_sa": 2.91
}
```

So per target you get **N molecules × {SMILES + 2D structure + 10 drug-likeness properties + 11 ADMET predictions + Vina docking + synthesis routes for the top survivors + composite ranking}**, plus **aggregate statistics on the cohort**. Multiply by however many generator backends ran on that target (up to 3: RDKit + TargetDiff + Pocket2Mol).

A 50-target batch with `mode=targetdiff num_samples=30` produces ~1,500 ranked molecules in the dashboard, each with the schema above.

#### Idempotency, cost, and failure handling

**Idempotency.** Any target whose generation stage in the requested `mode` succeeded in the past 24 h is silently skipped — repeated workflow invocations during one work session don't re-pay for finished work. Override with `force: true`.

**Cost guard.** Before any pod is provisioned, the dispatcher queries the RunPod balance and refuses to start if it's less than the worst-case spend (`parallelism × waves × 1.5 h × $0.50 × 1.5`). Worst-case assumes every pod runs to its 90-min ceiling; actual cost is usually 3-4× lower because pods finish well before the fuse. A 30-target batch on 5 pods historically lands at ~$3-4 of RunPod credit.

**Failure handling.** Three retries per target. A pod reclaimed by RunPod's preemption logic is re-provisioned automatically. A pod whose pipeline truly fails (bad PDB, no detectable pocket, OOM, etc.) writes a `.failed` sentinel — the dispatcher records it and the batch continues. The summary in the Action's job log lists every target's outcome at the end, and the same data lives in the `batch-summary` artifact.

#### Debugging a run that didn't produce a dashboard update

1. Open the failed Action run. Check the "Run the batch" step's tail for the per-target outcome summary.
2. Look for `.failed` sentinels in R2: `rclone cat r2:<bucket>/sentinels/<sentinel-key>.failed` — they carry the exit code and the last 50 log lines of the pod's pipeline.
3. If a target consistently fails the first stage (P2Rank), the PDB likely has no detectable pocket — try a co-crystal structure instead of an apo form.
4. If a target consistently times out, your `mode` is probably wrong for the GPU (TargetDiff on a non-CUDA-11.7 GPU will fall back to CPU and exceed the fuse) — pick a different RunPod GPU type via the `RUNPOD_GPU_TYPE` repo variable, or switch the target to `mode=rdkit`.
5. If the workflow itself errors out before pods even provision, the most common cause is a missing or mis-named secret — re-check the six values listed in [One-time setup](#one-time-setup).

### Full pipeline, local GPU (no Docker)

TargetDiff diffusion can also run on a local NVIDIA GPU without Docker or RunPod — the orchestrator dispatches each stage straight into conda environments. This is the path to use when the cu117 Docker image cannot target your GPU, notably **NVIDIA Blackwell (RTX 50-series, `sm_120`)**.

Pick the TargetDiff environment by GPU generation — both install under the conda env name `targetdiff_env`, only the file differs:

| Local GPU | Env spec | PyTorch / CUDA |
|---|---|---|
| Ampere / Ada / Turing (RTX 20/30/40, A-series) | `envs/env_targetdiff.yml` | 1.13 / cu117 |
| Blackwell (RTX 50-series, `sm_120`) | `envs/env_targetdiff_blackwell.yml` | 2.8 / cu128 |

```bash
# One-time setup
conda env create -f autonomous_drug_discovery/envs/env_orchestrator.yml          # base
conda env create -f autonomous_drug_discovery/envs/env_targetdiff_blackwell.yml  # or env_targetdiff.yml
scripts/apply_targetdiff_patches.sh        # PyTorch 2.x + NumPy submodule patches

# Run the full pipeline on the local GPU
conda run -n base python autonomous_drug_discovery/orchestrator.py run \
    autonomous_drug_discovery/data/processed/<TARGET>.pdb \
    --mode targetdiff --device cuda --num_samples 30
```

`--device` accepts `auto` (default — detects a GPU), `cuda`, or `cpu`. P2Rank must also be installed locally (see Stage 1). The Docker image is unaffected; it keeps its cu117 environment for the cloud path.

### Syncing campaign results across machines

```bash
make pull       # download data/ from R2
make push       # upload data/ to R2
```

R2 is the canonical store. Treat your local `data/` as a cache.

## Running Individual Stages

Each stage can be run independently against an existing manifest, useful for re-running with different parameters or debugging:

```bash
# Stage 1 — pocket detection
docker run --rm -v $(pwd)/data:/app/data ghcr.io/vladsft/autonomous-drug-discovery \
  orchestrator.py ingest data/processed/1M17.pdb

# Stage 2 — generation
docker run --rm -v $(pwd)/data:/app/data ghcr.io/vladsft/autonomous-drug-discovery \
  orchestrator.py generate data/processed/1M17_manifest.json --mode rdkit

# Stage 3 — screening
docker run --rm -v $(pwd)/data:/app/data ghcr.io/vladsft/autonomous-drug-discovery \
  orchestrator.py screen data/candidates/generated_molecules.sdf

# Stage 4 — docking
docker run --rm -v $(pwd)/data:/app/data ghcr.io/vladsft/autonomous-drug-discovery \
  orchestrator.py dock data/processed/1M17_manifest.json --mode production

# Stage 5 — ranking
docker run --rm -v $(pwd)/data:/app/data ghcr.io/vladsft/autonomous-drug-discovery \
  orchestrator.py rank data/results/docking_results.csv --screening_json data/screened/screening_report.json
```

These work the same way inside `make cloud-run` — the entrypoint is identical.

## Stage Reference

### Stage 1: Pocket Detection

**Input:** `.pdb` file
**Output:** `data/processed/{stem}_manifest.json` with pocket location, score, probability, and pocket atom PDB path.
**Backend:** P2Rank (ML-based, default). fpocket is a geometry-based fallback selectable via `--backend fpocket`.

### Stage 2: Molecule Generation

**Input:** `manifest.json` from Stage 1
**Output:** `generated_molecules.sdf` with 3D coordinates

**Backends:**

| Mode | Paradigm | Speed (per mol) | Pocket-aware | GPU needed? |
|---|---|---|---|---|
| `rdkit` | Fragment combinatorial | ~10 ms | Size-only | No |
| `targetdiff` | E(3)-equivariant diffusion | ~30 s GPU / ~12 min CPU | Yes (3D) | Yes, in practice |
| `pocket2mol` | Autoregressive GNN | ~7 s GPU | Yes (3D) | Yes, in practice |
| `simulation` | Stub (single benzene) | instant | No | No |

**Pocket2Mol is currently deferred.** Its pretrained checkpoint was hosted on a Google Drive folder that has since been deleted, and we never recovered a local copy. The mode is wired into the orchestrator but will fail without the checkpoint. See `plan.md` for the recovery plan.

**Parameters** (edit in `modules/02_generation/run_generation.py`'s `_DEFAULT_PARAMS_BY_MODE`, or pass `--num_samples` to the orchestrator):
- `num_samples`: number of molecules to generate — pass `--num_samples` to the orchestrator
- `device`: compute device for TargetDiff / Pocket2Mol — pass `--device auto|cpu|cuda` (`auto` detects a GPU)
- `seed`: random seed for reproducibility
- `sampling_steps`: denoising steps for TargetDiff (default 1000)

### Stage 3: Screening

**Input:** `.sdf` file from Stage 2
**Output:** `screened_molecules.sdf` (survivors), `screening_report.json` (per-molecule properties + ADMET annotations + attrition).

**Filters applied** (configurable in `modules/03_screening/default_scoring_config.json`):

| Filter | Threshold | Rationale |
|---|---|---|
| Molecular weight | ≤ 500 | Lipinski Rule of Five |
| LogP | ≤ 5 | Lipinski Rule of Five |
| H-bond donors | ≤ 5 | Lipinski Rule of Five |
| H-bond acceptors | ≤ 10 | Lipinski Rule of Five |
| SA Score | ≤ 5.0 | Synthetic accessibility |
| QED | ≥ 0.3 | Drug-likeness |
| PAINS | = 0 | No promiscuous substructures |

ADMET-AI enrichment (104 properties) runs automatically on surviving molecules.

**Note on current thresholds.** Survival rates are currently 73-98%, which is too loose for a useful candidate shortlist. Phase 1 of the plan tightens these to target 40-60% survival; Phase 2 calibrates them against expert input.

### Stage 4: Docking

**Input:** `manifest.json` (for receptor + pocket centre) and the screened SDF directory
**Output:** `docking_results.csv` (sorted by affinity), `docked_*.pdbqt` (best pose per ligand)

**Modes:**
- `production` — full AutoDock Vina pipeline (PDB to PDBQT, grid maps, exhaustive docking)
- `triage` — fast SMILES-based docking via TDC Oracle (requires `pytdc`)
- `simulation` — dummy scores for testing

**Parameters** (edit in `modules/04_docking/run_docking.py`):
```python
DEFAULT_PARAMS = {
    "exhaustiveness": 8,    # higher = more thorough, slower
    "num_modes": 9,         # binding poses to generate
    "energy_range": 3,      # kcal/mol range for pose clustering
    "box_size": [20, 20, 20],  # docking grid size in Ångstroms
}
```

### Stage 5: Multi-criteria Ranking

**Input:** `docking_results.csv` + `screening_report.json`
**Output:** `ranked_candidates.json`

Combines docking score, ADMET profile, and (optionally) AiZynthFinder retrosynthetic feasibility into a composite score: `0.5 × docking + 0.3 × ADMET + 0.2 × synthesis`. Pass `--aizynth_config <config.yml>` to score the top-N candidates with AiZynthFinder, which runs in the separate `aizynth_env` (`envs/env_aizynth.yml`); without it the synthesis term stays neutral and contributes a constant.

## Agent Planner (LLM-driven, experimental)

`agent_planner.py` wraps the pipeline in an LLM loop that can inspect outputs between stages and adapt strategy:

```bash
# Deterministic mode (no API key needed)
docker run --rm -v $(pwd)/data:/app/data ghcr.io/vladsft/autonomous-drug-discovery \
  agent_planner.py --target data/processed/1M17.pdb --mode production

# With LLM adaptation
docker run --rm -v $(pwd)/data:/app/data -e DISCOVERY_LLM_API_KEY ghcr.io/vladsft/autonomous-drug-discovery \
  agent_planner.py --target data/processed/1M17.pdb --mode production --max_iterations 5
```

Without an API key the planner falls back to the deterministic pipeline. The full agent design (cascade + Bayesian recommender + Sonnet-in-the-loop + Obsidian) is Phase 3 of the plan.

## Validation and Benchmarking

Three cancer targets with crystallographic ground truth ship in `data/processed/` (1M17, 2HYY, 6P3D). The benchmark script compares all completed campaigns:

```bash
docker run --rm -v $(pwd)/data:/app/data ghcr.io/vladsft/autonomous-drug-discovery benchmark.py
```

For detailed validation results, see [testing-guide.md](testing-guide.md).

## Output Directory Structure

```
data/
  processed/                     # Shared across campaigns (keyed by PDB stem)
    1M17_manifest.json
    1M17_p2rank/
  campaign_fd4fad48/             # Per-campaign isolation
    candidates/
      generated_molecules.sdf
      run_metadata.json
    screened/
      screened_molecules.sdf
      screening_report.json
      run_metadata.json
    results/
      docking_results.csv
      docked_mol_0000.pdbqt
      run_metadata.json
  telemetry.db                   # Shared telemetry database
```

Per-stage invocations write to `data/candidates/`, `data/screened/`, `data/results/` instead (no campaign isolation). Use the full `make run` / `make cloud-run` path for proper isolation.

## Adding a New Target

1. Obtain a PDB file from <https://www.rcsb.org/>.
2. Place it in `data/processed/`.
3. Run the pipeline:
   ```bash
   make run TARGET=<stem> MODE=production
   ```

For best results the PDB should be well-resolved (<3.0 Å) in a relevant conformation. Co-crystal structures with a known ligand are ideal for validation.
