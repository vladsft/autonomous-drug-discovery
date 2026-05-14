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
docker run --rm -v $(pwd)/data:/app/data ghcr.io/<you>/agent-harness \
  orchestrator.py ingest data/processed/1M17.pdb

# Stage 2 — generation
docker run --rm -v $(pwd)/data:/app/data ghcr.io/<you>/agent-harness \
  orchestrator.py generate data/processed/1M17_manifest.json --mode rdkit

# Stage 3 — screening
docker run --rm -v $(pwd)/data:/app/data ghcr.io/<you>/agent-harness \
  orchestrator.py screen data/candidates/generated_molecules.sdf

# Stage 4 — docking
docker run --rm -v $(pwd)/data:/app/data ghcr.io/<you>/agent-harness \
  orchestrator.py dock data/processed/1M17_manifest.json --mode production

# Stage 5 — ranking
docker run --rm -v $(pwd)/data:/app/data ghcr.io/<you>/agent-harness \
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
- `num_samples`: number of molecules to generate (default varies per backend)
- `seed`: random seed for reproducibility
- `device`: `cpu` or `cuda` (used by TargetDiff / Pocket2Mol)
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

Combines docking score, ADMET profile, and (optionally) AiZynthFinder retrosynthetic feasibility into a composite score: `0.5 × docking + 0.3 × ADMET + 0.2 × synthesis`. AiZynth is currently a stubbed slot; populated on top-N candidates only.

## Agent Planner (LLM-driven, experimental)

`agent_planner.py` wraps the pipeline in an LLM loop that can inspect outputs between stages and adapt strategy:

```bash
# Deterministic mode (no API key needed)
docker run --rm -v $(pwd)/data:/app/data ghcr.io/<you>/agent-harness \
  agent_planner.py --target data/processed/1M17.pdb --mode production

# With LLM adaptation
docker run --rm -v $(pwd)/data:/app/data -e DISCOVERY_LLM_API_KEY ghcr.io/<you>/agent-harness \
  agent_planner.py --target data/processed/1M17.pdb --mode production --max_iterations 5
```

Without an API key the planner falls back to the deterministic pipeline. The full agent design (cascade + Bayesian recommender + Sonnet-in-the-loop + Obsidian) is Phase 3 of the plan.

## Validation and Benchmarking

Three cancer targets with crystallographic ground truth ship in `data/processed/` (1M17, 2HYY, 6P3D). The benchmark script compares all completed campaigns:

```bash
docker run --rm -v $(pwd)/data:/app/data ghcr.io/<you>/agent-harness benchmark.py
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
