# Pipeline Usage Guide

## Overview

The pipeline takes a protein structure (PDB file) and produces a ranked list of candidate drug molecules through 4 stages:

```
PDB file  -->  [1. Pocket Detection]  -->  [2. Molecule Generation]  -->  [3. Screening]  -->  [4. Docking]  -->  Ranked candidates
```

Each stage reads the output of the previous stage and writes structured output (JSON manifests, SDF files, CSV results). Everything is logged to a SQLite telemetry database.

## Running the Full Pipeline

```bash
cd autonomous_drug_discovery

# Production mode: real generation, screening, and docking
conda run -n base python orchestrator.py run data/processed/1M17.pdb --mode production

# Simulation mode: stub data, tests pipeline plumbing only
conda run -n base python orchestrator.py run data/processed/1M17.pdb --mode simulation

# TargetDiff mode: diffusion-based generation (requires targetdiff_env)
conda run -n base python orchestrator.py run data/processed/1M17.pdb --mode targetdiff
```

Each full pipeline run creates a unique campaign ID (e.g. `campaign_fd4fad48`) and writes outputs to `data/<campaign_id>/`.

## Running Individual Stages

You can run each stage independently. This is useful for re-running a single stage with different parameters or debugging.

### Stage 1: Pocket Detection

Identifies binding pockets on the protein surface.

```bash
conda run -n base python orchestrator.py ingest data/processed/1M17.pdb
```

**Input:** `.pdb` file

**Output:** `data/processed/{stem}_manifest.json` containing:
- Pocket location (center coordinates)
- Pocket score and probability
- Path to pocket atom PDB file

**Backend selection:** P2Rank is the default. To use fpocket, edit the `run_ingestion` call or run the module directly:

```bash
conda run -n base python modules/01_ingestion/run_pocket.py \
  --pdb data/processed/1M17.pdb \
  --output_dir data/processed \
  --backend fpocket
```

### Stage 2: Molecule Generation

Generates candidate molecules that fit the detected pocket.

```bash
# RDKit fragment-based (default for production)
conda run -n base python orchestrator.py generate data/processed/1M17_manifest.json --mode rdkit

# TargetDiff diffusion model
conda run -n base python orchestrator.py generate data/processed/1M17_manifest.json --mode targetdiff

# Simulation stub
conda run -n base python orchestrator.py generate data/processed/1M17_manifest.json --mode simulation
```

**Input:** `manifest.json` from Stage 1

**Output:** `generated_molecules.sdf` containing 100 molecules (default) with 3D coordinates

**Parameters** (edit in `modules/02_generation/run_generation.py`):
- `num_samples`: number of molecules to generate (default: 100)
- `seed`: random seed for reproducibility (default: 42)

### Stage 3: Screening

Filters molecules for drug-likeness, synthetic accessibility, and toxicity.

```bash
conda run -n base python orchestrator.py screen data/candidates/generated_molecules.sdf
```

**Input:** `.sdf` file from Stage 2

**Output:**
- `screened_molecules.sdf` — molecules that passed all filters
- `screening_report.json` — per-molecule properties, pass/fail status, ADMET annotations
- `run_metadata.json` — execution metadata

**Filters applied** (configurable in `modules/03_screening/default_scoring_config.json`):

| Filter | Threshold | Rationale |
|---|---|---|
| Molecular weight | <= 500 | Lipinski Rule of Five |
| LogP | <= 5 | Lipinski Rule of Five |
| H-bond donors | <= 5 | Lipinski Rule of Five |
| H-bond acceptors | <= 10 | Lipinski Rule of Five |
| SA Score | <= 5.0 | Synthetic accessibility |
| QED | >= 0.3 | Drug-likeness |
| PAINS | = 0 | No promiscuous substructures |

ADMET-AI enrichment (104 properties) runs automatically on surviving molecules if installed.

**Customizing filters:** edit `default_scoring_config.json`. Each entry in `filter_thresholds` maps a property name to `{"max": N}`, `{"min": N}`, or `{"equals": N}`. No code change required.

### Stage 4: Docking

Estimates binding affinity of each molecule in the protein pocket.

```bash
conda run -n base python orchestrator.py dock data/processed/1M17_manifest.json --mode production
```

**Input:** `manifest.json` (for receptor + pocket center) and the screened SDF directory

**Output:**
- `docking_results.csv` — columns: `ligand_id, smiles, affinity` sorted by score
- `docked_*.pdbqt` — best pose for each ligand

**Modes:**
- `production` — full AutoDock Vina pipeline (PDB to PDBQT conversion, grid maps, exhaustive docking)
- `triage` — fast SMILES-based docking via TDC Oracle
- `simulation` — dummy scores for testing

**Parameters** (edit in `modules/04_docking/run_docking.py`):

```python
DEFAULT_PARAMS = {
    "exhaustiveness": 8,    # higher = more thorough, slower
    "num_modes": 9,         # binding poses to generate
    "energy_range": 3,      # kcal/mol range for pose clustering
    "box_size": [20, 20, 20],  # docking grid size in Angstroms
}
```

## Agent Planner (LLM-Driven)

The agent planner wraps the pipeline in an LLM loop that can inspect outputs between stages and adapt strategy:

```bash
# Deterministic mode (no API key needed)
conda run -n base python agent_planner.py --target data/processed/1M17.pdb --mode production

# With LLM adaptation (requires Gemini API key)
DISCOVERY_LLM_API_KEY=your_key conda run -n base python agent_planner.py \
  --target data/processed/1M17.pdb \
  --mode production \
  --max_iterations 5
```

Without an API key, the agent planner falls back to running the deterministic pipeline in sequence.

## Validation and Benchmarking

Three cancer targets with crystallographic ground truth are included in `data/processed/` (1M17, 2HYY, 6P3D). Run the benchmark comparison across all completed campaigns:

```bash
conda run -n base python benchmark.py
```

For detailed validation results, experimental methodology, and what "good" looks like, see [testing-guide.md](testing-guide.md).

## Output Directory Structure

When running the full pipeline (`orchestrator.py run`), outputs are organized per-campaign:

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

When running individual stages via `orchestrator.py ingest/generate/screen/dock`, outputs go to the default shared directories (`data/candidates/`, `data/screened/`, `data/results/`).

## Adding a New Target

1. Obtain a PDB file for your protein of interest from the [RCSB PDB](https://www.rcsb.org/)
2. Place it in `data/processed/`
3. Run the pipeline:

```bash
conda run -n base python orchestrator.py run data/processed/YOUR_TARGET.pdb --mode production
```

For best results, the PDB should contain a well-resolved structure (resolution < 3.0 A) of the protein in a relevant conformation. Co-crystal structures with a known ligand are ideal for validation.
