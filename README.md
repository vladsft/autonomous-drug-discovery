# Autonomous Drug Discovery Pipeline

An end-to-end computational drug discovery system that takes a protein structure (PDB file) and produces a ranked list of candidate drug molecules. The pipeline automates pocket detection, molecule generation, drug-likeness screening, and molecular docking — with every step logged to a telemetry database for reproducibility and analysis.

## What it does

Given a protein target, the pipeline:

1. **Detects binding pockets** on the protein surface using fpocket
2. **Generates candidate molecules** using fragment-based combinatorial chemistry (RDKit), sized to fit the detected pocket
3. **Screens candidates** against drug-likeness filters (Lipinski Rule of Five, QED, synthetic accessibility, PAINS toxicity alerts)
4. **Docks survivors** into the binding pocket using AutoDock Vina to estimate binding affinity
5. **Ranks and reports** the results, with all intermediate data logged to SQLite

The system has been validated against EGFR (PDB: 1M17), where it correctly identifies the erlotinib binding site and produces molecules with physically meaningful docking scores (-3.9 to -9.1 kcal/mol).

## Quick start

### Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- [fpocket](https://github.com/Discngine/fpocket) binary installed at `/home/<user>/fpocket/bin/fpocket` (or update the path in `modules/01_ingestion/run_pocket.py`)

### Install dependencies

```bash
conda install -n base -c conda-forge rdkit vina -y
pip install meeko gemmi
```

### Run the full pipeline

```bash
cd autonomous_drug_discovery

# Production mode: real generation, screening, and docking
conda run -n base python orchestrator.py run data/processed/1M17.pdb --mode production

# Simulation mode: stub data, for testing the pipeline plumbing
conda run -n base python orchestrator.py run data/processed/1UYD.pdb --mode simulation
```

### Run individual stages

```bash
conda run -n base python orchestrator.py ingest data/processed/1M17.pdb
conda run -n base python orchestrator.py generate data/processed/1M17_manifest.json --mode rdkit
conda run -n base python orchestrator.py screen data/candidates/generated_molecules.sdf
conda run -n base python orchestrator.py dock data/processed/1M17_manifest.json --mode production
```

## Repository structure

```
.
├── autonomous_drug_discovery/       # Main application
│   ├── orchestrator.py              # CLI entrypoint — runs stages individually or as full pipeline
│   ├── agent_planner.py             # LLM-driven adaptive orchestration (optional, requires API key)
│   ├── telemetry.py                 # SQLite telemetry database (runs + molecule_scores tables)
│   ├── plan.md                      # Architectural plan and roadmap
│   │
│   ├── modules/
│   │   ├── 01_ingestion/
│   │   │   └── run_pocket.py        # fpocket wrapper — detects binding pockets on a PDB file
│   │   ├── 02_generation/
│   │   │   └── run_generation.py    # Molecule generator — fragment-based (RDKit) or TargetDiff
│   │   ├── 03_screening/
│   │   │   ├── run_screening.py     # Drug-likeness filters — Lipinski, QED, SA, PAINS
│   │   │   └── default_scoring_config.json  # Filter thresholds (editable)
│   │   └── 04_docking/
│   │       └── run_docking.py       # AutoDock Vina docking — simulation, or production mode
│   │
│   ├── data/
│   │   └── processed/
│   │       ├── 1M17.pdb             # EGFR kinase (erlotinib co-crystal) — validation target
│   │       └── 1UYD.pdb             # UDP-GlcNAc epimerase — test target
│   │
│   ├── envs/                        # Conda environment specs (for reference)
│   │   ├── env_orchestrator.yml
│   │   ├── env_docking.yml
│   │   └── env_targetdiff.yml
│   │
│   └── tests/
│       ├── test_screening.py
│       └── test_telemetry.py
│
├── data/                            # Reference data (CIF dictionaries, utilities)
│   ├── keepResidues.txt
│   └── extractModifiedResidueCodes.py
│
└── reports/
    ├── north_star.md
    └── config.js
```

## Pipeline stages in detail

### Stage 1: Pocket detection (`01_ingestion/run_pocket.py`)

**Tool:** fpocket (Voronoi tessellation + alpha sphere clustering)

**Input:** A `.pdb` file containing a protein structure.

**Output:**
- `{stem}_manifest.json` — lists all detected pockets, ranked by druggability score. The best pocket (pocket1) is selected automatically.
- `{stem}_out/pockets/` — individual pocket PDB files (`pocket1_atm.pdb` = protein atoms lining the pocket, `pocket1_vert.pqr` = cavity geometry and scoring metadata).

**Key metrics in the `.pqr` headers:** Drug Score (0-1), Pocket Score, Volume, Hydrophobicity, Polarity.

**Validation:** On EGFR (1M17), pocket1 captures 50% of the known erlotinib binding residues and the docking box center falls 6.3 Angstroms from the co-crystallized ligand position.

### Stage 2: Molecule generation (`02_generation/run_generation.py`)

**Modes:**
- `simulation` — writes a single stub molecule (benzene) for pipeline testing.
- `rdkit` — fragment-based combinatorial generation using RDKit. Assembles drug-like scaffolds (indole, quinazoline, piperidine, etc.) with functional group substituents and linkers. Uses BRICS decomposition for additional diversity. Molecules are sized to fit the pocket (estimated from pocket radius). 3D conformers are generated and MMFF-optimized.
- `targetdiff` — wraps the [TargetDiff](https://github.com/guanjq/targetdiff) diffusion model. Requires a separate conda environment and pretrained checkpoint (not included). Reserved for future use.

**Input:** `manifest.json` from Stage 1 (specifically the `best_pocket` path).

**Output:** `generated_molecules.sdf` — 100 molecules by default, each with 3D coordinates, SMILES, and a unique `molecule_id`.

**Configuration:** `num_samples` (default 100), `seed` (default 42).

### Stage 3: Screening (`03_screening/run_screening.py`)

**Tool:** RDKit (with optional MolScore backend if installed)

**Input:** An `.sdf` file of candidate molecules.

**Output:**
- `screened_molecules.sdf` — only molecules that passed all filters.
- `screening_report.json` — per-molecule properties, pass/fail status, and attrition summary.

**Filters applied (configurable in `default_scoring_config.json`):**

| Filter | Threshold | Rationale |
|--------|-----------|-----------|
| Molecular weight | <= 500 | Lipinski Rule of Five |
| LogP | <= 5 | Lipinski Rule of Five |
| H-bond donors | <= 5 | Lipinski Rule of Five |
| H-bond acceptors | <= 10 | Lipinski Rule of Five |
| SA Score | <= 5.0 | Synthetic accessibility |
| QED | >= 0.3 | Drug-likeness |
| PAINS | = 0 | No toxic substructure alerts |

**Adding a new filter:** Edit `default_scoring_config.json`. Each entry in `filter_thresholds` maps a property name to `{"max": N}`, `{"min": N}`, or `{"equals": N}`. No code change required.

### Stage 4: Docking (`04_docking/run_docking.py`)

**Tool:** AutoDock Vina (Python API) + Meeko (ligand PDBQT preparation)

**Modes:**
- `simulation` — returns hardcoded dummy scores.
- `production` — full Vina docking pipeline: converts receptor PDB to PDBQT with proper AutoDock atom typing (C, A, N, NA, OA, SA, HD), prepares each ligand via Meeko, computes Vina grid maps centered on the pocket centroid, docks with exhaustiveness=8 and 9 poses per ligand, reports the best binding affinity.

**Input:** `manifest.json` (for receptor PDB and pocket centroid) + candidates directory containing an `.sdf` file.

**Output:**
- `docking_results.csv` — columns: `ligand_id, smiles, affinity` (kcal/mol, more negative = stronger binding). Sorted best-to-worst.
- Docked pose files (`docked_mol_XXXX.pdbqt`) for the top pose of each ligand.

**Typical score ranges:** Drug-like molecules against real targets score between -4 and -11 kcal/mol. Scores near 0 indicate a problem with receptor preparation.

## Telemetry database

All pipeline runs are logged to `data/telemetry.db` (SQLite). Two tables:

**`runs`** — one row per module execution:
- `run_id` (UUID), `campaign_id`, `module_name`, `started_at`, `completed_at`, `status`
- `input_path`, `output_path`, `parameters` (JSON), `error_trace`, `git_commit`

**`molecule_scores`** — one row per molecule per stage:
- `molecule_id`, `smiles`, `qed`, `sa_score`, `logp`, `mol_weight`, `docking_score`
- `passed_triage` (1/0), `stage_eliminated` (reason string)

Query examples:
```bash
# All campaigns and their status
sqlite3 data/telemetry.db "SELECT campaign_id, module_name, status FROM runs ORDER BY started_at;"

# Top molecules across all campaigns
sqlite3 data/telemetry.db "SELECT smiles, docking_score, qed FROM molecule_scores WHERE docking_score IS NOT NULL ORDER BY docking_score LIMIT 10;"

# Attrition by stage
sqlite3 data/telemetry.db "SELECT stage_eliminated, COUNT(*) FROM molecule_scores WHERE passed_triage = 0 GROUP BY stage_eliminated;"
```

## Agent planner

`agent_planner.py` wraps the pipeline in an LLM-driven loop that can inspect outputs between stages and adapt strategy. Currently supports Google Gemini (requires `DISCOVERY_LLM_API_KEY` environment variable). Without an API key, it falls back to the deterministic pipeline.

```bash
# Deterministic mode (no API key needed)
conda run -n base python agent_planner.py --target data/processed/1M17.pdb --mode production

# With LLM adaptation
DISCOVERY_LLM_API_KEY=your_key conda run -n base python agent_planner.py --target data/processed/1M17.pdb --mode production --max_iterations 5
```

The agent planner is a **recommendation engine**, not an autonomous decision-maker. See `plan.md` Layer 4 for the design philosophy.

## Configuration

### Screening thresholds

Edit `modules/03_screening/default_scoring_config.json`. Example — loosening the QED filter:

```json
"desc_QED": {"min": 0.2, "reason": "QED Drug-likeness >= 0.2"}
```

### Docking parameters

In `modules/04_docking/run_docking.py`, `DEFAULT_PARAMS`:

```python
DEFAULT_PARAMS = {
    "exhaustiveness": 8,    # higher = more thorough but slower
    "num_modes": 9,         # number of binding poses to generate
    "energy_range": 3,      # kcal/mol range for pose clustering
    "box_size": [20, 20, 20],  # docking grid size in Angstroms
}
```

### Generation parameters

In `modules/02_generation/run_generation.py`, `DEFAULT_PARAMS`:

```python
DEFAULT_PARAMS = {
    "num_samples": 100,  # number of molecules to generate
}
```

The fragment library (scaffolds, substituents, linkers) is defined at the top of the same file and can be extended.

## Current status and roadmap

### Working now
- Pocket detection (fpocket) — real, validated
- Molecule generation (RDKit fragment-based) — real, pocket-aware
- Screening (RDKit) — real, Lipinski + QED + SA + PAINS
- Docking (AutoDock Vina) — real, production-grade PDBQT preparation
- Campaign telemetry — full logging of all stages
- Validated against EGFR (1M17) with erlotinib binding site recovery

### Planned next (see `plan.md` for full roadmap)
- Multi-target validation benchmark (EGFR, BCR-ABL, BRAF V600E)
- TargetDiff integration for structure-aware diffusion-based generation
- ADMET prediction (absorption, metabolism, toxicity)
- GNINA CNN-based rescoring
- Retrosynthetic accessibility analysis (AiZynthFinder)
- Agent planner with empirical strategy selection

## Dependencies

| Package | Purpose | Install |
|---------|---------|---------|
| RDKit | Molecular property calculation, fragment generation, PAINS filters | `conda install -c conda-forge rdkit` |
| AutoDock Vina | Molecular docking (binding affinity scoring) | `conda install -c conda-forge vina` |
| Meeko | Ligand PDBQT preparation for Vina | `pip install meeko` |
| gemmi | Receptor PDB parsing and atom typing | `pip install gemmi` |
| fpocket | Binding pocket detection | [Build from source](https://github.com/Discngine/fpocket) |

## License

See repository license file.
