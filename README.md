# Autonomous Drug Discovery Pipeline

An end-to-end computational drug discovery system that takes a protein structure (PDB file) and produces a ranked list of candidate drug molecules. The pipeline automates pocket detection, molecule generation, drug-likeness screening, and molecular docking — with every step logged to a telemetry database for reproducibility and analysis.

## What it does

Given a protein target, the pipeline:

1. **Detects binding pockets** on the protein surface using P2Rank (ML-based, default) or fpocket (fallback)
2. **Generates candidate molecules** using one of three backends:
   - RDKit fragment-based combinatorial chemistry (default, fast, CPU)
   - Pocket2Mol autoregressive generator (pocket-conditioned, ~7 s/molecule on GPU)
   - TargetDiff E(3)-equivariant diffusion (pocket-conditioned, slower, highest-fidelity 3D)
3. **Screens candidates** against drug-likeness filters (Lipinski, QED, SA, PAINS) and annotates with ADMET-AI predictions (104 properties: toxicity, absorption, metabolism, etc.)
4. **Docks survivors** into the binding pocket using AutoDock Vina to estimate binding affinity
5. **Ranks and reports** the results, with all intermediate data logged to SQLite

Validated against three cancer targets with crystallographic ground truth:

| Target | Disease | Best Dock Score | Pocket Accuracy |
|--------|---------|----------------|-----------------|
| EGFR (1M17) | Lung cancer | -9.32 kcal/mol | 2.7 A from erlotinib, 82% residue overlap |
| BCR-ABL (2HYY) | Leukemia | -12.56 kcal/mol | 2.7 A from imatinib, 92% residue overlap |
| BRAF V600E (6P3D) | Melanoma | -10.40 kcal/mol | 3.1 A from ponatinib, 89% residue overlap |

## Quick start

### Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- [P2Rank](https://github.com/rdk/p2rank) (default pocket detection) — requires Java 17+
- [fpocket](https://github.com/Discngine/fpocket) (fallback pocket detection)

### Install dependencies

```bash
# Step 1 — conda packages (compiled extensions)
conda install -n base -c conda-forge rdkit vina openjdk=17 -y

# Step 2 — pip packages
pip install -r autonomous_drug_discovery/requirements.txt
```

Download P2Rank:
```bash
cd ~
wget https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz
tar -xzf p2rank_2.5.1.tar.gz && rm p2rank_2.5.1.tar.gz
```

### Run the full pipeline

```bash
cd autonomous_drug_discovery

# Production mode: real generation (RDKit), screening, and docking
conda run -n base python orchestrator.py run data/processed/1M17.pdb --mode production

# Pocket2Mol mode: autoregressive pocket-conditioned generation (requires pocket2mol_env)
conda run -n base python orchestrator.py run data/processed/6P3D.pdb --mode pocket2mol

# TargetDiff mode: diffusion-based generation (requires targetdiff_env)
conda run -n base python orchestrator.py run data/processed/6P3D.pdb --mode targetdiff

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

## Documentation

| Doc | Question it answers |
|-----|---------------------|
| This README | What is this repo? Quick start, structure, dependencies |
| [docs/north-star.md](docs/north-star.md) | What are we building and why? Vision, market, strategy, roadmap |
| [docs/testing-guide.md](docs/testing-guide.md) | How do I know the science works? Plain-language explanations, validation experiments, glossary |
| [docs/pipeline-guide.md](docs/pipeline-guide.md) | How do I use the pipeline? Commands, modes, parameters, output structure |
| [docs/installation.md](docs/installation.md) | How do I set this up? Prerequisites, install steps, troubleshooting |
| [docs/telemetry-guide.md](docs/telemetry-guide.md) | How do I query the data? DB schema, SQL queries, Python API |
| [docs/targetdiff-setup.md](docs/targetdiff-setup.md) | How do I set up TargetDiff? Separate env, standalone testing, performance |
| [docs/pocket2mol-setup.md](docs/pocket2mol-setup.md) | How do I set up Pocket2Mol? Separate env, checkpoint download, troubleshooting |
| [autonomous_drug_discovery/plan.md](autonomous_drug_discovery/plan.md) | What is the architectural design? Layer model, agent planner philosophy, future direction |

## Repository structure

```
.
├── autonomous_drug_discovery/       # Main application
│   ├── orchestrator.py              # CLI entrypoint — runs stages individually or as full pipeline
│   ├── agent_planner.py             # LLM-driven adaptive orchestration (optional, requires API key)
│   ├── telemetry.py                 # SQLite telemetry database (runs + molecule_scores tables)
│   │
│   ├── modules/
│   │   ├── 01_ingestion/
│   │   │   └── run_pocket.py        # P2Rank / fpocket — detects binding pockets on a PDB file
│   │   ├── 02_generation/
│   │   │   └── run_generation.py    # Molecule generator — fragment-based (RDKit) or TargetDiff
│   │   ├── 03_screening/
│   │   │   ├── run_screening.py     # Drug-likeness filters — Lipinski, QED, SA, PAINS, ADMET-AI
│   │   │   └── default_scoring_config.json  # Filter thresholds (editable)
│   │   ├── 04_docking/
│   │   │   └── run_docking.py       # AutoDock Vina docking — simulation, triage, or production mode
│   │   └── 05_ranking/              # Reserved for the planned multi-criteria ranker (empty placeholder)
│   │
│   ├── data/
│   │   └── processed/               # PDB files and manifests
│   │
│   ├── envs/                        # Conda environment specs
│   └── tests/
│       ├── test_screening.py
│       └── test_telemetry.py
│
├── docs/                            # Documentation (see table above)
├── reports/                         # Visualizations (TargetDiff molecule comparisons)
└── data/                            # Reference data (CIF dictionaries, utilities)
```

## Pipeline stages in detail

### Stage 1: Pocket detection (`01_ingestion/run_pocket.py`)

**Tool:** P2Rank (default, ML-based) or fpocket (fallback, geometry-based). Selectable via `--backend p2rank|fpocket`.

**Input:** A `.pdb` file containing a protein structure.

**Output:**
- `{stem}_manifest.json` — best pocket location, score, probability, and pre-computed center coordinates.
- `{stem}_p2rank/` for P2Rank — contains `{stem}.pdb_predictions.csv`, `{stem}.pdb_residues.csv`, and one `{stem}_pocket{N}_atm.pdb` per ranked pocket (e.g. `1M17_pocket1_atm.pdb`).
- `{stem}_out/` for fpocket — contains `pockets/pocket{N}_atm.pdb` plus fpocket's per-pocket info files.

**P2Rank advantages over fpocket:** 10-20 percentage point better recall on standard benchmarks. On EGFR (1M17), P2Rank places the pocket 2.7 A from the known drug (vs fpocket's 6.3 A) with 82% residue overlap (vs 53%).

### Stage 2: Molecule generation (`02_generation/run_generation.py`)

**Modes:**
- `simulation` — writes a single stub molecule (benzene) for pipeline testing.
- `rdkit` — fragment-based combinatorial generation using RDKit. Assembles drug-like scaffolds (indole, quinazoline, piperidine, etc.) with functional group substituents and linkers. Uses BRICS decomposition for additional diversity. Molecules are sized to fit the pocket (estimated from pocket radius). 3D conformers are generated and MMFF-optimized.
- `pocket2mol` — wraps the [Pocket2Mol](https://github.com/pengxingang/Pocket2Mol) autoregressive generator. Builds molecules atom-by-atom inside the pocket using a graph neural network. ~11× faster than TargetDiff. Requires the `pocket2mol_env` conda environment and a pretrained checkpoint. See [docs/pocket2mol-setup.md](docs/pocket2mol-setup.md).
- `targetdiff` — wraps the [TargetDiff](https://github.com/guanjq/targetdiff) E(3)-equivariant diffusion model. Denoises from random noise over 1000 steps, conditioned on the pocket shape. Highest-fidelity 3D generation but slowest. Requires the `targetdiff_env` conda environment. See [docs/targetdiff-setup.md](docs/targetdiff-setup.md).

**Input:** `manifest.json` from Stage 1 (specifically the `best_pocket` path).

**Output:** `generated_molecules.sdf` — 100 molecules by default, each with 3D coordinates, SMILES, and a unique `molecule_id`.

**Configuration:** `num_samples` (default 100), `seed` (default 42).

### Stage 3: Screening (`03_screening/run_screening.py`)

**Tool:** MolScore (preferred) or RDKit fallback + ADMET-AI (104-property toxicity/absorption/metabolism prediction)

The screening module uses MolScore as the primary backend for descriptor calculation and PAINS filtering. If MolScore is not installed, it falls back to hand-rolled RDKit filters. ADMET-AI enrichment runs on survivors regardless of backend.

**Input:** An `.sdf` file of candidate molecules.

**Output:**
- `screened_molecules.sdf` — only molecules that passed all filters.
- `screening_report.json` — per-molecule properties, pass/fail status, ADMET annotations, and attrition summary.
- `run_metadata.json` — module execution metadata (backend, counts, status).

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

**Adding a new filter:** Edit `default_scoring_config.json`. Each entry in `filter_thresholds` maps a property name (prefixed `desc_` for descriptors, `filter_` for substructure filters) to `{"max": N}`, `{"min": N}`, or `{"equals": N}`. No code change required.

### Stage 4: Docking (`04_docking/run_docking.py`)

**Tool:** AutoDock Vina (Python API) + Meeko (ligand PDBQT preparation)

**Modes:**
- `simulation` — returns hardcoded dummy scores for pipeline testing.
- `triage` — fast SMILES-based docking via the TDC Oracle (Vina under the hood). Lower setup overhead than production; requires `pytdc`. Falls back to simulation if unavailable.
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

### Working now (M1 + M2 complete)
- Pocket detection (P2Rank, default) — ML-based, validated against crystallography
- Molecule generation (RDKit fragment-based) — pocket-aware sizing
- Screening (MolScore + ADMET-AI) — Lipinski + QED + SA + PAINS + 104 ADMET properties
- Docking (AutoDock Vina) — production-grade, uses P2Rank pocket centers
- Campaign telemetry — full logging of all stages
- Validated against 3 cancer targets (EGFR, BCR-ABL, BRAF V600E) with crystallographic ground truth
- Benchmark comparison script (`benchmark.py`)

### Recently added (M2.6)
- Pocket2Mol autoregressive generation wired into the orchestrator (`--mode pocket2mol`)
- TargetDiff diffusion generation wired into the orchestrator (`--mode targetdiff`)

### Planned next
- Empirical comparison: Pocket2Mol vs TargetDiff vs RDKit on the validated cancer targets
- GNINA CNN-based rescoring
- AiZynthFinder retrosynthetic feasibility
- Multi-criteria ranker (`modules/05_ranking/`, currently an empty placeholder) — combines docking, ADMET, and synthesis-feasibility into a single ranked output
- Domain expert review (M3)
- Agent planner with empirical strategy selection

## Dependencies

| Package | Purpose | Install |
|---------|---------|---------|
| P2Rank | ML-based pocket detection (default) | [Download binary](https://github.com/rdk/p2rank/releases) |
| Java 17+ | Required by P2Rank | `conda install -c conda-forge openjdk=17` |
| RDKit | Molecular property calculation, fragment generation, PAINS filters | `conda install -c conda-forge rdkit` |
| MolScore | Primary screening backend — descriptor calculation and PAINS filtering | `pip install molscore` |
| AutoDock Vina | Molecular docking (binding affinity scoring) | `conda install -c conda-forge vina` |
| Meeko | Ligand PDBQT preparation for Vina | `pip install meeko` |
| gemmi | Receptor PDB parsing and atom typing | `pip install gemmi` |
| ADMET-AI | 104-property ADMET prediction | `pip install admet-ai` |
| fpocket | Pocket detection (fallback) | [Build from source](https://github.com/Discngine/fpocket) |

## License

See repository license file.
