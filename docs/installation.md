# Installation Guide

## Prerequisites

- Linux or macOS (tested on Ubuntu 22.04+)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- Java 17+ (for P2Rank)
- ~4 GB disk space (dependencies + pretrained models)

## Step 1: Conda Packages

RDKit, Vina, and Java need conda because they ship compiled extensions:

```bash
conda install -n base -c conda-forge rdkit vina openjdk=17 -y
```

## Step 2: Pip Packages

```bash
cd autonomous_drug_discovery
pip install -r requirements.txt
```

This installs:

| Package | Purpose |
|---|---|
| molscore | Primary screening backend (descriptor calculation, PAINS filters) |
| admet-ai | 104-property ADMET prediction (toxicity, absorption, metabolism) |
| meeko | Ligand PDBQT preparation for AutoDock Vina |
| gemmi | Receptor PDB parsing and AutoDock atom typing |
| pytdc | TDC Oracle for triage docking mode |
| numpy, pandas, scipy | Numerical utilities |
| pyyaml | YAML config generation for TargetDiff |

## Step 3: P2Rank (Pocket Detection)

P2Rank is the default pocket detection backend. It's a standalone Java binary:

```bash
cd ~
wget https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz
tar -xzf p2rank_2.5.1.tar.gz && rm p2rank_2.5.1.tar.gz
```

The pipeline expects P2Rank at `~/p2rank_2.5.1/prank`. To use a different location, set the environment variable:

```bash
export P2RANK_BIN=/your/path/to/prank
```

Verify it works:

```bash
~/p2rank_2.5.1/prank predict -f autonomous_drug_discovery/data/processed/1M17.pdb -o /tmp/p2rank_test
```

## Step 4: fpocket (Optional Fallback)

fpocket is a geometry-based pocket detector, used as a fallback when P2Rank is unavailable:

```bash
cd ~
git clone https://github.com/Discngine/fpocket.git
cd fpocket
make
```

The pipeline expects fpocket at `~/fpocket/bin/fpocket`. To use a different location:

```bash
export FPOCKET_BIN=/your/path/to/fpocket
```

## Step 5: TargetDiff (Optional, GPU Recommended)

TargetDiff is an E(3)-equivariant diffusion model for structure-based drug design. It requires a separate conda environment due to PyTorch version constraints.

See [targetdiff-setup.md](targetdiff-setup.md) for full instructions.

This step is optional. The default RDKit fragment-based generator works without TargetDiff.

## Verify Installation

Run the pipeline in simulation mode (no GPU, no external tools needed beyond P2Rank):

```bash
cd autonomous_drug_discovery
conda run -n base python orchestrator.py run data/processed/1M17.pdb --mode simulation
```

Expected output: campaign completes with all 4 stages showing `success`.

Run in production mode to verify the full stack:

```bash
conda run -n base python orchestrator.py run data/processed/1M17.pdb --mode production
```

This runs real molecule generation, screening with ADMET-AI, and full Vina docking. Takes 2-5 minutes depending on hardware.

## Troubleshooting

### `P2Rank binary not found`

Set the environment variable to point at your P2Rank installation:

```bash
export P2RANK_BIN=/path/to/p2rank_2.5.1/prank
```

### `Java not found` or `Unsupported class file major version`

P2Rank requires Java 17+:

```bash
java -version  # should show 17+
conda install -c conda-forge openjdk=17 -y
```

### `ModuleNotFoundError: No module named 'vina'`

Vina must be installed via conda, not pip:

```bash
conda install -c conda-forge vina -y
```

### `ModuleNotFoundError: No module named 'rdkit'`

RDKit must be installed via conda:

```bash
conda install -c conda-forge rdkit -y
```

### ADMET-AI warnings about `num_workers`

This is a PyTorch Lightning info message, not an error. It can be safely ignored. To suppress:

```bash
export PYTHONWARNINGS="ignore::UserWarning"
```

### Screening shows 0% survival rate

Check your `default_scoring_config.json` thresholds. If you modified them too aggressively, all molecules will be filtered. Reset to defaults:

```bash
git checkout modules/03_screening/default_scoring_config.json
```
