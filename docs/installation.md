# Installation Guide

## Prerequisites

- Linux or macOS (tested on Ubuntu 22.04+)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- Java 17+ (for P2Rank, installed via conda below)
- ~8 GB disk space (dependencies + pretrained models, including optional deep-learning envs)

### Installing Miniconda (if you don't have it)

```bash
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -p ~/miniconda3
~/miniconda3/bin/conda init bash
# Restart your shell, or source ~/.bashrc

# Accept Anaconda Terms of Service for default channels
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
```

## Step 1: Conda Packages

RDKit, Vina, and Java need conda because they ship compiled extensions:

```bash
conda install -n base -c conda-forge rdkit vina openjdk=17 -y
```

## Step 2: Pip Packages

```bash
cd autonomous_drug_discovery
conda run -n base pip install -r requirements.txt
```

This installs:

| Package | Purpose |
|---|---|
| molscore | Primary screening backend (descriptor calculation, PAINS filters) |
| admet-ai | 104-property ADMET prediction (toxicity, absorption, metabolism) |
| meeko | Ligand PDBQT preparation for AutoDock Vina |
| gemmi | Receptor PDB parsing and AutoDock atom typing |
| numpy, pandas, scipy | Numerical utilities |
| pyyaml | YAML config generation for TargetDiff / Pocket2Mol |

**Note on `pytdc`:** The TDC Oracle (used by `--mode triage` in docking) pins `scikit-learn==1.2.2`, which does not compile on Python 3.13. It is therefore commented out in `requirements.txt`. Triage docking mode is unavailable on Python 3.13+, but production Vina docking is unaffected. To enable triage mode, install the base env on Python ≤3.11 and `pip install PyTDC` manually.

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

## Step 5: Deep-Learning Generation Backends (Optional)

The pipeline supports two optional deep-learning molecule generators, each in its own conda environment. Both are optional — the default RDKit fragment-based generator works without either.

### Pocket2Mol (fast autoregressive generator, recommended)

Pocket2Mol is ~11× faster than TargetDiff and easier to run on modest hardware. See [pocket2mol-setup.md](pocket2mol-setup.md) for full instructions. Summary:

```bash
conda env create -f autonomous_drug_discovery/envs/env_pocket2mol.yml
git clone https://github.com/pengxingang/Pocket2Mol.git \
    autonomous_drug_discovery/modules/02_generation/pocket2mol
# Download pretrained_Pocket2Mol.pt (44.9 MB) to pocket2mol/ckpt/
```

### TargetDiff (diffusion generator, highest-fidelity 3D)

TargetDiff uses diffusion sampling (1000 denoising steps), producing higher-fidelity 3D geometry but much slower. See [targetdiff-setup.md](targetdiff-setup.md) for full instructions. Summary:

```bash
conda env create -f autonomous_drug_discovery/envs/env_targetdiff.yml
# Clone TargetDiff + download checkpoint — see targetdiff-setup.md
```

## Step 6: LLM Agent Planner (Optional)

`agent_planner.py` adds an LLM-driven adaptive loop on top of the deterministic pipeline. It is optional — without these dependencies (or without `DISCOVERY_LLM_API_KEY`), the planner falls back to a deterministic run.

```bash
conda run -n base pip install langchain langchain-core langchain-google-genai
export DISCOVERY_LLM_API_KEY=your_gemini_api_key
# Optional: export DISCOVERY_LLM_PROVIDER=google   # default; "openai" / "anthropic" are reserved
```

These packages are intentionally excluded from `requirements.txt` so the core pipeline stays slim.

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

### `CondaToSNonInteractiveError: Terms of Service have not been accepted`

Happens on fresh Miniconda installs before any environment is created. Accept the ToS for the default channels:

```bash
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
```

### `conda: command not found` in subprocess / orchestrator

The orchestrator uses `subprocess.check_call` to invoke conda, which does not inherit shell aliases. If conda is not on the system PATH (common when conda is only initialized in `~/.bashrc` interactively), set `CONDA_EXE`:

```bash
export CONDA_EXE=$HOME/miniconda3/bin/conda
```

The orchestrator reads this variable and falls back to the bare `conda` command if unset.
