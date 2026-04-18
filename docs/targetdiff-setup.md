# TargetDiff Setup Guide

TargetDiff is an E(3)-equivariant diffusion model for structure-based drug design. It generates molecules conditioned on the 3D shape of a protein binding pocket, producing candidates that are geometrically complementary to the target.

TargetDiff is optional. The default RDKit fragment-based generator works without it.

## Why a Separate Environment

TargetDiff requires Python 3.8 and PyTorch 1.13, which conflict with the main pipeline's dependencies. It runs in its own conda environment (`targetdiff_env`) and is invoked as a subprocess.

## Step 1: Create the Conda Environment

```bash
cd autonomous_drug_discovery
conda env create -f envs/env_targetdiff.yml
```

This creates the `targetdiff_env` environment with:
- Python 3.8
- PyTorch 1.13 + CUDA 11.7
- PyTorch Geometric (pyg, pytorch-cluster, pytorch-scatter, pytorch-sparse)
- OpenBabel, RDKit, BioPython

## Step 2: Verify the Checkpoint

The pretrained diffusion model checkpoint should already be present at:

```
modules/02_generation/targetdiff/pretrained_models/pretrained_diffusion.pt
```

If missing, download it from the [TargetDiff repository](https://github.com/guanjq/targetdiff) releases.

A second checkpoint for property prediction is also included:

```
modules/02_generation/targetdiff/pretrained_models/egnn_pdbbind_v2016.pt
```

## Step 3: Test Standalone

Verify TargetDiff works independently before using it through the pipeline:

```bash
conda run -n targetdiff_env python modules/02_generation/targetdiff/scripts/sample_for_pocket.py \
  modules/02_generation/targetdiff/configs/sampling.yml \
  --pdb_path data/processed/6P3D_p2rank/6P3D_pocket1_atm.pdb \
  --result_path /tmp/targetdiff_test \
  --device cpu \
  --num_samples 2
```

This should produce:
- `/tmp/targetdiff_test/sample.pt` — raw tensor output
- `/tmp/targetdiff_test/sdf/` — individual SDF files for each generated molecule

On CPU, expect ~12 minutes per molecule (1000 denoising steps). On a GPU, this drops to ~1-2 minutes.

## Step 4: Run Through the Pipeline

Once standalone works, use it through the orchestrator:

```bash
# Single stage
conda run -n base python orchestrator.py generate data/processed/6P3D_manifest.json --mode targetdiff

# Full pipeline
conda run -n base python orchestrator.py run data/processed/6P3D.pdb --mode targetdiff
```

The orchestrator automatically:
1. Generates a `sampling_config.yml` with the correct checkpoint path
2. Invokes `sample_for_pocket.py` in the `targetdiff_env` conda environment
3. Consolidates individual SDF outputs into a single `generated_molecules.sdf`
4. Assigns pipeline-compatible molecule IDs and SMILES properties

## Configuration

The default sampling config (`modules/02_generation/targetdiff/configs/sampling.yml`):

```yaml
model:
  checkpoint: ./pretrained_models/pretrained_diffusion.pt

sample:
  seed: 2021
  num_samples: 100
  num_steps: 1000       # denoising steps (more = higher quality, slower)
  pos_only: False
  center_pos_mode: protein
  sample_num_atoms: prior
```

When run through the pipeline, a temporary config is generated that overrides `num_samples` and uses the absolute checkpoint path.

## Performance

| Setting | Time per Molecule | Total (100 molecules) |
|---|---|---|
| CPU (single core) | ~12 min | ~20 hours |
| GPU (NVIDIA V100) | ~1-2 min | ~2-3 hours |
| GPU (NVIDIA A100) | ~30 sec | ~50 min |

For initial testing, generate 2-5 molecules on CPU. For production campaigns, a GPU is strongly recommended.

## Proof-of-Concept Results

TargetDiff was validated standalone on the BRAF V600E pocket (6P3D), generating molecules from random noise:

| Molecule | MW | QED | Dock Score | Ligand Efficiency | Tanimoto to Ponatinib |
|---|---|---|---|---|---|
| `C1=CN=CC=C(c2cccc(Nc3ccncc3)c2)C1` | 261 | 0.899 | -7.59 kcal/mol | 0.38 | 0.186 |
| `COc1cnc(C(=O)NCc2cccc(C)n2)cn1` | 258 | 0.889 | -7.38 kcal/mol | 0.39 | 0.154 |

Both molecules are drug-like, compact, and structurally novel (low Tanimoto similarity to known drugs). Docking scores are moderate but within the "worth investigating" range. Ligand efficiency (0.38-0.39) is excellent.

## Troubleshooting

### `conda: command not found` inside subprocess

The pipeline invokes TargetDiff via `conda run -n targetdiff_env python ...`. If conda is not on the PATH in subprocess contexts, add it:

```bash
export PATH="$HOME/miniconda3/bin:$PATH"
```

### `CUDA out of memory`

Reduce `batch_size` in the generation parameters. The orchestrator passes `--batch_size 16` by default. Try 4 or 8:

Edit `DEFAULT_PARAMS` in `modules/02_generation/run_generation.py`:

```python
DEFAULT_PARAMS = {
    "batch_size": 8,   # reduce if GPU memory is limited
    ...
}
```

### `TargetDiff produced no SDF files`

This usually means molecule reconstruction failed for all samples. Check:
1. The pocket PDB actually contains protein atoms (not empty)
2. The pocket is reasonable size (P2Rank probability > 0.5)
3. Try increasing `num_samples` — TargetDiff has a ~50-70% reconstruction success rate

### Environment conflicts

If the `targetdiff_env` fails to create, try installing packages incrementally:

```bash
conda create -n targetdiff_env python=3.8 -y
conda activate targetdiff_env
conda install pytorch=1.13.0 pytorch-cuda=11.7 -c pytorch -c nvidia -y
conda install pyg pytorch-cluster pytorch-scatter pytorch-sparse -c pyg -y
conda install rdkit openbabel biopython -c conda-forge -y
pip install lmdb easydict
```
