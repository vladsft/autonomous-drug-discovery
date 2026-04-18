# Pocket2Mol Setup Guide

Pocket2Mol is an autoregressive, pocket-conditioned molecular generator. Unlike diffusion models (which denoise from random noise over many steps), Pocket2Mol builds molecules atom-by-atom inside the binding pocket using a graph neural network. It is roughly **11× faster than TargetDiff** (~7 s/molecule on GPU vs ~78 s for TargetDiff) while producing molecules of comparable quality.

Pocket2Mol is optional. The default RDKit fragment-based generator works without it, and TargetDiff remains available for diffusion-based generation.

## Why a Separate Environment

Pocket2Mol pins PyTorch 1.10.1 + CUDA 11.3 + RDKit 2022.03, which conflict with both the base pipeline (Python 3.10+, latest RDKit) and the TargetDiff environment (PyTorch 1.13 + CUDA 11.7). It runs in its own conda environment (`pocket2mol_env`) and is invoked as a subprocess.

## Step 1: Create the Conda Environment

```bash
cd autonomous_drug_discovery
conda env create -f envs/env_pocket2mol.yml
```

This creates the `pocket2mol_env` environment with:
- Python 3.8
- PyTorch 1.10.1 + CUDA 11.3
- PyTorch Geometric (pyg, pytorch-cluster, pytorch-scatter, pytorch-sparse)
- RDKit 2022.03, BioPython, PyYAML
- lmdb, easydict (via pip)

The environment will run on CPU even without a CUDA 11.3 runtime; it's just much slower.

## Step 2: Clone the Pocket2Mol Repository

```bash
cd autonomous_drug_discovery/modules/02_generation
git clone https://github.com/pengxingang/Pocket2Mol.git pocket2mol
```

The orchestrator expects the repo at `modules/02_generation/pocket2mol/`.

## Step 3: Download the Pretrained Checkpoint

The checkpoint is hosted on Google Drive. Use `gdown`:

```bash
pip install gdown   # (add --break-system-packages on PEP 668 systems)
gdown --folder "https://drive.google.com/drive/folders/1KfdOczjUPITPhIvCuBmnj4xFTV-iI2xB" \
      -O modules/02_generation/pocket2mol/ckpt/
```

`gdown` nests the downloaded folder one level deep. Flatten it:

```bash
mv modules/02_generation/pocket2mol/ckpt/ckpt/pretrained_Pocket2Mol.pt \
   modules/02_generation/pocket2mol/ckpt/pretrained_Pocket2Mol.pt
rmdir modules/02_generation/pocket2mol/ckpt/ckpt
```

Verify the checkpoint is 44.9 MB and at the right path:

```bash
ls -lh modules/02_generation/pocket2mol/ckpt/pretrained_Pocket2Mol.pt
# -rw-r--r-- 1 user user 45M ... pretrained_Pocket2Mol.pt
```

## Step 4: Test Standalone

Before running through the pipeline, verify Pocket2Mol works on its own:

```bash
cd autonomous_drug_discovery/modules/02_generation/pocket2mol

conda run -n pocket2mol_env python sample_for_pdb.py \
  --pdb_path ../../data/processed/6P3D.pdb \
  --center "33.0,30.0,24.0" \
  --outdir /tmp/pocket2mol_test \
  --device cpu
```

(The center coordinates above are placeholders — use real values from a P2Rank manifest.)

Output lands in `/tmp/pocket2mol_test/{timestamp}/SDF/*.sdf`, one `.sdf` file per successfully generated molecule, plus a `SMILES.txt`.

## Step 5: Run Through the Pipeline

```bash
# Single stage
conda run -n base python orchestrator.py generate data/processed/6P3D_manifest.json --mode pocket2mol

# Full pipeline
conda run -n base python orchestrator.py run data/processed/6P3D.pdb --mode pocket2mol
```

The orchestrator automatically:
1. Reads the pocket center from the P2Rank manifest (or computes the centroid from a pocket PDB if using fpocket).
2. Writes a per-run `pocket2mol_run_config.yml` with the absolute checkpoint path and requested sample count.
3. Invokes `sample_for_pdb.py` in the `pocket2mol_env` environment.
4. Collects all SDFs from `{outdir}/{timestamp}/SDF/*.sdf` and consolidates them into a single `generated_molecules.sdf`.
5. Assigns pipeline-compatible molecule IDs, SMILES properties, and a `generator=pocket2mol` tag.

## Configuration

The repo ships a default sampling config at `pocket2mol/configs/sample_for_pdb.yml`:

```yaml
model:
  checkpoint: ./ckpt/pretrained_Pocket2Mol.pt

sample:
  seed: 2020
  num_samples: 100
  beam_size: 300
  max_steps: 50
  threshold:
    focal_threshold: 0.5
    pos_threshold: 0.25
    element_threshold: 0.3
    hasatom_threshold: 0.6
    bond_threshold: 0.4
```

When run through the orchestrator, a per-run copy is written to the campaign output directory with an absolute checkpoint path and overridden `num_samples`. The bounding box size is fixed at 23 Å (the value Pocket2Mol was trained on); adjust via `--bbox_size` if your pocket is unusually large or small.

## Performance

| Setting | Time per Molecule | Total (100 molecules) |
|---|---|---|
| GPU (NVIDIA V100/A100) | ~7 sec | ~12 min |
| CPU (multi-core) | ~1–3 min | 2–5 hours |

Compared to TargetDiff:

| Generator | GPU time / mol | Generation paradigm | 3D pocket-aware |
|---|---|---|---|
| Pocket2Mol | ~7 s | Autoregressive | Yes |
| TargetDiff | ~78 s | Diffusion (1000 steps) | Yes |
| RDKit | ~10 ms | Fragment combinatorial | No (size-only) |

## Comparison with TargetDiff

Both are structure-based generative models conditioned on a 3D pocket, but they differ in:

| Feature | Pocket2Mol | TargetDiff |
|---|---|---|
| Paradigm | Autoregressive (atom-by-atom) | Diffusion (denoise from noise) |
| Speed | ~11× faster | Slower (1000 denoising steps) |
| Input | Full protein PDB + center coords | Pocket-only PDB |
| Use case | Breadth / fast exploration | Depth / refinement |

The pipeline's `plan.md` recommends Pocket2Mol for broad exploration when entering new chemical space, and TargetDiff for refinement around a known scaffold region.

## Troubleshooting

### `conda: command not found` inside subprocess

The orchestrator respects the `CONDA_EXE` environment variable and falls back to `conda`. If conda is not on the PATH in subprocess contexts (common on non-interactive shells), either export it:

```bash
export PATH="$HOME/miniconda3/bin:$PATH"
```

Or set `CONDA_EXE` explicitly:

```bash
export CONDA_EXE=$HOME/miniconda3/bin/conda
```

### `Pocket2Mol produced no SDF files`

Known "initialization failed" behavior on some targets (Pocket2Mol GitHub issue #38). Try:
1. Confirming the pocket center coordinates are reasonable (inside the protein).
2. Using a different pocket from the P2Rank predictions CSV.
3. Increasing `--bbox_size` if the pocket is large.
4. Re-running — initialization has a random component.

### `Initialization failed` in logs

Pocket2Mol samples initial atoms from the pocket; if none of the beam candidates survive early filtering, nothing is written. This is different from the TargetDiff "50% reconstruction failure" — it's a beam search collapse. Same mitigations as above.

### Environment creation fails

If `conda env create` times out or fails:

```bash
conda create -n pocket2mol_env python=3.8 -y
conda activate pocket2mol_env
conda install pytorch=1.10.1 cudatoolkit=11.3 -c pytorch -c nvidia -y
conda install pyg pytorch-cluster pytorch-scatter pytorch-sparse -c pyg -y
conda install rdkit=2022.03 biopython pyyaml -c conda-forge -y
pip install lmdb easydict
```

### CUDA version mismatch

Pocket2Mol officially targets CUDA 11.3. If your system has a newer CUDA driver (e.g. 11.7 or 12.x), PyTorch's bundled CUDA runtime should still work — but if you see `CUDA error: no kernel image is available`, fall back to `--device cpu` or rebuild PyTorch against your system CUDA.

## Sources

- Repo: <https://github.com/pengxingang/Pocket2Mol>
- Paper: Peng et al., "Pocket2Mol: Efficient Molecular Sampling Based on 3D Protein Pockets" (ICML 2022)
- Checkpoint folder: <https://drive.google.com/drive/folders/1KfdOczjUPITPhIvCuBmnj4xFTV-iI2xB>
