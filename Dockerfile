# syntax=docker/dockerfile:1.7
#
# Agent Harness — Autonomous Drug Discovery pipeline.
#
# One image holds everything a run needs: the conda environments (solved once,
# here, not on every contributor's laptop), the P2Rank binary, and the
# TargetDiff diffusion weights. The same image runs on a laptop and on a
# RunPod GPU pod — the only difference at run time is `--gpus all`.
#
# Phase 1 ships TargetDiff only. Pocket2Mol is intentionally excluded: its
# environment and checkpoint are not part of this build (see plan.md).
#
# Layer order is deliberate — environments and weights sit *above* the code
# COPY, so editing a Python file rebuilds in seconds, not in a conda solve.

FROM nvidia/cuda:11.7.1-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive \
    CONDA_DIR=/opt/conda \
    PATH=/opt/conda/bin:/opt/conda/condabin:$PATH \
    P2RANK_BIN=/opt/p2rank_2.5.1/prank \
    PYTHONUNBUFFERED=1

# --- OS packages -------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        git curl wget ca-certificates bzip2 procps unzip build-essential \
    && rm -rf /var/lib/apt/lists/*

# rclone — used by cloud runs to sync campaign outputs to Cloudflare R2.
RUN curl -fsSL https://rclone.org/install.sh | bash

# --- Miniforge (conda + mamba, conda-forge as the default channel) -----------
RUN curl -fsSL -o /tmp/miniforge.sh \
        https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    && bash /tmp/miniforge.sh -b -p ${CONDA_DIR} \
    && rm /tmp/miniforge.sh \
    && conda config --set always_yes yes \
    && conda clean -afy

# --- Conda environments ------------------------------------------------------
# base           — orchestrator + every CPU stage: ingestion (P2Rank), RDKit
#                  generation, screening, AutoDock Vina docking, ranking.
# targetdiff_env — TargetDiff diffusion generation (PyTorch 1.13 / CUDA 11.7).
#
# No docking_env: the orchestrator runs Stage 4 inside `base` (which carries
# `vina` + `meeko`), so a separate docking environment would never be used.
COPY autonomous_drug_discovery/envs/env_orchestrator.yml /tmp/envs/
COPY autonomous_drug_discovery/envs/env_targetdiff.yml   /tmp/envs/
RUN mamba env update -n base -f /tmp/envs/env_orchestrator.yml \
    && mamba env create     -f /tmp/envs/env_targetdiff.yml \
    && conda clean -afy \
    && rm -rf /tmp/envs

# --- P2Rank (pocket detection, Stage 1) --------------------------------------
# Java is supplied by `openjdk` inside the `base` environment, which is where
# Stage 1 runs — so `prank` finds a JVM on PATH at run time.
RUN curl -fsSL -o /tmp/p2rank.tar.gz \
        https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz \
    && tar -xzf /tmp/p2rank.tar.gz -C /opt \
    && rm /tmp/p2rank.tar.gz \
    && test -x ${P2RANK_BIN}

# --- TargetDiff weights ------------------------------------------------------
# Mirrored on the Hugging Face Hub (the upstream Google Drive folders are dead).
# Fetched here, above the code COPY, so a source edit does not re-download them.
# The weights are excluded from the build context by .dockerignore — this layer
# is their only path into the image.
ARG HF_WEIGHTS=https://huggingface.co/vladsft/agent-harness-weights/resolve/main
ARG TD_WEIGHTS_DIR=/app/autonomous_drug_discovery/modules/02_generation/targetdiff/pretrained_models
RUN mkdir -p ${TD_WEIGHTS_DIR} \
    && curl -fsSL -o ${TD_WEIGHTS_DIR}/pretrained_diffusion.pt  ${HF_WEIGHTS}/pretrained_diffusion.pt \
    && curl -fsSL -o ${TD_WEIGHTS_DIR}/egnn_pdbbind_v2016.pt    ${HF_WEIGHTS}/egnn_pdbbind_v2016.pt

# --- Application code --------------------------------------------------------
COPY . /app
WORKDIR /app

# Re-apply the NumPy-deprecation patches to the pinned TargetDiff submodule
# (idempotent — a no-op if the build context was already patched), then prove
# the checkpoints deserialise inside the environment that will load them.
RUN bash scripts/apply_targetdiff_patches.sh \
    && conda run -n targetdiff_env python scripts/verify_targetdiff_weights.py

RUN install -m 0755 scripts/docker-entrypoint.sh /usr/local/bin/entrypoint.sh

# Runs land here; the data/ volume is mounted at ./data relative to this dir.
WORKDIR /app/autonomous_drug_discovery
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["orchestrator.py", "--help"]
