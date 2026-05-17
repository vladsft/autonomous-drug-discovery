# Adaptive Discovery Orchestrator — Architectural Plan

This is the canonical plan. The [README](../README.md) and other docs under [`docs/`](../docs/) reference this file; if they disagree, this file wins. Last major revision: 2026-05-11.

## Where we are

**M1 (Working Pipeline) and M2 (Validation) are complete.** The 5-stage deterministic pipeline (P2Rank → generation → screening + ADMET → Vina docking → multi-criteria ranking) runs end-to-end on real proteins, with full SQLite telemetry. It is validated against three cancer targets with crystallographic ground truth:

| Target | Disease | Best Dock Score | Pocket Distance | Residue Overlap |
|---|---|---|---|---|
| EGFR (1M17) | Lung cancer | -9.32 kcal/mol | 2.7 Å | 82% |
| BCR-ABL (2HYY) | CML | -12.59 kcal/mol | 2.7 Å | 92% |
| BRAF V600E (6P3D) | Melanoma | -11.20 kcal/mol | 3.1 Å | 89% |

**Current friction (the reason this plan exists).** The pipeline runs on the dev machine, but is not portable:
- TargetDiff and Pocket2Mol checkpoints live on **Google Drive folders that have been deleted** (confirmed 2026-05-11). The dev box still has `pretrained_diffusion.pt` (33 MB, Mar 2023) from before the takedown. Pocket2Mol's `pretrained_Pocket2Mol.pt` was never recovered locally.
- Four conda envs (`base`, `targetdiff_env`, `pocket2mol_env`, `docking_env`) are recipes-only. A new machine spends 15-30 minutes resolving each, and Pocket2Mol's CUDA 11.3 pin doesn't run on Blackwell.
- Generated `data/campaign_*/` directories are git-ignored and trapped on whichever machine produced them.
- The dev box has no GPU; lid-close suspends every CPU job; one TargetDiff molecule is ~12 CPU-minutes (vs ~30 s on an A100).
- A 2026-05-11 attempt to run TargetDiff on 2HYY also surfaced a `PYTHONPATH` bug in the wrapper that meant TargetDiff had never worked end-to-end through the orchestrator from a fresh shell (now fixed in `run_generation.py`).

**Demo forcing function.** The professor presentation is the next milestone. The deliverable is **20 quality candidate molecules** rendered through the dashboard, end-to-end through the pipeline, on five targets. Everything in Phase 1 below is scoped to that deliverable; nothing extra.

## Architecture (target state)

```
                      ┌────────────────────┐
                      │  GitHub (code +    │
                      │  Actions CI build) │
                      └─────────┬──────────┘
                                │
                  ┌─────────────┴──────────────┐
                  │                            │
                  ▼                            ▼
        ┌───────────────────┐         ┌──────────────────┐
        │ GHCR              │         │ Hugging Face Hub │
        │ ghcr.io/.../      │         │ <org>/td-weights │
        │ agent-harness:    │         │ <org>/p2m-weights│
        │ <git-sha>         │         │ (108 MB total)   │
        └─────────┬─────────┘         └────────┬─────────┘
                  │ docker pull                │ baked into image at build time
       ┌──────────┴───────────┐                │
       │                      │                │
       ▼                      ▼                │
┌─────────────┐     ┌──────────────────┐       │
│ Laptop A    │     │ RunPod GPU pod   │       │
│ (you)       │     │ (3090 / A40,     │       │
│ Docker CLI  │     │  ephemeral)      │       │
└──────┬──────┘     └────────┬─────────┘       │
       │                     │                 │
       └─────────┬───────────┘                 │
                 │                             │
                 ▼                             │
        ┌────────────────────┐                 │
        │ Cloudflare R2      │  ←──────────────┘
        │ (campaigns/,       │     (image pulls weights once at build,
        │  telemetry.db,     │      so cloud runs are immediately usable)
        │  large artefacts)  │
        └────────┬───────────┘
                 │ static export
                 ▼
        ┌────────────────────┐
        │ GitHub Pages       │     ← professor opens this URL
        │ (dashboard/*)      │
        └────────────────────┘
```

### Operating principles (load-bearing)

1. **Reproducibility is an artefact, not a recipe.** Don't ship "run `conda env create -f x.yml` and hope." Ship a binary image that already has the env solved. Pay the cost once in CI, not on every contributor's laptop.
2. **Where the work runs and where you sit are different questions.** Editing happens on the laptop. Heavy generation runs on cloud GPU. State lives in a third place. Coupling them was the bug.
3. **Container parity.** Local laptop == cloud GPU == CI runner. Same Docker image, same entrypoint, only difference is `--gpus all`. Eliminates an entire class of "works on mine, not yours" failures.
4. **Make the lid irrelevant.** If closing a laptop can lose work, the architecture is wrong, not the laptop. Long jobs go off the laptop.
5. **Sequence by demo-criticality.** One week to professor day. Every step delivers a working improvement; if Day 4 slips, Days 1-3 still ship a usable system.
6. **Containers, not abstractions, are the line we hold.** No Kubernetes. No Modal refactor yet. No Nix. Conda inside Docker. This is enough.

### Component decisions

#### Docker, not bootstrap.sh, is the install path
**Decision:** Multi-stage `Dockerfile` from `nvidia/cuda:11.7.1-runtime-ubuntu22.04`. Three conda envs (`orchestrator`, `targetdiff`, `docking`) installed at build time. P2Rank tarred in. Model weights pulled from HF Hub during build via `huggingface-cli download`. Final image ~7-8 GB, pushed to GHCR by GitHub Actions on every `main` commit.

**Why:** Solves env reproducibility, weight distribution, P2Rank placement, and CUDA pinning at once. New contributor flow is `docker pull && docker run --gpus all -v $(pwd)/data:/app/data ghcr.io/<you>/agent-harness orchestrator.py run ...`.

**Rejected:** *bootstrap.sh + conda envs from yml.* Still costs 30+ min per machine, still drifts. *Nix / pixi.* Better in theory; worse in practice for a 1-week sprint with this ML stack. *Pure pip / uv.* Conda is how RDKit + Vina + PyG + the diffusion envs ship. Fighting that is gold-plating.

**Concession:** `targetdiff_env` is pinned to CUDA 11.7 — fine on any Ampere/Turing GPU. When Blackwell becomes the demo machine, a second image variant (`:cu121`) gets added; that's a Phase 3+ chore, not now.

#### Cloud GPU on RunPod, wrapped in a script
**Decision:** `scripts/cloud_run.sh` calls the RunPod REST API (or `runpodctl`) to:
1. Spin up an RTX 3090 / A40 pod using the published GHCR image as its container
2. Mount a persistent network volume that holds `data/`
3. Run `orchestrator.py run <target> --mode targetdiff` inside it
4. Stream logs back to your terminal
5. On exit, sync `data/campaign_*/` to R2
6. Terminate the pod

**Why:** Removes manual pod-launching tedium. Removes "did I forget to shut down the GPU?" cost overruns. Removes the lid-close failure mode entirely (the pod is independent of your laptop).

**Rejected:** *Manual web UI workflow.* What we have today. Fine for one-off; bad for repeating 5 times. *Modal.* Strictly better long-term answer for GPU functions, but requires turning `subprocess.check_call(conda run …)` into `@modal.function`. ~1 day refactor. Deferred to Phase 3+. *AWS Batch / Vertex AI.* Enterprise overkill at two-laptop scale.

#### State lives in object storage
**Decision:** Cloudflare R2 bucket holds `data/campaign_*/`, `data/telemetry.db`, and large artefacts. `make sync-up` / `make sync-down` wrap `rclone bisync`. The Docker image ships pre-configured with `rclone`; credentials come from a mounted `.env`.

**Why:** R2 is S3-compatible, ~$0.015/GB-month, zero egress. For ~5 GB of campaigns, ~$0.08/month. Single source of truth across laptops + cloud GPU.

**Rejected:** *Tailscale + rsync between laptops.* Free and simpler, but requires both machines online simultaneously. R2 lets either of us pick up where the other left off without coordination. *Git LFS.* Wrong tool for binary, large, not-version-relevant artefacts. *SQLite as the only shared state.* R2 doesn't mount as a filesystem, so the DB lives locally; at end of each run we sync the file. Acceptable at 1-2 contributors.

#### Dashboard via GitHub Pages, auto-deployed
**Decision:** GitHub Action on `main` commits runs `scripts/regenerate_dashboard.py` (reads telemetry + ranking outputs from R2, rewrites `dashboard/professor_demo.js`), then pushes `dashboard/` to the `gh-pages` branch. Pages serves it.

**Why:** Professor gets a URL — `https://<you>.github.io/agent-harness/` — that always reflects latest state. No "send me the latest file." Free.

**Rejected:** *Cloudflare Pages.* Equally good; chose GitHub Pages because everything else is already on GitHub. *A live FastAPI server on a small VPS.* Future work; for the demo, static is correct.

#### Makefile is the user interface
**Decision:** Top-level `Makefile` with these targets, and nothing else:

```
make bootstrap        # one-time: install Docker, configure rclone, pull image
make pull             # pull latest data from R2
make push             # push local data to R2
make run TARGET=2HYY MODE=targetdiff NUM=30   # local docker run
make cloud-run TARGET=2HYY MODE=targetdiff NUM=30   # RunPod provision + run
make dashboard        # regenerate dashboard JSON locally
make deploy           # push dashboard to GitHub Pages
make logs CAMPAIGN=campaign_xxx     # tail logs from R2
make clean            # nuke local cache
```

**Why:** A Makefile is a contract. Tab-completable, language-agnostic, low-magic, every contributor immediately knows what's possible.

**Rejected:** *Justfile.* Nicer but requires installing `just`. Make is everywhere. *A Python CLI wrapper.* Adds a maintenance surface; Make + bash is sufficient.

#### CI is mandatory infrastructure, not polish
**Decision:** GitHub Actions on every `main` push:
1. Lints (ruff)
2. Runs the test suite inside the Docker image
3. Builds the image with cache layers
4. Pushes to GHCR with `:latest` and `:<git-sha>` tags
5. Deploys the dashboard if `dashboard/` or `data/<latest-campaigns>` changed

**Why:** Without this, the Docker image drifts from the code. With it, "what's running in cloud GPU" = "what's tagged `:latest`" = "what's on `main`."

## Roadmap (4 phases)

### Phase 1 — Containerised cloud-GPU run on 5 targets, 20 quality results

**Goal:** The architecture above shipped end-to-end. Five targets (1M17, 2HYY, 6P3D, KRAS G12C, JAK2) × ~30 TargetDiff molecules each, full pipeline through screening + docking + ranking + dashboard. Output is the 20-quality-result corpus the professor sees.

**Success criteria:**
- `docker pull` + `make run TARGET=… MODE=simulation` works on both contributors' laptops with no manual intervention
- `make cloud-run` provisions a RunPod GPU, runs TargetDiff, syncs results to R2, tears down the pod, all from a single command
- Five targets × ~30 TargetDiff molecules generated, screened, docked, ranked, logged to telemetry
- ≥ 20 candidates labelled "quality" (passed tightened screening + docking better than median of generation set; ideally within 1 kcal/mol of known-drug baseline on at least three targets)
- Dashboard loads at `https://<you>.github.io/agent-harness/`; auto-deploys on every commit; no manual JSON editing
- Total cloud spend < $5

**Why this is Phase 1:** The professor presentation is the next milestone, and everything past Phase 1 is built on having a credible, reproducible result set in hand. Anything else (Pocket2Mol, Bayesian recommender, Sonnet agent) is downstream of "do we have 20 results we believe in?"

**Pocket2Mol is explicitly deferred.** Its checkpoint isn't recovered, its env doesn't run on Blackwell, and TargetDiff alone covers the deliverable. We revisit Pocket2Mol post-demo only if TargetDiff results disappoint.

### Phase 2 — Show the professor; capture qualitative feedback

**Goal:** Walk a medicinal chemistry advisor through the dashboard for the Phase 1 targets. Capture which ranked candidates they think are sensible, which look like generative junk, what they'd modify, what's missing. The conversation includes the forward-looking architecture (Layer 4 below) so we shape Phase 3 with their input.

**Success criteria:**
- 1-2 hour review session with at least one medicinal chemistry expert (Imperial Chemistry / Cancer Research / etc.)
- Per-candidate annotations captured: `promising | borderline | reject` plus free-text reason
- Tightened screening thresholds re-calibrated against expert judgment (Phase 1 used a first-pass guess)
- Decision: which targets are worth running the agent loop on in Phase 3

**Why this is Phase 2:** M3 (expert review) is the gate north-star calls non-negotiable. Without expert validation that current output is *interpretable*, no amount of agent intelligence makes it useful.

### Phase 3 — Adaptive Layer (Sonnet + Bayesian + Obsidian)

**Goal:** Wire the four-part adaptive layer (Layer 4 below) into the orchestrator. Three sub-deliverables, shippable independently:

- **3a. Obsidian campaign emitter** — every `orchestrator.py run` writes `campaigns/{date}-{target}-{id}/{summary.md, top10/{mol_xxx.md}, route_trees/...}` with `[[wikilinks]]` and embedded structure SVGs. Independent of agent; ships first because it's the highest-leverage / lowest-risk piece.
- **3b. Bayesian strategy recommender** — Thompson sampling over (backend, num_samples, beam_size, screening thresholds) conditional on pocket descriptors. Cold-started from Phase 1 telemetry; updates after each new campaign. Outputs a recommendation + confidence band; chemist can override.
- **3c. Sonnet agent shell** — wraps `agent_planner.py` with retrieval over telemetry + open-source corpora (CrossDocked2020, BindingDB, ChEMBL), structured recommendations citing evidence, an attrition-funnel watcher for mid-campaign sanity checks, and a chemist-facing report writer. Every recommendation is logged, never auto-applied.

**Success criteria:**
- All three sub-deliverables shipped end-to-end on at least one validated target
- Agent recommendations logged with rationale; baseline campaigns also retained for comparison
- Phase 2 advisor can answer "did the agent's recommendations match yours, and where did it diverge?"

### Phase 4 — Bayesian evaluation; investor / university narrative

**Goal:** Quantify whether the Phase 3 agent loop produces measurably better candidates than the Phase 1 deterministic baseline. Report posteriors with credible intervals — never point estimates — on every lift metric: composite-score uplift, recall against expert-`promising` annotations, calibration of the recommender's posterior. Honest, even if the lift is null.

**Success criteria:**
- Quantitative comparison on a held-out set of targets with Bayesian credible intervals
- Calibration plot for the recommender: are 90% credible intervals actually 90% empirically?
- Honest write-up of where the agent helps, where it doesn't, and where the system still needs the chemist

**Why Bayesian here:** deep-learning chemistry is a graveyard of point-estimate hype. Reporting a single number ("agent improves Vina by 0.4 kcal/mol!") invites the BenevolentAI critique. Posteriors signal the kind of scientific rigour that survives peer review.

## Layer 4 — Adaptive Orchestration (architecture for Phase 3+)

Only after Phases 1-2 are solid does adaptive orchestration become meaningful. This layer has four cooperating sub-components, each independently shippable and each contributing to the *compounding* differentiation (no single tool here is novel; the value is the loop).

**4a. Cascaded generation pipeline (no agent required).** Pocket2Mol for breadth (~7 s/mol on GPU) → triage by docking + ADMET → TargetDiff to refine the top 10-20 (denoise around them as 3D seeds) → AiZynthFinder synthesis check → multi-criteria rank. The cascade pattern captures the "best of both" without any agent: speed where it matters, fidelity where it pays off. Gated on Pocket2Mol weight recovery.

**4b. Bayesian strategy selection.** For each new target, the question "which generator + parameters?" is a multi-armed bandit. Maintain a posterior over each backend's expected composite-score distribution conditional on pocket descriptors. Use Thompson sampling: occasionally explore weaker arms, mostly exploit the leading one, posterior tightens with each campaign. The *honest* version of "the system learns" — no fine-tuning, no RL on noisy proxies; principled exploration/exploitation over a small set of well-defined choices. Multi-criteria rewards blunt the ICLR 2025 "Vina-bigger-is-better" reward-hacking trap.

**4c. Sonnet-in-the-loop (LLM agent).** A Claude Sonnet agent sits next to the orchestrator with three concrete jobs:
- *Retrieval over campaigns.* Surface similar past campaigns from telemetry + open-source corpora; recommend a starting strategy with cited evidence.
- *Mid-campaign sanity checks.* Watch attrition funnels and ADMET distributions; flag pathological runs ("95% Lipinski violations", "all mols < 200 Da") and suggest a parameter adjustment with rationale.
- *Reporting.* Translate the multi-criteria ranked output into a chemist-readable campaign report.

Crucially: the agent never decides autonomously. It surfaces structured recommendations the chemist can accept, modify, or reject. Every recommendation is logged so we can later measure the agent's hit rate.

**4d. Obsidian knowledge graph (the moat).** Each campaign emits a folder of Markdown notes with embedded SMILES, scorecards, 2D structure SVGs, and `[[wikilinks]]` to other campaigns sharing scaffolds, pockets, or expert annotations. Tags surface common patterns (`#kinase`, `#hERG-flag`, `#3-step-synth`). Over a year, this becomes a queryable knowledge graph of *what worked, what didn't, and why* — sitting in the chemist's existing tool, not behind a custom UI.

**What this layer is NOT.** Not a generative-model fine-tuning loop (needs wet-lab data we don't have). Not a fully autonomous agent (overpromising and unsafe). Not a replacement for the chemist's judgment.

## Supporting technical work (referenced from the phases above)

The numbered items below are lower-level technical tasks that the phases pull from. Many can be done in parallel within a phase.

### Tighten screening thresholds (Phase 1)
Current survival rates are 73-98% (target: 40-60%). The generator is over-producing junk and/or filters are too loose. Tighten conservatively as a config-only change to `modules/03_screening/default_scoring_config.json`. Phase 1 picks a first-pass guess (Lipinski-strict, QED ≥ 0.5, SA ≤ 4.5, MW ≤ 450); Phase 2 re-calibrates with expert input.

### `--num_samples` CLI flag (Phase 1)
`run_generation.py` currently only respects the value baked into `_DEFAULT_PARAMS_BY_MODE`. Add a `--num_samples` argument so campaigns can be sized from the Makefile / `make` target without editing source. Single-line change.

### Pocket2Mol checkpoint recovery (Phase 3 prep)
The checkpoint is the one piece of state we cannot self-mirror — we never had a local copy. Path: open a polite issue on `pengxingang/Pocket2Mol` asking for a re-host; in parallel, ask anyone in our academic network who's run Pocket2Mol recently. Until recovered, Pocket2Mol stays out of the build.

### Pocket2Mol + TargetDiff env rebuild for Blackwell
Both upstream envs (Pocket2Mol: pytorch 1.10 / cu113, TargetDiff: pytorch 1.13 / cu117) target pre-Blackwell architectures and cannot execute on an RTX 50-series GPU (sm_120) — CUDA 11.x ships no kernels for it.

**TargetDiff — done.** `envs/env_targetdiff_blackwell.yml` rebuilds `targetdiff_env` on PyTorch 2.8 / CUDA 12.8 with prebuilt PyG wheels (`torch-scatter/cluster/sparse`) from <https://data.pyg.org/whl/> — no source compile needed. The submodule's PyTorch 2.x break (`torch.load` now defaults to `weights_only=True`) is handled by `targetdiff_patches/02-torch2-compat.patch`, and `run_generation.py` / `orchestrator.py` now take a `--device` flag that auto-detects the GPU. This is the **local-GPU** counterpart to the cloud path: the Docker image and RunPod keep the cu117 env (it runs fine on RunPod's Ampere GPUs); the cu128 env is the no-Docker path for a local Blackwell box. See `docs/pipeline-guide.md` → "Full pipeline, local GPU".

**Pocket2Mol — still pending.** Its cu113 env has not been rebuilt; deferred with the rest of the Pocket2Mol work. Diagnostic that a rebuilt env works: `python -c "import torch; x=torch.randn(2000,2000).cuda(); [x@x.T for _ in range(10)]; torch.cuda.synchronize(); print('ok')"` should return in <1 s.

### GNINA rescoring (post-Phase 2)
Add GNINA CNN-based rescoring alongside Vina. Needs GPU; download binary from <https://github.com/gnina/gnina/releases>. Only worthwhile after screening is properly calibrated and we have evidence Vina's r=0.4-0.6 correlation is the limiting factor.

### AiZynthFinder retrosynthetic feasibility on all candidates (post-Phase 2)
Currently a stubbed slot in the Stage 5 ranker. Apply to top-ranked candidates only — slow and heavyweight. Provides route-level synthesis feasibility, much more informative than SA score alone.

### Agent Planner (Phase 3)
By this point we have real telemetry data, validated benchmarks, and expert-informed parameter ranges. Now we can build the Sonnet agent shell with actual knowledge to encode, not guesses.

## Immediate Action Plan (week of 2026-05-11)

Strict sequence. Each day's output is the next day's input. If anything cracks, stop and stabilise before continuing.

### Day 1 — Foundations (4-5 hours)
- [x] **Pin TargetDiff as a proper git submodule.** ✅ Done 2026-05-11 (commit on `main`). Submodule at `guanjq/targetdiff` SHA `142f1eb…`. Local NumPy-deprecation patches stored under `autonomous_drug_discovery/modules/02_generation/targetdiff_patches/` and re-applied by `scripts/apply_targetdiff_patches.sh` (idempotent). Fresh clones need `git submodule update --init --recursive && scripts/apply_targetdiff_patches.sh`.
- [x] **HF Hub weights mirror.** ✅ Done 2026-05-11. Public model repo: <https://huggingface.co/vladsft/agent-harness-weights>. Holds `pretrained_diffusion.pt`, `egnn_pdbbind_v2016.pt`, and `data/telemetry.db` (M1/M2 campaign history). Fetch with `hf download vladsft/agent-harness-weights <file> --local-dir <dest>`. This is the canonical weights source now that the upstream Google Drive folders are dead — the Dockerfile and any fresh setup pull from here.
- [ ] **Skeleton `Dockerfile`** at repo root. Multi-stage from `nvidia/cuda:11.7.1-runtime-ubuntu22.04`. Miniconda, three envs (`env_orchestrator.yml`, `env_targetdiff.yml`, `env_docking.yml`), P2Rank tarball, weights from `huggingface-cli download`, then `scripts/apply_targetdiff_patches.sh` at the end of the build. `COPY . /app`, entrypoint. Local `docker build -t agent-harness:dev .` succeeds. ~3 hours.
- [ ] **Smoke test:** `docker run agent-harness:dev orchestrator.py run data/processed/1M17.pdb --mode simulation` completes successfully. ~30 min.

### Day 2 — CI + GHCR + R2 (5-6 hours)
- [ ] **GitHub Actions workflow** at `.github/workflows/build.yml`. Build on push to `main`, push to GHCR with `:latest` and `:<sha>`. Layer caching via the `actions/cache` + `docker/build-push-action`. ~2 hours.
- [ ] **Make GHCR image public.** Repository → Packages → settings. ~5 min.
- [ ] **Top-level `Makefile`** with `bootstrap`, `pull`, `push`, `run`, `clean` targets. ~1 hour.
- [ ] **Cloudflare R2 bucket** + API token + `rclone.conf` template + `docs/development.md` section on credentials. ~1 hour.
- [ ] **End-to-end on laptop:** `make pull && make run TARGET=1M17 MODE=rdkit` produces results locally; `make push` puts them in R2; `make pull` from the other laptop retrieves them. ~1 hour.

### Day 3 — Cloud GPU wiring (5-6 hours)
- [ ] **RunPod account + persistent network volume** (~50 GB) configured to mount at `/app/data`. ~30 min.
- [ ] **`scripts/cloud_run.sh`.** Provisions a pod (`runpodctl` or REST), runs the container, syncs to R2 on exit, terminates. ~3 hours including API debugging.
- [ ] **`make cloud-run TARGET=2HYY MODE=targetdiff NUM=30`** end-to-end on one target. Verify output in R2 and via `make pull` locally. ~2 hours.

### Day 4 — Demo corpus generation (~3 hours active + ~3 hours waiting)
- [ ] **Tighten screening thresholds** in `default_scoring_config.json` to ~50% survival target. ~30 min.
- [ ] **Add `--num_samples` CLI flag** to `run_generation.py`. ~15 min.
- [ ] **`make cloud-run` × 5 targets** (1M17, 2HYY, 6P3D, KRAS G12C, JAK2 — fetch PDBs for the latter two). NUM=30 each. ~3 hours wallclock batched. Total spend ~$2.
- [ ] **Spot-check the SDFs** in PyMOL or via `Chem.MolToImage`. ~1 hour.

### Day 5 — Dashboard polish + deploy (5-6 hours)
- [ ] **`scripts/regenerate_dashboard.py`.** Reads `telemetry.db` + ranking outputs from R2, rewrites `dashboard/professor_demo.js` with per-target tabs, top-4 cards per target (SVG + dock score + QED + SA + ADMET flags), attrition funnel, known-drug baseline row. ~3 hours.
- [ ] **GitHub Pages workflow.** Auto-deploy `dashboard/` on `main` push. ~1 hour.
- [ ] **Share the URL** with the professor. Done.

### Day 6-7 — Buffer
Anything that overran. Last-mile polish. Demo dry run. Notes for the meeting.

### Cost envelope

| Item | Cost |
|---|---|
| GHCR (public image) | $0 |
| Hugging Face Hub (public weights) | $0 |
| GitHub Pages | $0 |
| Cloudflare R2 (~5 GB, no egress) | ~$0.10/month |
| RunPod RTX 3090, ~3 hours total | ~$1 |
| GitHub Actions CI minutes | $0 (public repo, free tier) |
| **One-time demo total** | **~$1.10** |
| **Steady state monthly** | **~$0.10** |

If we keep a GPU pod running ~8 hours/week post-demo for iteration: + ~$8/month. Still trivial.

### What I'm explicitly *not* doing this week

| Skipped | Reason |
|---|---|
| Pocket2Mol env rebuild | Checkpoint not recovered. Deferred to post-demo. |
| GNINA rescoring | Post-Phase 2, after threshold calibration with expert input. |
| AiZynth on every candidate | Slow. Top-N only, post-Phase 2. |
| Bayesian recommender | Phase 3. Cold-starts from Phase 1 telemetry. |
| Sonnet-in-the-loop | Phase 3. Needs telemetry corpus. |
| Obsidian emitter | Phase 3. Independent of agent, but downstream of dashboard. |
| Modal / Replicate migration | Phase 3+ refactor. RunPod via Docker is sufficient now. |
| Full VM / disk snapshot | Docker gives 95% of the value at 10% of the size. |
| Kubernetes / k3s | Two laptops + one GPU pod. K8s pays off at ~10+ workloads. |
| Live FastAPI dashboard | Static is enough for the demo; live UI is Phase 3+. |
| Custom Python CLI | Makefile is sufficient. |

## What we learned from the 2026-05-11 2HYY attempt

The failed run on 2HYY surfaced four issues, three of which are now reflected in this plan and one of which is a code fix already applied:

1. **`ModuleNotFoundError: utils`** (code fix applied). `scripts/sample_for_pocket.py` does `import utils.misc as misc` relative to the TargetDiff repo root, but Python sets `sys.path[0]` to the script's *own* directory (`scripts/`), not the cwd. Setting `cwd=TARGETDIFF_REPO` wasn't enough. Fixed in `run_generation.py` by injecting `PYTHONPATH=TARGETDIFF_REPO` into the subprocess environment.
2. **Google Drive checkpoint folder is dead.** Confirmed permanent. The dev box has a Mar-2023 local copy of the TargetDiff diffusion weights; Pocket2Mol's were never recovered. → Phase 1 mirrors these to HF Hub during the Docker build; Pocket2Mol is deferred.
3. **CPU inference is impractically slow.** 12 min/mol on a single-core baseline becomes much slower on a 2-thread laptop. 60 minutes of denoising produced zero usable molecules before the run was killed. → Phase 1 moves all heavy generation to a RunPod GPU.
4. **Lid-close and an ill-considered `systemctl restart systemd-logind` both kill local runs.** → Architectural fix: long jobs go off the laptop entirely (Step 12 / `cloud_run.sh`).

These are infrastructure problems, not pipeline problems. The pipeline itself ran exactly as designed up to the points where the environment failed it.

## Principles

**Ship the deterministic pipeline before the adaptive one.** A reliable fixed pipeline that produces good results is more valuable than an adaptive one that produces unreliable results intelligently.

**Every claim must be verifiable.** If the system says a molecule has a docking score of -8.2, a human must be able to rerun that docking and get the same number. If the system ranks molecule A above molecule B, the reasoning must be traceable through logged data.

**The tools are not the product. The orchestration and the judgment layer are the product.** Any research group can install fpocket and Vina. What they cannot easily do is wire them into a reproducible, logged, quality-controlled workflow that enforces scientific discipline at every stage.

**Do not automate judgment you do not yet have.** The agent planner is the long-term differentiator, but it must be built on empirical knowledge, not assumed intelligence. Until you have run enough real campaigns to know what good looks like, the system should execute, log, and present — not decide.

**Reproducibility is an artefact.** Recipes drift. Artefacts don't. Whenever there is a choice between "document the steps" and "produce the binary," choose the binary.
