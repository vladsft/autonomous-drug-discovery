# Adaptive Discovery Orchestrator — Architectural Plan

This is the canonical plan. The [README](../README.md) and other docs under [`docs/`](../docs/) reference this file; if they disagree, this file wins. Last major revision: 2026-05-21 (Phase 1.5 — fire-and-forget batch driver — added).

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
                          ┌────────────────────────────┐
                          │  GitHub                    │
                          │  ├─ build.yml   (CI build) │
                          │  ├─ pages.yml   (deploy)   │
                          │  └─ batch.yml   (dispatch) │◄── workflow_dispatch
                          └────────────┬───────────────┘    (you, once, from
                                       │                     the Actions UI)
              ┌────────────────────────┼──────────────────┐
              │                        │                  │
              ▼                        ▼                  ▼
   ┌──────────────────┐      ┌────────────────┐    ┌────────────────────┐
   │ Hugging Face Hub │      │ GHCR           │    │ GitHub Action       │
   │ <org>/td-weights │      │ ghcr.io/.../   │    │ (batch_cloud_run)   │
   │ (~108 MB)        │      │ harness:<sha>  │    │ Python orchestrator │
   └──────────────────┘      └────────┬───────┘    │ — semaphore=N pods  │
              ▲                       │            │ — sentinel polling  │
              │ baked in              │            │ — cost guard        │
              │ at image-build        │            │ — telemetry merge   │
              │                       ▼            │ — dashboard regen   │
              │             ┌──────────────────┐   │ — auto-commit       │
              │             │ Laptop A         │   └─────────┬───────────┘
              │             │ (interactive,    │             │
              │             │  per-target run) │             │ provisions N pods,
              │             └────────┬─────────┘             │ each running the
              │                      │                       │ same image
              │                      │            ┌──────────┼──────────┐
              │                      │            ▼          ▼          ▼
              │                      │     ┌───────────┐┌───────────┐┌───────────┐
              │                      │     │ RunPod #1 ││ RunPod #2 ││ RunPod #k │
              │                      │     │ TARGET=A  ││ TARGET=B  ││ TARGET=…  │
              │                      │     │ ephemeral ││ ephemeral ││ ephemeral │
              │                      │     └─────┬─────┘└─────┬─────┘└─────┬─────┘
              │                      │           │            │            │
              │                      └──────┬────┴────────────┴────────────┘
              │                             ▼
              │                  ┌────────────────────────┐
              │                  │ Cloudflare R2          │
              └──── weights ──── │ campaigns/<id>/        │
                                 │ sentinels/<id>.done    │
                                 │ telemetry.db (merged)  │
                                 │ dashboard/*.json       │
                                 └───────────┬────────────┘
                                             │ batch action pushes
                                             │ regenerated dashboard/
                                             ▼
                                 ┌────────────────────────┐
                                 │ GitHub Pages           │ ← professor URL
                                 │ multi-target selector  │
                                 └────────────────────────┘
```

Two paths into the same architecture. **Per-target, interactive (Phase 1):** `make cloud-run TARGET=X` from a laptop provisions one pod. **Batched, hands-free (Phase 1.5):** one `workflow_dispatch` of `batch.yml` provisions N pods concurrently, merges telemetry, regenerates the dashboard, deploys. The pods run the same image either way — only the orchestrator above them changes.

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

#### Batch dispatcher: a GitHub Action fans out N pods
**Decision:** A second GitHub Actions workflow, `.github/workflows/batch.yml`, triggered exclusively by `workflow_dispatch`, runs a Python orchestrator (`scripts/batch_cloud_run.py`) that:

1. Reads inputs: `targets` (space-separated list), `mode` (rdkit / targetdiff / pocket2mol), `num_samples`, `parallelism` (default 5), `force` (default false — skip targets with a fresh successful campaign).
2. Pre-flights: checks the RunPod balance via GraphQL; refuses to start if balance < `parallelism × ~$0.30 × estimated_hours × 1.5`. Fetches any missing PDBs from RCSB into R2.
3. Provisions pods in a worker pool (`asyncio.Semaphore(parallelism)`) — *not* an Actions matrix. Pods are independent, Actions runners are scarce; one runner provisioning N pods is the right shape.
4. Each pod runs the existing `pod_campaign.sh` end-to-end and, on success/failure, writes a sentinel file at `r2:bucket/sentinels/<campaign_id>.done` (or `.failed`) before exiting.
5. The orchestrator polls R2 every 60 s for sentinels; collects outcomes; replaces dead pods (RunPod is preemptible) up to a retry ceiling.
6. After all pods drain: `scripts/merge_telemetry.py` consolidates per-campaign telemetry DBs into the canonical R2 `telemetry.db`. `scripts/regenerate_dashboard.py --all-targets` rebuilds the dashboard JSON. The action commits `dashboard/` and pushes to `main`. The Pages workflow picks the change up and redeploys.

**Why a Python orchestrator inside the Action, not Actions matrix:** GitHub free-tier minutes are budgeted; spinning up N Actions runners for ~3 hours of *wall clock* burns ~3N × 60 = a lot. One runner spending those minutes *waiting on R2* burns much less. The runner is doing I/O and ten lines of GraphQL — fits comfortably in a free-tier job.

**Why R2 sentinels, not RunPod status polling:** RunPod's pod-status API is eventually consistent and lies about exit codes when the container terminates ungracefully. A sentinel file written by the pod *just before exit* is authoritative: present ⇒ the pod's payload completed; absent + pod gone ⇒ failure. The orchestrator never has to trust RunPod's view of "did the work finish?"

**Why all of this is the right layer, not Modal / Vertex AI / Kubernetes:** Modal and Vertex are better fits *once we have a stable batch surface to migrate*. Today we don't. The Action + bash + RunPod path is what we already half-built — finishing it costs 2-3 days, and the result is a strict superset of the per-target flow (the per-target path remains unchanged). Modal becomes a Phase 3+ refactor once we know what the surface should look like.

**Concession:** RunPod GPUs can be reclaimed; a long batch may lose a pod mid-run. The orchestrator handles this with at-most-three retries per target. If a target genuinely refuses to run (PDB malformed, weights missing), it gets logged as `.failed` and the batch continues — one bad target must not poison the rest.

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

## Roadmap (5 phases — Phase 1.5 added 2026-05-21)

The shape is: prove the primitive (1), scale it hands-free (1.5), get it in front of an expert (2), make it adaptive (3), and measure honestly (4). Each phase produces something demonstrable before the next one is allowed to start.

### Phase 1 — Containerised cloud-GPU run on a target set, 20 quality results

**Goal:** The architecture above shipped end-to-end. A target set × ~30 TargetDiff molecules each, full pipeline through screening + docking + ranking + dashboard. Output is the 20-quality-result corpus the professor sees.

**Target set (actual, as of 2026-05-21):**

| Target | Disease | Status |
|---|---|---|
| 1M17 (EGFR) | Lung cancer | ✅ Phase-1 campaigns complete (RDKit + Pocket2Mol) |
| 2HYY (BCR-ABL) | CML | ✅ Phase-1 campaigns complete |
| 6P3D (BRAF V600E) | Melanoma | ✅ Phase-1 campaigns complete |
| 8P1L | Internal validation | ✅ TargetDiff campaign complete on local Blackwell GPU |

**Deferred from this phase:** KRAS G12C and JAK2 were in the original plan but not run yet — they're follow-up work, not a blocker for the professor demo. The deliverable (20 quality candidates with the new infrastructure) is met by the four targets above. 1UYD / 3KRR / 6OIM PDBs are checked into `data/processed/` for downstream use but not part of the Phase 1 campaign corpus.

**Success criteria:**
- `docker pull` + `make run TARGET=… MODE=simulation` works on both contributors' laptops with no manual intervention
- `make cloud-run` provisions a RunPod GPU, runs TargetDiff, syncs results to R2, tears down the pod, all from a single command
- The target set × ~30 molecules each generated, screened, docked, ranked, logged to telemetry
- ≥ 20 candidates labelled "quality" (passed tightened screening + docking better than median of generation set; ideally within 1 kcal/mol of known-drug baseline on at least three targets)
- Dashboard loads at `https://<you>.github.io/agent-harness/`; auto-deploys on every commit; no manual JSON editing
- Total cloud spend < $5

**Why this is Phase 1:** The professor presentation is the next milestone, and everything past Phase 1 is built on having a credible, reproducible result set in hand. Anything else (Pocket2Mol, Bayesian recommender, Sonnet agent) is downstream of "do we have 20 results we believe in?"

**Pocket2Mol is explicitly deferred.** Its checkpoint isn't recovered, its env doesn't run on Blackwell, and TargetDiff alone covers the deliverable. We revisit Pocket2Mol post-demo only if TargetDiff results disappoint.

### Phase 1.5 — Fire-and-forget batch over 30-50 targets

**Goal:** One button-press triggers a campaign over 30-50 enzymes. No laptop involvement after kickoff. ~3 hours later, the professor URL shows a multi-target dashboard with all the results.

This is the load-bearing piece between Phase 1 (one-target primitive proved working) and Phase 2 (professor demo). Without it, "scale the corpus" remains a manual chore — 50 invocations of `make cloud-run`, 50 chances to forget one, no laptop allowed to sleep for ~25 hours. With it, the professor demo becomes a function of *target list curation*, not of human babysitting.

**Success criteria:**
- One `workflow_dispatch` of `batch.yml`, given a target list of 30-50 PDB codes, runs the full pipeline on every target on RunPod GPUs.
- The orchestrator runs `parallelism` (default 5) pods concurrently; total wallclock ≈ ceil(N / parallelism) × ~30 min.
- Per-target failure is recorded as a `.failed` sentinel and reported in the Action's summary; the batch never aborts on a single target.
- Telemetry from every pod is merged into the canonical R2 `telemetry.db`.
- The dashboard is regenerated for **all** completed targets and auto-deployed to GitHub Pages without any local commit.
- One notification at the end: GitHub's built-in Action-completed email is enough.
- Total cloud spend < $20 for a 50-target batch (5 concurrent × ~30 min/target × ~$0.30/h = $5 best case; budget 4× for retries + provisioning delay).

**Detailed flow (T+ from clicking "Run workflow"):**

| Time | Where it runs | What happens |
|---|---|---|
| T+0  | You, in browser | Open `Actions → batch_cloud_run → Run workflow`. Fill `targets="1M17 2HYY 6P3D 8P1L …"`, `mode=targetdiff`, `num_samples=30`. Click run. |
| T+5s | Action runner | `scripts/batch_cloud_run.py` parses inputs, validates target codes (4-character alphanumeric), de-dupes. |
| T+10s | Action runner | Query R2 `telemetry.db`: for each target, if a successful campaign exists newer than 24 h (and `force=false`), drop it from the work-list. Print the skip-list. |
| T+15s | Action runner | Cost guard: query RunPod's GraphQL `myself { credits }`. Compute `estimated_spend = parallelism × $0.30 × ceil(len(work_list) / parallelism) × 1.5`. If credits insufficient, fail immediately. |
| T+30s | Action runner | For each missing PDB: `curl https://files.rcsb.org/download/<TARGET>.pdb` → `rclone copy` to `r2:bucket/processed/`. |
| T+1m  | Action runner | Spin up the first `parallelism` pods via RunPod GraphQL `podFindAndDeployOnDemand`. Each pod runs `pod_campaign.sh`. Record `pod_id → target` mapping locally. |
| T+1m  | Pods (parallel) | Each pod: `rclone copy r2:bucket data --progress` → `python orchestrator.py run data/processed/<TARGET>.pdb --mode … --num_samples …` → `rclone copy data r2:bucket --progress` → write `r2:bucket/sentinels/<campaign_id>.done` → exit 0. On error: write `.failed` + error trace, exit non-zero. |
| ~T+30m | Action runner | Polling loop, 60 s tick: list `r2:bucket/sentinels/`. For each new `.done`, mark the target complete; cycle a new pod in from the work-list. For each `.failed`, decrement that target's retry budget and re-queue if budget > 0. |
| T+~3h | Action runner | All sentinels accounted for. Print per-target outcomes. |
| T+~3h+2m | Action runner | `rclone copy r2:bucket/campaigns/ ./data/` (only the new ones, by sentinel timestamp). Run `scripts/merge_telemetry.py` against the canonical R2 `telemetry.db` to fold in per-campaign DBs. Push the merged DB back to R2. |
| T+~3h+3m | Action runner | `python scripts/regenerate_dashboard.py --all-targets` reads the merged telemetry and emits a multi-target `dashboard/professor_demo.{js,json}`. |
| T+~3h+4m | Action runner | `git add dashboard/`, commit with `[skip ci]` on the build workflow but *not* on Pages, `git push` via the workflow's `GITHUB_TOKEN`. |
| T+~3h+5m | Pages workflow | Auto-fires on `paths: ['dashboard/**']` change. Deploys to Pages. |
| T+~3h+7m | You | Email: "batch_cloud_run completed". Open the Pages URL → multi-target dashboard live. |

**Failure modes the design handles:**
- **Pod reclaimed mid-run.** RunPod can terminate a preemptible pod. Detected by `pod_id` disappearing from the status endpoint without a `.done` sentinel; the orchestrator retries (up to 3 times) on a fresh pod.
- **One target's pipeline fails (bad PDB / no pocket).** Pod writes `.failed`, exits non-zero. Orchestrator records it, moves on. Final report flags it; the rest of the batch is unaffected.
- **RunPod outage / GraphQL flakiness.** GraphQL calls wrapped in `tenacity`-style retry-with-backoff (or bare `requests` retry, no need to add the dep). After 5 minutes of failure to provision, the action errors with a clean message instead of silently hanging.
- **R2 outage.** The pod will fail to push results; sentinel won't appear; orchestrator retries. If R2 stays down > 10 min, the orchestrator aborts and reports.
- **Action runner timeout (6h limit on free tier).** For 50 targets at parallelism=5, wallclock is ~5 hours — comfortably under. For 100 targets, bump parallelism to 10 (free-tier RunPod balance permitting) or split into two batch runs.
- **A pod hangs forever (e.g., P2Rank deadlocks on a malformed PDB).** Each pod sets a 90-minute wallclock fuse inside `pod_campaign.sh`; the orchestrator additionally times pods out at 75 min and explicitly terminates them via GraphQL.

**Failure modes the design does NOT handle (by intent, deferred to Phase 3+):**
- **Cross-batch deduplication of identical campaigns** (same target, same mode, same num_samples, recently). The 24h skip-list is the only check; bumping `num_samples` or `mode` produces a new campaign. That's the right behaviour for an exploration phase.
- **Live progress UI.** The Action's stdout is the only progress view. Building a real dashboard for "what's running right now" is Phase 3+ work (this is what an Obsidian campaign emitter would surface).
- **Cost amortisation across batches.** Each batch provisions its own pods. Keeping one pod warm for serial use would save $1-2 per batch but adds complexity (and the lid-close problem reappears in a different form).

**Multi-target dashboard schema (extension of the current single-target shape):**

```jsonc
{
  "targets": {
    "1M17": { "pdb": "1M17", "name": "EGFR (…)", "disease": "NSCLC",
              "known_drug": "Erlotinib", "backends": { … same shape … } },
    "2HYY": { "pdb": "2HYY", "name": "BCR-ABL", "disease": "CML",
              "known_drug": "Imatinib", "backends": { … } },
    "…": { … }
  },
  "default_target": "1M17",          // the first key with a complete campaign
  "generated_at": "2026-05-…",
  "pipeline": { … unchanged … },
  "admet_pass_rules": { … unchanged … },
  "generator_descriptions": { … unchanged … }
}
```

The dashboard JS picks up a top-level *target selector* (dropdown in the header); each target keeps the existing backend tabs / molecule cards / attrition funnel verbatim. No structural rewrite — one element added, one nested level of indexing.

**Components to build (in dependency order):**

1. **Pod-side sentinel.** Append three lines to `scripts/pod_campaign.sh` so the pod writes `sentinels/<campaign_id>.done` (or `.failed` with the exit code and last 50 log lines) right before `exit`. ~15 min.
2. **`scripts/fetch_pdb.py`.** `curl rcsb.org/download/<TARGET>.pdb` → R2. Idempotent; skips PDBs already in R2. ~30 min.
3. **`scripts/batch_cloud_run.py`.** The Python orchestrator described above. Async pool, sentinel polling, RunPod GraphQL with retries, cost guard, retry budget. Reuses `cloud_run.sh`'s GraphQL mutations (copy them into Python; don't shell out — easier to retry and error-handle). ~6-8 hours.
4. **Multi-target `regenerate_dashboard.py`.** Add `--all-targets` flag: enumerate distinct targets in telemetry, run the existing per-target assembly for each, splice into the new top-level schema. ~2 hours including a new test.
5. **Multi-target dashboard JS.** Add a `<select>` element + `setActiveTarget(pdb)` function; everything else stays. ~2 hours including manual browser test on the existing committed dataset.
6. **`.github/workflows/batch.yml`.** `workflow_dispatch` with the inputs above. Single job: install Python deps (`pip install requests`), run `python scripts/batch_cloud_run.py`, configure git for the commit + push step at the end. ~1 hour.
7. **Secrets.** Repo Settings → Secrets and variables → Actions: `RUNPOD_API_KEY`, `RUNPOD_NETWORK_VOLUME_ID`, `R2_ACCESS_KEY_ID`, `R2_SECRET_ACCESS_KEY`, `R2_ENDPOINT`, `R2_BUCKET`. The Action exports these as env vars exactly the way `.env` already shapes them, so `cloud_run.sh` and `batch_cloud_run.py` consume identical config. ~15 min.
8. **Docs update.** New section in `docs/pipeline-guide.md` titled "Batch runs (cloud, hands-free)" with a single example. Update `README.md` Quick Start to mention the batch path alongside `make cloud-run`. ~30 min.

**Total estimated work: ~2 days end-to-end, including a real 5-target dry-run on RunPod to shake out timing and idempotency bugs.**

**Cost envelope (running once, 50 targets, parallelism=5):**

| Item | Cost |
|---|---|
| RunPod 3090 × 5 pods × ~3h wallclock × $0.30/h | ~$4.50 best case |
| Retry budget (assume ~15% pod-loss rate) | +$0.70 |
| R2 storage delta (~5 GB extra per 50-target batch) | ~$0.08/mo |
| GitHub Actions minutes (single ~3h job, free-tier ceiling 2000 min/month) | $0 |
| **One-shot total** | **~$5-6** |

**Why this is Phase 1.5, not Phase 2:** Phase 2 is *show the professor*. Showing them 4 targets manually-curated is qualitatively different from showing them 40 targets the system found on its own. The pitch lands harder with the second story. Phase 1.5 is the difference between a demo and a product.

### Phase 2 — Show the professor; capture qualitative feedback

**Goal:** Walk a medicinal chemistry advisor through the **Phase 1.5 dashboard** — the 30-50-target multi-selector view, not the 4-target Phase 1 dashboard. Capture which ranked candidates they think are sensible, which look like generative junk, what they'd modify, what's missing. The conversation includes the forward-looking architecture (Layer 4 below) so we shape Phase 3 with their input.

**Why the dashboard width matters here:** Showing a chemist 4 targets invites the response "interesting, run more." Showing them 40 targets invites the response "let's talk about which 5 are actually promising and why" — which is the conversation we want. The shift from "demo" to "tool the chemist already wants to use" happens between Phase 1.5 and Phase 2.

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
- [x] **Tighten screening thresholds** in `default_scoring_config.json`. Done — lead-like window (MW 250-450, heavy atoms 18-35, QED ≥ 0.5, SA ≤ 4.5) + chemical-sanity filter.
- [x] **Add `--num_samples` CLI flag** to `run_generation.py`. Done.
- [x] **Generate the target-set corpus.** Done — RDKit + Pocket2Mol on 1M17/2HYY/6P3D + TargetDiff on 8P1L (local Blackwell, since `:cu117` can't target sm_120). KRAS G12C / JAK2 deferred.
- [ ] **`make cloud-run` smoke against one target** to prove the RunPod path. Pending — script exists, end-to-end verification still owed.
- [ ] **Spot-check the SDFs** in PyMOL or via `Chem.MolToImage`. ~1 hour.

### Day 5 — Dashboard polish + deploy (5-6 hours)
- [x] **`scripts/regenerate_dashboard.py`.** Auto-discovers the latest successful campaign per (target, backend) from `telemetry.db` and emits `dashboard/professor_demo.js`. `make dashboard` is now a no-arg call.
- [x] **GitHub Pages workflow.** `.github/workflows/pages.yml` deploys `dashboard/` on every `main` push. Activates the moment PR #1 merges.
- [ ] **Share the URL** with the professor. Pending until PR #1 merges.

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

## Immediate Action Plan (week of 2026-05-21) — Phase 1.5 batch driver

Sequenced so each day's output is the next day's input. Roughly 2 working days end-to-end if there are no RunPod-side surprises.

### Day 1 — Pod-side sentinels + PDB fetch + Python orchestrator skeleton (~5-6 hours)
- [x] **Append sentinel writes to `scripts/pod_campaign.sh`.** Three lines: on success, `rclone touch r2:${R2_BUCKET}/sentinels/${CAMPAIGN_ID}.done`. On non-zero exit, `rclone copy <(echo "$ERR_TAIL") r2:.../${CAMPAIGN_ID}.failed`. ~30 min.
- [x] **`scripts/fetch_pdb.py`.** `argparse --targets`, `requests.get(https://files.rcsb.org/download/<T>.pdb)`, write to a temp dir, `rclone copy` to `r2:bucket/processed/`. Skips PDBs already in R2. Unit-tested against `https://files.rcsb.org/download/1M17.pdb` over the wire (cheap, public). ~45 min.
- [x] **`scripts/batch_cloud_run.py` — first cut.** Stub the orchestrator: parse args, load targets, call `fetch_pdb.py`, print a work-list, exit. No RunPod calls yet. ~1 hour.
- [x] **`scripts/batch_cloud_run.py` — RunPod pool.** Ported the GraphQL mutations to Python with `requests`. Implemented as `ThreadPoolExecutor(max_workers=parallelism)` instead of `asyncio.Semaphore` — same parallelism shape, no asyncio surface area, no new deps. ~3 hours.
- [x] **Local dry-run.** `python scripts/batch_cloud_run.py --targets 1M17 2HYY --parallelism 2 --dry-run` prints the dispatch plan without provisioning pods. ~30 min.

### Day 2 — Polling loop + workflow + smoke test (~5-7 hours)
- [x] **Sentinel-polling loop.** Every 60 s, list `r2:bucket/sentinels/`, mark targets `done` or `failed`, cycle a new pod from the work-list, enforce retry budget (3 retries per target). Outer pod timeout enforces a 90-min ceiling per attempt.
- [x] **Cost guard.** `myself { clientBalance }` query, refuse to start if balance < worst-case estimate. Best-effort: if the field schema shifts, we log and proceed rather than block the batch.
- [x] **Telemetry merge + dashboard regen at the end.** `regenerate_dashboard.py --all-targets` is invoked by the workflow's tail. `merge_telemetry.py` is wired in as the integration point — today's pods write directly to the canonical DB so it's a no-op, marked for activation when per-campaign DBs land.
- [x] **`.github/workflows/batch.yml`.** `workflow_dispatch` with inputs, single ubuntu-latest job (350-min ceiling), env vars wired from repo Secrets, runs `python scripts/batch_cloud_run.py`. Concurrency-grouped so two batches can't race.
- [ ] **First real batch — pending user.** Add repo secrets (RUNPOD_API_KEY, RUNPOD_NETWORK_VOLUME_ID, R2_*), then click Run workflow with `targets="1M17 2HYY"`, `parallelism=2`, `num_samples=10`. Verify R2 sentinels appear and no orphaned pods remain.

### Day 3 — Multi-target dashboard + production batch (~4-5 hours)
- [x] **Multi-target schema in `regenerate_dashboard.py`.** `--all-targets` enumerates distinct targets; new top-level `targets` map; single-target path uses the same schema for consistency. 6 unit tests pass.
- [x] **Multi-target dashboard JS.** `<select id="target-select">` in the header + `setActiveTarget(pdb)` plumbing. Backend pills, summary stats, molecule table, detail panel all re-render on target change. Dock-slider min is the union over all targets × backends so the range is stable.
- [ ] **Production batch — pending user.** Curate the 30-50 target list. Click Run workflow. Walk away. Come back to the Pages URL.

### Known issues (batch driver)

#### Bot dashboard commits have no image — FIXED 2026-05-29

**Symptom (was):** `batch.yml` pinned the pod image to `:${github.sha}` (so RunPod couldn't reuse a stale `:latest`). But the tail of every *successful* batch commits the regenerated `dashboard/` back to `main` as a bot commit via the default `GITHUB_TOKEN`. GitHub does **not** fire workflows from `GITHUB_TOKEN` pushes (loop prevention) — so that commit never triggered `build.yml` and no `:<sha>` image was ever built for it. Any batch dispatched while `HEAD` was that bot commit pinned to a tag GHCR didn't have; pods retried the pull until the 40-min `STARTUP_GRACE_S` aborted them. Confirmed live 2026-05-27 (run #26539113240 on commit `179b92c`, dashboard commit from the prior batch — pod stuck `RUNNING` with `runtime: null`, RunPod logs showing `manifest unknown`).

**Fix:** `batch_cloud_run.py` now runs an **image preflight** before provisioning any pod (`preflight_image` / `image_manifest_exists`). It queries the GitHub Packages API via `gh api /users/<owner>/packages/container/<name>/versions` (also tries `/orgs/…` as fallback) to confirm the requested tag exists; if not, it falls back to `:latest` and mutates `env["IMAGE"]` so every subsequent pod picks up the resolved tag. `:latest` is always safe here because dashboard commits change only `dashboard/`, never anything baked into the image. The workflow also gained `permissions: packages: read` and a `GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}` env on the run step so `gh` is authenticated for the (private) package. Non-ghcr.io references and any probe failure return True (graceful — a flaky probe must never block a real batch; the pod surfaces the real pull error if it actually can't pull).

**Verification:** the new probe correctly returns False for the exact broken tag `:179b92c3155726ceb1fa356e8ef5647bc50dd294`, True for `:latest` and the known-good sha `:ee3e04fd…`, False for garbage tags, and True (graceful) for non-ghcr.io references.

#### Cloud TargetDiff env was unpinned → `libc10_cuda.so` import crash — FIXED 2026-05-29

**Symptom (was):** the first cloud batch to reach TargetDiff generation (run #26663134404, both 1M17 and 2HYY, all 3 retries each) crashed at *import* time — before generating a single molecule — with `OSError: libc10_cuda.so: cannot open shared object file`, raised while `torch_geometric` eagerly imported `torch_cluster`. Root cause: `env_targetdiff.yml` pinned `pytorch=1.13.0` but left the entire PyG stack (`pyg`, `pytorch-cluster/scatter/sparse`) unpinned on the conda `pyg` channel. The image rebuild triggered by the preflight commit (68d9da7) re-solved that env from scratch and pulled a torch build without `libc10_cuda.so` paired with CUDA-built PyG extensions that demand it — a mismatch. The same broken solve also made `torch.cuda.is_available()` return False in `targetdiff_env`, so `--device auto` resolved to `cpu` even on an RTX 3090 pod (a second symptom of the one root cause). The preflight fix worked perfectly here (`[preflight] image OK: …:68d9da7`); this was a separate latent bug it merely exposed by forcing a rebuild.

**Fix:** rewrote `env_targetdiff.yml` to mirror the proven `env_targetdiff_blackwell.yml` structure — conda-forge for the non-torch deps, and a pinned pip stack for torch + PyG from the prebuilt cu117 wheels: `torch==1.13.0+cu117` (from the PyTorch cu117 index, which bundles `libc10_cuda.so`), `torch-scatter==2.1.1+pt113cu117`, `torch-cluster==1.6.1+pt113cu117`, `torch-sparse==0.6.17+pt113cu117` (from `data.pyg.org`), and `torch-geometric==2.3.1`. This makes the build reproducible and, because torch is now a real CUDA build, fixes the device detection too (`--device auto` → cuda on a GPU pod). Needs an image rebuild to take effect; `Dockerfile` installs the env unchanged via `mamba env create`.

**Verified live (cloud batch 26667598702):** the rebuilt image cleared the `libc10_cuda.so` crash and `--device auto` correctly resolved to `cuda` on the RTX 3090 — but surfaced a *third* domino: `lmdb` (in the pip section) fell back to its cffi backend, which JIT-compiles a C module via gcc at first import, and that compile failed in the cu117/py38 image (`cffi.VerificationError: CompileError: gcc failed`). **Follow-on fix:** moved `lmdb` from pip → conda-forge `python-lmdb` (ships a prebuilt binding + bundled liblmdb, no runtime compiler). One more rebuild needed.

### What I'm explicitly *not* doing in this batch (defer to Phase 3+)

| Skipped | Reason |
|---|---|
| Replacing RunPod with Modal | Migrating off the GraphQL API once the batch shape is proven is a clean Phase 3 refactor, not a Phase 1.5 detour. |
| Live progress dashboard | Action stdout + email is enough for hands-free. Live UI is Obsidian/Phase 3. |
| Smart target-list generation | The target list is curated by a human until the Bayesian recommender ships (Phase 3b). |
| Cross-campaign rerank | Each batch is independent; the multi-target dashboard does not rerank across batches. |

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
