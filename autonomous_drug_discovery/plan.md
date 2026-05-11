# Adaptive Discovery Orchestrator — Architectural Plan

## Current State

**M1 (Working Pipeline) and M2 (Validation) are complete.**

The pipeline runs end-to-end with real tools on real proteins. All five stages are production-grade:

- **Pocket detection**: P2Rank (ML-based, default) and fpocket (geometry-based fallback) both integrated. P2Rank is the default.
- **Molecule generation**: RDKit fragment-based generation integrated. TargetDiff (diffusion) and Pocket2Mol (autoregressive) are both wired into the orchestrator as selectable backends via `--mode targetdiff` / `--mode pocket2mol`. Empirical comparison across validated targets is the next step.
- **Screening**: MolScore (primary backend) or RDKit fallback, with ADMET-AI enrichment (104 properties) on survivors.
- **Docking**: AutoDock Vina (production), TDC Oracle (triage), and simulation stubs all supported.
- **Ranking**: Multi-criteria final ranker (`05_ranking/run_ranking.py`) blends docking + ADMET into a composite score and writes `ranked_candidates.json`. AiZynthFinder retrosynthetic feasibility is wired as an optional add-on for top-N candidates.
- **Telemetry**: Full SQLite logging across all stages.

Validated against three cancer targets with crystallographic ground truth:

| Target | Best Dock | Pocket Distance | Residue Overlap |
|--------|-----------|-----------------|-----------------|
| EGFR (1M17) | -9.32 kcal/mol | 2.7 A | 82% |
| BCR-ABL (2HYY) | -12.59 kcal/mol | 2.7 A | 92% |
| BRAF V600E (6P3D) | -11.20 kcal/mol | 3.1 A | 89% |

The immediate priorities are: M3 domain expert review (non-negotiable), empirical comparison of the generation backends (RDKit / Pocket2Mol / TargetDiff) on validated targets, and tightening screening thresholds with expert input.

## Demo Constraint — 20 quality results for the professor (Active)

The professor presentation is the next forcing function. The deliverable is **20 quality candidate molecules** surfaced through the dashboard, end-to-end through the deterministic pipeline. "Quality" here means: passed tightened screening, docked credibly (better than known-drug baseline or comparable, depending on target), and accompanied by a chemist-readable rationale row in the dashboard. Everything in the near-term roadmap is filtered by whether it serves this deliverable.

**Headline math.** 20 quality results = ~5 targets × ~4 surviving top-ranked candidates. With tightened screening at ~50% survival (Step 4) and docking yielding a workable top decile, that requires ~30-40 generated molecules per target. On a rented GPU at ~30-60 s/molecule for TargetDiff, this is 30-60 GPU-minutes per target, ~3-4 GPU-hours total, < $2 on RunPod / Vast.ai. Wallclock end-to-end (including screening + docking + ranking on CPU) is ~half a day.

**Targets for the demo.** Reuse the three already validated (1M17, 2HYY, 6P3D) and add two more from the existing roadmap (KRAS G12C, JAK2). That gives five targets × four quality molecules = the 20-result floor.

**No gold-plating.** The pipeline already works end-to-end (M1) and validates against crystal ground truth (M2). The presentation does *not* need: a new generator, GNINA rescoring (Step 5), AiZynth on every molecule (Step 6), the Bayesian recommender (3b), Sonnet-in-the-loop (3c), or the Obsidian emitter (3a). It needs the existing pipeline run reproducibly on five targets and rendered cleanly in the existing dashboard. The four blockers are reproducibility (Step 11), cloud GPU (Step 12), checkpoint mirror (Step 1), and dashboard polish (Step 10 scope-cut, see below).

**Easy setup for both contributors.** "Both of us cloning fresh and running on a new box" is in scope. "Anyone in the world cloning and running" is out of scope. A `bash scripts/bootstrap.sh` that produces a working pipeline within one coffee break on a clean Ubuntu/macOS box is the bar.

## Architecture

The system has four layers, listed in order of implementation priority.

### Layer 1 — The Pipeline (Build First)

A deterministic, reproducible workflow that takes a protein structure and produces ranked candidate molecules. No adaptive logic, no LLM reasoning — just tools chained together correctly.

```
PDB file
  → Structure preparation (BioPython)
  → Pocket detection (fpocket, optionally P2Rank)
  → Pocket ranking and selection
  → Molecule generation (TargetDiff or Pocket2Mol)
  → Fast screening (RDKit: QED, SA score, Lipinski, basic toxicity flags)
  → Docking (AutoDock Vina or GNINA)
  → Ranking and report generation
```

Each stage reads from and writes to a standardised intermediate format. Every stage can be run independently and inspected. The pipeline is configured via a YAML file specifying which tools to use at each stage and what thresholds to apply. No dynamic decisions — the configuration is set before the run and does not change during it.

This layer is the foundation. Nothing above it matters if this layer does not produce scientifically credible output.

### Layer 2 — The Telemetry Database (Build Alongside Layer 1)

A structured record of everything the pipeline does. This is not optional infrastructure — it is how you debug, validate, and eventually learn.

Schema covers:

- **Campaigns**: target ID, PDB source, date, pipeline configuration used, who reviewed the output.
- **Pockets**: per-pocket geometry, druggability score, volume, polarity descriptors, whether it was selected for generation.
- **Candidates**: SMILES, generation method, which pocket it was designed for, all computed properties (QED, SA, Lipinski violations, ADMET flags), docking score, final rank.
- **Attrition log**: at each filtering stage, how many molecules entered and how many survived, with the reason for each rejection.
- **Annotations** (added later, by humans): whether a candidate was considered promising by a domain expert, whether it was synthesised, and if so, what happened.

Use SQLite initially. Migrate to PostgreSQL only when you need concurrent access or the dataset outgrows a single file. Do not over-engineer the database before you have data in it.

### Layer 3 — Validation Framework (Build After Layer 1 Works)

The pipeline is worthless until you can demonstrate that it recovers known results. This is not a nice-to-have — it is the scientific credibility gate.

Select 3-5 well-characterised protein targets where the following are publicly known:
- Crystal structure with a bound ligand (so you know the correct binding site).
- Multiple known active compounds (so you can measure recall).
- Published SAR data (so a chemist can evaluate whether generated molecules make sense).

Suggested targets: EGFR (non-small cell lung cancer), BCR-ABL (chronic myeloid leukaemia), BRAF V600E (melanoma). These are textbook oncology targets with extensive public data.

For each target, run the full pipeline and measure:
- Does pocket detection identify the known binding site?
- Do generated molecules occupy the same pharmacophoric space as known actives?
- Do known actives score well in the docking stage? (If your scoring function ranks known drugs poorly, the scoring function is broken, not the drugs.)
- What is the attrition rate at each stage, and is it reasonable?

Document the results honestly. Where the pipeline fails, diagnose why. This produces the evidence needed to justify the next layers.

### Layer 4 — Adaptive Orchestration: Cascade + Sonnet-in-the-loop + Knowledge Graph

Only after Layers 1-3 are solid does adaptive orchestration become meaningful. This layer has four cooperating sub-components, each independently shippable and each contributing to the *compounding* differentiation the project depends on (no single tool here is novel; the value is the loop).

**4a. Cascaded generation pipeline (no agent required).** Use Pocket2Mol for breadth (~7 sec/mol on GPU, novel scaffold exploration) → triage by docking + ADMET → TargetDiff to refine the top 10-20 (denoise around them as 3D seeds) → AiZynthFinder synthesis check → multi-criteria rank. The cascade pattern captures the "best of both" without any agent: speed where it matters, fidelity where it pays off. This is the lowest-risk, highest-value architectural change.

**4b. Bayesian strategy selection.** For each new target, the question "which generator + parameters?" is a multi-armed bandit. Maintain a posterior over each backend's expected composite-score distribution conditional on pocket descriptors (volume, polarity, residue composition, druggability score). Use **Thompson sampling** to choose: occasionally explore weaker arms, mostly exploit the leading one, and the posterior tightens with each campaign. This is the *honest* version of "the system learns" — no fine-tuning, no RL on noisy proxies; a principled exploration/exploitation trade-off over a small set of well-defined choices. Multi-criteria rewards (composite score from Stage 5) blunt the ICLR 2025 "Vina-bigger-is-better" reward-hacking trap.

**4c. Sonnet-in-the-loop (LLM agent).** A Claude Sonnet agent sits next to the orchestrator with three concrete jobs:
- *Retrieval over campaigns.* Given a new target, surface similar past campaigns from the telemetry DB and the open-source corpus (CrossDocked2020, BindingDB, ChEMBL) and recommend a starting strategy with cited evidence.
- *Mid-campaign sanity checks.* Watch attrition funnels and ADMET distributions; flag pathological generation runs ("95% Lipinski violations", "all mols < 200 Da", "no aromatic rings against an aromatic-rich pocket") and suggest a parameter adjustment with rationale.
- *Reporting.* Translate the multi-criteria ranked output into a chemist-readable campaign report — three top candidates, why they ranked there, what's uncertain, what to test next.

Crucially: the agent never decides autonomously. It surfaces structured recommendations the chemist can accept, modify, or reject. Every recommendation is logged so we can later measure the agent's hit rate.

**4d. Obsidian knowledge graph (the moat).** Each campaign emits a folder of Markdown notes with embedded SMILES, scorecards, 2D structure SVGs, and `[[wikilinks]]` to other campaigns sharing scaffolds, pockets, or expert annotations. Tags surface common patterns (`#kinase`, `#hERG-flag`, `#3-step-synth`). Over a year, this becomes a queryable knowledge graph of *what worked, what didn't, and why* — sitting in the chemist's existing tool, not behind a custom UI. This is what makes the data advantage compound: every campaign makes the next campaign cheaper to scope and easier to interpret, and expert annotations live next to the data instead of in lab notebooks.

**What this layer is NOT.** Not a generative-model fine-tuning loop (needs wet-lab data we don't have). Not a fully autonomous agent (overpromising and unsafe). Not a replacement for the chemist's judgment (their domain expertise is the user's input, not the target).

## Tool Selection

### Pocket Detection
- **P2Rank**: ML-based, default backend. 10-20 percentage points better recall than fpocket on standard benchmarks; on EGFR (1M17) places the pocket 2.7 A from the known drug vs fpocket's 6.1 A.
- **fpocket**: Geometry-based fallback. Comparable to P2Rank on large/well-defined pockets (BCR-ABL, BRAF), but markedly worse on smaller/cryptic ones. Kept as a backup when Java is unavailable or for cross-checking.

### Molecule Generation
- **Pocket2Mol**: Fast autoregressive sampling, E(3)-equivariant. Good for generating diverse scaffolds quickly when exploring new chemical space. Lower per-molecule compute cost. Use for breadth.
- **TargetDiff / DiffSBDD**: Diffusion-based, higher fidelity 3D generation. Better for refining specific binding interactions once you have a promising scaffold region. Higher compute cost. Use for depth.

Both should be available as interchangeable modules behind a common interface. The pipeline configuration or (later) the agent selects which to use.

### Screening
- **RDKit**: Core dependency. QED, SA score, Lipinski/Veber rules, PAINS filters, basic physicochemical property calculation. Non-negotiable — install and integrate first.
- **ADMET predictors** (ADMETlab, pkCSM, or similar): Add after RDKit basics are working. These provide heuristic estimates of absorption, metabolism, toxicity. Useful for filtering, but treat as approximate signals, not ground truth.

### Docking and Scoring
- **AutoDock Vina**: Well-established, fast, good enough for initial ranking. Start here.
- **GNINA** (optional upgrade): CNN-based rescoring on top of Vina poses. Better discrimination than Vina alone, but adds complexity. Add after Vina is working and you can measure whether GNINA rescoring actually changes your rankings meaningfully.

### Synthetic Accessibility
- **RDKit SA score**: Basic but fast. Include in the screening stage.
- **Retrosynthetic analysis** (AiZynthFinder or similar): Much more informative than SA score — actually proposes synthesis routes. Add as a later enhancement for top-ranked candidates only (it is slow and heavyweight).

## Roadmap (4 phases)

The roadmap is organised as four phases, in order. Earlier phases are prerequisites for later ones — don't skip ahead. Each phase has concrete success criteria; if they aren't met, stop and diagnose before moving on.

### Phase 1 — Reproducible cloud-GPU run on 5 targets, 20 quality results for the professor

**Goal:** Five validated targets (1M17, 2HYY, 6P3D + KRAS G12C + JAK2), ~30 TargetDiff molecules each, full pipeline through to the dashboard. Output is the 20-quality-result corpus the professor sees. Pocket2Mol stays optional — TargetDiff alone covers the deliverable; Pocket2Mol's env-rebuild work (Step 9) is deferred until after the demo unless TargetDiff results disappoint.

**Critical path (in order):**
1. **Step 11 — Codebase unification.** Submodule TargetDiff, mirror the two `.pt` checkpoints to Hugging Face (Drive is gone, see Step 1), write `scripts/bootstrap.sh`. Both contributors can clone-and-run after this.
2. **Step 12 — Cloud GPU spin-up.** Stand up a RunPod RTX 3090 (~$0.25/hr) with the existing `targetdiff_env` (CUDA 11.7 — no env rebuild needed for the demo). Run `bootstrap.sh` on it.
3. **Step 4 — Tighten screening thresholds** in `default_scoring_config.json` to bring survival from 73-98% down to ~50% so the dashboard isn't drowning in junk.
4. **Generate + dock + rank, five targets × 30 molecules.** Single CLI loop over targets, logged to telemetry. Wallclock ~half a day.
5. **Step 10 (scope-cut, "quality dashboard") — re-render `professor_demo.js`** with the new corpus: per-target tab, top-4 ranked molecule cards (SVG + dock score + QED + SA + key ADMET flags), attrition funnel, an explicit comparison row vs known-drug baseline. Nothing more.

**Success criteria:**
- `bash scripts/bootstrap.sh` produces a working pipeline on a fresh box in < 30 minutes (both of us verify on our own machines)
- Five targets × ~30 TargetDiff molecules generated, screened, docked, ranked, logged
- ≥ 20 candidates labelled "quality" (tightened-screen + docking better than median of generation set; ideally within 1 kcal/mol of known-drug baseline on at least three targets)
- Dashboard loads in a browser, shows the 20 with structures and scores; no manual JSON editing required to rerun
- Total cloud spend < $5

Why this is Phase 1: the professor presentation is the next milestone, and everything past Phase 1 is built on having a credible, reproducible result set in hand. Anything else (Pocket2Mol, Bayesian recommender, Sonnet agent) is downstream of "do we have 20 results we believe in?"

### Phase 2 — Show the professor; capture qualitative feedback
**Goal:** Walk a medicinal chemistry advisor through the dashboard for the targets from Phase 1. Capture which ranked candidates they think are sensible, which look like generative junk, what they'd modify, what's missing. The conversation includes our forward-looking architecture (cascade, Bayesian recommender, Sonnet-in-the-loop, Obsidian knowledge graph) so we can shape Phase 3 with their input.

**Success criteria:**
- 1–2 hour review session with at least one medicinal chemistry expert (Imperial Chemistry / Cancer Research / etc.)
- Per-candidate annotations captured: `(promising | borderline | reject)` + free-text reason
- Tightened screening thresholds (current 73-98% survival → target 40-60%) calibrated against expert judgment
- Decision: which targets are worth running the agent loop on in Phase 3

Why this is Phase 2: M3 is the gate north-star.md calls non-negotiable. Without expert validation that current output is *interpretable*, no amount of agent intelligence makes it useful. The expert is also the one who tells us where the Bayesian reward signal should point.

### Phase 3 — Sonnet-in-the-loop: agent learns from past campaigns + open-source data
**Goal:** Wire a Claude Sonnet agent into the pipeline. It retrieves over the local telemetry DB and the open-source SBDD corpus (CrossDocked2020, BindingDB, ChEMBL), recommends generation strategy and parameters per new target, watches mid-campaign attrition, and writes the chemist-facing report. Bayesian strategy selection (Thompson sampling over backend × parameter combos) sits underneath the agent and provides the principled exploration/exploitation control. Obsidian campaign emitter ships in this phase too — every campaign produces a folder of cross-linked MD notes.

See Layer 4 above for the architectural detail. Three sub-deliverables:
- **3a. Obsidian campaign emitter** — every `orchestrator.py run` writes `campaigns/{date}-{target}-{id}/{summary.md, top10/{mol_xxx.md}, route_trees/...}` with `[[wikilinks]]` and embedded structure SVGs. Independent of agent — ships first.
- **3b. Bayesian recommender** — Thompson sampling over (backend, num_samples, beam_size, screening thresholds) given pocket descriptors. Cold-started from the Phase 1 telemetry; updates after each new campaign. Outputs a recommendation + a confidence band; chemist can override.
- **3c. Sonnet agent shell** — wraps `agent_planner.py` with retrieval (telemetry + open-source), structured recommendations citing evidence, attrition-funnel watcher, and report writer. Logs every recommendation so we can later measure hit rate.

**Success criteria:**
- All three sub-deliverables shipped end-to-end on at least one validated target
- Agent recommendations logged with rationale; baseline campaigns also retained for comparison
- Chemist (Phase 2 advisor) can answer "did the agent's recommendations match yours, and where did it diverge?"

### Phase 4 — Bayesian evaluation of progress; investor / university narrative
**Goal:** Quantify whether the agent loop produces measurably better candidates than the deterministic baseline. Use a Bayesian evaluation framework — credible intervals on the lift in composite score, multi-criteria recall against expert-annotated `promising` molecules, calibration of the Bayesian recommender's posterior. Present the result honestly (positive, negative, or null) with credible intervals, not point estimates.

**Success criteria:**
- Quantitative comparison: agent-loop campaigns vs deterministic-pipeline campaigns on a held-out set of targets, with Bayesian credible intervals on every lift metric
- Calibration plot for the recommender: are 90% credible intervals actually 90% empirically?
- Honest write-up of where the agent helps, where it doesn't, and where the system still needs the chemist
- If lift is real and credible: this is the substantive case for a grant / seed funding (academic labs paying for hosted access, or non-dilutive grants to scale targets and agent quality)
- If lift is null: that's also publishable and tells us where to invest next (better retrieval? better reward? more wet-lab grounding?)

Why Bayesian here: deep-learning chemistry is a graveyard of point-estimate hype. Reporting a single number ("agent improves Vina by 0.4 kcal/mol!") invites the BenevolentAI critique. Reporting a *posterior over lift* with explicit prior, sample size, and credible interval signals the kind of scientific rigour that convinces sceptical reviewers and serious investors.

---

## Supporting technical work (referenced from the phases above)

The numbered steps below are the lower-level technical tasks that the four phases pull from. Many can be done in parallel within a phase.

### Step 1 — Mirror the TargetDiff and Pocket2Mol checkpoints (Phase 1 blocker)
Google Drive folder `1-ftaIrTXjWFhw3-0Twkrs5m0yX6CNarz` returns 404 (confirmed 2026-05-10, re-confirmed 2026-05-11 — treat as permanent). The dev machine still has a copy of `pretrained_diffusion.pt` (33 MB, dated Mar 2023) and `egnn_pdbbind_v2016.pt` (30 MB) from before the takedown. `pretrained_Pocket2Mol.pt` (44.9 MB) is **not** on this machine and would need to come from a collaborator or a paper-author request.

Mirror plan, in this order:
1. **Hugging Face Hub** under our own org (`huggingface-cli upload`). Two model repos: `<org>/targetdiff-pretrained` (holding both `pretrained_diffusion.pt` and `egnn_pdbbind_v2016.pt`) and `<org>/pocket2mol-pretrained` (held empty until Pocket2Mol checkpoint is recovered). HF Hub gives versioned URLs and is free for public repos. ~5 min once `huggingface_hub` is installed.
2. **Fallback: GitHub Release on a fork.** If HF is for any reason off-limits, attach the two `.pt` files to a GitHub release on our internal fork of the agent-harness repo — also free, supports single files up to 2 GB, fetchable via `gh release download` or a plain URL.
3. **Update `docs/targetdiff-setup.md` Step 2 and `docs/pocket2mol-setup.md` Step 3** to point at the mirror, and have `scripts/bootstrap.sh` (Step 11) do the fetch automatically.
4. **For Pocket2Mol weights specifically:** open a polite issue on `pengxingang/Pocket2Mol` asking if they can re-host, and in parallel ask anyone in our academic network who's run Pocket2Mol recently. This is the only checkpoint we cannot mirror from local copies.

Until the mirror exists, every new contributor (and every cloud-GPU rental) is blocked. Until the Pocket2Mol checkpoint is recovered, Pocket2Mol stays out of Phase 1 — it's not on the critical path for the 20-result demo, so this is acceptable.

### Step 2 — Empirical Backend Comparison (TargetDiff vs Pocket2Mol vs RDKit)
TargetDiff and Pocket2Mol are wired into the orchestrator (`--mode targetdiff` / `--mode pocket2mol`), but as of 2026-05-10 neither produces output end-to-end on this machine:

- **TargetDiff** — checkpoint cannot be downloaded (Google Drive folder returns 404; no HuggingFace/Zenodo mirror found). The earlier 2026-04-18 standalone POC (2 BRAF molecules, ~7 kcal/mol) is therefore not reproducible without a backup checkpoint.
- **Pocket2Mol** — install + checkpoint work, but `sample_for_pdb.py`'s InitSample beam search is single-threaded CPU-bound and stalls for 10+ minutes even at `beam_size=100` / `num_samples=3` on EGFR, BRAF, and the repo's own `4yhj.pdb` example. GPU sits at 4–38% util; the bottleneck is the Python controller. Reproduced cleanly across three targets, so this is a setup limitation, not a per-pocket bug.

To unblock Step 2: either (a) source a TargetDiff checkpoint via paper authors / a backup, or (b) patch Pocket2Mol's beam-search loop to use multi-process or torch-vectorised candidate evaluation. Until one of these is done, the comparison can only run on RDKit.

### Step 3 — Per-Campaign Output Directories (DONE)
Campaigns now write to `data/{campaign_id}/{candidates,screened,results}/` when using the full `run` command. Independent stage invocations still use shared `data/candidates/`, `data/screened/`, `data/results/` directories.

### Step 4 — Tighten Screening Thresholds
Current survival rates are 73-98% (target: 40-60%). This means the generator is over-producing junk or filters are too loose. Tighten based on Step 1 expert input. This is a config-only change (`default_scoring_config.json`).

### Step 5 — GNINA Rescoring (After Steps 1-4)
Add GNINA CNN-based rescoring alongside Vina. Needs GPU; download binary from github.com/gnina/gnina/releases. Only worthwhile after screening is properly calibrated.

### Step 6 — AiZynthFinder Retrosynthetic Feasibility
`pip install aizynthfinder` (needs policy models). Apply to top-ranked candidates only — it is slow and heavyweight. Provides route-level synthesis feasibility, far more informative than SA score alone.

### Step 7 — Agent Planner (Only After M3 Feedback)
By this point you have real telemetry data, validated benchmarks, and expert-informed parameter ranges. Now you can build the agent planner with actual knowledge to encode, not guesses.

### Step 8 — Reproducible GPU Workflow for Diffusion Backends
The dev machine now has an NVIDIA RTX 5060 (8 GB VRAM, CUDA 13.2), so TargetDiff and Pocket2Mol can run at full speed locally. The remaining work is making this reproducible for collaborators who may be on different hardware:

- **Local (now):** verify both `targetdiff_env` (CUDA 11.7) and `pocket2mol_env` (CUDA 11.3) work against the system CUDA 13.2 driver via PyTorch's bundled runtime. Document the verification command for each.
- **Portable (campaigns at scale):** standardise on a single template, check it into the repo as a setup script, and document the spin-up sequence in `docs/installation.md` so any contributor can reproduce a GPU run from a fresh box. Cloud target: RunPod / Vast.ai RTX 3090 or A10 (~$0.25-0.40/hr) covers a TargetDiff run on all three validated targets with `num_samples=100` in ~1 hour.
- **Free-tier escape hatch:** Google Colab (T4) for one-off interactive demos. Reproducing the conda envs in Colab is awkward but doable with `pip` equivalents.

Why this is a step and not infrastructure: the empirical comparison in Step 2 is gated on it producing reproducible results across machines, and the agent planner in Step 7 needs the comparison data to make backend selection decisions.

### Step 9 — Rebuild deep-learning conda envs for Blackwell+ GPUs
Confirmed 2026-05-10: both `pocket2mol_env` (pytorch 1.10 + cu113) and `targetdiff_env` (pytorch 1.13 + cu117) target pre-Blackwell architectures. The RTX 5060 (sm_120) on this machine triggers PTX-JIT in the driver's forward-compatibility path, which is so slow per kernel that Pocket2Mol's autoregressive loop never finishes. CPU mode is the only functional workaround today.

Plan:
1. **Pocket2Mol**: new `pocket2mol_env_v2` with `pytorch=2.4.*` + `pytorch-cuda=12.1`, latest `torch_geometric` + matching `torch-cluster`/`torch-scatter`/`torch-sparse` wheels from the PyG wheel index (https://data.pyg.org/whl/). Run `sample_for_pdb.py` against the existing `pretrained_Pocket2Mol.pt` checkpoint — checkpoints are weight files only, they don't pin PyTorch versions. Expect 0–3 small import patches (PyG API churn around `MessagePassing`).
2. **TargetDiff**: same approach for `targetdiff_env_v2` (pytorch 2.4 + cu121 + matching PyG). Gated on first recovering `pretrained_diffusion.pt` (Google Drive 404, paper authors / backup needed).
3. Verify with the trivial diagnostic: `python -c "import torch; x=torch.randn(2000,2000).cuda(); [x@x.T for _ in range(10)]; torch.cuda.synchronize(); print('ok')"` should return in <1 sec, not 2+ minutes.
4. Re-run the empirical backend comparison (Step 2). With a working GPU, Pocket2Mol's quoted ~7 sec/molecule and TargetDiff's ~30 sec/molecule become realistic for full 100-molecule campaigns.

Estimated effort: 30–60 min per env, mostly conda dependency babysitting. The model code itself uses standard PyG ops (no custom CUDA kernels), so the upgrade risk is low.

### Step 10 — Interactive Pipeline Visualisation Frontend
Build a web UI that surfaces the pipeline's work as it happens — campaigns in flight, attrition funnel, generated molecules, docking poses, telemetry queries — instead of leaving everything in CSV/SDF/SQLite. The visual goal is the kind of "live agentic workflow" panel popularised by recent demos (e.g. Chris Yoo's OpenAI hackathon project): clear stage-by-stage status, intermediate artefacts visible, expert-friendly review surface.

Mockups under `reports/ui_mockups/` capture the intended look and feel. This is downstream of Steps 1-7 — the UI should visualise a system that already produces trustworthy results.

**Demo scope-cut (Phase 1).** For the professor presentation, the existing static dashboard (`dashboard/index.html` + `dashboard/professor_demo.js`) is the deliverable, not a live web app. The work is regenerating `professor_demo.js` from the new five-target corpus, adding a per-target tab, surfacing the top-4 ranked molecules per target as cards (SVG + dock score + QED + SA + key ADMET flags), and showing the attrition funnel and a known-drug-baseline comparison row. The live "agentic workflow" panel is a Phase 3+ deliverable; the static dashboard is enough for the 20-result demo. A one-line `scripts/regenerate_dashboard.py` that reads telemetry and rewrites the JSON is the only new code needed for the demo.

### Step 11 — Codebase unification (Phase 1 critical path)
"Cloning fresh and running" must be one command. Today it is six commands across three docs, plus two checkpoints from a now-dead Drive folder, plus a TargetDiff repo that `.gitignore` tells you to "add as submodule" but is not in `.gitmodules`. Fix:

1. **Pin TargetDiff as a proper git submodule.** Remove the `autonomous_drug_discovery/modules/02_generation/targetdiff/` line from the top-level `.gitignore`. Add a `[submodule "…/targetdiff"]` block to `.gitmodules` pointing at `https://github.com/guanjq/targetdiff.git` at the commit currently checked out on the dev machine (`git -C autonomous_drug_discovery/modules/02_generation/targetdiff/ rev-parse HEAD` is the pin). Now `git submodule update --init --recursive` brings both TargetDiff and Pocket2Mol in deterministically.
2. **Mirror the checkpoints (Step 1) and fetch them automatically.** `scripts/bootstrap.sh` calls `huggingface-cli download <org>/targetdiff-pretrained pretrained_diffusion.pt egnn_pdbbind_v2016.pt --local-dir modules/02_generation/targetdiff/pretrained_models/`.
3. **One bootstrap script does everything.** `scripts/bootstrap.sh` performs, idempotently:
   - `git submodule update --init --recursive`
   - `conda env create -f envs/env_orchestrator.yml` (creates / updates the `base`-side env — see point 4 below)
   - `conda env create -f envs/env_targetdiff.yml`
   - `conda env create -f envs/env_docking.yml`
   - `conda env create -f envs/env_pocket2mol.yml` (optional; gated on Pocket2Mol checkpoint availability)
   - Fetch P2Rank to `~/p2rank_2.5.1/` via `wget + tar` if not already present
   - Fetch checkpoints from the HF mirror
   - Run `python orchestrator.py run data/processed/1M17.pdb --mode simulation` as a smoke test
4. **Consolidate the base-env story.** `docs/installation.md` currently says `conda install -n base -c conda-forge rdkit vina openjdk=17 -y` but `envs/env_orchestrator.yml` exists in the repo and is the proper source of truth. Delete the manual `conda install` from `installation.md` Step 1; have it (and the bootstrap script) call `conda env create -f envs/env_orchestrator.yml` instead. Same for `env_docking.yml` — it exists but isn't mentioned in `installation.md`.
5. **Sanity-check on a fresh checkout.** After the above, the test is: blow away `~/miniconda3/envs/targetdiff_env`, re-clone the repo into a scratch directory, run `bash scripts/bootstrap.sh`, run `orchestrator.py run … --mode simulation`. Should succeed without manual intervention. Both contributors do this independently on their own machines.

Effort: ~half a day total (the actual diff is small; the work is making sure the bootstrap is genuinely reproducible).

### Step 12 — Cloud GPU workflow (Phase 1 critical path)
The dev box has no GPU and the targetdiff_env is CPU-only there; one TargetDiff molecule costs ~12 CPU-minutes, so a 30-mol × 5-target campaign on CPU is two days of wallclock — not viable for the demo.

**Provider choice for the demo: RunPod, RTX 3090 or A40, ~$0.25-0.40/hr.** Cheap, persistent network volumes, ssh access, and the `targetdiff_env` (PyTorch 1.13 + CUDA 11.7) runs unmodified on Ampere GPUs — no env rebuild required (Step 9 stays deferred until after the demo). Modal is a stronger long-term answer but requires non-trivial refactor of `run_generation_targetdiff` from `subprocess.check_call` to a `@modal.function`; deferred to Phase 3+.

**Workflow:**
1. Launch RunPod pod (RTX 3090 or A40, Ubuntu 22.04, ~50 GB disk).
2. `git clone <our-repo>` + `bash scripts/bootstrap.sh` (Step 11). Wait one coffee.
3. `python orchestrator.py run data/processed/<target>.pdb --mode targetdiff --num_samples 30` for each of the five targets (or a thin wrapper script that loops). Add a `--num_samples` CLI flag to `run_generation.py` while we're at it — currently the only way to override is to edit `_DEFAULT_PARAMS_BY_MODE`. Single-line change.
4. `rsync` the resulting `data/campaign_*/` directories back to the dev box for dashboard regeneration.
5. Tear down the pod. Cost ≈ ($0.30/hr × ~3 hours) ≈ $1.

**Why not Modal yet:** Modal billing-per-second and Python-native env builds would be ideal, but it requires turning the subprocess-based `run_generation_targetdiff` into a Modal function. Worth doing once the pipeline stabilises; not worth doing for one demo.

**Why not pin a specific provider in the script:** the cloud step is the only one a human definitely does interactively (launching a pod, picking a region, generating an SSH key). Encoding it in shell automation is overkill for this scale. A 20-line README section in `docs/installation.md` is sufficient.

## Principles

**Ship the deterministic pipeline before the adaptive one.** A reliable fixed pipeline that produces good results is more valuable than an adaptive one that produces unreliable results intelligently.

**Every claim must be verifiable.** If the system says a molecule has a docking score of -8.2, a human must be able to rerun that docking and get the same number. If the system ranks molecule A above molecule B, the reasoning must be traceable through logged data.

**The tools are not the product. The orchestration and the judgment layer are the product.** Any research group can install fpocket and Vina. What they cannot easily do is wire them into a reproducible, logged, quality-controlled workflow that enforces scientific discipline at every stage.

**Do not automate judgment you do not yet have.** The agent planner is the long-term differentiator, but it must be built on empirical knowledge, not assumed intelligence. Until you have run enough real campaigns to know what good looks like, the system should execute, log, and present — not decide.

---

## Immediate actions (week of 2026-05-11)

Ordered. Each step's output is the next step's input. None depends on a tool we don't already have.

1. **Mirror checkpoints to Hugging Face** (Step 1). Upload `pretrained_diffusion.pt` and `egnn_pdbbind_v2016.pt` from the dev machine to `<org>/targetdiff-pretrained`. Effort: 30 min. Unblocks every fresh clone.
2. **Pin TargetDiff as a submodule** (Step 11.1). Remove the `.gitignore` line, add to `.gitmodules`, commit. Effort: 15 min.
3. **Write `scripts/bootstrap.sh`** (Step 11.3). Submodule init, conda env create × 4, P2Rank fetch, checkpoint fetch from HF, smoke-test. Effort: ~3 hours including testing on a clean directory.
4. **Add a `--num_samples` flag** to `run_generation.py` so we don't have to keep editing `_DEFAULT_PARAMS_BY_MODE` to control campaign size (Step 12.3). Effort: 15 min.
5. **Both contributors run `bash scripts/bootstrap.sh` on their own machines** and confirm a simulation-mode pipeline run completes. Effort: 15 min each.
6. **Tighten screening thresholds** in `default_scoring_config.json` (Step 4). Pick conservative defaults (Lipinski-strict, QED ≥ 0.5, SA ≤ 4.5, MW ≤ 450). Effort: 15 min — calibrate against the existing telemetry.
7. **Spin up a RunPod RTX 3090** (Step 12). Clone, bootstrap, run TargetDiff on the five targets with `num_samples=30` each. Effort: ~3 hours wallclock, mostly waiting; ~$1.
8. **Regenerate the dashboard** (Step 10 scope-cut). Write `scripts/regenerate_dashboard.py` that reads telemetry + ranking outputs and rewrites `dashboard/professor_demo.js` with per-target tabs and top-4 cards. Effort: half a day.
9. **Sanity-review the 20 results** — spot-check the SDFs in a chemistry viewer (PyMOL or just RDKit MolToImage) before showing the professor. Effort: ~1 hour.

Total: under one focused working week, ~$3 of cloud spend. Anything not in this list is post-demo work, including the Pocket2Mol env rebuild (Step 9), GNINA rescoring (Step 5), AiZynth on every molecule (Step 6), Bayesian recommender (Phase 3b), Sonnet-in-the-loop (Phase 3c), and the Obsidian emitter (Phase 3a).