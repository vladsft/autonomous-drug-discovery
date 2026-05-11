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

### Phase 1 — Get TargetDiff working on cloud GPU; pull more validated targets
**Goal:** Reproducible TargetDiff + Pocket2Mol runs on a rented GPU instance, applied across more than the current three kinase targets (e.g. add KRAS G12C, JAK2, CDK4/6, ER-alpha, AR — all well-characterised oncology targets with crystal structures). No agent yet, no Sonnet — just the deterministic cascade running cleanly at scale.

**Success criteria:**
- TargetDiff `pretrained_diffusion.pt` recovered (paper authors / colleague / mirror search) — see Step 1 below for the blocker
- `pocket2mol_env_v2` and `targetdiff_env_v2` rebuilt on PyTorch 2.4 / cu121 (Step 9 below) so we get the quoted ~7 sec/mol and ~30 sec/mol on a real GPU
- 5–10 targets, 100 mols each, full pipeline (P2Rank → cascade gen → ADMET screen → Vina dock → AiZynth on top 10 → Stage 5 rank), all logged to telemetry
- Cost estimate: $5-10 on RunPod / Vast.ai for an RTX 3090 or A10 (~$0.25-0.40/hr × ~10 hrs total)

Why this is Phase 1: without diverse, reproducible per-target results we have nothing to show the professor. Anything more sophisticated downstream is built on this data.

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

### Step 1 — Recover TargetDiff checkpoint (Phase 1 blocker)
Google Drive folder `1-ftaIrTXjWFhw3-0Twkrs5m0yX6CNarz` returns 404 (confirmed 2026-05-10). Need to source `pretrained_diffusion.pt` from one of: paper authors (open an issue on the GitHub repo), an academic colleague who has a copy, or an alternative mirror. Until this is recovered, Phase 1 cannot ship TargetDiff results — only Pocket2Mol + RDKit.

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

## Principles

**Ship the deterministic pipeline before the adaptive one.** A reliable fixed pipeline that produces good results is more valuable than an adaptive one that produces unreliable results intelligently.

**Every claim must be verifiable.** If the system says a molecule has a docking score of -8.2, a human must be able to rerun that docking and get the same number. If the system ranks molecule A above molecule B, the reasoning must be traceable through logged data.

**The tools are not the product. The orchestration and the judgment layer are the product.** Any research group can install fpocket and Vina. What they cannot easily do is wire them into a reproducible, logged, quality-controlled workflow that enforces scientific discipline at every stage.

**Do not automate judgment you do not yet have.** The agent planner is the long-term differentiator, but it must be built on empirical knowledge, not assumed intelligence. Until you have run enough real campaigns to know what good looks like, the system should execute, log, and present — not decide.