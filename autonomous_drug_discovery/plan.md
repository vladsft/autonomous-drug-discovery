# Adaptive Discovery Orchestrator — Architectural Plan

## Current State

**M1 (Working Pipeline) and M2 (Validation) are complete.**

The pipeline runs end-to-end with real tools on real proteins. All four stages are production-grade:

- **Pocket detection**: P2Rank (ML-based, default) and fpocket (geometry-based fallback) both integrated. P2Rank is the default.
- **Molecule generation**: RDKit fragment-based generation integrated. TargetDiff (diffusion) and Pocket2Mol (autoregressive) are both wired into the orchestrator as selectable backends via `--mode targetdiff` / `--mode pocket2mol`. Empirical comparison across validated targets is the next step.
- **Screening**: MolScore (primary backend) or RDKit fallback, with ADMET-AI enrichment (104 properties) on survivors.
- **Docking**: AutoDock Vina (production), TDC Oracle (triage), and simulation stubs all supported.
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

### Layer 4 — The Agent Planner (Build Last)

Only after Layers 1-3 are solid does adaptive orchestration become meaningful. The agent planner requires two things that do not exist yet: real telemetry data from real campaigns, and empirical knowledge of which parameter adjustments improve outcomes for which pocket types. Both come from running the pipeline repeatedly and analysing the results.

When the time comes, the agent's scope is:

**Strategy selection.** Given a new target and its pocket profile, query the telemetry database for similar past campaigns. Recommend which generation method to use (Pocket2Mol for broad exploration of novel space, TargetDiff for refinement around known scaffolds), what filter thresholds to start with, and how many candidates to generate.

**Mid-campaign adjustment.** Monitor attrition rates during a run. If screening is killing more than 95% of candidates, the generation parameters may be miscalibrated — flag this and suggest loosening or tightening constraints. This does not require an LLM. A rule-based system with well-chosen thresholds, informed by data from Layer 3, is more reliable and more interpretable than an LLM making judgment calls about chemistry.

**Reporting.** Summarise campaign results in natural language. Explain why candidates were ranked as they were, what assumptions were made, and where uncertainty is highest. This is a genuine LLM strength — translating structured data into readable narrative — and is the most defensible use of an LLM in this system.

Do not build the agent planner as an autonomous decision-maker. Build it as a recommendation engine that a human chemist can override at every step. Autonomy is earned through demonstrated reliability, not assumed at launch.

## Tool Selection

### Pocket Detection
- **fpocket**: Already integrated. Fast, reliable, well-validated. Keep as primary.
- **P2Rank** (optional addition): ML-based, sometimes catches pockets fpocket misses. Worth adding as a second opinion, not a replacement.

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

## Immediate Next Steps (In Order)

### Step 1 — M3: Domain Expert Review (Highest Priority)
Engage a medicinal chemistry or cancer research collaborator to review pipeline output on the three validated targets. Ask: are the generated molecules sensible, are the rankings meaningful, are the failure modes expected? This feedback shapes filter thresholds, scoring weights, and generation constraints. Without medicinal chemistry judgment, the pipeline produces numbers without meaning.

### Step 2 — Empirical Backend Comparison (TargetDiff vs Pocket2Mol vs RDKit)
TargetDiff and Pocket2Mol are now wired into the orchestrator (`--mode targetdiff` / `--mode pocket2mol`). Run each on the three validated cancer targets (EGFR, BCR-ABL, BRAF) and compare: docking score distribution, QED/SA/ADMET distributions, novelty (Tanimoto to known drugs), and runtime. Document which backend wins under which conditions — the answer will inform the agent planner's strategy selection in Step 7.

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

## Principles

**Ship the deterministic pipeline before the adaptive one.** A reliable fixed pipeline that produces good results is more valuable than an adaptive one that produces unreliable results intelligently.

**Every claim must be verifiable.** If the system says a molecule has a docking score of -8.2, a human must be able to rerun that docking and get the same number. If the system ranks molecule A above molecule B, the reasoning must be traceable through logged data.

**The tools are not the product. The orchestration and the judgment layer are the product.** Any research group can install fpocket and Vina. What they cannot easily do is wire them into a reproducible, logged, quality-controlled workflow that enforces scientific discipline at every stage.

**Do not automate judgment you do not yet have.** The agent planner is the long-term differentiator, but it must be built on empirical knowledge, not assumed intelligence. Until you have run enough real campaigns to know what good looks like, the system should execute, log, and present — not decide.