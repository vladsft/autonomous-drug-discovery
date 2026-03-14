# North Star: Adaptive Discovery Orchestrator

## What This Is

A computational drug discovery pipeline that takes a disease-relevant protein target and produces a ranked shortlist of candidate molecules worth synthesising — with full traceability of every decision made along the way.

It is not a replacement for medicinal chemistry expertise. It is infrastructure that automates the tedious, error-prone, repetitive parts of early-stage computational screening so that a trained chemist can focus on the decisions that actually require human judgment.

## The Problem

Early-stage drug discovery involves a repetitive workflow: identify a target protein, find plausible binding sites, generate or screen candidate molecules, filter for drug-likeness and toxicity, score for binding affinity, and decide what to test next. Today, most academic and small-biotech teams do this by manually chaining together a patchwork of open-source tools — fpocket, AutoDock Vina, RDKit, various generative models — with custom scripts, inconsistent file formats, and no systematic record of what was tried, what failed, and why.

This means three things go wrong consistently:

1. **Wasted cycles.** Teams repeat experiments that have already been tried (by them or by others using similar targets) because there is no structured memory of past campaigns.
2. **Silent attrition.** Molecules fail downstream — in synthesis, in assay, in ADMET profiling — for reasons that were detectable computationally but weren't checked because the pipeline didn't enforce it.
3. **Expert bottleneck.** Computational chemists spend a disproportionate amount of time on pipeline plumbing (format conversion, job management, results aggregation) rather than on the scientific reasoning that actually moves the project forward.

## What We Are Building

A modular, orchestrated pipeline with five stages:

### Stage 1 — Target Ingestion and Context

Input: a PDB structure file for a protein of interest.

The system prepares the structure for downstream analysis and gathers available context: known ligands from public databases (ChEMBL, PDB ligand records), relevant literature on the target's biology, existing drug programmes against the same target or family. The output is a structured target profile that informs every subsequent stage.

This stage is important because generative models and docking tools are only as good as the context they operate within. A molecule designed without awareness of known SAR (structure-activity relationships) for the target is a molecule designed in ignorance.

### Stage 2 — Structural Analysis and Pocket Detection

Using fpocket and related tools, the system identifies candidate binding pockets on the protein surface. Each pocket is characterised by volume, druggability score, and geometric properties.

Where possible, the system cross-references detected pockets against known binding sites in the literature. The goal is not just to find pockets, but to prioritise them — which sites are most likely to be pharmacologically relevant, and which are structural artefacts.

The output is one or more ranked binding sites with associated geometric data, ready to condition the generative stage.

### Stage 3 — Molecule Generation

Given a prioritised binding pocket, the system generates candidate molecules designed to complement the pocket geometry. This uses structure-conditioned generative models (currently TargetDiff; the architecture supports swapping in alternatives as the field evolves).

Critical constraint: generated molecules must be filtered immediately for synthetic accessibility. A molecule that cannot be reasonably synthesised in a lab is not a candidate — it is noise. This stage must enforce that boundary aggressively, because generative models have no intrinsic incentive to produce makeable molecules.

The output is a set of structurally novel candidate molecules, each associated with the pocket it was designed for.

### Stage 4 — Scoring and Filtering

Every candidate passes through a series of computational checks:

- **Drug-likeness.** Lipinski and related heuristics via RDKit. Not predictive on their own, but effective at eliminating obvious failures.
- **ADMET estimation.** Absorption, distribution, metabolism, excretion, and toxicity heuristics. These are approximate, not definitive — the system must communicate confidence levels, not binary pass/fail.
- **Binding affinity scoring.** Docking against the target pocket using established tools. Docking scores are noisy and weakly correlated with true affinity — the system treats them as a ranking signal, not a measurement.
- **Synthetic accessibility.** SA scores and retrosynthetic feasibility checks where possible.

No single filter is reliable alone. The value is in the combination — molecules that survive all filters are meaningfully more likely to be worth a chemist's attention than unfiltered generative output.

The output is a ranked candidate list with per-molecule scorecards showing every metric computed and the assumptions behind each.

### Stage 5 — Reporting and Handoff

The system produces a structured research artefact for human review:

- Ranked candidates with full rationale (why this molecule, for this pocket, with these scores).
- Explicit assumptions and limitations (where the scoring is weak, where the model is extrapolating, where literature support is thin).
- Suggested next steps for experimental validation, ordered by information value — what experiment would most efficiently confirm or eliminate the top candidates.

This is the deliverable. Not a cure, not a clinical candidate, not a paper. A prioritised, transparent, auditable set of hypotheses that helps a researcher decide what to make and test next.

## The Orchestration Layer

The five stages above are not novel individually. What is novel is wiring them into a single reproducible workflow with:

- **Full provenance tracking.** Every decision, parameter, score, and intermediate output is logged in a structured database. Any result can be traced back to the exact inputs and configuration that produced it.
- **Campaign-level memory.** When multiple campaigns are run against related targets, the system retains what was tried and what failed. Over time, this builds an institutional memory that no individual researcher carries.
- **Adaptive strategy (future).** An AI planning layer that adjusts pipeline parameters mid-campaign based on observed attrition rates. If a generative model is producing molecules that consistently fail the ADMET filter, the planner tightens generation constraints or switches models. This is not implemented yet and depends on having sufficient real campaign data to learn from. We do not pretend it works today.

## What This Is Not

**It is not a drug.** Nothing this system produces is a therapeutic. It produces hypotheses for experimental validation.

**It is not a replacement for domain expertise.** The system's output requires evaluation by someone who understands medicinal chemistry, structural biology, and the specific disease context. The system accelerates their work — it does not replace their judgment.

**It is not validated yet.** The pipeline exists as a working skeleton. The individual tools are established, but the integrated pipeline has not yet been benchmarked against known outcomes. Validation — running well-studied targets where the right answers are known and checking whether the system recovers them — is the immediate next milestone.

**It is not the adaptive agent yet.** The learning-across-campaigns vision is the long-term goal, not the current capability. The current system runs a fixed pipeline with configurable parameters. Adaptive planning requires real experimental feedback data that we do not yet have.

## Milestones

### M1 — Working Pipeline (Current)
Pipeline runs end-to-end on a real protein target with real tools (not simulation stubs). RDKit screening functional. Real docking scores. Output reviewed by a domain expert.

### M2 — Validation Against Known Targets
Run the pipeline against 3-5 well-characterised targets (e.g., EGFR, BCR-ABL, BRAF V600E) where known active compounds exist. Measure whether the system recovers known ligands or proposes structurally similar candidates. This is the scientific credibility gate — if the system cannot recover what is already known, it cannot be trusted to propose what is unknown.

### M3 — Domain Expert Integration
A computational chemistry or cancer research collaborator actively reviews pipeline output and feeds back on quality: are the generated molecules sensible, are the rankings meaningful, are the failure modes expected or surprising. This feedback shapes filter thresholds, scoring weights, and generation constraints.

### M4 — First Novel Campaign
Run the pipeline against a target with genuine unmet need — where the answer is not known — and produce a candidate set that a research group considers worth synthesising. This is the point at which the system produces real scientific value.

### M5 — Adaptive Planning
With sufficient campaign history (M4 repeated across multiple targets), implement the adaptive planning layer. The system begins adjusting strategy based on observed attrition patterns. This is the long-term differentiator, but it earns its existence only after M1-M4 are solid.

## Why We Might Succeed

The individual tools are mature and freely available. The integration — turning a manual, ad hoc workflow into a reproducible, logged, orchestrated pipeline — is genuine engineering value that most small research teams cannot build themselves. The team combines software engineering capability with chemistry training, and has access to domain expertise through academic collaborators.

The competitive landscape is dominated by well-funded companies building proprietary end-to-end platforms. We are not competing with them. We are building an open, modular, transparent alternative for academic labs and small biotechs who cannot afford Schrödinger licenses or Recursion partnerships but who still need to move from target to candidate efficiently.

## Why We Might Fail

The generative models may not produce molecules that real chemists take seriously. The scoring layer may not be discriminative enough to meaningfully rank candidates. The adaptive planning vision may require more data than we can realistically accumulate. And the team may not secure the sustained domain expertise needed to validate the pipeline's scientific output.

These are not reasons to stop. They are the specific risks we are working to retire, in order, starting with M1.