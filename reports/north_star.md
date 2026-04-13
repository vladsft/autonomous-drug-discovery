# North Star: Adaptive Discovery Orchestrator

## What This Is

A computational drug discovery pipeline that takes a disease-relevant protein target and produces a ranked shortlist of candidate molecules worth synthesising — with full traceability of every decision made along the way.

It is not a replacement for medicinal chemistry expertise. It is infrastructure that automates the tedious, error-prone, repetitive parts of early-stage computational screening so that a trained chemist can focus on the decisions that actually require human judgment.

## The Problem

Early-stage drug discovery involves a repetitive workflow: identify a target protein, find plausible binding sites, generate or screen candidate molecules, filter for drug-likeness and toxicity, score for binding affinity, and decide what to test next. Today, most academic and small-biotech teams do this by manually chaining together a patchwork of open-source tools — P2Rank, AutoDock Vina, RDKit, various generative models — with custom scripts, inconsistent file formats, and no systematic record of what was tried, what failed, and why.

This means three things go wrong consistently:

1. **Wasted cycles.** Teams repeat experiments that have already been tried (by them or by others using similar targets) because there is no structured memory of past campaigns.
2. **Silent attrition.** Molecules fail downstream — in synthesis, in assay, in ADMET profiling — for reasons that were detectable computationally but weren't checked because the pipeline didn't enforce it.
3. **Expert bottleneck.** Computational chemists spend a disproportionate amount of time on pipeline plumbing (format conversion, job management, results aggregation) rather than on the scientific reasoning that actually moves the project forward.

No single open-source pipeline integrates the full workflow from target input through pocket detection, generative design, docking validation, ADMET filtering, and retrosynthetic feasibility assessment. This integration gap is our primary opportunity.

## The Market

The AI drug discovery market sits at approximately **$2.5–5B in 2025**, projected to reach **$8–16B by 2030** at 20–30% CAGR. Oncology dominates at 73% of all AI drug discovery studies. Over 200 AI-influenced drugs are in clinical development, with no FDA approval yet — the first is anticipated in 2026–2027.

The competitive landscape is dominated by well-funded platforms:

- **Schrödinger** ($180–200M annual revenue, ~1,750 customers, $50K–500K+/year licenses) — most commercially mature, zasocitinib in Phase 3
- **Recursion-Exscientia** ($688M merger) — largest public AI drug discovery company
- **Isomorphic Labs** ($600M Series A, $3B in pharma partnerships) — Alphabet/DeepMind
- **Insilico Medicine** — fully AI-native drug (rentosertib) completed Phase IIa, target to candidate in 18 months at ~$150K cost

A critical financial reality: the **50:1 ratio** between announced "biobucks" and actual upfront payments reveals deep industry caution. Big pharma is hedging — paying options, not commitments.

**We are not competing with these companies.** We are building an open, modular, transparent alternative for academic labs and small biotechs who cannot afford Schrödinger licenses but still need to move from target to candidate efficiently. This academic/research segment is growing at the fastest CAGR in the market and is currently underserved by integrated tooling.

## What We Build

A modular, orchestrated pipeline with five stages:

### Stage 1 — Pocket Detection

Input: a PDB structure file for a protein of interest.

The system identifies candidate binding pockets on the protein surface using P2Rank (ML-based, default) or fpocket (geometry-based fallback). Each pocket is ranked by predicted ligandability with pre-computed center coordinates for downstream docking.

P2Rank outperforms fpocket by 10–20 percentage points on standard benchmarks. Our validation confirms this: on EGFR (1M17), P2Rank places the pocket 2.7 A from the known drug binding site with 82% residue overlap, versus fpocket's 6.3 A and 53%.

### Stage 2 — Molecule Generation

Given a prioritised binding pocket, the system generates candidate molecules. Two backends:

- **RDKit fragment-based** (current default): Assembles drug-like scaffolds with functional groups and linkers, sized to fit the detected pocket. Fast, no GPU required, produces valid chemistry.
- **TargetDiff diffusion** (ready, not yet in pipeline): E(3)-equivariant diffusion model that designs molecules conditioned on the 3D pocket shape. Should produce better-fitting candidates but requires hours on CPU (minutes on GPU).

Critical caveat from the literature: an ICLR 2025 paper demonstrated that SBDD models routinely generate molecules with better Vina docking scores than known ligands — but this improvement is largely an artifact of generating larger molecules, not better binders. We must always evaluate molecular weight alongside docking score.

### Stage 3 — Scoring and Filtering

Every candidate passes through multiple computational checks:

- **Drug-likeness.** Lipinski rules, QED, and related heuristics via MolScore (preferred) or RDKit fallback. Effective at eliminating obvious failures.
- **ADMET estimation.** ADMET-AI provides 104 predictions covering absorption, metabolism, toxicity, CYP450 inhibition, hERG liability, and more. AUROC >0.85 for 20 of 31 classification tasks on the TDC benchmark. Performance degrades on novel scaffolds — which is exactly where AI-generated molecules live.
- **Synthetic accessibility.** SA scores flag molecules that would be impractical to make. Retrosynthetic feasibility via AiZynthFinder (planned) will provide route-level assessment.
- **PAINS filters.** Removes molecules with known promiscuous substructures.

No single filter is reliable alone. The value is in the combination.

### Stage 4 — Docking

Binding affinity estimation via AutoDock Vina, using the pocket center from Stage 1 to place the docking box. Each molecule is prepared via Meeko, docked with exhaustiveness=8, and scored.

Honest limitation: the Pearson correlation between docking scores and experimental binding affinity is only **r = 0.4–0.6**. Docking scores are a ranking signal, not a measurement. Virtual screening hit rates from docking run 1–5% traditionally, improving to 10–30% with ML-enhanced approaches. We treat docking as one signal among many, not the final answer.

### Stage 5 — Reporting and Handoff

The system produces a structured research artefact for human review:

- Ranked candidates with per-molecule scorecards showing every metric computed.
- ADMET profiles highlighting toxicity risks and metabolic liabilities.
- Explicit assumptions and limitations.
- All results logged to a SQLite telemetry database for full provenance tracking.

This is the deliverable. Not a cure, not a clinical candidate, not a paper. A prioritised, transparent, auditable set of hypotheses that helps a researcher decide what to make and test next.

## Where We Are

### Validated Results (M1 + M2 complete)

The pipeline runs end-to-end with real tools on real proteins. Validated against three well-characterised cancer targets where the answers are known:

| Target | Disease | Known Drug | Best Dock Score | Avg Dock | Pocket Distance | Residue Overlap |
|--------|---------|------------|----------------|----------|-----------------|-----------------|
| EGFR (1M17) | Lung cancer | Erlotinib | -9.32 kcal/mol | -6.58 | 2.7 A | 82% |
| BCR-ABL (2HYY) | Leukemia | Imatinib | -12.59 kcal/mol | -9.25 | 2.7 A | 92% |
| BRAF V600E (6P3D) | Melanoma | Ponatinib | -11.20 kcal/mol | -8.39 | 3.1 A | 89% |

The pipeline independently finds molecules scoring in the same range as known approved drugs, without any knowledge of those drugs. Pocket detection places the docking box within 3 A of the crystallographic drug position for all three targets. This is the scientific credibility gate — the system recovers what is already known.

P2Rank outperformed fpocket in a head-to-head comparison on all three targets: EGFR pocket placement improved from 6.1 A to 2.7 A, and residue overlap improved from 53% to 82%. BCR-ABL and BRAF showed comparable geometry (both methods placed pockets well) with P2Rank maintaining a slight edge in residue overlap.

### TargetDiff Proof-of-Concept (M2.5)

The TargetDiff E(3)-equivariant diffusion model was run standalone on the BRAF V600E pocket, generating molecules from random noise conditioned on the 3D pocket shape (1000 denoising steps, CPU inference, ~12 min/molecule):

| Molecule | MW | QED | Dock Score | Ligand Efficiency | Tanimoto to Ponatinib |
|----------|----|-----|-----------|-------------------|----------------------|
| `C1=CN=CC=C(c2cccc(Nc3ccncc3)c2)C1` | 261 | 0.899 | -7.59 kcal/mol | 0.38 | 0.186 |
| `COc1cnc(C(=O)NCc2cccc(C)n2)cn1` | 258 | 0.889 | -7.38 kcal/mol | 0.39 | 0.154 |

Both molecules are drug-like, compact, and structurally novel (low similarity to known drugs). Docking scores are moderate (-7.4 to -7.6), weaker than the best RDKit-generated candidates (-11.2) but within the "worth investigating" range. Ligand efficiency (0.38-0.39) is excellent. TargetDiff is not yet integrated into the orchestrator pipeline.

### Current Tool Stack

| Component | Tool | Status |
|-----------|------|--------|
| Pocket detection | P2Rank (ML-based) | Integrated, default |
| Pocket detection | fpocket (geometry-based) | Integrated, fallback |
| Molecule generation | RDKit fragment-based | Integrated |
| Molecule generation | TargetDiff diffusion | Env ready, standalone POC done, not in orchestrator |
| Screening | MolScore (primary) + RDKit fallback + ADMET-AI (104 properties) | Integrated |
| Docking | AutoDock Vina + Meeko | Integrated |
| Telemetry | SQLite | Integrated |
| Benchmarking | benchmark.py | Integrated |

## Where We Are Going

### M3 — Domain Expert Integration (next)

A computational chemistry or cancer research collaborator reviews pipeline output and feeds back: are the generated molecules sensible, are the rankings meaningful, are the failure modes expected or surprising. This feedback shapes filter thresholds, scoring weights, and generation constraints.

**This is non-negotiable.** Without medicinal chemistry judgment, the pipeline produces numbers without meaning. Imperial's Chemistry department or Cancer Research center should be engaged now.

### M4 — First Novel Campaign

Run the pipeline against a target with genuine unmet need — where the answer is not known — and produce a candidate set that a research group considers worth synthesising. This is the point at which the system produces real scientific value.

Oncology kinase targets are the most computationally tractable starting point: abundant crystal structures, well-defined binding pockets, extensive SAR data. KRAS G12D, novel EGFR resistance mutations, or emerging kinase targets would be strong M4 candidates.

### M5 — Adaptive Planning

With sufficient campaign history, implement an AI planning layer that adjusts pipeline parameters mid-campaign based on observed attrition. This is the long-term differentiator, but it earns its existence only after M3-M4 are solid.

### Technical Roadmap

Near-term (M3 timeline):
- Wire TargetDiff into orchestrator pipeline (standalone POC done, needs config YAML integration)
- Per-campaign output directories (prevent file collisions between concurrent runs)
- Add GNINA CNN-based rescoring alongside Vina (needs GPU)
- Add AiZynthFinder retrosynthetic feasibility
- Tighten screening thresholds based on expert feedback (current survival: 73-98%, target: 40-60%)

Medium-term:
- AlphaFold/ESMFold integration for targets without crystal structures
- Multi-pocket docking (top 3 pockets, not just best)
- Docker containerization for reproducible deployment
- Publishable methods paper (Journal of Cheminformatics or JCIM)

## Strategic Positioning

**What we are:** An open, modular, well-engineered orchestration layer that makes state-of-the-art computational drug discovery accessible to researchers who lack the software engineering capacity to build pipelines themselves.

**Who we serve:** Academic labs, small biotechs, and chemistry groups at institutions who have domain expertise but not the engineering resources to chain P2Rank, TargetDiff, RDKit, Vina, ADMET-AI, and AiZynthFinder into a reproducible, documented workflow.

**How we differentiate:** Not by algorithmic novelty (the individual tools are freely available), but by integration quality, reproducibility, transparency, and honest communication of limitations. The COVID Moonshot project proved that open-source, community-driven drug discovery can produce real candidates. The hunger for accessible, integrated tools is genuine.

**Business model path:** Freemium open-source — release the pipeline openly to build adoption and credibility, then monetize through cloud-hosted premium features, consulting/customization, or partnership deals based on demonstrated utility.

## What This Is Not

**It is not a drug.** Nothing this system produces is a therapeutic. It produces hypotheses for experimental validation. The wet lab bottleneck is real: synthesising a single compound costs $500–5,000+ and takes weeks.

**It is not a replacement for domain expertise.** The system's output requires evaluation by someone who understands medicinal chemistry, structural biology, and the specific disease context. The system accelerates their work — it does not replace their judgment.

**It is not competing with Schrödinger or Recursion.** A 2–3 person team cannot compete on pipeline depth, proprietary data, or clinical validation. We compete on accessibility, transparency, and community adoption in an underserved market segment.

## Honest Risk Assessment

**The scoring function gap undermines all downstream claims.** If docking scores correlate with experimental binding at only r = 0.4–0.6, then ranking candidates by docking score is ranking by a noisy proxy. Every result we present must acknowledge this limitation.

**The commoditization clock is ticking.** AlphaFold 3 was published in May 2024; Chai Discovery open-sourced a comparable model by September 2024. Over 200 foundation models for drug discovery have been published since 2022. Generic generative models are being rapidly commoditized. Our value must come from integration quality, not algorithmic uniqueness.

**ADMET prediction degrades on novel chemistry.** ADMET-AI performs well on known scaffolds but poorly on out-of-distribution structures — which is exactly what generative models produce. We must communicate confidence levels, not binary pass/fail.

**Defensibility is limited for pure software.** Durable moats come from proprietary experimental data, integrated wet-dry lab platforms, and clinical pipeline IP — not from software that chains freely available tools. Long-term defensibility requires community adoption and real-world validation data.

**The BenevolentAI cautionary tale.** Founded 2013, raised hundreds of millions, went public at ~$2B valuation, then watched their lead clinical asset fail Phase 2a with no benefit over placebo. The target and mechanism were already well-known — AI's contribution was questionable. Stock collapsed from >$9 to <$2. Lesson: AI is not a substitute for good target selection, and overpromising destroys credibility.

## Operating Principles

1. **Don't oversell.** Position honestly: we generate ranked hypotheses for experimental validation, not "discover drugs." Intellectual honesty, rare in this space, is itself a differentiator.
2. **Validate before claiming.** Every capability must be benchmarked against known outcomes before being used on novel targets.
3. **Show your work.** Full provenance tracking, explicit limitations, auditable decisions. If a result can't be traced back to its inputs, it doesn't count.
4. **The domain expert is the customer.** Build for the computational chemist who needs to decide what to synthesise next, not for the investor who wants to hear about AI.
5. **Publish early and openly.** A well-documented, benchmarked, open-source pipeline paper establishes credibility and attracts the community contributions that any future commercial model requires.
