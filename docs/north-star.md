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

The AI drug discovery market sits at approximately **$2.5-5B in 2025**, projected to reach **$8-16B by 2030** at 20-30% CAGR. Oncology dominates at 73% of all AI drug discovery studies. Over 200 AI-influenced drugs are in clinical development, with no FDA approval yet — the first is anticipated in 2026-2027.

The competitive landscape is dominated by well-funded platforms:

- **Schrodinger** ($180-200M annual revenue, ~1,750 customers, $50K-500K+/year licenses) — most commercially mature, zasocitinib in Phase 3
- **Recursion-Exscientia** ($688M merger) — largest public AI drug discovery company
- **Isomorphic Labs** ($600M Series A, $3B in pharma partnerships) — Alphabet/DeepMind
- **Insilico Medicine** — fully AI-native drug (rentosertib) completed Phase IIa, target to candidate in 18 months at ~$150K cost

A critical financial reality: the **50:1 ratio** between announced "biobucks" and actual upfront payments reveals deep industry caution. Big pharma is hedging — paying options, not commitments.

**We are not competing with these companies.** We are building an open, modular, transparent alternative for academic labs and small biotechs who cannot afford Schrodinger licenses but still need to move from target to candidate efficiently. This academic/research segment is growing at the fastest CAGR in the market and is currently underserved by integrated tooling.

## What We Build

A modular, orchestrated pipeline with six stages (a synthesizability gate was added at Stage 2.5 in 2026-05; see the caveat below):

### Stage 1 — Pocket Detection

Input: a PDB structure file for a protein of interest. The system identifies candidate binding pockets using P2Rank (ML-based, default) or fpocket (geometry-based fallback). Each pocket is ranked by predicted ligandability with pre-computed center coordinates for downstream docking.

### Stage 2 — Molecule Generation

Given a prioritised binding pocket, the system generates candidate molecules. The recommended mode is **cascade** — RDKit fragment-based (CPU, the *makeable* workhorse) + TargetDiff diffusion (E(3)-equivariant, highest-fidelity 3D, GPU) merged into one pool. **Pocket2Mol was dropped** (Blackwell-incompatible, CPU-only there, little unique makeable matter). TargetDiff is now understood as a *binding-mode proposer*, not a finished-candidate producer.

Two critical caveats from the literature, both confirmed empirically here:
1. **Docking-score inflation** (ICLR 2025): SBDD models routinely beat known ligands on Vina largely by generating *larger* molecules, not better binders — so we evaluate molecular weight / ligand efficiency alongside docking score, and our composite ranker down-weights size.
2. **Synthesizability collapse** (GenBench3D; measured here 2026-05): 3D generators produce mostly *unmakeable* molecules — **0 of 28** top kinase candidates had a retrosynthetic route; TargetDiff scored 0% makeable vs RDKit's 30%. This is exactly the "silent attrition" pain point above, and it is why synthesizability is now an enforced **gate** (Stage 2.5), not a hope.

### Stage 2.5 — Synthesizability Gate

Before docking, candidates pass through an AiZynthFinder retrosynthetic search (with an optional fast RAScore pre-filter); only molecules with a real route to purchasable building blocks proceed. This turns "we should check synthesizability eventually" into "the pipeline refuses to advance molecules that can't be made" — directly closing the silent-attrition gap that kills computational candidates downstream.

### Stage 3 — Scoring and Filtering

Every candidate passes through drug-likeness (Lipinski, QED), ADMET estimation (ADMET-AI, 104 predictions), synthetic accessibility, and PAINS filters. No single filter is reliable alone. The value is in the combination.

### Stage 4 — Docking

Binding affinity estimation via AutoDock Vina. Honest limitation: the Pearson correlation between docking scores and experimental binding affinity is only **r = 0.4-0.6**. We treat docking as one signal among many, not the final answer.

### Stage 5 — Reporting and Handoff

Ranked candidates with per-molecule scorecards, ADMET profiles, explicit assumptions and limitations. All results logged to a SQLite telemetry database for full provenance tracking. The chemist-facing dashboard (multi-backend toggle, AiZynth route trees, ADMET pass/fail badges) is the primary review surface; in parallel, each campaign emits an Obsidian-friendly Markdown bundle so per-target knowledge accumulates in the chemist's own tool over time (see Adaptive Layer below).

This is the deliverable. Not a cure, not a clinical candidate, not a paper. A prioritised, transparent, auditable set of hypotheses that helps a researcher decide what to make and test next.

### The Adaptive Layer (in development)

The deterministic five-stage pipeline above is the foundation. On top of it sits a four-part adaptive layer that turns the system from a one-shot pipeline into one that learns from its own campaigns:

- **Cascade generation (shipped, 2026-05)** — the original "Pocket2Mol for breadth, TargetDiff for refinement" idea evolved, on contact with data, into something sharper: **RDKit for makeable breadth + TargetDiff for novel binding modes**, unified by a **pharmacophore bridge** that uses TargetDiff's docked poses as binding-mode hypotheses and ranks *makeable* candidates by how faithfully they reproduce them, with AiZynthFinder enforcing synthesizability throughout. This is the realized "best of both" — TargetDiff's pocket insight, expressed in chemistry you can actually make. (Pocket2Mol was dropped; see [plan.md](../autonomous_drug_discovery/plan.md) "Pipeline v2".)
- **Bayesian strategy selection** — Thompson sampling over (backend, parameters) given pocket descriptors. The system maintains a posterior over which configuration works for which pocket family, exploring less as the posterior tightens. Multi-criteria reward (composite score) blunts the known Vina-bigger-is-better reward-hacking trap.
- **Sonnet-in-the-loop** — A Claude Sonnet agent retrieves over local telemetry + open-source corpora (CrossDocked2020, BindingDB, ChEMBL), surfaces structured recommendations with cited evidence, watches mid-campaign attrition for pathological runs, and writes the chemist-facing campaign report. Recommendations are never auto-applied; the chemist accepts or overrides each one, and every decision is logged.
- **Obsidian knowledge graph** — Each campaign emits a folder of cross-linked Markdown notes that becomes a queryable record of what worked, what failed, and what the expert thought. The graph compounds with each campaign — this is what makes the data advantage durable.

What the adaptive layer is *not*: it is not a generative-model fine-tuning loop (no wet-lab labels), not a fully autonomous agent (overpromising and unsafe), and not a replacement for the chemist's judgment.

## Where We Are

### Validated (M1 + M2 complete)

The pipeline runs end-to-end with real tools on real proteins. Validated against three well-characterised cancer targets where the answers are known:

| Target | Disease | Best Dock Score | Pocket Distance | Residue Overlap |
|--------|---------|-----------------|-----------------|-----------------|
| EGFR (1M17) | Lung cancer | -9.32 kcal/mol | 2.7 A | 82% |
| BCR-ABL (2HYY) | Leukemia | -12.59 kcal/mol | 2.7 A | 92% |
| BRAF V600E (6P3D) | Melanoma | -11.20 kcal/mol | 3.1 A | 89% |

The pipeline independently finds molecules scoring in the same range as known approved drugs, without any knowledge of those drugs. For the full experiment log with per-test breakdowns (fpocket vs P2Rank, TargetDiff POC), see [testing-guide.md](testing-guide.md).

## Where We Are Going

The roadmap is four phases. Each is a prerequisite for the next; we don't skip ahead. See [`autonomous_drug_discovery/plan.md`](../autonomous_drug_discovery/plan.md) for the detailed technical decomposition under each phase.

### Phase 1 — Containerise + cloud GPU + 20 quality results

Three intertwined deliverables: (a) ship the whole stack as a Docker image (envs + weights + binaries baked in, pushed to GHCR by CI); (b) wrap RunPod GPU rental in `make cloud-run` so heavy generation runs off-laptop; (c) generate, screen, dock, and rank ~30 TargetDiff molecules across five oncology kinases (1M17, 2HYY, 6P3D, KRAS G12C, JAK2), surfacing the top 20 in a GitHub-Pages-hosted dashboard. **No agent yet. No Sonnet. Just the deterministic five-stage pipeline running cleanly at scale, reproducibly, on either contributor's laptop.** TargetDiff and Pocket2Mol weights are mirrored to Hugging Face Hub since the upstream Google Drive folders are dead. Pocket2Mol stays deferred until its checkpoint is recovered. See [`autonomous_drug_discovery/plan.md`](../autonomous_drug_discovery/plan.md) for the full architecture and 7-day action plan.

### Phase 2 — Show the professor (M3, non-negotiable)

Walk a medicinal chemistry advisor through the dashboard for the Phase 1 targets. Capture per-candidate annotations (`promising` / `borderline` / `reject`) plus free-text reasoning. Calibrate filter thresholds. Get explicit feedback on the proposed adaptive layer (cascade + Bayesian + Sonnet + Obsidian) so we shape Phase 3 around what they need.

**Without medicinal chemistry judgment, the pipeline produces numbers without meaning.** Imperial Chemistry or Cancer Research collaborators should be engaged once Phase 1 has data to show.

### Phase 3 — Sonnet-in-the-loop + Bayesian recommender + Obsidian

Wire the four-part Adaptive Layer described above into the orchestrator. Three sub-deliverables, in order:

- **Obsidian campaign emitter** — every `orchestrator.py run` writes a folder of cross-linked MD notes (per-mol scorecards, route trees, embedded structure SVGs, tags). Independent of the agent; ships first because it's the highest-leverage / lowest-risk piece. The accumulating knowledge graph is the durable moat.
- **Bayesian strategy recommender** — Thompson sampling over (backend, num_samples, beam_size, screening thresholds) conditional on pocket descriptors. Cold-started from Phase 1 telemetry; updates after each new campaign. Honest framing of "the system learns": no fine-tuning, no RL on noisy proxies, just a principled exploration/exploitation trade-off over a small set of well-defined choices.
- **Sonnet agent shell** — wraps the existing `agent_planner.py` with retrieval over telemetry + open-source corpora (CrossDocked2020, BindingDB, ChEMBL), structured recommendations citing evidence, an attrition-funnel watcher for mid-campaign sanity checks, and a chemist-facing report writer. Every recommendation is logged, never auto-applied.

### Phase 4 — Bayesian evaluation; investor / university narrative

Quantify whether the Phase 3 agent loop produces measurably better candidates than the Phase 1 deterministic baseline. Report posteriors with credible intervals — never point estimates — on every lift metric: composite-score uplift, recall against expert-`promising`-annotated molecules, calibration of the recommender's posterior. Honest, even if the lift is null.

If the lift is real and credible: this is the defensible case for non-dilutive grants or seed funding to scale targets, agent quality, and (eventually) wet-lab validation partnerships. If the lift is null: that result is also publishable and tells us where to invest next (better retrieval? better reward signal? more wet-lab grounding?). Reporting credible intervals — not point estimates — is the antidote to the BenevolentAI failure mode of investor-friendly numbers that don't survive peer review.

### Recently completed (snapshot, 2026-05-11)

- Wired TargetDiff and Pocket2Mol into the orchestrator (`--mode targetdiff` / `--mode pocket2mol`)
- Per-campaign output directories (prevent file collisions)
- AiZynthFinder retrosynthetic feasibility integrated and used in the multi-criteria ranker (Stage 5)
- Stage 5 (multi-criteria ranker) shipped — composite score = 0.5·docking + 0.3·ADMET + 0.2·synthesis
- Multi-backend chemist dashboard (`dashboard/index.html`) with backend toggle, ADMET badges, and AiZynth route-tree viewer
- Pocket2Mol patched (`models/maskfill.py`) so CPU mode works end-to-end on Blackwell hardware where the pinned PyTorch can't run on GPU
- Fixed a long-standing `PYTHONPATH` bug in the TargetDiff wrapper that prevented end-to-end runs from a fresh shell (the `utils.misc` import failed unless `cwd` and `PYTHONPATH` both pointed at the TargetDiff repo)
- Confirmed the TargetDiff and Pocket2Mol Google Drive checkpoint folders are permanently dead; Phase 1 mirrors them to Hugging Face Hub from the dev-box copy (only TargetDiff weights survived; Pocket2Mol weights need to be sourced from a collaborator)

### Medium-term (post-Phase 4)

- AlphaFold / ESMFold integration for targets without crystal structures
- Multi-pocket docking (top 3 pockets, not just best)
- GNINA CNN-based rescoring alongside Vina
- Modal / Replicate migration so generation runs as a serverless function instead of a subprocess (Docker image stays the artefact; only the launcher changes)
- Publishable methods paper (Journal of Cheminformatics or JCIM) — most defensible after Phase 4 has Bayesian evaluation data

## Strategic Positioning

**What we are:** An open, modular, well-engineered orchestration layer that makes state-of-the-art computational drug discovery accessible to researchers who lack the software engineering capacity to build pipelines themselves.

**Who we serve:** Academic labs, small biotechs, and chemistry groups at institutions who have domain expertise but not the engineering resources to chain P2Rank, TargetDiff, RDKit, Vina, ADMET-AI, and AiZynthFinder into a reproducible, documented workflow.

**How we differentiate:** Not by algorithmic novelty (the individual tools are freely available), but by:
1. **Integration quality, reproducibility, transparency, and honest communication of limitations.** The COVID Moonshot project proved that open-source, community-driven drug discovery can produce real candidates. The hunger for accessible, integrated tools is genuine.
2. **A compounding knowledge graph.** Every campaign produces structured, cross-linked Obsidian notes that accumulate per-target SAR, expert annotations, and what-worked/what-didn't. Over time this is the durable moat — much more than the orchestration layer itself.
3. **Bayesian rigour over computational bravado.** The adaptive layer is a Thompson-sampling recommender, not a black-box RL agent; reported lifts come with credible intervals, not point estimates. This signals the kind of scientific honesty that a sceptical reviewer or grant committee actually trusts.

**Business model path:** Freemium open-source — release the pipeline openly to build adoption and credibility, then monetize through cloud-hosted premium features, consulting/customization, or partnership deals based on demonstrated utility.

## What This Is Not

**It is not a drug.** Nothing this system produces is a therapeutic. It produces hypotheses for experimental validation. The wet lab bottleneck is real: synthesising a single compound costs $500-5,000+ and takes weeks.

**It is not a replacement for domain expertise.** The system's output requires evaluation by someone who understands medicinal chemistry, structural biology, and the specific disease context. The system accelerates their work — it does not replace their judgment.

**It is not competing with Schrodinger or Recursion.** A 2-3 person team cannot compete on pipeline depth, proprietary data, or clinical validation. We compete on accessibility, transparency, and community adoption in an underserved market segment.

## Honest Risk Assessment

**The scoring function gap undermines all downstream claims.** If docking scores correlate with experimental binding at only r = 0.4-0.6, then ranking candidates by docking score is ranking by a noisy proxy. Every result we present must acknowledge this limitation.

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
6. **Posteriors, not point estimates.** Every empirical claim downstream of a campaign — "the agent improves rank quality by X" — must be reported with a credible interval, the prior used, and the sample size. Single-number results invite the BenevolentAI failure mode (investor-friendly headline, doesn't survive peer review). Bayesian framing is both more honest and more durable.
