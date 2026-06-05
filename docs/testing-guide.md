# How This Pipeline Works and How to Validate It

A plain-language guide to the science behind the pipeline and every validation experiment run to date. For operational usage (commands, parameters, modes), see [pipeline-guide.md](pipeline-guide.md). For querying results, see [telemetry-guide.md](telemetry-guide.md).

## The Big Picture

Imagine you're trying to find a key that fits a specific lock. The "lock" is a protein in the human body that causes disease. The "key" is a small molecule (a potential drug) that fits into a pocket on that protein and blocks it from working.

This pipeline automates the process of:
1. Finding the keyhole (pocket detection)
2. Making candidate keys (molecule generation)
3. Throwing away obviously bad keys (screening)
4. Testing which keys fit (docking)

It does all of this computationally — no lab, no chemicals, no test tubes. The output is a ranked list of "keys worth trying in a real lab."

---

## The Stages, Explained

### Stage 1: Pocket Detection (P2Rank or fpocket)

**What it does:** Takes the 3D shape of a disease protein and finds indentations ("pockets") on its surface where a drug molecule could physically sit.

**Real-world analogy:** Imagine a golf ball. The dimples are like pockets. Some dimples are deeper and more sheltered — those are better drug targets. This stage finds all the dimples and ranks them by how "druggable" they are.

**Two backends available:**
- **P2Rank (default, recommended):** ML-based pocket detection. Outperforms fpocket by 10-20 percentage points on standard benchmarks. Outputs a probability score (0-1) and pre-computed pocket center coordinates used directly for docking.
- **fpocket (fallback):** Geometry-based pocket detection. Uses Druggability Score (0-1) for ranking.

**What to look for:**
- **Pocket probability/Druggability Score** (0 to 1): How likely this pocket is to bind a drug molecule. Above 0.5 is promising. Above 0.8 is excellent.
- **Pocket distance to known drug site**: For validated targets, P2Rank places the pocket within 2.7-3.1 Angstroms of the crystallographic drug position.
- The pipeline picks the top-ranked pocket automatically.

### Stage 2: Molecule Generation (cascade: RDKit + TargetDiff)

**What it does:** Creates candidate molecules for the pocket. The recommended `cascade` mode runs two generators and merges them: **RDKit** combines chemical building blocks (rings, chains, functional groups) like LEGO bricks — fast, CPU, and crucially the *makeable* one; **TargetDiff** is a diffusion model that grows a 3D molecule directly inside the pocket — novel and high-fidelity, but its molecules are usually hard to synthesise. (Pocket2Mol, a third backend, was dropped — it can't run on the newest NVIDIA GPUs.)

**Real-world analogy:** RDKit is a chef combining familiar ingredients into plausible dishes; TargetDiff is an avant-garde chef inventing wholly new dishes that look stunning but may use ingredients no shop sells. You want both at the table — then a hard rule (the synthesizability gate, next) that throws out anything you can't actually cook.

**What to look for:**
- How many molecules each backend generated (typically ~50 each in cascade)
- TargetDiff output is potent-looking but often unmakeable — that's expected, and the gate handles it
- The generation uses the pocket size to decide how big the molecules should be

### Stage 3: Screening (MolScore + ADMET-AI)

**What it does:** Filters out molecules that would be bad drugs — too big, too oily, too hard to make, or known to cause problems in the body. Then runs surviving molecules through ADMET-AI for 104-property toxicity and pharmacokinetic profiling.

**Backend:** MolScore is the primary backend for descriptor calculation and PAINS filtering. If MolScore is not installed, it falls back to hand-rolled RDKit filters that compute the same properties. ADMET-AI enrichment runs on survivors regardless of backend.

**Real-world analogy:** Like quality control at a factory. Before expensive testing, throw away anything that's obviously defective. A molecule that's too large won't be absorbed by the gut. A molecule that's too oily will stick to everything. A molecule that's toxic will hurt the patient. The ADMET-AI step is like running blood tests on a job candidate — checking liver toxicity, heart safety (hERG), cancer-causing potential (AMES), and whether the molecule can actually reach its target in the body.

**Key filters applied (RDKit):**
- **Molecular weight** (must be <500): Bigger molecules can't get into cells easily
- **LogP** (must be <5): Measures oiliness. Too oily = won't dissolve in blood
- **QED** (0-1, higher is better): A combined "drug-likeness" score. Above 0.5 is decent
- **SA Score** (1-10, lower is better): How hard it would be to actually make this molecule in a lab. Below 4 is practical
- **PAINS filters**: Catches molecules with chemical patterns known to give false positive results in experiments

**ADMET-AI properties (key ones from 104 total):**
- **hERG**: Heart toxicity risk (lower = safer)
- **AMES**: Mutagenicity/cancer-causing potential (lower = safer)
- **DILI**: Drug-induced liver injury risk (lower = safer)
- **CYP inhibition**: Whether it interferes with liver enzymes that metabolize other drugs
- **Caco2/HIA**: Can it be absorbed through the gut?
- **BBB**: Can it cross the blood-brain barrier?
- **Clearance**: How quickly the body removes it
- **LD50**: Lethal dose estimate (higher = safer)

**What to look for:**
- **Survival rate**: What percentage passes all filters. 40-60% is healthy. If it's very low, the generator is making junk; if very high, the filters are too loose. After the thresholds were tightened (lead-like window: MW 250-450, QED ≥ 0.5, SA ≤ 4.5), real runs land in a sensible band — e.g. RDKit ~46%, TargetDiff ~23% on 8c4v — replacing the earlier too-loose 73-98%.

### Stage 2.5: Synthesizability Gate (AiZynthFinder)

**What it does:** Before spending docking effort, the gate asks a blunt question of each molecule: *can a chemist actually make this, from chemicals you can buy?* It runs a retrosynthetic tree search (AiZynthFinder, optionally pre-filtered by the fast RAScore model) and keeps only molecules with a complete route to purchasable building blocks.

**Real-world analogy:** A recipe that calls for "one unicorn egg" scores great on paper and is useless in the kitchen. The gate is the line cook who refuses to plate anything that can't be sourced.

**Why it exists:** When measured honestly, the 3D generators fail here badly — **0 of 28** top kinase candidates had any synthetic route, and TargetDiff scored **0% makeable** vs RDKit's **30%**. Without this gate the pipeline happily ranks beautiful, potent, *uncreatable* molecules. With it, synthesizability is enforced, not hoped for.

**What to look for:**
- **Makeable fraction**: how many survivors have a route. If a backend contributes ~0% (TargetDiff often does), that's the signal to lean on RDKit and treat the diffusion output as binding-mode inspiration only.
- **Route length**: 1-3 steps is genuinely tractable; 5-6 steps is a real project.

### Stage 4: Docking (AutoDock Vina)

**What it does:** Computationally "tries" each surviving molecule inside the pocket, calculating how strongly it would bind.

**Real-world analogy:** Like trying different keys in a lock and measuring how snugly each one fits. The software rotates and positions each molecule millions of times to find the best fit.

**What to look for:**
- **Docking score** (in kcal/mol, more negative = stronger binding):
  - **Better than -9**: Excellent — very strong predicted binding
  - **-7 to -9**: Strong — worth investigating further
  - **-5 to -7**: Moderate — might work but not a top candidate
  - **Worse than -5**: Weak — probably won't bind effectively
- Real approved drugs typically score between -7 and -11 kcal/mol

### Stage 5 + 5.5: Ranking and the Pharmacophore Bridge

**What it does:** Stage 5 combines docking, ADMET, and synthesizability into one composite score. In `cascade` mode a **pharmacophore bridge** (Stage 5.5) then does something clever with TargetDiff's output: instead of throwing away its (usually unmakeable) molecules, it treats their docked poses as a *map of how to bind the pocket* — which groups should point where — and re-ranks the **makeable** candidates by how well they reproduce that map.

**Real-world analogy:** TargetDiff is an architect who sketches a beautiful building out of materials that don't exist. Rather than bin the sketch, we keep the *floor plan* (where the load-bearing walls go) and ask: which buildable design best matches it? That match is rewarded in the final ranking.

**What to look for:**
- **`pharmacophore_match`** (0-1) on each ranked candidate: how faithfully a makeable molecule reproduces TargetDiff's binding mode. High match + makeable + good dock = the molecules you actually want.
- If *no* makeable molecule scores a high match, that's an honest finding: TargetDiff's binding mode isn't reachable with synthesizable chemistry on this target.

---

## How to Validate: "Is This Pipeline Working?"

The key question is: **if we run this on a protein where we already know the answer, does it find molecules similar to known drugs?**

### What we're testing (M2 Validation)

We picked three well-studied protein targets where drugs already exist:

| Target | Disease | Known Drug | What good looks like |
|--------|---------|------------|---------------------|
| EGFR (1M17) | Lung cancer | Erlotinib | Docking scores around -7 to -9.5 |
| BCR-ABL (2HYY) | Leukemia | Imatinib | Docking scores around -7.5 to -10 |
| BRAF V600E (6P3D) | Melanoma | Ponatinib | Docking scores around -8 to -10.5 |

### What "passing" looks like

For each target, we want to see:

1. **Pocket detection finds the right pocket.** The known drug binding site should be the top-ranked pocket (or close to it). We can check this by comparing the detected pocket's location to the known binding site from X-ray crystallography.

2. **Generated molecules dock reasonably.** We're NOT expecting the pipeline to re-discover erlotinib. That would be extraordinary. We're looking for molecules that score in the same ballpark — say within 2-3 kcal/mol of the known drug's score. This means the pipeline can independently find chemistry that "fits" the target.

3. **Screening is calibrated.** 40-60% survival rate means the generator and screener are balanced. If almost everything passes or almost everything fails, something is misconfigured.

4. **Results are consistent across targets.** If the pipeline works on EGFR but completely fails on BCR-ABL, that's a red flag. Different targets will have different scores, but the pipeline should produce reasonable results across all three.

### What "failing" looks like

- Docking scores are all worse than -5 kcal/mol (weak binding everywhere)
- The pocket detector picks a pocket far from the known drug binding site
- Zero molecules survive screening
- Only 1-2 molecules dock successfully (out of 50+)

---

## Validation Tests Run

This section documents every validation experiment that has been completed, with results.

### Test 1: Full Pipeline with fpocket (2026-03-21)

**What:** End-to-end pipeline (fpocket pocket detection -> RDKit generation -> RDKit screening -> Vina docking) on 3 cancer targets with known drugs.

**Setup:** fpocket backend, 100 molecules generated per target, exhaustiveness=8 for Vina.

| Target | PDB | Disease | Known Drug | Molecules Docked | Best Score | Avg Score |
|--------|-----|---------|------------|-----------------|------------|-----------|
| BRAF V600E | 6P3D | Melanoma | Ponatinib (-10.4) | 39 | -10.40 kcal/mol | -8.76 |
| BCR-ABL | 2HYY | Leukemia | Imatinib (-9.1) | 39 | -12.56 kcal/mol | -9.63 |
| EGFR | 1M17 | Lung cancer | Erlotinib (-8.5) | 98 | -8.18 kcal/mol | -5.70 |

**Verdict:** Pipeline generates molecules scoring in the same range as approved drugs. BCR-ABL results are particularly strong. EGFR had lower survival through screening (fewer molecules docked) but still produced competitive scores.

### Test 2: Full Pipeline with P2Rank + ADMET-AI (2026-03-21)

**What:** Same 3 targets, but using P2Rank (ML-based) for pocket detection and ADMET-AI for 104-property ADMET profiling after screening.

**Setup:** P2Rank backend, 100 molecules generated per target, ADMET-AI annotation on survivors, Vina docking with exhaustiveness=8.

| Target | PDB | Molecules Docked | Best Score | Avg Score | Pocket Distance | Residue Overlap |
|--------|-----|-----------------|------------|-----------|-----------------|-----------------|
| EGFR | 1M17 | 93 | -9.32 kcal/mol | -6.58 | 2.7 A | 82% |
| BRAF V600E | 6P3D | 84 | -11.20 kcal/mol | -8.39 | 3.1 A | 89% |
| BCR-ABL | 2HYY | 84 | -12.59 kcal/mol | -9.25 | 2.7 A | 92% |

**Verdict:** P2Rank improved results across all targets. EGFR best score improved from -8.18 to -9.32 (14% better). BRAF improved from -10.40 to -11.20. More molecules survived to docking (84-93 vs 39). Pocket placement within 3.1 A of crystallographic drug position for all targets.

### Test 3: Crystallographic Validation — P2Rank vs fpocket (2026-03-21)

**What:** Compared detected pocket locations against X-ray crystallography data (the "ground truth" of where drugs actually bind).

**Method:** Measured distance from detected pocket center to the centroid of the co-crystallized drug in each PDB structure, and computed residue overlap between detected pocket and known binding site residues.

| Target | P2Rank Distance | fpocket Distance | P2Rank Residue Overlap | fpocket Residue Overlap |
|--------|----------------|-----------------|----------------------|----------------------|
| EGFR (1M17) | 2.7 A | 6.1 A | 82% | 53% |
| BCR-ABL (2HYY) | 2.7 A | 2.7 A | 92% | ~90% |
| BRAF V600E (6P3D) | 3.1 A | 3.1 A | 89% | ~85% |

**Verdict:** P2Rank places the docking box closer to the real drug binding site, especially for EGFR where it is 3.4 A closer than fpocket. Residue overlap is consistently higher with P2Rank. Both methods perform well on BCR-ABL and BRAF (large, well-defined pockets).

### Test 4: TargetDiff Diffusion Generation (2026-03-21 to 2026-04-02; subsequent retries failed)

**What:** Generated molecules from noise using the TargetDiff E(3)-equivariant diffusion model, conditioned on the BRAF V600E (6P3D) binding pocket shape. Run standalone (not through the orchestrator pipeline).

**Status (2026-05-11):** This result is *not* currently reproducible through the orchestrator on a clean machine. The upstream Google Drive checkpoint folder is gone; a separate `PYTHONPATH` bug in the wrapper (since fixed) had prevented orchestrator-driven runs from succeeding; and CPU-only inference is too slow to be practical. Phase 1 of [`autonomous_drug_discovery/plan.md`](../autonomous_drug_discovery/plan.md) mirrors the surviving local checkpoint to Hugging Face and moves TargetDiff runs to a RunPod GPU.

**Setup:** TargetDiff checkpoint (pretrained), 1000 denoising steps, CPU inference (~12 min/molecule). Two separate runs producing one molecule each.

| Molecule | SMILES | MW | QED | Dock Score | Ligand Eff. | Tanimoto to Ponatinib |
|----------|--------|----|-----|-----------|-------------|----------------------|
| Mol 1 | `C1=CN=CC=C(c2cccc(Nc3ccncc3)c2)C1` | 261 | 0.899 | -7.59 kcal/mol | 0.38 | 0.186 |
| Mol 2 | `COc1cnc(C(=O)NCc2cccc(C)n2)cn1` | 258 | 0.889 | -7.38 kcal/mol | 0.39 | 0.154 |

**Interpretation:**
- Both molecules are drug-like (QED >0.88), synthesisable, and dock competitively to BRAF.
- Docking scores (-7.4 to -7.6) are moderate — weaker than the pipeline's best RDKit-generated molecules (-11.2) but still in the "worth investigating" range.
- Low Tanimoto similarity to Ponatinib (0.15-0.19) means these are structurally novel — the model did not copy known drugs.
- Ligand efficiency (0.38-0.39) is excellent — compact molecules with good binding per atom.
- 50% reconstruction failure rate (1 of 2 molecules failed in a 2-molecule batch) is expected for diffusion models.
- Visualizations available in `reports/` (mol1_vs_ponatinib.png, mol2_vs_ponatinib.png, targetdiff_all_vs_ponatinib.png).

### Test 5: Docking Protocol Validation — imatinib redock into 1IEP (2026-05-30)

**What:** Retrospective ("self-docking") validation of the docking protocol. Took the ABL kinase crystal structure 1IEP, removed its co-crystallised imatinib (Gleevec), and asked Vina to dock it back from scratch. Compares the docked pose against the true crystal pose by RMSD.

**Result:** Pose **RMSD 0.93 Å** (threshold for "correct" is < 2.0 Å) — PASS, and tighter than typical. Affinity −12.17 kcal/mol.

**Why it matters:** This calibrates the instrument. It means docking scores *on this target* are trustworthy — the 1IEP generated candidates' numbers sit on a validated protocol, unlike targets with no co-crystal ligand (e.g. 8c4v, where docking is unvalidatable and every score carries an asterisk). Surfaced and fixed two real bugs in the validation path: a multi-copy-ligand extraction error (1IEP has imatinib in two chains) and a bad RDKit import.

### Test 6: Synthesizability Audit — the finding that reshaped the pipeline (2026-05-30)

**What:** Ran a real AiZynthFinder retrosynthetic search (full USPTO templates + 17M-compound ZINC stock) over the top candidates from 10 kinase targets, and over all screened survivors of the three 8c4v generators.

**Result:**
- **0 of 28** top kinase candidates had a complete synthetic route — the 3D generators produce molecules that score well but can't be made.
- Makeability by backend (8c4v): **RDKit 30%** (1–3 step routes), **Pocket2Mol 18%** (4-step), **TargetDiff 0%**.
- The strongest binders were also the least makeable (flat fused-aromatic intercalators with high mutagenicity flags).

**Why it matters:** This is *the* finding behind Pipeline v2. It promoted synthesizability from an after-the-fact score to an enforced **Stage 2.5 gate**, motivated the **cascade** generator mix (RDKit for makeable matter, TargetDiff for binding-mode novelty), and justified **dropping Pocket2Mol**. Honest negative results like this — measured, not assumed — are exactly what the project's "every claim verifiable" principle is for.

### Summary of All Campaigns in Telemetry

Total campaigns: 7 (6 successful end-to-end, 1 failed at ingestion)

| Campaign | Backend | Target | Stages Completed | Date |
|----------|---------|--------|-----------------|------|
| campaign_955423eb | fpocket | BRAF V600E (6P3D) | All 4 | 2026-03-21 |
| campaign_fd4fad48 | fpocket | BCR-ABL (2HYY) | All 4 | 2026-03-21 |
| campaign_3a003826 | fpocket | EGFR (1M17) | All 4 | 2026-03-21 |
| campaign_92df6fd3 | P2Rank | EGFR (1M17) | Ingestion only (failed) | 2026-03-21 |
| campaign_c0f9df2f | P2Rank | EGFR (1M17) | All 4 | 2026-03-21 |
| campaign_ea9d4f1c | P2Rank | BRAF V600E (6P3D) | All 4 | 2026-03-21 |
| campaign_c6122296 | P2Rank | BCR-ABL (2HYY) | All 4 | 2026-03-21 |

TargetDiff experiments (standalone, not in telemetry): 2 molecules generated and docked for BRAF V600E.

---

## Common-Sense Sanity Checks

Here are checks anyone can do, no chemistry background required:

### 1. "Did anything actually run?"
Check the telemetry database for successful stages. See [telemetry-guide.md](telemetry-guide.md) for query examples. You should see multiple "success" entries. "failed" entries are fine — they happen when runs are interrupted.

### 2. "Are the numbers physically meaningful?"
- Docking scores should be negative (positive scores = molecule is repelled, not attracted)
- Molecular weights should be 100-600 (smaller = too simple, bigger = can't enter cells)
- QED should be 0-1 (anything else is a bug)
- SA scores should be 1-10

### 3. "Is the pipeline consistent?"
Run the same target twice. The molecule generation uses a fixed random seed, so results should be identical. If they differ, something is non-deterministic (not necessarily bad, but worth knowing).

### 4. "How does it compare to known drugs?"
Look up the known drug for each target (Erlotinib, Imatinib, Ponatinib) and compare: are our best docking scores in the same range? Are our molecules of similar size (molecular weight)? The benchmark script does this comparison automatically — see [pipeline-guide.md](pipeline-guide.md).

### 5. "Is the pocket the right one?"
For EGFR (1M17), the known erlotinib binding site is near residues around the ATP-binding cleft. P2Rank and fpocket both output which amino acid residues are in each pocket. A domain expert can verify whether these match the known binding site. For a non-expert: the pocket with the highest score should be picked. P2Rank's pocket center should be within ~3 Angstroms of the known drug binding site for well-characterized targets.

---

## Glossary

| Term | Plain English |
|------|--------------|
| **PDB file** | A 3D map of a protein's shape, from X-ray experiments |
| **Pocket** | An indentation on a protein where a drug can sit |
| **Druggability score** | How likely a pocket is to bind a drug (0-1) |
| **SMILES** | A text code that describes a molecule's structure (like a chemical formula but more detailed) |
| **Docking** | Computationally fitting a molecule into a pocket |
| **kcal/mol** | Unit for binding energy. More negative = stronger binding |
| **QED** | Drug-likeness score (0-1, higher = more drug-like) |
| **SA Score** | How hard it is to make a molecule in a lab (1-10, lower = easier) |
| **LogP** | How oily/water-repelling a molecule is |
| **PAINS** | Chemical patterns known to give false positives in experiments |
| **Campaign** | One complete run of the pipeline on one target |
| **SDF file** | A file containing 3D structures of multiple molecules |
| **PDBQT** | A modified PDB format with charge info, needed by the docking software |
| **Vina** | AutoDock Vina — the docking software that predicts binding strength |
| **RDKit** | An open-source chemistry toolkit for molecular calculations |
| **fpocket** | Geometry-based software that finds pockets on protein surfaces |
| **P2Rank** | ML-based pocket detection tool — more accurate than fpocket, now the default |
| **ADMET-AI** | AI model predicting 104 drug safety/pharmacokinetic properties |
| **ADMET** | Absorption, Distribution, Metabolism, Excretion, Toxicity — what happens to a drug in the body |
| **hERG** | Heart toxicity test — drugs that block this ion channel can cause fatal heart arrhythmias |
| **TargetDiff** | Diffusion model that generates molecules from noise, conditioned on a 3D pocket shape |
| **Ligand efficiency** | Binding strength per atom — better for comparing molecules of different sizes |
| **Tanimoto similarity** | How structurally similar two molecules are (0=nothing in common, 1=identical) |
| **Angstrom (A)** | Unit of distance = 0.1 nanometers. Atoms are ~1-2 A across |
| **Residue overlap** | What percentage of known binding site amino acids the detected pocket captures |
