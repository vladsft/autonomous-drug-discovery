# How This Drug Discovery Pipeline Works — A Plain-Language Guide

## The Big Picture

Imagine you're trying to find a key that fits a specific lock. The "lock" is a protein in the human body that causes disease. The "key" is a small molecule (a potential drug) that fits into a pocket on that protein and blocks it from working.

This pipeline automates the process of:
1. Finding the keyhole (pocket detection)
2. Making candidate keys (molecule generation)
3. Throwing away obviously bad keys (screening)
4. Testing which keys fit (docking)

It does all of this computationally — no lab, no chemicals, no test tubes. The output is a ranked list of "keys worth trying in a real lab."

---

## The Four Stages, Explained

### Stage 1: Pocket Detection (fpocket)

**What it does:** Takes the 3D shape of a disease protein and finds indentations ("pockets") on its surface where a drug molecule could physically sit.

**Real-world analogy:** Imagine a golf ball. The dimples are like pockets. Some dimples are deeper and more sheltered — those are better drug targets. This stage finds all the dimples and ranks them by how "druggable" they are.

**What to look for:**
- **Druggability Score** (0 to 1): How likely this pocket is to bind a drug molecule. Above 0.5 is promising. Above 0.8 is excellent.
- **Volume**: How big the pocket is. Bigger pockets can accommodate larger, more complex molecules.
- The pipeline picks the pocket with the highest Druggability Score automatically.

### Stage 2: Molecule Generation (RDKit)

**What it does:** Creates ~100 random drug-like molecules from chemical building blocks (rings, chains, functional groups), combining them like LEGO bricks.

**Real-world analogy:** Like a chef combining ingredients (salt, sugar, spices, different proteins) into dishes. Most random combinations taste bad, but a trained chef knows which building blocks work well together. Our "chef" uses rules from medicinal chemistry.

**What to look for:**
- How many molecules were generated (typically 100)
- The generation uses the pocket size to decide how big the molecules should be

### Stage 3: Screening (RDKit)

**What it does:** Filters out molecules that would be bad drugs — too big, too oily, too hard to make, or known to cause problems in the body.

**Real-world analogy:** Like quality control at a factory. Before expensive testing, throw away anything that's obviously defective. A molecule that's too large won't be absorbed by the gut. A molecule that's too oily will stick to everything. A molecule that's toxic will hurt the patient.

**Key filters applied:**
- **Molecular weight** (must be <500): Bigger molecules can't get into cells easily
- **LogP** (must be <5): Measures oiliness. Too oily = won't dissolve in blood
- **QED** (0-1, higher is better): A combined "drug-likeness" score. Above 0.5 is decent
- **SA Score** (1-10, lower is better): How hard it would be to actually make this molecule in a lab. Below 4 is practical
- **PAINS filters**: Catches molecules with chemical patterns known to give false positive results in experiments

**What to look for:**
- **Survival rate**: What percentage passes all filters. 40-60% is healthy. If it's very low, the generator is making junk. If it's very high, the filters may be too loose.

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

## How to Run It Yourself

### Prerequisites
- The conda environment is already set up with all dependencies
- You need a `.pdb` file — a 3D structure of a protein (downloadable free from rcsb.org)

### Running the full pipeline on a target

```bash
cd /home/vladsft/agent-harness/autonomous_drug_discovery

# Run the full pipeline on EGFR (already downloaded)
conda run -n base python orchestrator.py run data/processed/1M17.pdb --mode production

# Run on BCR-ABL
conda run -n base python orchestrator.py run data/processed/2HYY.pdb --mode production

# Run on BRAF V600E
conda run -n base python orchestrator.py run data/processed/6P3D.pdb --mode production
```

Each run takes 5-20 minutes depending on the target (docking is the slow part).

### Running the benchmark comparison

After at least one pipeline run completes:

```bash
conda run -n base python benchmark.py
```

This prints a formatted report comparing all completed targets.

### Running individual stages

If you want to understand each stage separately:

```bash
# Stage 1: Find pockets on a protein
conda run -n base python orchestrator.py ingest data/processed/1M17.pdb

# Stage 2: Generate molecules (needs manifest from stage 1)
conda run -n base python orchestrator.py generate data/processed/1M17_manifest.json --mode production

# Stage 3: Screen generated molecules
conda run -n base python orchestrator.py screen data/candidates/generated_molecules.sdf

# Stage 4: Dock screened molecules (needs manifest from stage 1)
conda run -n base python orchestrator.py dock data/processed/1M17_manifest.json --mode production
```

### Trying a new protein target

1. Go to https://www.rcsb.org and search for your protein of interest
2. Download the `.pdb` file
3. Place it in `data/processed/`
4. Run: `conda run -n base python orchestrator.py run data/processed/YOUR_FILE.pdb --mode production`

### Looking at results

**Docking results** are in CSV format:
```bash
# View top docked molecules (sorted by score, most negative first)
cat data/results/docking_results.csv
```

**Screening report** shows which molecules passed and why others were rejected:
```bash
cat data/screened/screening_report.json | python -m json.tool | head -50
```

**Telemetry database** has everything in SQLite:
```bash
conda run -n base python -c "
import sqlite3
conn = sqlite3.connect('data/telemetry.db')

# See all campaigns
for row in conn.execute('SELECT campaign_id, module_name, status FROM runs ORDER BY started_at'):
    print(row)

# See best molecules across ALL runs
print('\nTop 10 molecules by docking score:')
for row in conn.execute('''
    SELECT smiles, docking_score, qed, sa_score
    FROM molecule_scores
    WHERE docking_score IS NOT NULL
    ORDER BY docking_score LIMIT 10
'''):
    print(f'  Score: {row[1]:.2f} kcal/mol | QED: {row[2]} | SA: {row[3]} | {row[0][:60]}')
"
```

---

## Common-Sense Sanity Checks

Here are checks anyone can do, no chemistry background required:

### 1. "Did anything actually run?"
```bash
conda run -n base python -c "
import sqlite3
conn = sqlite3.connect('data/telemetry.db')
for row in conn.execute('SELECT status, COUNT(*) FROM runs GROUP BY status'):
    print(f'{row[0]}: {row[1]} runs')
"
```
You should see multiple "success" entries. "failed" or "abandoned" entries are fine — they happen when runs are interrupted.

### 2. "Are the numbers physically meaningful?"
- Docking scores should be negative (positive scores = molecule is repelled, not attracted)
- Molecular weights should be 100-600 (smaller = too simple, bigger = can't enter cells)
- QED should be 0-1 (anything else is a bug)
- SA scores should be 1-10

### 3. "Is the pipeline consistent?"
Run the same target twice. The molecule generation uses a fixed random seed, so results should be identical. If they differ, something is non-deterministic (not necessarily bad, but worth knowing).

### 4. "How does it compare to known drugs?"
Look up the known drug for each target (Erlotinib, Imatinib, Ponatinib) and compare:
- Are our best docking scores in the same range?
- Are our molecules of similar size (molecular weight)?
- The benchmark script does this comparison automatically.

### 5. "Is the pocket the right one?"
For EGFR (1M17), the known erlotinib binding site is near residues around the ATP-binding cleft. The fpocket output lists which amino acid residues are in each pocket. A domain expert can verify whether these match the known binding site. For a non-expert: the pocket with the highest Druggability Score should be picked, and it should have a large volume (>500 cubic Angstroms for kinases).

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
| **fpocket** | Software that finds pockets on protein surfaces |
