# Agent Onboarding — Autonomous Drug Discovery Pipeline

## What this is
An end-to-end computational drug discovery pipeline. Input: a protein PDB file. Output: ranked candidate drug molecules. Four stages: pocket detection (fpocket) -> molecule generation (RDKit) -> drug-likeness screening (RDKit) -> molecular docking (AutoDock Vina). Everything logged to SQLite telemetry.

## Critical: Python environment
- **Always use:** `conda run -n base python <script>` — this has RDKit, Vina, Meeko, gemmi
- **Never use:** bare `python3` or `/usr/bin/python3` — system Python 3.12 lacks all chemistry deps
- Conda Python is at `/home/vladsft/miniconda3/bin/python` (3.13)
- There is also a `targetdiff_env` conda env — currently broken (PyTorch symbol error), needs reinstall

## Repository layout
```
autonomous_drug_discovery/
  orchestrator.py              # CLI entrypoint — runs full pipeline or individual stages
  telemetry.py                 # SQLite schema + logging (tables: runs, molecule_scores)
  benchmark.py                 # M2 validation report across targets
  plan.md                      # Architectural plan (Steps 1-7)
  agent_planner.py             # LLM-driven orchestration (needs API key, not yet used)

  modules/
    01_ingestion/run_pocket.py       # fpocket wrapper, pocket ranking by Druggability Score
    02_generation/run_generation.py  # Fragment-based (rdkit mode) or diffusion (targetdiff mode)
    02_generation/targetdiff/        # TargetDiff repo (cloned, but env broken + checkpoint missing)
    03_screening/run_screening.py    # RDKit filters: QED, SA, LogP, MW, PAINS
    03_screening/default_scoring_config.json  # Filter thresholds (editable)
    04_docking/run_docking.py        # AutoDock Vina via Python API + Meeko for ligand prep

  data/
    processed/        # Input PDBs + fpocket output (*_out/ dirs, *_manifest.json)
    candidates/       # Generated molecules (generated_molecules.sdf)
    screened/          # Filtered molecules (screened_molecules.sdf, screening_report.json)
    results/           # Docking output (docking_results.csv, docked_*.pdbqt)
    telemetry.db       # SQLite — source of truth for all campaign data

reports/
  north_star.md       # Project vision, milestones M1-M5
  testing_guide.md    # Plain-language guide for non-experts
```

## How to run
```bash
cd /home/vladsft/agent-harness/autonomous_drug_discovery

# Full pipeline on a target
conda run -n base python orchestrator.py run data/processed/1M17.pdb --mode production

# Benchmark comparison across all completed targets
conda run -n base python benchmark.py

# Individual stages
conda run -n base python orchestrator.py ingest <pdb_file>
conda run -n base python orchestrator.py generate <manifest.json> --mode production
conda run -n base python orchestrator.py screen <sdf_file>
conda run -n base python orchestrator.py dock <manifest.json> --mode production
```

## Milestone status (as of 2026-03-21)
- **M1 (Working Pipeline):** DONE
- **M2 (Validation):** DONE — all 3 targets validated
- **M3 (Expert Review):** NOT STARTED — need domain expert feedback
- **M4 (Novel Campaign):** NOT STARTED
- **M5 (Adaptive Planning):** NOT STARTED

## M2 validation results
| Target | PDB | Disease | Best Dock Score | Pocket-to-Drug Dist | Residue Overlap |
|--------|-----|---------|----------------|--------------------|-----------------|
| EGFR | 1M17 | Lung cancer | -8.18 kcal/mol | 6.3 A | 53% (9/17) |
| BCR-ABL | 2HYY | Leukemia | -12.56 kcal/mol | 2.7 A | 92% (23/25) |
| BRAF V600E | 6P3D | Melanoma | -10.40 kcal/mol | 3.1 A | 89% (25/28) |

BCR-ABL and BRAF pocket detection nearly perfectly matches crystallographic data. EGFR is weaker (shallower pocket, druggability 0.529 vs ~0.99).

## Key tools and binaries
- **fpocket:** `/home/vladsft/fpocket/bin/fpocket` — pocket detection
- **RDKit 2025.03.6:** molecule generation, screening, property calculation
- **AutoDock Vina:** docking (Python API via `from vina import Vina`)
- **Meeko:** ligand PDBQT preparation for Vina
- **gemmi:** receptor PDB parsing and atom typing for PDBQT conversion

## Known issues and limitations
1. **Concurrent runs collide** — all stages write to shared dirs with hardcoded filenames (generated_molecules.sdf, screened_molecules.sdf, docking_results.csv). Run targets sequentially, or fix orchestrator to use per-campaign output dirs.
2. **Screening filters too loose** — survival rates 73-98%, should be 40-60%. Tighten thresholds in default_scoring_config.json.
3. **EGFR pocket offset** — 6.3A from known binding site, 53% overlap. Pocket centroid calculation could be improved (weight by alpha sphere proximity to cavity interior).
4. **TargetDiff not functional** — env broken (PyTorch libtorch_cpu.so symbol error), pretrained weights missing from pretrained_models/. Needs: reinstall env, download checkpoint, test sample_for_pocket.py. Runs on CPU (slow, hours) or GPU (minutes).
5. **Docking box is fixed 20x20x20 A** — may need enlarging for bigger pockets or better centering.
6. **2HYY has 4 chains** — pocket1 maps to chain C's Imatinib site. Multi-chain PDBs need care.

## What to work on next
**Highest impact improvements (in order):**
1. Get TargetDiff diffusion generation working — molecules designed for pocket geometry will dock far better than random fragment assembly
2. Improve pocket centroid accuracy — use alpha sphere centers, not atom averages
3. Add per-campaign output directories to orchestrator — prevents file collisions
4. Tighten screening thresholds with domain expert input (M3)
5. Multi-pocket docking — try top 3 pockets instead of just best one

## Telemetry database
```sql
-- Completed campaigns with all stages
SELECT campaign_id, module_name, status FROM runs ORDER BY started_at;

-- Best molecules across all campaigns
SELECT smiles, docking_score, qed, sa_score FROM molecule_scores
WHERE docking_score IS NOT NULL ORDER BY docking_score LIMIT 20;
```

## Git / SSH
- Remote: `git@github.com:vladsft/autonomous-drug-discovery.git` (SSH configured)
- Branch: main (3 commits + uncommitted work from this session)
- Uncommitted: pocket ranking fix, collision fix, benchmark.py, testing_guide.md
