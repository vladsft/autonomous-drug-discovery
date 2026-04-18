# Telemetry Database Guide

Every pipeline run is logged to a SQLite database at `autonomous_drug_discovery/data/telemetry.db`. This document covers the schema, common queries, and how to use the data for analysis.

## Location

```
autonomous_drug_discovery/data/telemetry.db
```

The database is created automatically on the first pipeline run. You can specify a custom path with `--db_path` on any orchestrator command.

## Schema

### `runs` table

One row per module execution. Tracks what ran, when, whether it succeeded, and what parameters were used.

| Column | Type | Description |
|---|---|---|
| `run_id` | TEXT (PK) | UUID, unique per execution |
| `campaign_id` | TEXT | Groups all stages of a single pipeline run (e.g. `campaign_fd4fad48`) |
| `module_name` | TEXT | Which stage ran: `01_ingestion`, `02_generation`, `03_screening`, `04_docking`, `agent_planner` |
| `started_at` | TEXT | ISO 8601 UTC timestamp |
| `completed_at` | TEXT | ISO 8601 UTC timestamp, NULL if still running |
| `status` | TEXT | `running`, `success`, or `failed` |
| `input_hash` | TEXT | SHA-256 of the input file (for deduplication/cache) |
| `input_path` | TEXT | Absolute path to the primary input file |
| `output_path` | TEXT | Absolute path to the primary output file |
| `parameters` | TEXT | JSON blob of all parameters used for this run |
| `error_trace` | TEXT | Full Python traceback on failure, NULL on success |
| `git_commit` | TEXT | Git commit hash of the model repo (TargetDiff), NULL otherwise |
| `notes` | TEXT | Free-text annotations |

### `molecule_scores` table

One row per molecule per stage. Populated by screening (properties + pass/fail) and docking (binding affinity).

| Column | Type | Description |
|---|---|---|
| `score_id` | INTEGER (PK) | Auto-increment |
| `run_id` | TEXT (FK) | References `runs.run_id` |
| `molecule_id` | TEXT | Internal ID (e.g. `mol_0042`) |
| `smiles` | TEXT | Canonical SMILES string |
| `qed` | REAL | Quantitative Estimate of Drug-likeness (0-1, higher = more drug-like) |
| `sa_score` | REAL | Synthetic Accessibility score (1-10, lower = easier to synthesize) |
| `logp` | REAL | Octanol-water partition coefficient |
| `mol_weight` | REAL | Molecular weight in Daltons |
| `docking_score` | REAL | Binding affinity in kcal/mol (more negative = stronger binding) |
| `passed_triage` | INTEGER | 1 = passed screening, 0 = rejected |
| `stage_eliminated` | TEXT | Reason for rejection (e.g. `desc_MolLogP_exceeded`), NULL if passed |

### Indexes

```
idx_runs_campaign    ON runs(campaign_id)
idx_runs_module      ON runs(module_name)
idx_runs_status      ON runs(status)
idx_molscores_run    ON molecule_scores(run_id)
idx_molscores_triage ON molecule_scores(passed_triage)
```

## Common Queries

All examples use `sqlite3` from the `autonomous_drug_discovery/` directory.

### List all campaigns

```bash
sqlite3 -header -column data/telemetry.db "
  SELECT campaign_id, module_name, status, started_at
  FROM runs
  ORDER BY started_at;
"
```

### Campaign summary (stages completed, pass/fail)

```bash
sqlite3 -header -column data/telemetry.db "
  SELECT
    campaign_id,
    COUNT(*) as total_stages,
    SUM(CASE WHEN status='success' THEN 1 ELSE 0 END) as succeeded,
    SUM(CASE WHEN status='failed' THEN 1 ELSE 0 END) as failed,
    MIN(started_at) as started
  FROM runs
  GROUP BY campaign_id
  ORDER BY started;
"
```

### Identify the target protein for each campaign

```bash
sqlite3 -header -column data/telemetry.db "
  SELECT
    campaign_id,
    REPLACE(input_path, RTRIM(input_path, REPLACE(input_path, '/', '')), '') as target_file,
    status
  FROM runs
  WHERE module_name = '01_ingestion'
  ORDER BY started_at;
"
```

### Top molecules by docking score

```bash
sqlite3 -header -column data/telemetry.db "
  SELECT smiles, docking_score, qed, sa_score, mol_weight
  FROM molecule_scores
  WHERE docking_score IS NOT NULL
  ORDER BY docking_score
  LIMIT 20;
"
```

### Screening attrition breakdown

Shows why molecules were rejected, sorted by frequency:

```bash
sqlite3 -header -column data/telemetry.db "
  SELECT stage_eliminated, COUNT(*) as count
  FROM molecule_scores
  WHERE passed_triage = 0
  GROUP BY stage_eliminated
  ORDER BY count DESC;
"
```

### Survivors with full properties

```bash
sqlite3 -header -column data/telemetry.db "
  SELECT smiles, qed, sa_score, logp, mol_weight
  FROM molecule_scores
  WHERE passed_triage = 1
  ORDER BY qed DESC
  LIMIT 20;
"
```

### Pass/fail statistics across all campaigns

```bash
sqlite3 -header -column data/telemetry.db "
  SELECT
    passed_triage,
    COUNT(*) as count,
    ROUND(AVG(qed), 3) as avg_qed,
    ROUND(AVG(sa_score), 2) as avg_sa,
    ROUND(AVG(mol_weight), 1) as avg_mw,
    ROUND(AVG(logp), 2) as avg_logp
  FROM molecule_scores
  WHERE qed IS NOT NULL
  GROUP BY passed_triage;
"
```

### Molecules for a specific campaign

```bash
sqlite3 -header -column data/telemetry.db "
  SELECT ms.molecule_id, ms.smiles, ms.docking_score, ms.qed, ms.passed_triage
  FROM molecule_scores ms
  JOIN runs r ON ms.run_id = r.run_id
  WHERE r.campaign_id = 'campaign_c0f9df2f'
  ORDER BY ms.docking_score;
"
```

### Failed runs with error traces

```bash
sqlite3 -header -column data/telemetry.db "
  SELECT campaign_id, module_name, started_at, error_trace
  FROM runs
  WHERE status = 'failed';
"
```

### Parameters used for a run

The `parameters` column is a JSON blob. Use SQLite's `json_extract` to pull specific values:

```bash
sqlite3 -header -column data/telemetry.db "
  SELECT
    campaign_id,
    module_name,
    json_extract(parameters, '$.mode') as mode,
    json_extract(parameters, '$.backend') as backend
  FROM runs
  ORDER BY started_at;
"
```

## Benchmark Script

For a formatted cross-target comparison report, use the built-in benchmark script:

```bash
cd autonomous_drug_discovery
conda run -n base python benchmark.py
```

This queries the telemetry DB, finds the most recent successful campaign for each validation target, and prints:

- Per-target screening survival rates
- Docking score distributions (best, median, worst)
- Strong binder counts
- Comparison against known drug affinity ranges
- Drug-likeness statistics (QED, SA score)

You can point it at a different database:

```bash
conda run -n base python benchmark.py --db_path /path/to/other/telemetry.db
```

## Python API

You can also query telemetry programmatically:

```python
from telemetry import TelemetryDB

db = TelemetryDB("data/telemetry.db")

# List all runs for a campaign
runs = db.query_runs(campaign_id="campaign_c0f9df2f")
for r in runs:
    print(f"{r['module_name']}: {r['status']}")

# Get molecules from a screening run
screening_runs = db.query_runs(campaign_id="campaign_c0f9df2f", module_name="03_screening")
molecules = db.query_molecules(screening_runs[0]["run_id"])
for m in molecules:
    print(f"{m['smiles']}: QED={m['qed']}, passed={m['passed_triage']}")

# Campaign-level summary
summary = db.get_campaign_summary("campaign_c0f9df2f")
print(summary)

db.close()
```

## Interpreting Scores (Quick Reference)

Use this table when reading query results. For plain-language explanations of what each metric means and why it matters, see [testing-guide.md](testing-guide.md).

| Metric | Good Range | Warning | Notes |
|---|---|---|---|
| Docking score | < -7 kcal/mol | > -5 | More negative = stronger binding. Correlation with experimental affinity is only r=0.4-0.6 |
| QED | 0.6 - 1.0 | < 0.3 | Quantitative drug-likeness. Known drugs typically score 0.5-0.9 |
| SA score | 1 - 3 | > 6 | Synthetic accessibility. 1 = trivial to make, 10 = extremely difficult |
| LogP | 0 - 5 | > 5 | Lipinski Rule of Five cutoff. High LogP = poor solubility |
| Mol weight | 150 - 500 | > 500 | Lipinski cutoff at 500 Da |
| Survival rate | 40 - 60% | > 80% or < 20% | Too high = filters too loose, too low = generation producing junk |

## Database Maintenance

The database uses WAL (Write-Ahead Logging) mode for concurrent read/write access. No manual maintenance is needed for normal use.

To reset and start fresh:

```bash
rm autonomous_drug_discovery/data/telemetry.db
```

A new database will be created automatically on the next pipeline run.

To export data for external analysis:

```bash
# CSV export
sqlite3 -header -csv data/telemetry.db "SELECT * FROM molecule_scores;" > molecules.csv
sqlite3 -header -csv data/telemetry.db "SELECT * FROM runs;" > runs.csv
```
