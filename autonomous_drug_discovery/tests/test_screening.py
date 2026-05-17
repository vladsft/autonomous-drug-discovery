"""
Unit tests for the Screening Module (Module 03: Fast Triage).

Run: python -m pytest tests/test_screening.py -v

Requires RDKit to be installed.
"""

import os
import sys
import json
import tempfile
import pytest
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.resolve()))

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

pytestmark = pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed")

# Import after path setup — directory starts with digit, use importlib
import importlib.util

screening_module = None
if HAS_RDKIT:
    _spec = importlib.util.spec_from_file_location(
        "run_screening",
        str(Path(__file__).parent.parent / "modules" / "03_screening" / "run_screening.py"),
    )
    if _spec and _spec.loader:
        screening_module = importlib.util.module_from_spec(_spec)
        _spec.loader.exec_module(screening_module)


# ---------------------------------------------------------------------------
# Test helpers: Create SDF files from SMILES
# ---------------------------------------------------------------------------

def _smiles_to_sdf(smiles_list, output_path):
    """Write a list of SMILES to an SDF file."""
    writer = Chem.SDWriter(str(output_path))
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            mol.SetProp("molecule_id", f"mol_{i:04d}")
            writer.write(mol)
    writer.close()


# ---------------------------------------------------------------------------
# Known test molecules
# ---------------------------------------------------------------------------

# Aspirin — a small approved drug (MW 180). Note: it is intentionally *below*
# the lead-like MW/heavy-atom floor in the default config, so it no longer
# passes the default triage (the floor exists to reject fragment-sized output).
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"

# Caffeine (small, MW 194)
CAFFEINE = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"

# Ibuprofen (small, MW 206)
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"

# Lead-like approved drugs — sized within the default config's 250-450 MW /
# 18-35 heavy-atom window, so they should survive the default triage.
WARFARIN = "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"
DIAZEPAM = "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"
PROPRANOLOL = "CC(C)NCC(O)COc1cccc2ccccc12"

# Very large molecule (should fail Lipinski MW > 500)
LARGE_MOL = "CC(=O)Oc1ccc(cc1)C(=O)NCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"

# Highly lipophilic (should fail LogP > 5)
LIPOPHILIC = "CCCCCCCCCCCCCCCCCCCC"  # Eicosane


@pytest.fixture
def work_dir(tmp_path):
    return tmp_path


class TestScreeningBackend:
    """Test that a screening backend is available and selected."""

    def test_backend_selected(self):
        assert screening_module.BACKEND in ("molscore", "rdkit")


class TestMolScoreScreening:
    """Test the MolScore/RDKit screening functions directly."""

    def test_valid_molecule_gets_properties(self):
        """A valid drug molecule should have all expected properties computed."""
        mol = Chem.MolFromSmiles(ASPIRIN)
        mols = [mol]
        smiles_list = [ASPIRIN]
        config = {"filter_thresholds": {}}

        results = screening_module._screen_molecules(
            mols, smiles_list, config, screening_module.BACKEND)

        assert len(results) == 1
        r = results[0]
        assert r["passed"] is True
        props = r["properties"]
        assert "desc_MolWt" in props
        assert "desc_MolLogP" in props
        assert "desc_QED" in props
        assert props["desc_MolWt"] > 0
        assert 0 <= props["desc_QED"] <= 1

    def test_none_mol_rejected(self):
        """A None molecule should be rejected as invalid_chemistry."""
        mols = [None]
        smiles_list = [None]
        config = {"filter_thresholds": {}}

        results = screening_module._screen_molecules(
            mols, smiles_list, config, screening_module.BACKEND)

        assert len(results) == 1
        assert results[0]["passed"] is False
        assert results[0]["eliminated_by"] == "invalid_chemistry"

    def test_mw_filter_rejects_large_mol(self):
        """A molecule exceeding MW 500 should be filtered out."""
        mol = Chem.MolFromSmiles(LARGE_MOL)
        if mol is None:
            pytest.skip("LARGE_MOL SMILES did not parse")
        mw = Descriptors.MolWt(mol)
        if mw <= 500:
            pytest.skip(f"LARGE_MOL MW={mw} does not exceed 500")

        mols = [mol]
        smiles_list = [LARGE_MOL]
        config = {"filter_thresholds": {"desc_MolWt": {"max": 500}}}

        results = screening_module._screen_molecules(
            mols, smiles_list, config, screening_module.BACKEND)

        assert results[0]["passed"] is False
        assert "desc_MolWt" in results[0]["eliminated_by"]

    def test_logp_filter_rejects_lipophilic(self):
        """Eicosane (very high LogP) should be filtered out."""
        mol = Chem.MolFromSmiles(LIPOPHILIC)
        if mol is None:
            pytest.skip("LIPOPHILIC SMILES did not parse")

        mols = [mol]
        smiles_list = [LIPOPHILIC]
        config = {"filter_thresholds": {"desc_MolLogP": {"max": 5}}}

        results = screening_module._screen_molecules(
            mols, smiles_list, config, screening_module.BACKEND)

        assert results[0]["passed"] is False
        assert "desc_MolLogP" in results[0]["eliminated_by"]

    def test_lead_like_drug_passes_default_config(self):
        """A lead-sized approved drug should pass the default scoring config."""
        mol = Chem.MolFromSmiles(DIAZEPAM)
        mols = [mol]
        smiles_list = [DIAZEPAM]

        cfg_path = screening_module.DEFAULT_CONFIG
        with open(cfg_path) as f:
            config = json.load(f)

        results = screening_module._screen_molecules(
            mols, smiles_list, config, screening_module.BACKEND)

        assert results[0]["passed"] is True

    def test_fragment_rejected_by_size_floor(self):
        """A fragment-sized drug (aspirin) is rejected by the lead-like floor."""
        mol = Chem.MolFromSmiles(ASPIRIN)
        cfg_path = screening_module.DEFAULT_CONFIG
        with open(cfg_path) as f:
            config = json.load(f)

        results = screening_module._screen_molecules(
            [mol], [ASPIRIN], config, screening_module.BACKEND)

        assert results[0]["passed"] is False
        # The lower MW / heavy-atom bound is what should fire.
        assert any("MolWt" in v or "NumHeavyAtoms" in v
                   for v in results[0]["violations"])


class TestFullScreeningPipeline:
    """Integration tests for the full screening pipeline."""

    def test_known_drugs_pass(self, work_dir):
        """Lead-sized known drugs should survive all default filters."""
        sdf_path = work_dir / "input.sdf"
        _smiles_to_sdf([WARFARIN, DIAZEPAM, PROPRANOLOL], sdf_path)

        output_dir = work_dir / "output"
        total_in, total_passed, report_path = screening_module.run_screening(
            str(sdf_path), str(output_dir)
        )

        assert total_in == 3
        assert total_passed >= 2  # At least 2 of 3 lead-like drugs should pass

        # Verify report
        with open(report_path) as f:
            report = json.load(f)
        assert report["total_input"] == 3
        assert report["total_passed"] >= 2
        assert "attrition_summary" in report
        assert "survivors" in report

    def test_screened_sdf_produced(self, work_dir):
        """Verify the output SDF file is created."""
        sdf_path = work_dir / "input.sdf"
        _smiles_to_sdf([ASPIRIN], sdf_path)

        output_dir = work_dir / "output"
        screening_module.run_screening(str(sdf_path), str(output_dir))

        output_sdf = output_dir / "screened_molecules.sdf"
        assert output_sdf.exists()

    def test_run_metadata_produced(self, work_dir):
        """Verify run_metadata.json is created."""
        sdf_path = work_dir / "input.sdf"
        _smiles_to_sdf([ASPIRIN], sdf_path)

        output_dir = work_dir / "output"
        screening_module.run_screening(str(sdf_path), str(output_dir))

        metadata_path = output_dir / "run_metadata.json"
        assert metadata_path.exists()
        with open(metadata_path) as f:
            metadata = json.load(f)
        assert metadata["status"] == "success"
        assert metadata["module"] == "03_screening"

    def test_attrition_counts(self, work_dir):
        """Verify attrition counts are accurate."""
        sdf_path = work_dir / "input.sdf"
        _smiles_to_sdf([ASPIRIN, LIPOPHILIC], sdf_path)

        output_dir = work_dir / "output"
        total_in, total_passed, report_path = screening_module.run_screening(
            str(sdf_path), str(output_dir)
        )

        with open(report_path) as f:
            report = json.load(f)
        assert report["total_input"] == 2
        assert report["total_passed"] + report["total_rejected"] == report["total_input"]

    def test_telemetry_integration(self, work_dir):
        """Verify telemetry DB is populated when db_path provided."""
        sdf_path = work_dir / "input.sdf"
        _smiles_to_sdf([ASPIRIN, CAFFEINE], sdf_path)

        db_path = str(work_dir / "test.db")
        output_dir = work_dir / "output"
        screening_module.run_screening(
            str(sdf_path), str(output_dir),
            db_path=db_path, campaign_id="test_campaign"
        )

        from telemetry import TelemetryDB
        db = TelemetryDB(db_path)
        runs = db.query_runs(campaign_id="test_campaign")
        assert len(runs) == 1
        assert runs[0]["status"] == "success"
        assert runs[0]["module_name"] == "03_screening"

        molecules = db.query_molecules(runs[0]["run_id"])
        assert len(molecules) == 2
        db.close()

    def test_empty_sdf(self, work_dir):
        """Empty SDF should either produce zero output or exit gracefully."""
        sdf_path = work_dir / "empty.sdf"
        with open(sdf_path, "w") as f:
            f.write("")  # Empty file

        output_dir = work_dir / "output"
        # The screening module calls sys.exit(1) on invalid/empty files
        # since RDKit's SDMolSupplier treats an empty file as invalid
        with pytest.raises(SystemExit) as exc_info:
            screening_module.run_screening(str(sdf_path), str(output_dir))
        assert exc_info.value.code == 1

        # Verify failure metadata was still written
        metadata_path = output_dir / "run_metadata.json"
        assert metadata_path.exists()
        with open(metadata_path) as f:
            metadata = json.load(f)
        assert metadata["status"] == "failed"
