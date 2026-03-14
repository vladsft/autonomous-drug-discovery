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
    from rdkit.Chem import AllChem
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
        else:
            # Write an invalid entry by writing raw text
            pass
    writer.close()


def _smiles_to_sdf_with_invalid(smiles_list, output_path):
    """Write SMILES to SDF, including an intentionally malformed entry."""
    with open(output_path, "w") as f:
        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                AllChem.EmbedMolecule(mol, randomSeed=42)
                f.write(Chem.MolToMolBlock(mol))
                f.write(f">  <molecule_id>\nmol_{i:04d}\n\n")
                f.write("$$$$\n")

        # Add an invalid mol block
        f.write("INVALID_MOL\n")
        f.write("M  END\n")
        f.write("$$$$\n")


# ---------------------------------------------------------------------------
# Known test molecules
# ---------------------------------------------------------------------------

# Aspirin (should pass all filters)
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"

# Caffeine (should pass)
CAFFEINE = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"

# Ibuprofen (should pass)
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"

# Very large molecule (should fail Lipinski MW > 500)
LARGE_MOL = "CC(=O)Oc1ccc(cc1)C(=O)NCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"

# Highly lipophilic (should fail LogP > 5)
LIPOPHILIC = "CCCCCCCCCCCCCCCCCCCC"  # Eicosane


@pytest.fixture
def work_dir(tmp_path):
    return tmp_path


class TestIndividualFilters:
    """Test each filter function independently."""

    def test_valid_chemistry_pass(self):
        mol = Chem.MolFromSmiles(ASPIRIN)
        passed, reason = screening_module.filter_valid_chemistry(mol)
        assert passed is True
        assert reason is None

    def test_valid_chemistry_fail_none(self):
        passed, reason = screening_module.filter_valid_chemistry(None)
        assert passed is False
        assert reason == "invalid_chemistry"

    def test_lipinski_pass(self):
        mol = Chem.MolFromSmiles(ASPIRIN)
        passed, reason = screening_module.filter_lipinski(mol)
        assert passed is True

    def test_lipinski_fail_mw(self):
        # Create a molecule with MW > 500
        mol = Chem.MolFromSmiles(LARGE_MOL)
        if mol is not None:
            from rdkit.Chem import Descriptors
            mw = Descriptors.MolWt(mol)
            if mw > 500:
                passed, reason = screening_module.filter_lipinski(mol)
                assert passed is False
                assert reason == "lipinski_violation"

    def test_qed_pass(self):
        mol = Chem.MolFromSmiles(ASPIRIN)
        passed, reason = screening_module.filter_qed(mol)
        assert passed is True

    def test_sa_score_pass(self):
        mol = Chem.MolFromSmiles(ASPIRIN)
        passed, reason = screening_module.filter_sa_score(mol)
        assert passed is True


class TestMoleculeProperties:
    """Test property computation."""

    def test_compute_properties(self):
        mol = Chem.MolFromSmiles(ASPIRIN)
        props = screening_module.compute_molecule_properties(mol)
        assert "smiles" in props
        assert "mol_weight" in props
        assert "logp" in props
        assert "qed" in props
        assert "sa_score" in props
        assert props["mol_weight"] > 0
        assert 0 <= props["qed"] <= 1


class TestFullScreeningPipeline:
    """Integration tests for the full screening pipeline."""

    def test_known_drugs_pass(self, work_dir):
        """Known drug molecules should survive all filters."""
        sdf_path = work_dir / "input.sdf"
        _smiles_to_sdf([ASPIRIN, CAFFEINE, IBUPROFEN], sdf_path)

        output_dir = work_dir / "output"
        total_in, total_passed, report_path = screening_module.run_screening(
            str(sdf_path), str(output_dir)
        )

        assert total_in == 3
        assert total_passed >= 2  # At least 2 of 3 known drugs should pass

        # Verify report
        with open(report_path) as f:
            report = json.load(f)
        assert report["total_input"] == 3
        assert report["total_passed"] >= 2
        assert "attrition_by_filter" in report
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
        # Sum of attrition + passed should equal total input
        total_attrition = sum(report["attrition_by_filter"].values())
        assert total_attrition + report["total_passed"] == report["total_input"]

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
        """Empty SDF should produce zero output gracefully."""
        sdf_path = work_dir / "empty.sdf"
        with open(sdf_path, "w") as f:
            f.write("")  # Empty file

        output_dir = work_dir / "output"
        total_in, total_passed, report_path = screening_module.run_screening(
            str(sdf_path), str(output_dir)
        )
        assert total_in == 0
        assert total_passed == 0
