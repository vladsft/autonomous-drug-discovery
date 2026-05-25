#!/usr/bin/env python
"""Build-time smoke test for the production docking stack (Stage 4).

Importing vina + meeko is necessary but NOT sufficient: a Meeko/RDKit version
mismatch imports cleanly and only fails when `prepare()` is actually called
(e.g. Meeko 0.6+ calls `mol.HasQuery()`, which this env's RDKit lacks). That
failure used to surface only on a GPU pod, an hour into a campaign. So here we
exercise the real path the docking module uses — receptor-free, CPU-only,
seconds — and fail the Docker build if any of it breaks.

Mirrors autonomous_drug_discovery/modules/04_docking/run_docking.py:
    MoleculePreparation().prepare(mol) -> PDBQTWriterLegacy.write_string(setup)
plus a Vina() instantiation. Keep this in lock-step with run_docking.py.
"""
import sys

import rdkit
import meeko
from rdkit import Chem
from rdkit.Chem import AllChem
from vina import Vina
from meeko import MoleculePreparation, PDBQTWriterLegacy

print(f"rdkit {rdkit.__version__} | meeko {meeko.__version__}")

# Instantiating Vina loads the compiled Boost/C++ extension — the layer that
# the libstdc++ CXXABI mismatch broke. A bare `import vina` doesn't always.
Vina(sf_name="vina")

# Prepare a real ligand to PDBQT exactly as run_docking._prepare_ligand_pdbqt does.
mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)
prep = MoleculePreparation()
prep.prepare(mol)
pdbqt = PDBQTWriterLegacy.write_string(prep.setup)[0]
if "ATOM" not in pdbqt and "HETATM" not in pdbqt:
    print("FAIL: ligand PDBQT has no atom records", file=sys.stderr)
    sys.exit(1)

print("docking stack OK in base: vina instantiated, ligand prepared to PDBQT")
