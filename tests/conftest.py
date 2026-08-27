#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Shared fixtures: CJSON input built the way Avogadro builds it."""

import json
import pathlib
import subprocess
import sys

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))


def cjson_from_smiles(smiles, seed=0xF00D):
    """The CJSON Avogadro would hand a plugin for this molecule."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    assert AllChem.EmbedMolecule(mol, params) == 0, smiles

    # Avogadro's CJSON carries Kekulé bond orders, never an aromatic order,
    # so the fixture has to kekulize or every ring bond would round to 2.
    Chem.Kekulize(mol, clearAromaticFlags=True)

    conformer = mol.GetConformer()
    coords = []
    for i in range(mol.GetNumAtoms()):
        position = conformer.GetAtomPosition(i)
        coords.extend([position.x, position.y, position.z])

    connections, orders = [], []
    for bond in mol.GetBonds():
        connections.extend([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
        orders.append(int(round(bond.GetBondTypeAsDouble())))

    return {
        "chemicalJson": 1,
        "atoms": {
            "elements": {"number": [a.GetAtomicNum() for a in mol.GetAtoms()]},
            "coords": {"3d": coords},
        },
        "bonds": {"connections": {"index": connections}, "order": orders},
    }


def avo_input(smiles, options=None, selected=()):
    return {
        "cjson": cjson_from_smiles(smiles),
        "selectedAtoms": list(selected),
        "charge": 0,
        "spin": 1,
        "options": options or {},
    }


def invoke(feature, payload):
    """Run the real entry point the way Avogadro does, and parse its JSON."""
    process = subprocess.run(
        [sys.executable, "-c", "from avogadro_stk import main; main()", feature],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        cwd=ROOT,
        env={**_env(), "PYTHONPATH": str(ROOT / "src")},
    )
    assert process.returncode == 0, process.stderr
    assert process.stdout.strip(), f"no output; stderr was:\n{process.stderr}"
    return json.loads(process.stdout)


def _env():
    import os
    return dict(os.environ)


@pytest.fixture
def tfb():
    """1,3,5-Triformylbenzene — the canonical tritopic aldehyde."""
    return "O=Cc1cc(C=O)cc(C=O)c1"


@pytest.fixture
def terephthalaldehyde():
    return "O=Cc1ccc(C=O)cc1"


@pytest.fixture
def dibromobenzene():
    return "Brc1ccc(Br)cc1"
