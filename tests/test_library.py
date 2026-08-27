#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""The library's declared topicity must be what stk actually finds.

An entry whose `topicity` is wrong produces a topology mismatch that looks
like the user's fault, so it is checked here rather than discovered in a
dialog.
"""

import json
import pathlib

import pytest
import stk
from rdkit import Chem

from avogadro_stk import blocks, library

ENTRIES = list(library.entries())
ROOT = pathlib.Path(__file__).resolve().parent.parent


@pytest.mark.parametrize("entry", ENTRIES, ids=lambda e: e["id"])
def test_smiles_parses(entry):
    assert Chem.MolFromSmiles(entry["smiles"]) is not None


@pytest.mark.parametrize(
    "entry", [e for e in ENTRIES if e.get("factory")], ids=lambda e: e["id"]
)
def test_factory_finds_declared_topicity(entry):
    block = stk.BuildingBlock(
        entry["smiles"], [getattr(stk, entry["factory"])()]
    )
    assert block.get_num_functional_groups() == entry["topicity"], (
        f"{entry['id']}: {entry['factory']} found "
        f"{block.get_num_functional_groups()}, declared {entry['topicity']}"
    )


@pytest.mark.parametrize("entry", ENTRIES, ids=lambda e: e["id"])
def test_factory_is_known(entry):
    factory = entry.get("factory")
    assert factory is None or factory in blocks.FACTORY_NAMES


def test_ids_are_unique():
    ids = [entry["id"] for entry in ENTRIES]
    assert len(ids) == len(set(ids))


def test_labels_are_unique():
    """Dropdown text is the only handle the dialog has on an entry."""
    labels = [library.label(entry) for entry in ENTRIES]
    assert len(labels) == len(set(labels))


@pytest.mark.parametrize("entry", ENTRIES, ids=lambda e: e["id"])
def test_label_round_trips(entry):
    assert library.lookup(library.label(entry)) is entry


def test_lookup_by_alias():
    assert library.lookup("TFB")["id"] == "tfb"
    assert library.lookup("tfb")["id"] == "tfb"


@pytest.mark.parametrize(
    "name", ["polymer_options", "macrocycle_options", "cage_options",
             "cof_options"]
)
def test_generated_options_match_library(name):
    """The dialogs are generated; drift means a dropdown entry cannot resolve."""
    path = ROOT / "src" / "avogadro_stk" / "customization" / f"{name}.json"
    options = json.loads(path.read_text(encoding="utf-8"))

    for key, option in options.items():
        if option.get("type") != "stringList":
            continue
        default = option.get("default", 0)
        assert 0 <= default < len(option["values"]), f"{name}.{key} default"
        for value in option["values"]:
            if value in (library.NONE_LABEL, library.CUSTOM_LABEL,
                         blocks.AUTO_LABEL):
                continue
            if key in ("fg_a", "fg_b"):
                assert blocks.factory_from_display(value) is not None
            elif key in ("partner", "monomer_b", "end_cap"):
                if value == "Keep functional groups":
                    continue
                assert library.lookup(value) is not None, f"{name}: {value}"
