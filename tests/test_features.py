#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""The JSON contract, exercised through the real entry point.

There is no GUI in CI, so each feature is fed the input JSON Avogadro would
send and its stdout is checked. Atom counts and formulas are asserted;
coordinates are not, because MCHammer only fixes its seed by default.
"""

from conftest import avo_input, cjson_from_smiles, invoke

TFB = "O=Cc1cc(C=O)cc(C=O)c1"
PDA = "Nc1ccc(N)cc1"
DIBROMOBENZENE = "Brc1ccc(Br)cc1"


def atom_count(result):
    return len(result["cjson"]["atoms"]["elements"]["number"])


def bond_count(result):
    return len(result["cjson"]["bonds"]["order"])


def assert_structure(result):
    """Invariants every structure-producing feature must satisfy."""
    assert "error" not in result, result.get("error")
    assert result["moleculeFormat"] == "cjson"
    # Distance-based re-perception would wreck a freshly collapsed structure.
    assert result["bond"] is False

    cjson = result["cjson"]
    numbers = cjson["atoms"]["elements"]["number"]
    coords = cjson["atoms"]["coords"]["3d"]
    assert len(coords) == 3 * len(numbers)

    index = cjson["bonds"]["connections"]["index"]
    orders = cjson["bonds"]["order"]
    assert len(index) == 2 * len(orders)
    assert all(0 <= i < len(numbers) for i in index)

    stk_block = cjson["properties"]["stk"]
    assert stk_block["topology"]
    assert stk_block["buildingBlocks"]
    assert stk_block["optimizer"] in ("MCHammer", "none")

    # The result must never let the user think this geometry is optimized.
    if stk_block["optimizer"] == "MCHammer":
        assert "MCHammer" in result["message"]


def test_analyze_reports_groups_and_topologies():
    result = invoke("analyze", avo_input(TFB))
    assert "error" not in result, result.get("error")
    # analyze must not touch the molecule
    assert "cjson" not in result

    message = result["message"]
    assert "Aldehyde" in message
    assert "Topicity: 3" in message
    assert "FourPlusSix" in message
    assert "Polymer: no" in message


def test_analyze_on_ditopic_offers_polymer():
    result = invoke("analyze", avo_input(DIBROMOBENZENE))
    message = result["message"]
    assert "Topicity: 2" in message
    assert "Polymer: yes" in message


def test_analyze_without_functional_groups():
    result = invoke("analyze", avo_input("c1ccccc1"))
    assert "No functional groups" in result["message"]


def test_homopolymer():
    result = invoke("polymer", avo_input(DIBROMOBENZENE, {
        "repeating_unit": "A", "num_repeats": 3, "cleanup": True,
    }))
    assert_structure(result)
    # 3 x C6H4Br2, two couplings each deleting two bromines
    assert atom_count(result) == 3 * 12 - 4
    assert "C18H12Br2" in result["message"]


def test_alternating_copolymer_via_library():
    result = invoke("polymer", avo_input("O=Cc1ccc(C=O)cc1", {
        "monomer_b": "p-Phenylenediamine (2-topic, Primary amine)",
        "fg_a": "Aldehyde",
        "repeating_unit": "AB", "num_repeats": 2, "cleanup": True,
    }))
    assert_structure(result)
    assert "imine" in result["message"]


def test_auto_picks_the_group_that_pairs_with_the_partner():
    """Auto must prefer the aldehyde over the alcohol, given an amine partner."""
    result = invoke("polymer", avo_input("O=Cc1cc(O)c(C=O)cc1O", {
        "monomer_b": "p-Phenylenediamine (2-topic, Primary amine)",
        "repeating_unit": "AB", "num_repeats": 2, "cleanup": True,
    }))
    assert_structure(result)
    assert "imine" in result["message"]
    recipe = result["cjson"]["properties"]["stk"]["buildingBlocks"]
    assert recipe[0]["factory"] == "AldehydeFactory"


def test_polymer_with_end_caps():
    result = invoke("polymer", avo_input(DIBROMOBENZENE, {
        "repeating_unit": "A", "num_repeats": 3,
        "end_cap": "Bromobenzene (1-topic, Bromide)", "cleanup": True,
    }))
    assert_structure(result)
    numbers = result["cjson"]["atoms"]["elements"]["number"]
    assert 35 not in numbers, "end caps should consume every bromine"


def test_macrocycle():
    result = invoke("macrocycle", avo_input(DIBROMOBENZENE, {
        "repeating_unit": "A", "num_repeats": 5, "cleanup": True,
    }))
    assert_structure(result)
    # 5 benzene rings plus the macrocyclic ring they close: the cyclomatic
    # number (bonds - atoms + 1) has to be 6, not the 5 an open chain gives.
    assert bond_count(result) - atom_count(result) + 1 == 6


def test_cage_four_plus_six():
    result = invoke("cage", avo_input(TFB, {
        "partner": "p-Phenylenediamine (2-topic, Primary amine)",
        "fg_a": "Aldehyde",
        "topology": "FourPlusSix", "cleanup": True,
    }))
    assert_structure(result)
    assert "FourPlusSix" in result["message"]
    assert "imine" in result["message"]
    assert result["cjson"]["properties"]["stk"]["topology"] == "cage.FourPlusSix"
    # 4 TFB (18 atoms) + 6 PDA (16 atoms) - 12 imine condensations (H2O each)
    assert atom_count(result) == 4 * 18 + 6 * 16 - 12 * 3


def test_cage_custom_smiles_partner():
    result = invoke("cage", avo_input(TFB, {
        "partner": "Custom SMILES…",
        "partner_smiles": "NCCN",
        "fg_a": "Aldehyde",
        "topology": "FourPlusSix", "cleanup": True,
    }))
    assert_structure(result)
    smiles = [b["smiles"] for b in result["cjson"]["properties"]["stk"]["buildingBlocks"]]
    assert "NCCN" in smiles


def test_cof_non_periodic_has_no_unit_cell():
    result = invoke("cof", avo_input(TFB, {
        "partner": "p-Phenylenediamine (2-topic, Primary amine)",
        "fg_a": "Aldehyde",
        "lattice": "Honeycomb", "periodic": False, "cleanup": True,
    }))
    assert_structure(result)
    assert "unitCell" not in result["cjson"]
    assert "finite fragment" in result["message"]


def test_cof_periodic_emits_unit_cell():
    result = invoke("cof", avo_input(TFB, {
        "partner": "p-Phenylenediamine (2-topic, Primary amine)",
        "fg_a": "Aldehyde",
        "lattice": "Honeycomb", "periodic": True, "cleanup": False,
    }))
    assert_structure(result)
    cell = result["cjson"]["unitCell"]
    assert cell["a"] > 0 and cell["b"] > 0 and cell["c"] > 0
    assert 0 < cell["gamma"] < 180
    assert "experimental" in result["message"]


def test_selection_restricts_reacting_groups():
    """A tritopic molecule used as a ditopic one by selecting two groups."""
    cjson = cjson_from_smiles(TFB)
    payload = {
        "cjson": cjson,
        # The three aldehyde carbons of triformylbenzene are atoms 1, 5 and 9;
        # selecting two of them should leave a ditopic block.
        "selectedAtoms": [0, 1, 4, 5],
        "options": {
            "fg_a": "Aldehyde", "use_selection": True,
            "repeating_unit": "A", "num_repeats": 2, "cleanup": True,
        },
    }
    result = invoke("polymer", payload)
    assert_structure(result)
    assert "selection" in result["message"].lower()


def test_provenance_records_both_blocks():
    result = invoke("cage", avo_input(TFB, {
        "partner": "p-Phenylenediamine (2-topic, Primary amine)",
        "fg_a": "Aldehyde", "topology": "FourPlusSix", "cleanup": True,
    }))
    recipe = result["cjson"]["properties"]["stk"]
    roles = {b["role"] for b in recipe["buildingBlocks"]}
    assert roles == {"A", "B"}
    assert recipe["buildingBlocks"][0]["topicity"] == 3
    assert recipe["buildingBlocks"][1]["topicity"] == 2
    assert recipe["pluginVersion"]
    assert recipe["version"]
