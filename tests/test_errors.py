#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Failures must arrive as an `error` naming the specific problem.

A traceback on stdout gives Avogadro nothing to display, and a vague message
gives the user nothing to do, so both are treated as bugs here.
"""

from conftest import avo_input, invoke

TFB = "O=Cc1cc(C=O)cc(C=O)c1"
PDA_LABEL = "p-Phenylenediamine (2-topic, Primary amine)"


def error_of(result):
    assert "error" in result, f"expected an error, got {result}"
    assert "cjson" not in result, "a failed run must not replace the molecule"
    message = result["error"]
    assert "Traceback" not in message
    return message


def test_topicity_mismatch_suggests_an_alternative():
    """Two ditopic blocks cannot make FourPlusSix; say what would."""
    message = error_of(invoke("cage", avo_input("O=Cc1ccc(C=O)cc1", {
        "partner": PDA_LABEL, "fg_a": "Aldehyde", "topology": "FourPlusSix",
    })))
    assert "FourPlusSix" in message
    assert "3-topic" in message
    assert "2" in message


def test_unparseable_smiles():
    message = error_of(invoke("cage", avo_input(TFB, {
        "partner": "Custom SMILES…", "partner_smiles": "this is not smiles",
        "fg_a": "Aldehyde", "topology": "FourPlusSix",
    })))
    assert "SMILES" in message


def test_empty_custom_smiles():
    message = error_of(invoke("cage", avo_input(TFB, {
        "partner": "Custom SMILES…", "partner_smiles": "  ",
        "fg_a": "Aldehyde", "topology": "FourPlusSix",
    })))
    assert "empty" in message.lower()


def test_no_molecule_open():
    result = invoke("analyze", {"cjson": {}, "options": {}})
    message = error_of(result)
    assert "no molecule" in message.lower()


def test_molecule_without_reactive_groups():
    message = error_of(invoke("polymer", avo_input("c1ccccc1", {
        "repeating_unit": "A", "num_repeats": 3,
    })))
    assert "functional group" in message.lower()


def test_oversized_request_is_refused():
    """A big lattice of big blocks must fail fast, before construction."""
    tetratopic = "O=Cc1ccc(cc1)C(c1ccc(C=O)cc1)(c1ccc(C=O)cc1)c1ccc(C=O)cc1"
    message = error_of(invoke("cof", avo_input(tetratopic, {
        "partner": PDA_LABEL, "fg_a": "Aldehyde",
        "lattice": "Kagome", "size_a": 6, "size_b": 6, "size_c": 6,
    })))
    assert "atoms" in message
    assert "smaller" in message.lower()


def test_tritopic_block_in_a_polymer_suggests_a_cage():
    """stk's own message here is unusable; ours must point somewhere."""
    message = error_of(invoke("polymer", avo_input(TFB, {
        "repeating_unit": "A", "num_repeats": 3,
    })))
    assert "3 reactive groups" in message
    assert "Cage" in message


def test_repeating_unit_names_the_offending_letter():
    message = error_of(invoke("polymer", avo_input("Brc1ccc(Br)cc1", {
        "repeating_unit": "AXB", "num_repeats": 2,
    })))
    assert "'X'" in message or "X" in message


def test_repeating_unit_referencing_absent_monomer():
    message = error_of(invoke("polymer", avo_input("Brc1ccc(Br)cc1", {
        "repeating_unit": "AB", "num_repeats": 2, "monomer_b": "(none)",
    })))
    assert "B" in message


def test_selection_required_when_option_is_on():
    message = error_of(invoke("polymer", avo_input("Brc1ccc(Br)cc1", {
        "use_selection": True, "repeating_unit": "A", "num_repeats": 2,
    })))
    assert "select" in message.lower()


def test_unknown_topology():
    message = error_of(invoke("cage", avo_input(TFB, {
        "partner": PDA_LABEL, "topology": "NotATopology",
    })))
    assert "NotATopology" in message


def test_metal_cage_is_refused_clearly():
    message = error_of(invoke("cage", avo_input(TFB, {
        "partner": PDA_LABEL, "topology": "M4L6Tetrahedron",
    })))
    assert "dative" in message.lower() or "metal" in message.lower()


def test_incompatible_functional_groups_are_refused():
    """A dibromide and a trialdehyde have the right topicities but no reaction.

    stk would build this: it deletes both groups' atoms and makes a C–C bond,
    handing back a hydrocarbon cage with no Br and no O in it. Reporting that
    as a result is worse than refusing.
    """
    message = error_of(invoke("cage", avo_input("BrCCCBr", {
        "partner": "Custom SMILES…", "partner_smiles": "O=CC(C=O)CC=O",
        "topology": "FourPlusSix", "cleanup": True,
    })))
    assert "Bromide" in message and "Aldehyde" in message
    assert "linkage" in message.lower()


def test_incompatible_self_reaction_is_refused():
    """An aldehyde does not homopolymerize; say so instead of building it."""
    message = error_of(invoke("polymer", avo_input("O=Cc1ccc(C=O)cc1", {
        "repeating_unit": "A", "num_repeats": 3, "cleanup": True,
    })))
    assert "Aldehyde" in message
    assert "primary amine" in message


def test_compatible_pairs_still_build():
    """The guard must not reject real chemistry: imine cage, bromide polymer."""
    for feature, payload in (
        ("cage", avo_input("NCCN", {
            "partner": "Custom SMILES…", "partner_smiles": "O=CC(C=O)CC=O",
            "topology": "FourPlusSix", "cleanup": True,
        })),
        ("polymer", avo_input("Brc1ccc(Br)cc1", {
            "repeating_unit": "A", "num_repeats": 3, "cleanup": True,
        })),
    ):
        result = invoke(feature, payload)
        assert "error" not in result, result.get("error")


def test_smiles_typed_without_switching_the_dropdown():
    """The obvious dialog mistake: fill in the SMILES, leave the list alone.

    Nobody types a SMILES and means "no partner", so the text wins over the
    untouched "(none)" rather than the run being refused.
    """
    result = invoke("cage", avo_input("NCCN", {
        "partner": "(none)", "partner_smiles": "O=CC(C=O)CC=O",
        "topology": "FourPlusSix", "cleanup": True,
    }))
    assert "error" not in result, result.get("error")
    assert "imine" in result["message"]


def test_no_partner_at_all_names_both_dialog_fields():
    message = error_of(invoke("cage", avo_input("NCCN", {
        "partner": "(none)", "partner_smiles": "", "topology": "FourPlusSix",
    })))
    assert "Partner Block" in message
    assert "Partner SMILES" in message
