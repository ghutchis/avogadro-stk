#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""The generated topology table is a contract with stk.

Every entry is built with dummy blocks of the recorded topicity. If stk
changes a topology's vertex functionality, this fails loudly instead of
producing a wrong structure or a confusing mismatch message.
"""

import pytest
import stk

from avogadro_stk import validate
from avogadro_stk.topologies import TOPOLOGIES

#: Simple blocks of each topicity we need to instantiate a topology.
DUMMIES = {
    2: ("BrCCBr", "BromoFactory"),
    3: ("BrCC(CBr)CBr", "BromoFactory"),
    4: ("BrCC(CBr)(CBr)CBr", "BromoFactory"),
    5: ("BrCc1c(CBr)c(CBr)c(CBr)c(CBr)c1", "BromoFactory"),
    6: ("BrCc1c(CBr)c(CBr)c(CBr)c(CBr)c1CBr", "BromoFactory"),
}

ORGANIC = {
    name: spec for name, spec in TOPOLOGIES.items() if not spec["metal"]
}


def dummy(topicity):
    smiles, factory = DUMMIES[topicity]
    return stk.BuildingBlock(smiles, [getattr(stk, factory)()])


def test_every_functionality_has_a_dummy():
    needed = {
        functionality
        for spec in ORGANIC.values()
        for functionality in spec["requirements"]
    }
    assert needed <= set(DUMMIES), f"no dummy block for {needed - set(DUMMIES)}"


@pytest.mark.parametrize("name", sorted(ORGANIC), ids=lambda n: n)
def test_topology_builds_with_recorded_topicity(name):
    spec = ORGANIC[name]
    kind, class_name = name.split(".", 1)
    cls = getattr(getattr(stk, kind), class_name)

    building_blocks = tuple(dummy(f) for f in sorted(spec["requirements"]))
    kwargs = {"building_blocks": building_blocks}
    if kind == "cof":
        kwargs["lattice_size"] = (1, 1, 1)

    constructed = stk.ConstructedMolecule(cls(**kwargs))
    assert constructed.get_num_atoms() > 0


@pytest.mark.parametrize("name", sorted(ORGANIC), ids=lambda n: n)
def test_requirements_match_vertex_count(name):
    spec = ORGANIC[name]
    assert sum(spec["requirements"].values()) == spec["vertices"]


def test_available_excludes_metal_topologies():
    for name in validate.available("cage"):
        assert not TOPOLOGIES[f"cage.{name}"]["metal"]


def test_describe_reads_naturally():
    text = validate.describe("cage", "FourPlusSix")
    assert text == "6 × 2-topic and 4 × 3-topic building blocks"


def test_check_accepts_matching_blocks():
    assert validate.check("cage", "FourPlusSix", 3, 2)


def test_check_names_the_missing_topicity():
    with pytest.raises(Exception) as excinfo:
        validate.check("cage", "FourPlusSix", 2, 2)
    message = str(excinfo.value)
    assert "FourPlusSix" in message
    # It must suggest something concrete, not just report incompatibility.
    assert "3-topic" in message
    assert any(word in message for word in ("Triformylbenzene", "TAPB", "Tris"))
