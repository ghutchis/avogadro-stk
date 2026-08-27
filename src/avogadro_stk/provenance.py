#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""The construction recipe, stored in the CJSON.

A constructed molecule is otherwise a dead end: the coordinates say nothing
about which blocks and topology produced them. Recording the recipe means the
structure can be rebuilt, and it is what a future `rebuild` feature reads to
re-open the dialog pre-filled.
"""

import stk

from . import blocks

try:  # the plugin is normally installed, but tests may run from a checkout
    from importlib.metadata import PackageNotFoundError, version
    try:
        PLUGIN_VERSION = version("avogadro-stk")
    except PackageNotFoundError:
        PLUGIN_VERSION = "0.0.0+dev"
except ImportError:  # pragma: no cover
    PLUGIN_VERSION = "0.0.0+dev"


def block_record(block, factory, role):
    return {
        "smiles": stk.Smiles().get_key(block),
        "factory": factory,
        "role": role,
        "topicity": block.get_num_functional_groups(),
    }


def record(cjson, topology, building_blocks, optimizer="MCHammer", **extra):
    """Attach the recipe to a CJSON object and return it.

    ``building_blocks`` is a sequence of ``block_record`` dicts.
    """
    entry = {
        "version": stk.__version__,
        "pluginVersion": PLUGIN_VERSION,
        "topology": topology,
        "buildingBlocks": list(building_blocks),
        "optimizer": optimizer,
    }
    entry.update({k: v for k, v in extra.items() if v is not None})

    properties = cjson.setdefault("properties", {})
    properties["stk"] = entry
    return cjson


def linkage(factory_a, factory_b):
    """Name of the bond being formed, for the result message."""
    return blocks.linkage_name(factory_a, factory_b)
