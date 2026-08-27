#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Regenerate the user-options JSON that Avogadro builds its dialogs from.

Avogadro reads the option lists from static files (a feature is never asked
for its options with the molecule in hand), so the dropdown entries have to be
baked out ahead of time. Generating them from ``linkers.json`` and the factory
list is what stops the dialogs drifting away from the library.

Run after editing the library::

    python tools/generate_options.py
"""

import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from avogadro_stk import blocks, library, validate  # noqa: E402
from avogadro_stk.features.chains import KEEP_GROUPS  # noqa: E402
from avogadro_stk.features.cof import DISPLAY as COF_LATTICES  # noqa: E402

OUT = ROOT / "src" / "avogadro_stk" / "customization"

GROUPS = [blocks.AUTO_LABEL] + blocks.display_names()


def fg_option(label, order, tip):
    return {
        "type": "stringList", "label": label, "default": 0, "order": order,
        "toolTip": tip, "values": GROUPS,
    }


def selection_option(order):
    return {
        "type": "boolean", "label": "Use Selected Atoms Only", "default": False,
        "order": order,
        "toolTip": "React only the functional groups that touch the current "
                   "selection, instead of every group in the molecule.",
    }


def cleanup_option(order):
    return {
        "type": "boolean", "label": "Clean Up Geometry", "default": True,
        "order": order,
        "toolTip": "Run stk's MCHammer to remove the worst of the idealised "
                   "vertex geometry. This is not a force field — optimize the "
                   "result afterwards.",
    }


def partner_values(tags):
    return [library.NONE_LABEL] + library.labels(tags=tags) + [library.CUSTOM_LABEL]


def smiles_option(order):
    return {
        "type": "string", "label": "Partner SMILES", "default": "", "order": order,
        "toolTip": f"Used only when the partner is '{library.CUSTOM_LABEL}'.",
    }


def polymer_options():
    caps = [KEEP_GROUPS] + library.labels(tags=["cap"])
    return {
        "fg_a": fg_option("Monomer A Groups", 0,
                          "Which groups of the current molecule react. Auto "
                          "picks the group that pairs with monomer B."),
        "use_selection": selection_option(1),
        "monomer_b": {
            "type": "stringList", "label": "Monomer B", "default": 0, "order": 2,
            "toolTip": "Second monomer. Leave as (none) for a homopolymer.",
            "values": partner_values(["polymer", "cap"]),
        },
        "fg_b": fg_option("Monomer B Groups", 3,
                          f"Only used with '{library.CUSTOM_LABEL}'; library "
                          f"entries carry their own group."),
        "b_smiles": smiles_option(4),
        "repeating_unit": {
            "type": "string", "label": "Repeating Unit", "default": "AB", "order": 5,
            "toolTip": "A is the current molecule, B is monomer B. 'AB' "
                       "alternates; 'A' is a homopolymer.",
        },
        "num_repeats": {
            "type": "integer", "label": "Number of Repeats", "default": 3,
            "minimum": 1, "maximum": 50, "order": 6,
            "toolTip": "How many times the repeating unit is repeated.",
        },
        "end_cap": {
            "type": "stringList", "label": "End Caps", "default": 0, "order": 7,
            "toolTip": "By default the chain ends keep their unreacted groups. "
                       "A monotopic block caps them instead.",
            "values": caps,
        },
        "cleanup": cleanup_option(8),
    }


def macrocycle_options():
    options = polymer_options()
    del options["end_cap"]
    options["repeating_unit"]["toolTip"] += " The chain is then closed into a ring."
    options["num_repeats"]["default"] = 4
    return options


def cage_options():
    topologies = validate.available("cage", periodic=False, metal=False)
    return {
        "fg_a": fg_option("Your Molecule's Groups", 0,
                          "Which groups of the current molecule react."),
        "use_selection": selection_option(1),
        "partner": {
            "type": "stringList", "label": "Partner Block", "default": 0, "order": 2,
            "toolTip": "The second building block. Its topicity must match "
                       "what the topology needs.",
            "values": partner_values(["cage"]),
        },
        "fg_b": fg_option("Partner Groups", 3,
                          f"Only used with '{library.CUSTOM_LABEL}'."),
        "partner_smiles": smiles_option(4),
        "topology": {
            "type": "stringList", "label": "Topology",
            "default": topologies.index("FourPlusSix"), "order": 5,
            "toolTip": "Cage topologies are named for their block counts: "
                       "FourPlusSix is 4 tritopic + 6 ditopic.",
            "values": topologies,
        },
        "isomer": {
            "type": "integer", "label": "Structural Isomer", "default": 0,
            "minimum": 0, "maximum": 20, "order": 6,
            "toolTip": "0 is stk's default vertex alignment. Higher values "
                       "rotate a vertex to give a different isomer.",
        },
        "cleanup": cleanup_option(7),
    }


def cof_options():
    lattices = list(COF_LATTICES.values())
    return {
        "fg_a": fg_option("Your Molecule's Groups", 0,
                          "Which groups of the current molecule react."),
        "use_selection": selection_option(1),
        "partner": {
            "type": "stringList", "label": "Partner Block", "default": 0, "order": 2,
            "toolTip": "The second building block.",
            "values": partner_values(["cof"]),
        },
        "fg_b": fg_option("Partner Groups", 3,
                          f"Only used with '{library.CUSTOM_LABEL}'."),
        "partner_smiles": smiles_option(4),
        "lattice": {
            "type": "stringList", "label": "Lattice",
            "default": lattices.index("Honeycomb"), "order": 5,
            "toolTip": "Honeycomb needs 2 tritopic + 3 ditopic; Square needs "
                       "1 tetratopic + 2 ditopic; Kagome 3 tetratopic + 6 "
                       "ditopic.",
            "values": lattices,
        },
        "size_a": {"type": "integer", "label": "Lattice Repeats a", "default": 1,
                   "minimum": 1, "maximum": 10, "order": 6,
                   "toolTip": "Unit cells along a. Atom counts grow fast."},
        "size_b": {"type": "integer", "label": "Lattice Repeats b", "default": 1,
                   "minimum": 1, "maximum": 10, "order": 7,
                   "toolTip": "Unit cells along b."},
        "size_c": {"type": "integer", "label": "Lattice Repeats c", "default": 1,
                   "minimum": 1, "maximum": 10, "order": 8,
                   "toolTip": "Unit cells along c. Layered COFs use 1."},
        "periodic": {
            "type": "boolean", "label": "Periodic", "default": False, "order": 9,
            "toolTip": "Emit a unit cell instead of a finite fragment. "
                       "Experimental — see the README.",
        },
        "cleanup": cleanup_option(10),
    }


BUILDERS = {
    "polymer": polymer_options,
    "macrocycle": macrocycle_options,
    "cage": cage_options,
    "cof": cof_options,
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, builder in BUILDERS.items():
        target = OUT / f"{name}_options.json"
        target.write_text(json.dumps(builder(), indent=2) + "\n", encoding="utf-8")
        print(f"wrote {target.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
