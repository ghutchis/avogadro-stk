#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Covalent organic frameworks — stk.cof.*.

Non-periodic output is the default: it is a finite molecular fragment of the
framework and behaves like any other molecule in Avogadro. Periodic output
carries a unit cell and bonds that wrap across the cell boundary, which
Avogadro handles less predictably, so it sits behind a checkbox.
"""

import stk

from .. import blocks, provenance, resolve, validate
from ..convert import ConversionError, constructed_to_cjson, formula
from ..result import molecule, summarize

#: The non-periodic lattices offered in the dialog.
LATTICES = ("Honeycomb", "Hexagonal", "Kagome", "Square", "LinkerlessHoneycomb")

DISPLAY = {
    "Honeycomb": "Honeycomb",
    "Hexagonal": "Hexagonal",
    "Kagome": "Kagome",
    "Square": "Square",
    "LinkerlessHoneycomb": "Linkerless Honeycomb",
}


def _lattice_name(text):
    wanted = (text or "Honeycomb").strip()
    for name, shown in DISPLAY.items():
        if wanted in (name, shown):
            return name
    raise ConversionError(
        f"Unknown lattice '{text}'. Choose one of: "
        f"{', '.join(DISPLAY.values())}."
    )


def run(avo_input, options):
    lattice = _lattice_name(options.get("lattice"))
    periodic = bool(options.get("periodic"))
    topology_name = f"Periodic{lattice}" if periodic else lattice
    spec = validate.spec("cof", topology_name)

    size = tuple(
        max(1, int(options.get(key, 1) or 1)) for key in ("size_a", "size_b", "size_c")
    )

    single_block = len(spec["requirements"]) == 1
    partner_block = resolve.partner(
        options.get("partner"),
        options.get("partner_smiles"),
        options.get("fg_b", blocks.AUTO_LABEL),
        minimum=2,
        role="partner block",
    )
    if partner_block is None and not single_block:
        raise ConversionError(
            f"{DISPLAY[lattice]} needs {validate.describe('cof', topology_name)}. "
            f"Choose a partner block from the library, or enter one as SMILES."
        )

    current_block = resolve.current(
        avo_input,
        options.get("fg_a", blocks.AUTO_LABEL),
        partner_block=partner_block,
        minimum=2,
        use_selection=bool(options.get("use_selection")),
    )

    validate.check(
        "cof", topology_name,
        current_block.topicity,
        partner_block.topicity if partner_block else None,
        label_a="your molecule",
        label_b=partner_block.name if partner_block else "the partner",
    )

    sizes = {current_block.topicity: current_block.block.get_num_atoms()}
    if partner_block is not None:
        sizes[partner_block.topicity] = partner_block.block.get_num_atoms()
    warning = validate.check_size(
        validate.estimate_atoms("cof", topology_name, sizes, size), "framework"
    )

    used = [current_block]
    if partner_block is not None:
        used.append(partner_block)

    optimizer = stk.MCHammer() if options.get("cleanup", True) else stk.NullOptimizer()
    kwargs = {
        "building_blocks": tuple(b.block for b in used),
        "lattice_size": size,
        "optimizer": optimizer,
    }
    if not periodic:
        # The non-periodic classes take an explicit flag; the Periodic* ones
        # do not accept it at all.
        kwargs["periodic"] = False

    try:
        topology = getattr(stk.cof, topology_name)(**kwargs)
        constructed = stk.ConstructedMolecule(topology)
        periodic_info = topology.construct().get_periodic_info() if periodic else None
    except Exception as exc:  # noqa: BLE001 - stk raises a family of types
        raise ConversionError(
            f"stk could not build the framework: {exc}"
        ) from None

    cjson = constructed_to_cjson(constructed, periodic_info)
    provenance.record(
        cjson,
        topology=f"cof.{topology_name}",
        building_blocks=[b.record(role) for b, role in zip(used, "AB")],
        optimizer="MCHammer" if options.get("cleanup", True) else "none",
        latticeSize=list(size),
        periodic=periodic,
    )

    warnings = [warning] if warning else []
    partner_factory = partner_block.factory if partner_block else current_block.factory
    linkage = blocks.linkage_name(current_block.factory, partner_factory)

    headline = (
        f"Built a {DISPLAY[lattice]} COF, {size[0]}×{size[1]}×{size[2]}: "
        f"{constructed.get_num_atoms()} atoms ({formula(constructed)}), "
        f"joined by {linkage} bonds."
    )

    if periodic:
        crossing = sum(
            1 for bond in constructed.get_bonds()
            if bond.get_periodicity() != (0, 0, 0)
        )
        warnings.append(
            f"Periodic output is experimental. The unit cell comes from stk's "
            f"idealised vertex positions, not a relaxed structure, and "
            f"{crossing} bond(s) wrap across the cell boundary — Avogadro "
            f"draws those as long bonds across the cell. Optimize with a "
            f"periodic method before trusting the cell parameters."
        )
        cell = cjson["unitCell"]
        headline += (
            f"\nCell: a={cell['a']:.2f} Å, b={cell['b']:.2f} Å, "
            f"c={cell['c']:.2f} Å, γ={cell['gamma']:.1f}°"
        )
    else:
        warnings.append(
            "This is a finite fragment of the framework, so the edges have "
            "unreacted groups. Tick 'Periodic' for a unit cell instead."
        )

    notes = [b.note for b in used if b.note]
    return molecule(
        cjson,
        summarize(
            [headline, *notes], warnings,
            cleanup=bool(options.get("cleanup", True)),
        ),
    )
