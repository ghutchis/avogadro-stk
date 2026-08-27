#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Polymer and macrocycle — the same construction, joined up or closed.

stk maps the letters of the repeating unit onto positions in the building
block tuple (A -> 0, B -> 1, ...), so the letters the user types and the
blocks that were actually resolved have to be kept in step; that bookkeeping
is most of what this module does.
"""

import stk

from .. import blocks, geometry, provenance, resolve, validate
from ..convert import ConversionError, constructed_to_cjson, formula
from ..result import molecule, summarize

LETTERS = "ABCDEFGH"

#: Beyond this many repeats the idealised ring geometry is bad enough to be
#: worth warning about.
MACROCYCLE_WARN = 20

#: Dialog wording for "leave the terminal groups alone", which is what stk
#: does by default — the chain ends keep their unreacted functional groups.
KEEP_GROUPS = "Keep functional groups"


def _repeating_unit(text, available):
    """Validate the user's repeating unit against the blocks that exist."""
    unit = (text or "").strip().upper().replace(" ", "")
    if not unit:
        raise ConversionError(
            "The repeating unit is empty. Enter something like 'AB' (monomer "
            "A then monomer B) or 'A' for a homopolymer."
        )

    allowed = LETTERS[:len(available)]
    for character in unit:
        if character not in allowed:
            if character in LETTERS:
                raise ConversionError(
                    f"The repeating unit '{unit}' refers to monomer "
                    f"'{character}', but only {' and '.join(allowed)} "
                    f"{'is' if len(allowed) == 1 else 'are'} defined. Choose a "
                    f"monomer for '{character}' or remove it from the unit."
                )
            raise ConversionError(
                f"'{character}' is not a valid monomer letter in the "
                f"repeating unit '{unit}'. Use letters from {allowed}."
            )

    for index, character in enumerate(unit):
        block = available[LETTERS.index(character)]
        terminal = index in (0, len(unit) - 1)
        if block.topicity < 2 and not terminal:
            raise ConversionError(
                f"Monomer {character} ({block.name}) has only "
                f"{block.topicity} reactive group, so it can only sit at the "
                f"start or end of the chain, not in the middle of '{unit}'."
            )
        if block.topicity > 2:
            # stk's own message here names internal objects and no way out.
            raise ConversionError(
                f"Monomer {character} ({block.name}) has {block.topicity} "
                f"reactive groups, but a chain needs exactly 2 — a third "
                f"group would branch it. Either choose a different functional "
                f"group, select just two groups with 'Use Selected Atoms "
                f"Only', or build a Cage or Covalent Framework instead, which "
                f"is what {block.topicity}-topic blocks are for."
            )
    return unit


def run(avo_input, options, cyclic=False):
    what = "macrocycle" if cyclic else "polymer"

    partner_block = resolve.partner(
        options.get("monomer_b"),
        options.get("b_smiles"),
        options.get("fg_b", blocks.AUTO_LABEL),
        minimum=1,
        role="monomer B",
    )
    current_block = resolve.current(
        avo_input,
        options.get("fg_a", blocks.AUTO_LABEL),
        partner_block=partner_block,
        minimum=2,
        use_selection=bool(options.get("use_selection")),
    )

    available = [current_block]
    if partner_block is not None:
        available.append(partner_block)

    unit = _repeating_unit(options.get("repeating_unit", "AB"), available)
    repeats = max(1, int(options.get("num_repeats", 3) or 1))

    blocks.check_pair(
        current_block.factory,
        partner_block.factory if partner_block else current_block.factory,
        label_a="your molecule",
        label_b=partner_block.name if partner_block else "your molecule",
    )

    cap_block = None
    if not cyclic:
        cap_choice = options.get("end_cap")
        if cap_choice and cap_choice != KEEP_GROUPS:
            cap_block = resolve.partner(
                cap_choice, None, blocks.AUTO_LABEL, minimum=1, role="end cap"
            )
            if cap_block.topicity != 1:
                raise ConversionError(
                    f"End caps must have exactly one reactive group, but "
                    f"{cap_block.name} has {cap_block.topicity}."
                )
            blocks.check_pair(
                current_block.factory, cap_block.factory,
                label_a="your molecule", label_b=cap_block.name,
            )

    # stk indexes the repeating unit into the building block tuple, so the cap
    # is appended and the unit is expanded here rather than repeated by stk.
    chain = available[:]
    if cap_block is not None:
        cap_letter = LETTERS[len(chain)]
        chain.append(cap_block)
        expanded, num_units = f"{cap_letter}{unit * repeats}{cap_letter}", 1
    else:
        expanded, num_units = unit, repeats

    counts = {letter: expanded.count(letter) for letter in set(expanded)}
    estimate = sum(
        chain[LETTERS.index(letter)].block.get_num_atoms() * count * num_units
        for letter, count in counts.items()
    )
    warning = validate.check_size(estimate, what)

    warnings = [warning] if warning else []
    if cyclic and repeats * len(unit) > MACROCYCLE_WARN:
        warnings.append(
            f"{repeats * len(unit)} units is a large ring; the starting "
            f"geometry will be poor even after cleanup."
        )

    plan = geometry.plan(options)
    factory = stk.macrocycle.Macrocycle if cyclic else stk.polymer.Linear
    try:
        topology = factory(
            building_blocks=tuple(b.block for b in chain),
            repeating_unit=expanded,
            num_repeating_units=num_units,
            optimizer=geometry.optimizer(plan),
        )
        constructed = stk.ConstructedMolecule(topology)
    except Exception as exc:  # noqa: BLE001 - stk raises a family of types
        raise ConversionError(
            f"stk could not build the {what}: {exc}"
        ) from None

    constructed, ff_warning = geometry.refine(constructed, plan)
    if ff_warning:
        warnings.append(ff_warning)

    cjson = constructed_to_cjson(constructed)
    provenance.record(
        cjson,
        topology=f"{'macrocycle.Macrocycle' if cyclic else 'polymer.Linear'}",
        building_blocks=[
            b.record(LETTERS[i]) for i, b in enumerate(chain)
        ],
        optimizer=plan.applied,
        repeatingUnit=expanded,
        numRepeatingUnits=num_units,
    )

    partner_factory = partner_block.factory if partner_block else current_block.factory
    linkage = blocks.linkage_name(current_block.factory, partner_factory)
    headline = (
        f"Built a {what} from {len(expanded) * num_units} monomer units: "
        f"{constructed.get_num_atoms()} atoms ({formula(constructed)}), "
        f"joined by {linkage} bonds."
    )
    notes = [b.note for b in chain if b.note]

    return molecule(
        cjson,
        summarize([headline, *notes], warnings,
                  geometry_note=geometry.note(plan)),
    )
