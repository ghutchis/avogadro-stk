#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Organic cages — stk.cage.*.

Metal-organic M<n>L<n> topologies are excluded: they need SingleAtom
functional groups on the metal centres and dative bond handling that Avogadro
does not yet represent well.
"""

import stk

from .. import blocks, geometry, provenance, resolve, validate
from ..convert import ConversionError, constructed_to_cjson, formula
from ..result import molecule, summarize

CURRENT_IS_A = "Building block A"


def run(avo_input, options):
    topology_name = (options.get("topology") or "FourPlusSix").strip()
    spec = validate.spec("cage", topology_name)
    if spec["metal"]:
        raise ConversionError(
            f"{topology_name} is a metal-organic cage. Those need dative "
            f"bonds to metal centres, which this plugin does not build yet."
        )

    partner_block = resolve.partner(
        options.get("partner"),
        options.get("partner_smiles"),
        options.get("fg_b", blocks.AUTO_LABEL),
        minimum=2,
        role="partner block",
    )
    if partner_block is None:
        raise ConversionError(
            "A cage needs two different building blocks, but no partner was "
            "given. Pick one from the 'Partner Block' list, or type a SMILES "
            "string into the 'Partner SMILES' field."
        )

    current_block = resolve.current(
        avo_input,
        options.get("fg_a", blocks.AUTO_LABEL),
        partner_block=partner_block,
        minimum=2,
        use_selection=bool(options.get("use_selection")),
    )

    validate.check(
        "cage", topology_name,
        current_block.topicity, partner_block.topicity,
        label_a="your molecule", label_b=partner_block.name,
    )
    blocks.check_pair(
        current_block.factory, partner_block.factory,
        label_a="your molecule", label_b=partner_block.name,
    )

    estimate = validate.estimate_atoms(
        "cage", topology_name,
        {
            current_block.topicity: current_block.block.get_num_atoms(),
            partner_block.topicity: partner_block.block.get_num_atoms(),
        },
    )
    warning = validate.check_size(estimate, "cage")

    alignments = _alignments(spec, options.get("isomer", 0))
    plan = geometry.plan(options)

    try:
        topology = getattr(stk.cage, topology_name)(
            building_blocks=(current_block.block, partner_block.block),
            vertex_alignments=alignments,
            optimizer=geometry.optimizer(plan),
        )
        constructed = stk.ConstructedMolecule(topology)
    except Exception as exc:  # noqa: BLE001 - stk raises a family of types
        raise ConversionError(f"stk could not build the cage: {exc}") from None

    constructed, ff_warning = geometry.refine(constructed, plan)

    cjson = constructed_to_cjson(constructed)
    provenance.record(
        cjson,
        topology=f"cage.{topology_name}",
        building_blocks=[current_block.record("A"), partner_block.record("B")],
        optimizer=plan.applied,
        vertexAlignments=alignments,
    )

    linkage = blocks.linkage_name(current_block.factory, partner_block.factory)
    headline = (
        f"Built a {topology_name} {linkage} cage: {spec['vertices']} building "
        f"blocks, {constructed.get_num_atoms()} atoms "
        f"({formula(constructed)}), {spec['edges']} new {linkage} bonds."
    )
    composition = (
        f"{_count(spec, current_block.topicity)} × your molecule "
        f"({current_block.topicity}-topic) and "
        f"{_count(spec, partner_block.topicity)} × {partner_block.name} "
        f"({partner_block.topicity}-topic)."
    )
    notes = [b.note for b in (current_block, partner_block) if b.note]

    warnings = [w for w in (warning, ff_warning) if w]
    return molecule(
        cjson,
        summarize([headline, composition, *notes], warnings,
                  geometry_note=geometry.note(plan)),
    )


def _count(spec, topicity):
    return spec["requirements"].get(topicity, 0)


def _alignments(spec, isomer):
    """Turn the isomer index into a vertex_alignments dict.

    stk enumerates structural isomers by rotating which edge each vertex
    aligns to. Index 0 is stk's default; higher indices walk the vertices in
    order, which is a cheap way to cycle through isomers without exposing the
    full combinatorial space.
    """
    try:
        isomer = int(isomer or 0)
    except (TypeError, ValueError):
        isomer = 0
    if isomer <= 0:
        return None

    # Only vertices that join three or more edges have a meaningful choice.
    tunable = sorted(
        functionality for functionality in spec["requirements"] if functionality > 2
    )
    if not tunable:
        return None

    vertices = spec["vertices"]
    return {(isomer - 1) % vertices: isomer % max(2, min(tunable))}
