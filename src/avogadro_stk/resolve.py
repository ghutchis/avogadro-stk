#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Turning dialog answers into stk building blocks.

Every feature asks the same two questions — what is the current molecule, and
what is it reacting with — so the answers are resolved in one place. The
current molecule is always resolved *after* the partner, because knowing the
partner's chemistry is what lets "Auto" pick the right functional group.
"""

from . import blocks, library
from .convert import ConversionError, cjson_to_rdkit, mol_from_smiles


class Block:
    """A building block plus the provenance needed to explain it later."""

    def __init__(self, block, factory, name, note=None):
        self.block = block
        self.factory = factory
        self.name = name
        self.note = note

    @property
    def topicity(self):
        return self.block.get_num_functional_groups()

    def record(self, role):
        from . import provenance
        return provenance.block_record(self.block, self.factory, role)


def partner(choice, smiles, group_choice, minimum=1, role="partner"):
    """Resolve the second building block from its three dialog options.

    Returns None when the user wants no partner at all.
    """
    text = (smiles or "").strip()

    if not choice or choice == library.NONE_LABEL:
        # Typing a SMILES without also switching the drop-down away from
        # "(none)" is the obvious mistake to make in this dialog, and it can
        # only mean one thing — nobody fills in a SMILES field and means "no
        # partner" — so honour the text rather than rejecting the run.
        if not text:
            return None
        choice = library.CUSTOM_LABEL

    if choice == library.CUSTOM_LABEL:
        if not text:
            raise ConversionError(
                f"'{library.CUSTOM_LABEL}' was selected but the SMILES field "
                f"is empty. Enter a SMILES string for the {role}."
            )
        factory = blocks.factory_from_display(group_choice)
        mol = mol_from_smiles(text, label=role)
        note = None
        if factory is None:
            factory, note = blocks.choose(blocks.detect(mol), minimum=minimum)
        block = blocks.from_smiles(text, factory, minimum=minimum, label=role)
        return Block(block, factory, text, note)

    entry = library.lookup(choice)
    if entry is None:
        raise ConversionError(
            f"'{choice}' is not in the building block library. Pick an entry "
            f"from the list or use {library.CUSTOM_LABEL}."
        )
    block = library.block(choice, minimum=minimum, role=role)
    return Block(block, entry["factory"], entry["name"])


def current(avo_input, group_choice, partner_block=None, minimum=1,
            use_selection=False):
    """Resolve the open molecule into a building block.

    ``use_selection`` restricts the reaction to functional groups touching the
    user's current selection.
    """
    mol = cjson_to_rdkit(avo_input.get("cjson", {}))
    factory = blocks.factory_from_display(group_choice)

    note = None
    if factory is None:
        factory, note = blocks.choose(
            blocks.detect(mol),
            partner_factory=partner_block.factory if partner_block else None,
            minimum=minimum,
        )

    selected = avo_input.get("selectedAtoms") or []
    if use_selection and not selected:
        raise ConversionError(
            "'Use selected atoms only' is on, but nothing is selected. Select "
            "the atoms of the groups that should react, or turn the option off."
        )

    block = blocks.build(
        mol,
        factory,
        selected_atoms=selected if use_selection else None,
        minimum=minimum,
    )
    if use_selection:
        note = (
            f"Restricted to the {block.get_num_functional_groups()} "
            f"{blocks.display(factory).lower()} group(s) touching the selection."
        )
    return Block(block, factory, "current molecule", note)
