#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Functional group detection and stk BuildingBlock construction.

stk needs to be told which functional groups react. Asking the user to name a
factory is a poor first experience, so every factory is run over the molecule
and the result is chosen automatically where the choice is unambiguous. The
`analyze` feature exists to show that reasoning when it is not.
"""

import stk

from .convert import ConversionError

#: Every functional group factory stk offers, minus the two abstract ones.
FACTORY_NAMES = (
    "AlcoholFactory",
    "AldehydeFactory",
    "AmideFactory",
    "BoronicAcidFactory",
    "BromoFactory",
    "CarboxylicAcidFactory",
    "DibromoFactory",
    "DifluoroFactory",
    "DiolFactory",
    "FluoroFactory",
    "IodoFactory",
    "PrimaryAminoFactory",
    "RingAmineFactory",
    "SecondaryAminoFactory",
    "TerminalAlkeneFactory",
    "TerminalAlkyneFactory",
    "ThioacidFactory",
    "ThiolFactory",
)

#: Factory pairs that form a recognised linkage. Used to break ties during
#: auto-selection and to name the bond in the result message. stk itself will
#: happily join any two functional groups; this table is about picking the
#: chemistry the user almost certainly meant.
REACTION_PAIRS = {
    frozenset({"AldehydeFactory", "PrimaryAminoFactory"}): "imine",
    frozenset({"AldehydeFactory", "RingAmineFactory"}): "imine",
    frozenset({"AldehydeFactory", "SecondaryAminoFactory"}): "enamine",
    frozenset({"CarboxylicAcidFactory", "PrimaryAminoFactory"}): "amide",
    frozenset({"CarboxylicAcidFactory", "SecondaryAminoFactory"}): "amide",
    frozenset({"CarboxylicAcidFactory", "AlcoholFactory"}): "ester",
    frozenset({"ThioacidFactory", "PrimaryAminoFactory"}): "thioamide",
    frozenset({"BoronicAcidFactory", "DiolFactory"}): "boronate ester",
    frozenset({"BromoFactory"}): "C–C",
    frozenset({"IodoFactory"}): "C–C",
    frozenset({"DibromoFactory"}): "C=C",
    frozenset({"TerminalAlkyneFactory"}): "diyne",
    frozenset({"TerminalAlkeneFactory"}): "C=C",
    frozenset({"ThiolFactory"}): "disulfide",
    frozenset({"AlcoholFactory"}): "ether",
    frozenset({"BoronicAcidFactory"}): "boroxine",
    frozenset({"BromoFactory", "IodoFactory"}): "C–C",
    frozenset({"ThiolFactory", "TerminalAlkeneFactory"}): "thioether",
}

#: Human-readable names, so messages do not read like source code.
DISPLAY_NAMES = {
    "AlcoholFactory": "Alcohol",
    "AldehydeFactory": "Aldehyde",
    "AmideFactory": "Amide",
    "BoronicAcidFactory": "Boronic acid",
    "BromoFactory": "Bromide",
    "CarboxylicAcidFactory": "Carboxylic acid",
    "DibromoFactory": "gem-Dibromide",
    "DifluoroFactory": "gem-Difluoride",
    "DiolFactory": "Diol",
    "FluoroFactory": "Fluoride",
    "IodoFactory": "Iodide",
    "PrimaryAminoFactory": "Primary amine",
    "RingAmineFactory": "Ring amine",
    "SecondaryAminoFactory": "Secondary amine",
    "TerminalAlkeneFactory": "Terminal alkene",
    "TerminalAlkyneFactory": "Terminal alkyne",
    "ThioacidFactory": "Thioacid",
    "ThiolFactory": "Thiol",
}


#: The dialogs offer "Auto" alongside the human-readable factory names.
AUTO_LABEL = "Auto"


def display(factory_name):
    return DISPLAY_NAMES.get(factory_name, factory_name)


def display_names():
    """Factory display names in dialog order."""
    return [DISPLAY_NAMES[name] for name in FACTORY_NAMES]


def factory_from_display(text):
    """Map a dialog string back to a factory name; None means 'Auto'."""
    if not text or text == AUTO_LABEL:
        return None
    for name, shown in DISPLAY_NAMES.items():
        if text == shown or text == name:
            return name
    raise ConversionError(
        f"Unknown functional group type '{text}'. Choose one of: "
        f"{', '.join(display_names())}."
    )


def linkage_name(factory_a, factory_b):
    """Name of the bond formed between two factories, or a generic fallback."""
    return REACTION_PAIRS.get(frozenset({factory_a, factory_b}), "new")


def reacts_with(factory_name):
    """Display names of the groups this one forms a recognised linkage with."""
    partners = []
    for pair in sorted(REACTION_PAIRS, key=lambda p: sorted(p)):
        if factory_name not in pair:
            continue
        other = next(iter(pair - {factory_name}), factory_name)
        partners.append(
            "themselves" if other == factory_name else display(other).lower()
        )
    return partners


def check_pair(factory_a, factory_b, label_a="your molecule",
               label_b="the partner"):
    """Refuse two functional groups that do not form a known linkage.

    stk joins any two functional groups it is given: it deletes each group's
    leaving atoms and bonds the remnants together. For a combination that is
    not real chemistry that silently yields a plausible-looking structure made
    of the wrong elements — a bromide reacted with an aldehyde loses every Br
    and every O and comes back a hydrocarbon — so the pairing is rejected here
    rather than built and shown to the user as a result.
    """
    if frozenset({factory_a, factory_b}) in REACTION_PAIRS:
        return

    options = reacts_with(factory_a)
    if options:
        hint = (
            f"{display(factory_a)} groups react with "
            f"{_listing(options)}. Choose a partner with one of those groups, "
            f"or pick different groups instead of Auto."
        )
    else:
        hint = (
            f"{display(factory_a)} has no partner group this plugin "
            f"recognises. Pick different groups instead of Auto."
        )

    same = factory_a == factory_b
    combination = (
        f"{display(factory_a)} groups ({label_a}) do not react with each other"
        if same else
        f"{display(factory_a)} ({label_a}) and {display(factory_b)} "
        f"({label_b}) do not form a known linkage"
    )
    raise ConversionError(
        f"{combination}. stk would still join them — it deletes the atoms of "
        f"both groups and bonds what is left — but the result would not be "
        f"the chemistry you asked for. {hint}"
    )


def _listing(items):
    """'a', 'a and b', 'a, b and c'."""
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


class FGMatch:
    """What one factory found in the molecule."""

    def __init__(self, factory, groups):
        self.factory = factory
        self.groups = list(groups)
        self.atom_ids = [tuple(sorted(g.get_atom_ids())) for g in self.groups]

    @property
    def count(self):
        return len(self.groups)

    @property
    def atom_set(self):
        return {atom for ids in self.atom_ids for atom in ids}


def detect(mol, factories=FACTORY_NAMES):
    """Run every factory over an RDKit molecule.

    Returns matches with at least one group, most groups first.
    """
    matches = []
    for name in factories:
        factory = getattr(stk, name, None)
        if factory is None:
            continue
        try:
            block = stk.BuildingBlock.init_from_rdkit_mol(
                mol, functional_groups=[factory()]
            )
        except Exception:  # noqa: BLE001 - a factory that cannot apply is not an error
            continue
        groups = list(block.get_functional_groups())
        if groups:
            matches.append(FGMatch(name, groups))
    matches.sort(key=lambda m: (-m.count, m.factory))
    return matches


def drop_subsumed(matches):
    """Remove factories whose matches are contained in a more specific one.

    A carboxylic acid also answers to ``AlcoholFactory`` and a gem-dibromide to
    ``BromoFactory``. Whenever every atom one factory claims is also claimed by
    another that matched at least as much of the molecule, the broader factory
    is noise and is dropped.
    """
    keep = []
    for match in matches:
        covered = False
        for other in matches:
            if other is match or other.count < match.count:
                continue
            if match.atom_set and match.atom_set < other.atom_set:
                covered = True
                break
        if not covered:
            keep.append(match)
    return keep


def choose(matches, partner_factory=None, minimum=2):
    """Pick the factory the user most likely meant.

    Returns ``(factory_name, note)``. Raises ConversionError when the choice is
    genuinely ambiguous, listing the candidates so the user can pick manually.
    """
    candidates = [m for m in drop_subsumed(matches) if m.count >= minimum]
    if not candidates:
        found = ", ".join(
            f"{display(m.factory)} ({m.count})" for m in matches
        ) or "none"
        raise ConversionError(
            f"No functional group appears at least {minimum} times in this "
            f"molecule, so it cannot be used as a building block. Groups "
            f"found: {found}. Add reactive groups (aldehydes, amines, "
            f"bromides, boronic acids…) or pick a factory manually."
        )

    if partner_factory:
        paired = [
            m for m in candidates
            if frozenset({m.factory, partner_factory}) in REACTION_PAIRS
        ]
        if len(paired) == 1:
            match = paired[0]
            return match.factory, (
                f"{display(match.factory)} chosen automatically: it forms "
                f"{linkage_name(match.factory, partner_factory)} bonds with "
                f"{display(partner_factory)}."
            )
        if len(paired) > 1:
            candidates = paired

    best = candidates[0]
    tied = [m for m in candidates if m.count == best.count]
    if len(tied) > 1:
        listing = ", ".join(f"{display(m.factory)} ({m.count})" for m in tied)
        raise ConversionError(
            f"Several functional groups are equally plausible here: {listing}. "
            f"Choose one explicitly instead of Auto."
        )
    return best.factory, (
        f"{display(best.factory)} chosen automatically ({best.count} groups)."
    )


def build(mol, factory_name, selected_atoms=None, minimum=1):
    """Make a BuildingBlock, optionally restricted to the user's selection.

    ``selected_atoms`` are Avogadro atom indices. When supplied, only the
    functional groups overlapping the selection react — this is how "join just
    the two bromines I highlighted" works, and it is the one thing a GUI can
    offer that a bare stk script cannot.
    """
    factory = getattr(stk, factory_name, None)
    if factory is None:
        raise ConversionError(f"Unknown functional group type '{factory_name}'.")

    block = stk.BuildingBlock.init_from_rdkit_mol(
        mol, functional_groups=[factory()]
    )
    groups = list(block.get_functional_groups())

    if selected_atoms:
        wanted = set(selected_atoms)
        groups = [g for g in groups if wanted & set(g.get_atom_ids())]
        if len(groups) < minimum:
            raise ConversionError(
                f"The selection touches {len(groups)} "
                f"{display(factory_name).lower()} group(s), but at least "
                f"{minimum} are needed. Select the atoms of the groups that "
                f"should react, or turn the selection option off."
            )
        block = stk.BuildingBlock.init_from_rdkit_mol(
            mol, functional_groups=groups
        )

    if block.get_num_functional_groups() < minimum:
        raise ConversionError(
            f"This molecule has {block.get_num_functional_groups()} "
            f"{display(factory_name).lower()} group(s); at least {minimum} "
            f"are needed."
        )
    return block


def from_smiles(smiles, factory_name, minimum=1, label="partner"):
    """Build a block from SMILES, embedding 3D coordinates as stk requires."""
    factory = getattr(stk, factory_name, None)
    if factory is None:
        raise ConversionError(f"Unknown functional group type '{factory_name}'.")
    try:
        block = stk.BuildingBlock(smiles, [factory()])
    except Exception as exc:  # noqa: BLE001 - stk/rdkit raise several types
        raise ConversionError(
            f"Could not build the {label} from SMILES '{smiles}' ({exc}). "
            f"Check the SMILES string."
        ) from None

    if block.get_num_functional_groups() < minimum:
        raise ConversionError(
            f"The {label} '{smiles}' has "
            f"{block.get_num_functional_groups()} "
            f"{display(factory_name).lower()} group(s); at least {minimum} "
            f"are needed."
        )
    return block
