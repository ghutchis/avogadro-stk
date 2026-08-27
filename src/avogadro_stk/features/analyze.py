#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Analyze Building Block — report what stk sees in the current molecule.

This is the feature that teaches the vocabulary: which groups were found, what
topicity that implies, and which topologies therefore become reachable. It
never modifies the molecule.
"""

from .. import blocks, validate
from ..convert import cjson_to_rdkit
from ..result import note
from ..topologies import TOPOLOGIES


def _group_table(matches, primary):
    lines = ["Detected functional groups:"]
    for match in matches:
        marker = "*" if match.factory == primary else " "
        sets = " / ".join(",".join(str(i) for i in ids) for ids in match.atom_ids[:6])
        if len(match.atom_ids) > 6:
            sets += ", …"
        lines.append(
            f" {marker} {blocks.display(match.factory):<18}"
            f"{match.count:>2} group{'s' if match.count != 1 else ' '}   "
            f"(atoms {sets})"
        )
    return lines


def _compatible(topicity):
    """Topologies this topicity can fill, grouped by kind."""
    found = {"cage": [], "cof": []}
    for key, entry in TOPOLOGIES.items():
        if entry["metal"] or entry["periodic"]:
            continue
        needs = set(entry["requirements"])
        if topicity not in needs:
            continue
        kind = entry["kind"]
        partner = sorted(needs - {topicity})
        name = key.split(".", 1)[1]
        found[kind].append((name, partner[0] if partner else topicity))
    return found


def _wrap(names, indent):
    """Wrap a comma-separated list to keep the message readable."""
    lines, current = [], ""
    for name in names:
        candidate = f"{current}, {name}" if current else name
        if len(candidate) > 62:
            lines.append(current)
            current = name
        else:
            current = candidate
    if current:
        lines.append(current)
    return [f"{indent}{line}" for line in lines]


def run(avo_input):
    mol = cjson_to_rdkit(avo_input.get("cjson", {}))
    matches = blocks.detect(mol)
    if not matches:
        return note(
            "No functional groups that stk recognises were found in this "
            "molecule.\n\nstk builds structures by joining reactive groups — "
            "aldehydes, amines, carboxylic acids, boronic acids, bromides, "
            "alkynes, thiols and so on. Add such groups, or use this molecule "
            "as a guest rather than a building block."
        )

    useful = blocks.drop_subsumed(matches)
    try:
        primary, reason = blocks.choose(matches)
        ambiguity = None
    except Exception as exc:  # noqa: BLE001 - ambiguity is normal here
        primary, reason, ambiguity = None, None, str(exc)

    lines = _group_table(matches, primary)

    dropped = [m for m in matches if m not in useful]
    if dropped:
        names = ", ".join(blocks.display(m.factory) for m in dropped)
        verb = "overlaps a more specific group and was" if len(dropped) == 1 \
            else "overlap more specific groups and were"
        lines.append(f"  ({names} {verb} set aside.)")

    if primary is None:
        lines += ["", "Topicity: ambiguous.", f"  {ambiguity}"]
        # Even without a single answer, the candidates all have a topicity, so
        # say what each one would make possible.
        for match in useful:
            if match.count >= 2:
                lines.append(
                    f"  As {'an' if blocks.display(match.factory)[0] in 'AEIOU' else 'a'} "
                    f"{blocks.display(match.factory).lower()} block it is "
                    f"{match.count}-topic."
                )
        return note("\n".join(lines))

    topicity = next(m.count for m in matches if m.factory == primary)
    lines += ["", f"Topicity: {topicity} (as "
                  f"{'an' if blocks.display(primary)[0] in 'AEIOU' else 'a'} "
                  f"{blocks.display(primary).lower()} building block)",
              f"  {reason}", "", "Compatible topologies:"]

    compatible = _compatible(topicity)
    for kind, heading in (("cage", "Cage"), ("cof", "COF")):
        items = compatible[kind]
        if not items:
            continue
        by_partner = {}
        for name, partner in items:
            by_partner.setdefault(partner, []).append(name)
        for partner, names in sorted(by_partner.items()):
            partner_text = (
                f"with a {partner}-topic partner" if partner != topicity
                else "on its own"
            )
            lines.append(f"  {heading}, {partner_text}:")
            lines += _wrap(sorted(names), "    ")

    if topicity == 2:
        lines.append("  Polymer: yes — this molecule is ditopic.")
        lines.append("  Macrocycle: yes.")
    else:
        lines.append(
            f"  Polymer: no (needs exactly 2 groups, this has {topicity})."
        )

    partners = validate.suggestions(3) if topicity == 2 else validate.suggestions(2)
    if partners:
        lines += ["", "Partner suggestions from the library: "
                      + ", ".join(partners) + "."]

    return note("\n".join(lines))
