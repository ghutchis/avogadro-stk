#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Topicity checks against the generated topology table.

A bare "incompatible topicity" is useless to someone meeting these topologies
for the first time, so every rejection here names the numbers involved and
points at something that would actually work.
"""

from . import library
from .convert import ConversionError
from .topologies import TOPOLOGIES

#: Refuse outright above this many atoms; warn above the soft limit.
HARD_ATOM_LIMIT = 10000
SOFT_ATOM_LIMIT = 2000


def available(kind, periodic=None, metal=False):
    """Topology names of one kind, in menu order."""
    names = []
    for name, spec in TOPOLOGIES.items():
        if spec["kind"] != kind or spec["metal"] != metal:
            continue
        if periodic is not None and spec["periodic"] != periodic:
            continue
        names.append(name.split(".", 1)[1])
    return sorted(names)


def spec(kind, name):
    key = f"{kind}.{name}"
    if key not in TOPOLOGIES:
        raise ConversionError(
            f"Unknown {kind} topology '{name}'. Available: "
            f"{', '.join(available(kind))}."
        )
    return TOPOLOGIES[key]


def describe(kind, name):
    """e.g. '4 × 3-topic and 6 × 2-topic building blocks'."""
    needs = spec(kind, name)["requirements"]
    parts = [f"{count} × {functionality}-topic"
             for functionality, count in sorted(needs.items())]
    return " and ".join(parts) + " building blocks"


def suggestions(topicity, exclude=()):
    """Library entries with a given topicity, for 'try one of these' messages."""
    names = [
        entry["name"] for entry in library.entries()
        if entry["topicity"] == topicity and entry["id"] not in exclude
    ]
    return names[:3]


def check(kind, name, topicity_a, topicity_b, label_a="your molecule",
          label_b="the partner"):
    """Confirm two blocks can fill a topology's vertices.

    Returns the requirement dict on success.
    """
    needs = spec(kind, name)["requirements"]
    wanted = set(needs)
    have = {topicity_a} if topicity_b is None else {topicity_a, topicity_b}

    if wanted == have:
        return needs

    required = describe(kind, name)
    detail = (
        f"{label_a} has {topicity_a} reactive group"
        f"{'' if topicity_a == 1 else 's'}"
    )
    if topicity_b is not None:
        detail += (
            f" and {label_b} has {topicity_b}"
        )

    hint = _hint(kind, name, wanted, have, topicity_a)
    raise ConversionError(f"{name} needs {required}, but {detail}. {hint}")


def _hint(kind, name, wanted, have, topicity_a):
    """The actionable half of a mismatch message."""
    missing = sorted(wanted - have)
    if missing:
        options = suggestions(missing[0])
        listing = ", ".join(options) if options else "a suitable linker"
        return (
            f"Nothing supplied is {missing[0]}-topic. Pick a "
            f"{missing[0]}-topic partner such as {listing}."
        )

    # Everything supplied is the wrong shape for this topology; point at one
    # that does fit what the user already has.
    fits = [
        key.split(".", 1)[1] for key, entry in TOPOLOGIES.items()
        if entry["kind"] == kind and not entry["metal"]
        and set(entry["requirements"]) == have
    ]
    if fits:
        return f"With these blocks try {', '.join(sorted(fits)[:4])} instead."
    if have == {2}:
        return "With two ditopic blocks try Macrocycle or Polymer instead."
    return "Choose a different topology or a different partner block."


def estimate_atoms(kind, name, atoms_by_topicity, lattice=(1, 1, 1)):
    """Rough atom count before construction, so huge jobs fail fast.

    Reactions delete atoms, so this overestimates a little; that is the safe
    direction for a guard.
    """
    needs = spec(kind, name)["requirements"]
    total = sum(
        count * atoms_by_topicity.get(functionality, 0)
        for functionality, count in needs.items()
    )
    if kind == "cof":
        total *= max(1, lattice[0]) * max(1, lattice[1]) * max(1, lattice[2])
    return total


def check_size(count, what="structure"):
    """Refuse absurd jobs; warn about merely large ones. Returns a warning."""
    if count > HARD_ATOM_LIMIT:
        raise ConversionError(
            f"This {what} would contain roughly {count} atoms, past the "
            f"{HARD_ATOM_LIMIT}-atom limit for this plugin. Use a smaller "
            f"lattice, fewer repeats, or a smaller topology."
        )
    if count > SOFT_ATOM_LIMIT:
        return (
            f"This is a large structure (~{count} atoms); construction and "
            f"rendering may be slow."
        )
    return None
