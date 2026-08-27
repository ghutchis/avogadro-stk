#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""The bundled building block library.

Avogadro's option dialogs are built from static JSON, so the dropdown entries
have to be plain strings. Those strings are generated from ``linkers.json`` by
``tools/generate_options.py`` and resolved back to library entries here;
``tests/test_library.py`` checks the two stay in step.
"""

import functools
import json
import pathlib

from . import blocks
from .convert import ConversionError

DATA = pathlib.Path(__file__).resolve().parent / "data" / "linkers.json"

NONE_LABEL = "(none)"
CUSTOM_LABEL = "Custom SMILES…"


@functools.lru_cache(maxsize=1)
def entries():
    """Every library entry, in file order."""
    with DATA.open(encoding="utf-8") as handle:
        return tuple(json.load(handle)["linkers"])


def label(entry):
    """Dropdown text: name, topicity and the group that reacts.

    Showing the topicity next to the name is what makes a topology mismatch
    predictable rather than a surprise at run time.
    """
    if not entry.get("factory"):
        return f"{entry['name']} (guest)"
    return (
        f"{entry['name']} ({entry['topicity']}-topic, "
        f"{blocks.display(entry['factory'])})"
    )


def labels(tags=None, topicity=None):
    """Dropdown values, optionally filtered by tag and topicity."""
    result = []
    for entry in entries():
        if tags and not (set(tags) & set(entry.get("tags", []))):
            continue
        if topicity is not None and entry["topicity"] != topicity:
            continue
        result.append(label(entry))
    return result


@functools.lru_cache(maxsize=1)
def _by_label():
    return {label(entry): entry for entry in entries()}


def lookup(text):
    """Resolve a dropdown string, an id or an alias to a library entry."""
    if not text:
        return None
    entry = _by_label().get(text)
    if entry is not None:
        return entry

    needle = text.strip().casefold()
    for entry in entries():
        names = [entry["id"], entry["name"], *entry.get("aliases", [])]
        if any(needle == name.casefold() for name in names):
            return entry
    return None


def block(text, minimum=1, role="partner"):
    """Build the stk BuildingBlock for a library entry named by dropdown text."""
    entry = lookup(text)
    if entry is None:
        raise ConversionError(
            f"'{text}' is not in the building block library. Choose an entry "
            f"from the list or use {CUSTOM_LABEL} to enter your own."
        )
    if not entry.get("factory"):
        raise ConversionError(
            f"'{entry['name']}' has no reactive groups — it is a guest "
            f"molecule, not a linker."
        )
    return blocks.from_smiles(
        entry["smiles"],
        entry["factory"],
        minimum=minimum,
        label=f"{role} ({entry['name']})",
    )
