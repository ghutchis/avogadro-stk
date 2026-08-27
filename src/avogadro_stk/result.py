#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Shaping the JSON Avogadro expects back on stdout."""

#: Every result says this. MCHammer removes the worst of the idealised vertex
#: geometry and nothing more, and a user who does not know that will take the
#: bond lengths seriously.
CLEANUP_NOTE = (
    "Geometry was cleaned up with MCHammer only — run a force field or "
    "semi-empirical optimization (Extensions ▸ Optimize Geometry) before "
    "using this structure."
)


def molecule(cjson, message=None):
    """A new structure replacing what is open.

    ``bond`` is stated explicitly: stk supplied exact connectivity, and
    distance-based re-perception on a freshly collapsed structure would both
    invent and destroy bonds.
    """
    payload = {
        "moleculeFormat": "cjson",
        "cjson": cjson,
        "bond": False,
    }
    if message:
        payload["message"] = message
    return payload


def note(message):
    """Text only — no change to the molecule that is open."""
    return {"message": message}


def failure(message):
    """Report a problem without disturbing the open molecule."""
    return {"error": message}


def summarize(lines, warnings=(), cleanup=True):
    """Assemble a result message from parts, dropping the empty ones."""
    parts = [line for line in lines if line]
    parts.extend(f"Warning: {w}" for w in warnings if w)
    if cleanup:
        parts.append(CLEANUP_NOTE)
    return "\n\n".join(parts)
