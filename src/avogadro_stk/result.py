#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Shaping the JSON Avogadro expects back on stdout."""



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


def summarize(lines, warnings=(), geometry_note=None):
    """Assemble a result message from parts, dropping the empty ones.

    ``geometry_note`` goes last and every result carries one: whatever was or
    was not run, a user who does not know how good this geometry is will take
    the bond lengths seriously. See ``geometry.note``.
    """
    parts = [line for line in lines if line]
    parts.extend(f"Warning: {w}" for w in warnings if w)
    if geometry_note:
        parts.append(geometry_note)
    return "\n\n".join(parts)
