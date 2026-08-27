#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Regenerate ``src/avogadro_stk/topologies.py`` by introspecting stk.

Every stk topology fixes how many functional groups each of its vertices
needs. Transcribing that by hand is how the table silently goes stale, so it
is derived here from the ``_vertex_prototypes`` / ``_edge_prototypes`` class
attributes and written out as a plain dict.

Run after upgrading stk::

    python tools/generate_topologies.py

``tests/test_topologies.py`` builds every entry with dummy blocks, so a stk
change that invalidates the table fails loudly rather than quietly.
"""

import collections
import inspect
import pathlib
import pprint

import stk

HEADER = '''#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Topology requirements, generated from stk {version} by tools/generate_topologies.py.

Do not edit by hand. Each entry records how many building-block vertices of
each functionality (number of functional groups) the topology needs, so a
mismatch can be reported before construction is attempted.
"""

STK_VERSION = "{version}"

# name -> dict(kind, vertices, edges, requirements, metal, periodic)
# `requirements` maps functionality -> number of vertices needing it.
TOPOLOGIES = {table}
'''


def requirements(cls):
    """Map functionality -> vertex count, from the topology's prototypes."""
    vertices = getattr(cls, "_vertex_prototypes", None)
    edges = getattr(cls, "_edge_prototypes", None)
    if not vertices or not edges:
        return None

    degree = collections.Counter()
    for edge in edges:
        degree[edge.get_vertex1_id()] += 1
        degree[edge.get_vertex2_id()] += 1

    counts = collections.Counter()
    for vertex in vertices:
        counts[degree.get(vertex.get_id(), 0)] += 1
    return dict(sorted(counts.items()))


def collect():
    table = {}
    for kind, module, base in (
        ("cage", stk.cage, stk.cage.Cage),
        ("cof", stk.cof, stk.cof.Cof),
    ):
        for name in sorted(dir(module)):
            if name.startswith("_"):
                continue
            cls = getattr(module, name)
            if not inspect.isclass(cls) or cls is base or not issubclass(cls, base):
                continue
            needs = requirements(cls)
            if not needs:
                continue
            # Metal-organic cages are named M<n>L<n>; they need dative bonds
            # and are not offered yet (see README).
            metal = name.startswith("M") and any(c.isdigit() for c in name)
            table[f"{kind}.{name}"] = {
                "kind": kind,
                "vertices": len(cls._vertex_prototypes),
                "edges": len(cls._edge_prototypes),
                "requirements": needs,
                "metal": metal,
                "periodic": name.startswith("Periodic"),
            }
    return table


def main():
    table = collect()
    target = (
        pathlib.Path(__file__).resolve().parent.parent
        / "src" / "avogadro_stk" / "topologies.py"
    )
    target.write_text(
        HEADER.format(
            version=stk.__version__,
            table=pprint.pformat(table, indent=0, width=79, sort_dicts=True),
        )
    )
    print(f"wrote {len(table)} topologies to {target}")


if __name__ == "__main__":
    main()
