#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""Conversion between Avogadro's CJSON, RDKit and stk.

Avogadro hands us the whole molecule as CJSON, which already carries exact
connectivity. That connectivity is rebuilt into an RDKit molecule rather than
re-perceived from distances, and the constructed molecule is written back as
CJSON the same way, so bonds survive the round trip untouched.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

# CJSON bond orders are integers; stk hands back floats.
_ORDER_TO_RDKIT = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
    4: Chem.BondType.QUADRUPLE,
}


class ConversionError(Exception):
    """A structure problem worth reporting to the user in plain language."""


def _flat_to_triples(flat):
    return [tuple(flat[i:i + 3]) for i in range(0, len(flat) - 2, 3)]


def cjson_to_rdkit(cjson):
    """Build an RDKit molecule from Avogadro's CJSON.

    Raises ConversionError with a message aimed at the user, never a traceback.
    """
    atoms = (cjson or {}).get("atoms", {}) or {}
    numbers = (atoms.get("elements", {}) or {}).get("number", []) or []
    coords = _flat_to_triples((atoms.get("coords", {}) or {}).get("3d", []) or [])

    if not numbers:
        raise ConversionError(
            "There is no molecule open. Draw or open a structure first, then "
            "run this command again."
        )
    if len(coords) < len(numbers):
        raise ConversionError(
            "The molecule has no 3D coordinates. Use Build ▸ Optimize Geometry "
            "to generate a 3D structure first."
        )

    mol = Chem.RWMol()
    charges = atoms.get("formalCharges", []) or []
    for i, number in enumerate(numbers):
        atom = Chem.Atom(int(number))
        if i < len(charges):
            atom.SetFormalCharge(int(charges[i]))
        # Avogadro's connectivity is authoritative; stop RDKit inventing H.
        atom.SetNoImplicit(True)
        mol.AddAtom(atom)

    bonds = (cjson or {}).get("bonds", {}) or {}
    index = (bonds.get("connections", {}) or {}).get("index", []) or []
    orders = bonds.get("order", []) or []
    for pair in range(len(index) // 2):
        begin, end = int(index[2 * pair]), int(index[2 * pair + 1])
        if begin == end or not (0 <= begin < len(numbers) and 0 <= end < len(numbers)):
            continue
        order = int(orders[pair]) if pair < len(orders) else 1
        mol.AddBond(begin, end, _ORDER_TO_RDKIT.get(order, Chem.BondType.SINGLE))

    conformer = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y, z) in enumerate(coords[:mol.GetNumAtoms()]):
        conformer.SetAtomPosition(i, (float(x), float(y), float(z)))
    mol.AddConformer(conformer, assignId=True)

    mol = mol.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # noqa: BLE001 - RDKit raises a family of types
        raise ConversionError(
            "RDKit could not make chemical sense of this structure "
            f"({exc}). Check for missing hydrogens (Build ▸ Add Hydrogens) or "
            "an atom with the wrong number of bonds."
        ) from None

    if _is_flat(coords):
        mol = _embed(mol)
    return mol


def _is_flat(coords, tolerance=1e-3):
    """True when every atom shares a z coordinate, i.e. a 2D sketch."""
    if len(coords) < 3:
        return False
    zs = [c[2] for c in coords]
    return max(zs) - min(zs) < tolerance


def _embed(mol):
    """Give a 2D sketch a 3D conformer; stk needs a real position matrix."""
    mol = Chem.AddHs(mol, addCoords=True)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ConversionError(
            "The molecule appears to be 2D and RDKit could not generate 3D "
            "coordinates for it. Try Build ▸ Optimize Geometry first."
        )
    return mol


def mol_from_smiles(smiles, label="molecule"):
    """A 3D RDKit molecule from SMILES, for functional group detection.

    Embedding once here is much cheaper than letting every candidate factory
    build its own stk BuildingBlock from the SMILES string.
    """
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        raise ConversionError(
            f"Could not parse the {label} SMILES '{smiles}'. Check the string."
        )
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ConversionError(
            f"RDKit could not generate 3D coordinates for the {label} "
            f"'{smiles}'."
        )
    return mol


def constructed_to_cjson(constructed, periodic_info=None):
    """Write an stk ConstructedMolecule out as a CJSON object.

    stk knows the exact connectivity of what it just built, so it is copied
    across verbatim. Avogadro is told not to re-perceive bonds (see the
    ``bond`` key in the result), because distance-based perception on a
    freshly collapsed structure invents and destroys bonds.
    """
    numbers = [atom.get_atomic_number() for atom in constructed.get_atoms()]
    positions = constructed.get_position_matrix()

    coords = []
    for row in positions:
        coords.extend(float(value) for value in row)

    connections, orders = [], []
    for bond in constructed.get_bonds():
        connections.append(bond.get_atom1().get_id())
        connections.append(bond.get_atom2().get_id())
        orders.append(int(round(bond.get_order())))

    cjson = {
        "chemicalJson": 1,
        "atoms": {
            "elements": {"number": numbers},
            "coords": {"3d": coords},
        },
        "bonds": {
            "connections": {"index": connections},
            "order": orders,
        },
    }

    if periodic_info is not None:
        cjson["unitCell"] = {
            "a": float(periodic_info.get_a()),
            "b": float(periodic_info.get_b()),
            "c": float(periodic_info.get_c()),
            "alpha": float(periodic_info.get_alpha()),
            "beta": float(periodic_info.get_beta()),
            "gamma": float(periodic_info.get_gamma()),
        }
    return cjson


def formula(constructed):
    """Hill-notation formula, for the result message."""
    table = Chem.GetPeriodicTable()
    counts = {}
    for atom in constructed.get_atoms():
        symbol = table.GetElementSymbol(atom.get_atomic_number())
        counts[symbol] = counts.get(symbol, 0) + 1

    def term(symbol, count):
        return f"{symbol}{count}" if count != 1 else symbol

    # Hill notation: carbon, then hydrogen, then everything else alphabetically.
    ordered = []
    for symbol in ("C", "H"):
        if symbol in counts:
            ordered.append(term(symbol, counts.pop(symbol)))
    ordered.extend(term(symbol, counts[symbol]) for symbol in sorted(counts))
    return "".join(ordered)
