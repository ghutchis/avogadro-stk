#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-stk/

"""What happens to the geometry after stk has assembled it.

stk places building blocks on idealised vertices and leaves them there.
MCHammer pulls that back to something with sane bond lengths across the new
bonds, but it is a collapser, not a force field: angles, torsions and any
non-bonded contact it did not have to fix are still wrong. Running RDKit's
MMFF94 over the result is cheap — RDKit is already a hard dependency — and
turns the output into something worth looking at.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

from . import validate

#: Dialog wording, in menu order. Also accepted as ``cleanup`` values.
NONE = "None"
COLLAPSE = "MCHammer"
REFINE = "MCHammer + MMFF94"
LEVELS = (NONE, COLLAPSE, REFINE)

#: Iteration cap. MMFF on a fresh cage converges well inside this; the cap is
#: here so a pathological structure cannot hang the Avogadro UI, which blocks
#: on this process.
MAX_ITERATIONS = 2000

#: Above this many atoms the force field is skipped rather than run. Avogadro
#: has no progress bar for this and its own optimizer is a better place to
#: spend the minutes.
ATOM_LIMIT = validate.SOFT_ATOM_LIMIT


class Plan:
    """The cleanup the user asked for, and what actually happened."""

    def __init__(self, level):
        self.level = level
        #: Filled in by refine() with the force field that really ran.
        self.applied = COLLAPSE if level != NONE else "none"

    @property
    def collapse(self):
        return self.level != NONE

    @property
    def wants_force_field(self):
        return self.level == REFINE


def plan(options):
    """Read the ``cleanup`` dialog answer.

    Older dialogs (and the option's own history) express this as a boolean, so
    both spellings are accepted.
    """
    value = options.get("cleanup", REFINE)
    if isinstance(value, bool):
        return Plan(COLLAPSE if value else NONE)
    text = (str(value) if value is not None else "").strip()
    for level in LEVELS:
        if text.lower() == level.lower():
            return Plan(level)
    # An unrecognised string is not worth failing a build over.
    return Plan(REFINE)


def optimizer(current):
    """The stk optimizer for a plan. Imported late to keep stk out of import
    order for callers that only want the option parsing."""
    import stk

    return stk.MCHammer() if current.collapse else stk.NullOptimizer()


def refine(constructed, current, periodic=False):
    """Optimize a constructed molecule with MMFF94, falling back to UFF.

    Returns ``(molecule, warning)``. The molecule is a clone carrying the new
    positions, or the original when nothing was run — stk molecules are
    immutable, so ``with_position_matrix`` is the only way to update them.
    A force field that cannot run is never fatal: the MCHammer geometry is
    still usable, so the reason is reported as a warning instead.
    """
    if not current.wants_force_field:
        return constructed, None

    if periodic:
        return constructed, (
            "The force field was skipped: RDKit optimizes an isolated "
            "molecule with no periodic boundary, which would pull the "
            "structure away from the unit cell. Use a periodic method to "
            "relax this framework."
        )

    count = constructed.get_num_atoms()
    if count > ATOM_LIMIT:
        return constructed, (
            f"The force field was skipped: {count} atoms is past the "
            f"{ATOM_LIMIT}-atom limit for optimizing inside this plugin. Run "
            f"Extensions ▸ Optimize Geometry on the result instead."
        )

    mol = constructed.to_rdkit_mol()
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # noqa: BLE001 - RDKit raises a family of types
        return constructed, (
            f"The force field was skipped: RDKit could not sanitize the "
            f"constructed molecule ({exc}). The geometry is MCHammer's."
        )

    field, converged = _optimize(mol)
    if field is None:
        return constructed, (
            "The force field was skipped: neither MMFF94 nor UFF has "
            "parameters for every element in this structure. The geometry is "
            "MCHammer's."
        )

    current.applied = f"MCHammer+{field}"
    warning = None
    if not converged:
        warning = (
            f"{field} did not converge in {MAX_ITERATIONS} steps; the "
            f"geometry is improved but not at a minimum."
        )
    return (
        constructed.with_position_matrix(
            position_matrix=mol.GetConformer().GetPositions(),
        ),
        warning,
    )


def _optimize(mol):
    """Run MMFF94 if it is parameterised, otherwise UFF.

    MMFF94 has no boron, which this plugin builds routinely — boronate esters
    and boroxine frameworks are two of the linkages in the table — so the UFF
    fallback is load-bearing rather than defensive.
    """
    if AllChem.MMFFHasAllMoleculeParams(mol):
        result = AllChem.MMFFOptimizeMolecule(mol, maxIters=MAX_ITERATIONS)
        return "MMFF94", result == 0
    if AllChem.UFFHasAllMoleculeParams(mol):
        result = AllChem.UFFOptimizeMolecule(mol, maxIters=MAX_ITERATIONS)
        return "UFF", result == 0
    return None, False


def note(current):
    """The standing caveat for the geometry a result carries."""
    if not current.collapse:
        return (
            "No cleanup was run, so this is stk's idealised vertex geometry — "
            "bond lengths across the new bonds are meaningless. Optimize "
            "before using this structure."
        )
    if current.applied.endswith("MMFF94") or current.applied.endswith("UFF"):
        field = current.applied.split("+", 1)[1]
        return (
            f"Geometry was collapsed with MCHammer and optimized with "
            f"{field}. That is a starting geometry, not a final one — use a "
            f"method suited to your chemistry before reporting results."
        )
    return (
        "Geometry was cleaned up with MCHammer only — run a force field or "
        "semi-empirical optimization (Extensions ▸ Optimize Geometry) before "
        "using this structure."
    )
