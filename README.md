# avogadro-stk

Build polymers, macrocycles, cages and covalent organic frameworks inside
[Avogadro 2](https://two.avogadro.cc/) using [`stk`](https://stk.readthedocs.io).

The molecule you have open becomes one building block; the other comes from a
bundled library of curated linkers, or from a SMILES string you type. stk works
out the connectivity exactly, so the result arrives in Avogadro with correct
bonds rather than bonds guessed from distances.

## Features

All under **Build ▸ Supramolecular**:

| Menu item | What it does |
| --- | --- |
| Analyze Building Block… | Reports the functional groups stk can see, the resulting topicity, and which topologies that makes reachable. Changes nothing. |
| Polymer… | `stk.polymer.Linear` — a chain from one or two monomers, optionally end-capped. |
| Macrocycle… | `stk.macrocycle.Macrocycle` — the same chain, closed into a ring. |
| Cage… | `stk.cage.*` — 18 organic cage topologies. |
| Covalent Framework… | `stk.cof.*` — honeycomb, hexagonal, kagome, square and linkerless honeycomb, finite or periodic. |

**Start with Analyze.** It names the vocabulary the other dialogs use
(topicity, factory, topology) against the molecule actually in front of you.

## Geometry: read this before using a structure

stk places building blocks on idealized vertices. The **Clean Up Geometry**
dropdown decides what happens next:

| Setting | What it does |
|---------|--------------|
| None | stk's raw vertex geometry. Bond lengths across the new bonds are meaningless. |
| MCHammer | `stk.MCHammer()` collapses the structure to sane bond lengths. Not a force field: angles, torsions and contacts are still wrong. |
| MCHammer + MMFF94 (default) | The above, then RDKit's MMFF94 over the whole assembly. |

MMFF94 has no parameters for boron, so boronate-ester and boroxine structures
fall back to **UFF** automatically; the result message says which ran. The
force field is skipped for periodic cells (RDKit has no periodic boundary and
would pull the structure out of the unit cell) and above 2000 atoms.

Even with MMFF94 this is a *starting* geometry. Run a method suited to your
chemistry (**Extensions ▸ Optimize Geometry**, or one of the xtb / MACE / ANI
plugins) before reporting anything. Every result message says how good the
geometry it carries actually is.

## Functional groups

stk joins molecules at *functional groups*, and needs to be told which ones.
Leaving the dropdown on **Auto** runs every stk factory over your molecule and
picks a group by:

1. discarding anything appearing fewer than twice;
2. discarding groups that overlap a more specific match (a carboxylic acid also
   answers to the alcohol factory);
3. preferring the group that forms a known linkage with the partner block —
   aldehyde + primary amine → imine, boronic acid + diol → boronate ester, and
   so on;
4. failing with the candidates listed rather than guessing, if that still ties.

**Use Selected Atoms Only** restricts the reaction to the groups touching your
current selection. This is how you use one of the three aldehydes on a tritopic
molecule, or join just the two bromines you highlighted — the one thing this
plugin can offer that a bare stk script cannot.

## Topicity has to match the topology

Every topology fixes how many blocks of each connectivity it needs;
`FourPlusSix` is 4 tritopic + 6 ditopic. Mismatches are caught before
construction and reported with the numbers involved and a concrete suggestion:

> `FourPlusSix` needs 6 × 2-topic and 4 × 3-topic building blocks, but your
> molecule has 2 reactive groups and p-Phenylenediamine has 2. Nothing supplied
> is 3-topic. Pick a 3-topic partner such as 1,3,5-Triformylbenzene, …

The library dropdown shows each entry's topicity for the same reason.

## Periodic frameworks are experimental

Non-periodic is the default and gives a finite fragment of the framework, with
unreacted groups around the edge — an ordinary molecule Avogadro handles
normally.

Ticking **Periodic** emits a `unitCell` instead. Two caveats, both reported in
the result message:

- The cell parameters come from stk's idealized vertex placement, not a relaxed
  structure, and are correspondingly generous.
- Bonds that wrap across the cell boundary are drawn by Avogadro as long bonds
  spanning the cell. This is cosmetic but looks wrong.

Optimize with a periodic method before trusting the cell.

## Not yet supported

- **Metal-organic cages** (`M4L6Tetrahedron`, `M12L24`, …) are recognised but
  refused: they need dative bonds to metal centres, which Avogadro does not
  represent well enough yet.
- **Host–guest complexes and rotaxanes** need a second *constructed* molecule,
  so they wait on a way to stash one.

## Install size

`stk` pulls in rdkit, scipy, pandas, matplotlib, seaborn and polars, so a fresh
environment is roughly **500 MB** — polars alone is ~180 MB and rdkit ~100 MB.
Most of that is stk's database and analysis tooling, which this plugin never
touches, but the dependencies are not optional upstream. `import stk` costs
about 1.5 s, so expect a short pause the first time a command runs.

## Development

```bash
pixi install
pixi run test          # pytest tests
pixi run regen         # regenerate topologies.py and the options JSON
```

Two files are generated and should not be hand-edited:

- `src/avogadro_stk/topologies.py` — every topology's vertex requirements,
  introspected from stk by `tools/generate_topologies.py`. `tests/test_topologies.py`
  rebuilds all 44 with dummy blocks, so an stk change that invalidates the
  table fails loudly.
- `src/avogadro_stk/customization/*_options.json` — the dialogs, generated from
  `data/linkers.json` by `tools/generate_options.py`. Avogadro reads option
  lists from static files, so they cannot be built at run time.

To add a linker, edit `src/avogadro_stk/data/linkers.json` and run
`pixi run regen`. `tests/test_library.py` checks every SMILES parses and that
the named factory finds exactly the declared number of groups.

## License

3-Clause BSD, matching Avogadro.
