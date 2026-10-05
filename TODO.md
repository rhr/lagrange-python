# TODO

Known problems found during the Python 3 port. These predate the port
(they affect the Python 2 code too) and were left unchanged.

## `lagrange.Nexus` does not work

`lagrange/nexus.py` is a trimmed 2005 copy of Biopython's NEXUS parser,
with its Biopython imports (`Seq`, `IUPAC`, `IUPACData`) commented out.
`Nexus.__init__` references `IUPAC`, so creating a `Nexus` object raises
`NameError`, and so does every character-matrix method that uses `Seq`.

Options:
- depend on Biopython and use the current `Bio.Nexus`, or
- cut `nexus.py` down to parsing only the taxa and trees blocks.

## Fractional dispersal constraints are truncated to 0

Commit `1d48c90` (2017-11-13, "make param arrays dtype=int") changed
`DECModel.Dmask` and `DECModel.Emask` to integer arrays.
`input.eval_decmodel` then does `model.Dmask[:] = dm`, which silently
truncates fractional values in `area_dispersal` (e.g. 0.5 becomes 0,
turning a reduced dispersal rate into no dispersal). Inputs where every
value is 0 or 1, like the Psychotria example, are unaffected.

If fractional dispersal matrices are supported, `Dmask` (and probably
`Emask`) should go back to `float`.

## Minor

- Debug `print` of `dp_array` in `DECModel.__init__`
  (`lagrange/decmodel_mp.py`) prints a line of zeros on every run.
- `lagrange/optimize_mp.py` is an older near-duplicate of
  `lagrange/optimize.py` and is not imported anywhere.
