# SIQS documentation graphs

These are the stable presentation assets used by the SIQS documentation
in `../README.md`, `../EVOLUTION.md`, and `../BENCHMARKS.md`.
They are generated benchmark output and should not be edited by hand.

The graphs come from two presentation workspaces, with copied plotting
code and frozen data inputs. The original benchmark collections and
generators remain unchanged.

| Graphs | Workspace relative to the repository root |
| :--- | :--- |
| Evolution and AMD comparisons | `misc/siqs-readme-20261020` |
| Apple single-thread and eight-thread comparisons | `misc/siqs-readme-20261027-upper` |

These development workspaces are not included in the repository or the
Perl distribution. The commands below require a checkout containing them
and should be run from the repository root.

To replot and publish just the four current Apple graphs:

```sh
python3 misc/siqs-readme-20261027-upper/render.py --publish
```

To replot and publish all nine graphs, run both commands in this order:

```sh
python3 misc/siqs-readme-20261020/refresh.py --publish
python3 misc/siqs-readme-20261027-upper/render.py --publish
```

The older publisher also copies its Apple graphs, so the newer Apple
publisher must run last to restore the current audited selection.
These commands replot the saved data; they do not collect new measurements.
FLINT and YAFU collection at 102 digits remains deferred.

Previous assets are retained in the presentation workspaces.
[provenance.json](provenance.json) records the source workspace, source
file, and SHA-256 of each published SVG.
