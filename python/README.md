# pyrisearch-tauso

A Python API over the RIsearch binary that
[risearch-tauso](https://pypi.org/project/risearch-tauso/) ships.

```bash
pip install pyrisearch-tauso
```

That pulls in `risearch-tauso` for the binary and `pyarrow` for reading what it
prints. To get only the binary and its command line, install `risearch-tauso` on
its own.

## One search, buffered

```python
import pyrisearch_tauso as ris

result = ris.run(["-q", "queries.fa", "-t", "targets.fa", "-s", "900", "-p2"])
for line in result.stdout.splitlines():
    query, query_start, query_end, target, target_start, target_end, score, energy = line.split("\t")
```

A run that fails raises `RIsearchError`, carrying what RIsearch wrote to stderr,
so an empty result means no hits rather than a search that never happened.

## One search, streamed

RIsearch prints a line per hit and that count grows with the product of the
inputs: sixteen 20-nt queries against 4000 records of 600 nt at `-s 900` come to
16.4 million lines. `stream()` reads them as they are produced.

```python
import pyarrow as pa
import pyrisearch_tauso as ris

parts = []
with ris.stream(["-q", "queries.fa", "-t", "targets.fa", "-s", "900"],
                columns=("query", "target", "energy")) as batches:
    for batch in batches:
        parts.append(
            pa.Table.from_batches([batch])
            .group_by(["query", "target"]).aggregate([("energy", "min")])
        )
```

Reducing each batch as it arrives bounds memory by `block_size`; collecting the
batches into one table bounds it by the number of hits instead. On the 16.4
million hit run above, that is 689 MB against 2106 MB for the buffered `run()`,
in the same wall time -- RIsearch's own work dominates all of them.

`stream()` reads the eight-column `-p2` table and passes `-p2` itself, so leave
`-p` out of `args`.

## Which matrix

RIsearch scores with one of five matrices, named on `Matrix` so they can be
reached by name rather than remembered:

```python
ris.Matrix.T99         # RNA-RNA
ris.Matrix.T04         # RNA-RNA, and what RIsearch uses when not told otherwise
ris.Matrix.SU95        # RNA-DNA
ris.Matrix.SU95_NO_GU  # RNA-DNA, wobble pairs left out
ris.Matrix.SLH04_NO_GU # DNA-DNA, wobble pairs left out
```

Each is a string, so `matrix=ris.Matrix.SU95_NO_GU` and `matrix="su95_noGU"`
are the same thing. Which interaction each one scores is RIsearch's own, from
the table in its `src/dsm.h`.

## One search, reduced

Reducing each batch as it arrives is the shape a large search wants, and it has
a trap in it: the smallest energy within a batch is the smallest of that batch,
and the smallest overall is only known once the partials are put together. A
`Reduction` names both stages, and `search_reduced` runs them.

```python
import pyarrow as pa
import pyrisearch_tauso as ris

def combine(batch):                       # one batch -> a partial
    return (pa.Table.from_batches([batch])
            .group_by(["query", "target"]).aggregate([("energy", "min")])
            .rename_columns(["query", "target", "energy"]))

def finalize(partials):                   # the partials -> the answer
    return (partials.group_by(["query", "target"]).aggregate([("energy", "min")])
            .rename_columns(["query", "target", "energy"]))

best = ris.search_reduced(
    queries=queries, targets=targets, min_score=900,
    reduction=ris.Reduction(
        columns=("query", "target", "energy"),
        combine=combine, finalize=finalize, empty={},
    ),
)
```

The reduction names the columns, so the search parses only those.

## Every hit at once

```python
table = ris.hits_table(queries=queries, targets=targets, min_score=900)
```

The whole result is held, so this is for a search small enough to look at. A
search that found nothing gives back a table with the hit schema and no rows,
so its columns read the same either way.

## Energy statistics at several cutoffs

`energy_stats` supplies a reduction for the Boltzmann sum and minimum energy.
Pass RT explicitly in kcal/mol; the library does not choose a temperature.

```python
cutoffs = [800, 1000, 1200]
stats = ris.search_reduced(
    queries=queries, targets=targets,
    min_score=min(cutoffs), neighborhood=0,
    reduction=ris.energy_stats(cutoffs, group_by=("query", "target"), rt=0.616),
)
# stats[cutoff][(query_id, target_id)].sum_exp and .min_energy
```

Each cutoff uses `score > cutoff`. Use `group_by=("query",)` to combine targets
and obtain query IDs as keys. The sum is `sum(exp(-energy / rt))`; it is not a
normalized occupancy probability. No hits returns `{}`. Sums retain the
existing batch-then-final aggregation order and floating-point behavior.
Partial tables remain in memory until finalization, so memory also depends on
how many groups each batch produces.

## Reusing targets across searches

```python
with ris.fasta_targets(target_sequences) as target_path:
    for queries in query_chunks:
        stats = ris.search_reduced(
            queries=queries, targets=target_path, min_score=800,
            reduction=ris.energy_stats([800], rt=0.616),
        )
```

The context writes the targets once and removes the temporary FASTA on exit,
including exceptions. An existing path is borrowed without deleting it. All
workers using the path must finish before the context exits.
