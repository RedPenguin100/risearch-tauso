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
    qname, qbeg, qend, tname, tbeg, tend, score, energy = line.split("\t")
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
                columns=("qname", "tname", "energy")) as batches:
    for batch in batches:
        parts.append(
            pa.Table.from_batches([batch])
            .group_by(["qname", "tname"]).aggregate([("energy", "min")])
        )
```

Reducing each batch as it arrives bounds memory by `block_size`; collecting the
batches into one table bounds it by the number of hits instead. On the 16.4
million hit run above, that is 689 MB against 2106 MB for the buffered `run()`,
in the same wall time -- RIsearch's own work dominates all of them.

`stream()` reads the eight-column `-p2` table and passes `-p2` itself, so leave
`-p` out of `args`.

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
            .group_by(["qname", "tname"]).aggregate([("energy", "min")])
            .rename_columns(["qname", "tname", "energy"]))

def finalize(partials):                   # the partials -> the answer
    return (partials.group_by(["qname", "tname"]).aggregate([("energy", "min")])
            .rename_columns(["qname", "tname", "energy"]))

best = ris.search_reduced(
    queries=queries, targets=targets, min_score=900,
    reduction=ris.Reduction(
        columns=("qname", "tname", "energy"),
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
