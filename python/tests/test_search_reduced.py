import pyarrow as pa
import pytest

import pyrisearch_tauso

QUERY_NAME = "q1"
TARGET_NAME = "t1"
QUERY = "GCTAGCTAGCTAGCTAGCTA"
TARGET = "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
QUERIES = {QUERY_NAME: QUERY}
TARGETS = {TARGET_NAME: TARGET}

MIN_SCORE = 500
NOTHING_REACHES_THIS = 100000


def best_energy_per_pair():
    """The smallest energy for each query and target, in two stages."""

    def combine(batch):
        return (
            pa.Table.from_batches([batch])
            .group_by(["query", "target"])
            .aggregate([("energy", "min")])
            .rename_columns(["query", "target", "energy"])
        )

    def finalize(table):
        merged = (
            table.group_by(["query", "target"])
            .aggregate([("energy", "min")])
            .rename_columns(["query", "target", "energy"])
        )
        return {
            (q, t): e
            for q, t, e in zip(
                merged.column("query").to_pylist(),
                merged.column("target").to_pylist(),
                merged.column("energy").to_pylist(),
            )
        }

    return pyrisearch_tauso.Reduction(
        columns=("query", "target", "energy"), combine=combine, finalize=finalize, empty={}
    )


def test_the_reduction_answers_over_the_whole_search():
    reduced = pyrisearch_tauso.search_reduced(
        queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE, reduction=best_energy_per_pair()
    )
    table = pyrisearch_tauso.hits_table(
        queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE
    )

    assert reduced == {(QUERY_NAME, TARGET_NAME): min(table.column("energy").to_pylist())}


def test_no_hits_gives_back_what_empty_says():
    reduced = pyrisearch_tauso.search_reduced(
        queries=QUERIES, targets=TARGETS, min_score=NOTHING_REACHES_THIS,
        reduction=best_energy_per_pair(),
    )

    assert reduced == {}


def test_the_reduction_names_the_columns():
    """Passing columns twice is a mistake worth refusing rather than resolving."""
    with pytest.raises(TypeError, match="off the reduction"):
        pyrisearch_tauso.search_reduced(
            queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE,
            reduction=best_energy_per_pair(), columns=("query",),
        )


def test_only_the_columns_a_reduction_reads_are_parsed():
    seen = []

    def combine(batch):
        seen.append(batch.schema.names)
        return pa.Table.from_batches([batch]).group_by(["query"]).aggregate([("energy", "min")])

    pyrisearch_tauso.search_reduced(
        queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE,
        reduction=pyrisearch_tauso.Reduction(
            columns=("query", "energy"), combine=combine, finalize=lambda t: t, empty=None
        ),
    )

    assert seen and all(names == ["query", "energy"] for names in seen)


def test_hits_table_holds_every_hit():
    table = pyrisearch_tauso.hits_table(queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE)

    assert table.num_rows > 0
    assert table.column_names == list(pyrisearch_tauso.HIT_COLUMNS)
    assert set(table.column("query").to_pylist()) == {QUERY_NAME}


def test_hits_table_keeps_its_columns_when_there_are_no_hits():
    """A caller reads the columns the same way whether or not anything was found."""
    table = pyrisearch_tauso.hits_table(
        queries=QUERIES, targets=TARGETS, min_score=NOTHING_REACHES_THIS
    )

    assert table.num_rows == 0
    assert table.column_names == list(pyrisearch_tauso.HIT_COLUMNS)


def test_hits_table_takes_a_column_subset():
    table = pyrisearch_tauso.hits_table(
        queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE, columns=("query", "energy")
    )

    assert table.column_names == ["query", "energy"]


def test_a_reduction_over_several_batches_matches_one_pass():
    """Small blocks force many batches, which is where a missing second stage shows."""
    small_blocks = pyrisearch_tauso.search_reduced(
        queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE,
        reduction=best_energy_per_pair(), block_size=64,
    )
    one_block = pyrisearch_tauso.search_reduced(
        queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE,
        reduction=best_energy_per_pair(),
    )

    assert small_blocks == one_block


def test_a_min_score_above_what_the_reduction_counts_is_refused():
    """A reduction that counts hits down to a score cannot be fed a search that
    stopped above it: the hits it never saw would be missing with nothing to say so."""
    reduction = pyrisearch_tauso.Reduction(
        columns=("query", "target", "energy"),
        combine=best_energy_per_pair().combine,
        finalize=best_energy_per_pair().finalize,
        empty={},
        min_score_at_most=MIN_SCORE,
    )
    with pytest.raises(ValueError, match="at most"):
        pyrisearch_tauso.search_reduced(
            queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE + 1, reduction=reduction
        )

    reduced = pyrisearch_tauso.search_reduced(
        queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE, reduction=reduction
    )
    assert (QUERY_NAME, TARGET_NAME) in reduced
