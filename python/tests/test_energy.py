import math
from pathlib import Path

import pyarrow as pa
import pytest
from pyrisearch_tauso import EnergyStats, energy_stats, fasta_targets, hits_table, search_reduced


@pytest.mark.parametrize("keys", [("query",), ("query", "target")])
def test_cutoffs_groups_and_batches(keys):
    # q1/t1 straddles batches; a hit exactly at 800 must be excluded from 800.
    table = pa.table(
        {
            "query": ["q1", "q2", "q1", "q1", "q1"],
            "target": ["t1", "t1", "t1", "t2", "t1"],
            "score": [801, 900, 1201, 1000, 800],
            "energy": [-1.0, -2.0, -3.0, -4.0, -100.0],
        }
    )
    reduction = energy_stats([1200, 800, 800, 1500], group_by=keys, rt=1.0)
    result = reduction.finalize(
        pa.concat_tables([reduction.combine(batch) for batch in table.to_batches(max_chunksize=2)])
    )
    q1 = "q1" if len(keys) == 1 else ("q1", "t1")
    assert result[1200] == {q1: EnergyStats(math.exp(3), -3.0)}
    assert result[1500] == {}
    assert result[800][q1].sum_exp == pytest.approx(math.exp(1) + math.exp(3) + (math.exp(4) if len(keys) == 1 else 0))
    assert result[800][q1].min_energy == (-4.0 if len(keys) == 1 else -3.0)


def test_empty_filtered_batch_has_mergeable_schema():
    reduction = energy_stats([1000], rt=0.616)
    batch = pa.table({"query": ["q"], "target": ["t"], "score": [1000], "energy": [-10.0]}).to_batches()[0]
    assert reduction.finalize(reduction.combine(batch)) == {1000: {}}


@pytest.mark.parametrize("rt", [0, -1, float("nan"), float("inf")])
def test_invalid_rt(rt):
    with pytest.raises(ValueError):
        energy_stats([800], rt=rt)


def test_reduction_matches_independent_searches():
    queries = {"q": "GCTAGCTAGCTAGCTAGCTA"}
    targets = {"t": "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"}
    cutoffs = [500, 800, 1000]
    result = search_reduced(
        queries,
        targets,
        min_score=min(cutoffs),
        neighborhood=0,
        reduction=energy_stats(cutoffs, rt=0.616),
        block_size=256,
    )
    for cutoff in cutoffs:
        energies = hits_table(queries, targets, min_score=cutoff, neighborhood=0)["energy"].to_pylist()
        assert energies
        stats = result[cutoff][("q", "t")]
        assert stats.sum_exp == pytest.approx(sum(math.exp(-e / 0.616) for e in energies), rel=1e-14)
        assert stats.min_energy == min(energies)
    assert search_reduced(queries, targets, min_score=100000, reduction=energy_stats([100000], rt=0.616)) == {}


def test_prepared_target_reused_and_cleaned_on_exception():
    with pytest.raises(RuntimeError):
        with fasta_targets({"t": "TAGCTAGCTAGCTAGCTAGC"}) as target:
            path = Path(target)
            assert path.exists()
            first = hits_table({"q": "GCTAGCTAGCTAGCTA"}, target, min_score=500)
            second = hits_table({"q": "GCTAGCTAGCTAGCTA"}, target, min_score=500)
            assert first.num_rows > 0
            assert first.equals(second)
            raise RuntimeError("caller failed")
    assert not path.exists()


def test_existing_target_is_borrowed(tmp_path):
    path = tmp_path / "target.fa"
    path.write_text(">t\nACGT\n")
    with fasta_targets(path) as target:
        assert Path(target) == path
    assert path.read_text() == ">t\nACGT\n"
