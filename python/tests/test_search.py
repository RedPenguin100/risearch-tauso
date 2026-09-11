import concurrent.futures
import os

import pyarrow as pa
import pytest

import pyrisearch_tauso
from pyrisearch_tauso import _scratch_dir, _search_args, _write_fasta

QUERIES = {"q1": "GCTAGCTAGCTAGCTAGCTA"}
TARGETS = {"t1": "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"}


def hits(**kwargs):
    kwargs.setdefault("min_score", 500)
    with pyrisearch_tauso.search(**kwargs) as batches:
        return pa.Table.from_batches(list(batches))


def test_sequences_given_as_mappings():
    table = hits(queries=QUERIES, targets=TARGETS)

    assert table.num_rows > 0
    assert table.column("qname")[0].as_py() == "q1"
    assert table.column("tname")[0].as_py() == "t1"


def test_sequences_given_as_pairs():
    pairs = hits(queries=[("q1", QUERIES["q1"])], targets=[("t1", TARGETS["t1"])])
    mappings = hits(queries=QUERIES, targets=TARGETS)

    assert pairs.num_rows == mappings.num_rows


def test_files_are_taken_where_they_lie(tmp_path):
    query = tmp_path / "q.fa"
    target = tmp_path / "t.fa"
    query.write_text(">q1\n" + QUERIES["q1"] + "\n")
    target.write_text(">t1\n" + TARGETS["t1"] + "\n")

    from_files = hits(queries=query, targets=target)
    from_sequences = hits(queries=QUERIES, targets=TARGETS)

    assert from_files.num_rows == from_sequences.num_rows


def test_one_side_may_be_a_file(tmp_path):
    target = tmp_path / "t.fa"
    target.write_text(">t1\n" + TARGETS["t1"] + "\n")

    assert hits(queries=QUERIES, targets=target).num_rows > 0


def test_sequences_are_written_as_they_are(tmp_path):
    """Nothing is reverse complemented on the way in: the caller decides that."""
    path = tmp_path / "written.fa"
    _write_fasta({"q1": "ACGGTTCAGGCATTACGAGT"}, path)

    assert path.read_text() == ">q1\nACGGTTCAGGCATTACGAGT\n"


def test_pairs_and_mappings_are_written_alike(tmp_path):
    from_mapping = tmp_path / "a.fa"
    from_pairs = tmp_path / "b.fa"
    _write_fasta({"q1": "ACGT", "q2": "TTTT"}, from_mapping)
    _write_fasta([("q1", "ACGT"), ("q2", "TTTT")], from_pairs)

    assert from_mapping.read_text() == from_pairs.read_text()


def test_what_is_written_is_cleaned_up():
    before = set(os.listdir(_scratch_dir())) if _scratch_dir().is_dir() else set()
    hits(queries=QUERIES, targets=TARGETS)
    after = set(os.listdir(_scratch_dir()))

    assert after <= before


def test_absent_options_stay_off_the_command_line():
    """RIsearch's own defaults govern what the caller did not ask for."""
    args = _search_args("q.fa", "t.fa", 900, None, None, None, False, ())

    assert args == ["-q", "q.fa", "-t", "t.fa", "-s", "900"]


def test_options_that_were_asked_for_are_passed():
    args = _search_args("q.fa", "t.fa", 900, "t04", 30, 5, True, ["-l", "20"])

    assert args == [
        "-q", "q.fa", "-t", "t.fa", "-s", "900",
        "-m", "t04", "-d", "30", "-n", "5", "-R", "-l", "20",
    ]


def test_searches_running_beside_each_other_keep_their_own_results():
    """Each search writes its own files, so threads cannot read each other's."""
    queries = [{f"q{i}": QUERIES["q1"]} for i in range(8)]

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        tables = list(pool.map(lambda q: hits(queries=q, targets=TARGETS), queries))

    for i, table in enumerate(tables):
        assert set(table.column("qname").to_pylist()) == {f"q{i}"}
    assert len({t.num_rows for t in tables}) == 1


def test_the_format_is_still_refused_through_search():
    with pytest.raises(NotImplementedError):
        with pyrisearch_tauso.search(
            queries=QUERIES, targets=TARGETS, min_score=500, extra_args=["-p3"]
        ):
            pass
