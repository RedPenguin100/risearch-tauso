import concurrent.futures
import subprocess
import os

import pyarrow as pa
import pytest

import pyrisearch_tauso
from pyrisearch_tauso import _scratch_dir, _search_args, _write_fasta

QUERY_NAME = "q1"
TARGET_NAME = "t1"
QUERY = "GCTAGCTAGCTAGCTAGCTA"
TARGET = "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
QUERIES = {QUERY_NAME: QUERY}
TARGETS = {TARGET_NAME: TARGET}

# Low enough that the pair above produces hits.
MIN_SCORE = 500


def collect(batches):
    return pa.Table.from_batches(list(batches))


def test_sequences_given_as_mappings():
    with pyrisearch_tauso.search(
        queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE
    ) as batches:
        table = collect(batches)

    assert table.num_rows > 0
    assert table.column("query")[0].as_py() == QUERY_NAME
    assert table.column("target")[0].as_py() == TARGET_NAME


def test_sequences_given_as_pairs():
    with pyrisearch_tauso.search(
        queries=[(QUERY_NAME, QUERY)], targets=[(TARGET_NAME, TARGET)], min_score=MIN_SCORE
    ) as batches:
        from_pairs = collect(batches)

    with pyrisearch_tauso.search(
        queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE
    ) as batches:
        from_mappings = collect(batches)

    assert from_pairs.num_rows == from_mappings.num_rows
    assert from_pairs.column("query")[0].as_py() == QUERY_NAME


def test_files_are_taken_where_they_lie(tmp_path):
    query_file = tmp_path / "q.fa"
    target_file = tmp_path / "t.fa"
    query_file.write_text(f">{QUERY_NAME}\n{QUERY}\n")
    target_file.write_text(f">{TARGET_NAME}\n{TARGET}\n")

    with pyrisearch_tauso.search(
        queries=query_file, targets=target_file, min_score=MIN_SCORE
    ) as batches:
        from_files = collect(batches)

    with pyrisearch_tauso.search(
        queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE
    ) as batches:
        from_sequences = collect(batches)

    assert from_files.num_rows == from_sequences.num_rows
    assert from_files.column("query")[0].as_py() == QUERY_NAME


def test_one_side_may_be_a_file(tmp_path):
    target_file = tmp_path / "t.fa"
    target_file.write_text(f">{TARGET_NAME}\n{TARGET}\n")

    with pyrisearch_tauso.search(
        queries=QUERIES, targets=target_file, min_score=MIN_SCORE
    ) as batches:
        table = collect(batches)

    assert table.num_rows > 0
    assert table.column("query")[0].as_py() == QUERY_NAME
    assert table.column("target")[0].as_py() == TARGET_NAME


def test_sequences_are_written_as_they_are(tmp_path):
    """Nothing is reverse complemented on the way in: the caller decides that."""
    path = tmp_path / "written.fa"
    _write_fasta({QUERY_NAME: QUERY}, path)

    assert path.read_text() == f">{QUERY_NAME}\n{QUERY}\n"


def test_pairs_and_mappings_are_written_alike(tmp_path):
    from_mapping = tmp_path / "a.fa"
    from_pairs = tmp_path / "b.fa"
    _write_fasta({QUERY_NAME: QUERY, TARGET_NAME: TARGET}, from_mapping)
    _write_fasta([(QUERY_NAME, QUERY), (TARGET_NAME, TARGET)], from_pairs)

    assert from_mapping.read_text() == from_pairs.read_text()


def test_what_is_written_is_cleaned_up():
    scratch = _scratch_dir()
    before = set(os.listdir(scratch)) if scratch.is_dir() else set()

    with pyrisearch_tauso.search(
        queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE
    ) as batches:
        collect(batches)

    assert set(os.listdir(scratch)) <= before


def test_absent_options_stay_off_the_command_line():
    """RIsearch's own defaults govern what the caller did not ask for."""
    args = _search_args(
        "q.fa", "t.fa", 900,
        matrix=None, extension_penalty=None, neighborhood=None,
        transpose=False, extra_args=(),
    )

    assert args == ["-q", "q.fa", "-t", "t.fa", "-s", "900"]


def test_options_that_were_asked_for_are_passed():
    args = _search_args(
        "q.fa", "t.fa", 900,
        matrix="t04", extension_penalty=30, neighborhood=5,
        transpose=True, extra_args=["-l", "20"],
    )

    assert args == [
        "-q", "q.fa", "-t", "t.fa", "-s", "900",
        "-m", "t04", "-d", "30", "-n", "5", "-R", "-l", "20",
    ]


def search_one_query(query_name):
    with pyrisearch_tauso.search(
        queries={query_name: QUERY}, targets=TARGETS, min_score=MIN_SCORE
    ) as batches:
        return collect(batches)


def test_searches_running_beside_each_other_keep_their_own_results():
    """Each search writes its own files, so threads cannot read each other's."""
    names = [f"q{i}" for i in range(8)]

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        tables = list(pool.map(search_one_query, names))

    for name, table in zip(names, tables):
        assert set(table.column("query").to_pylist()) == {name}
    assert len({table.num_rows for table in tables}) == 1


def test_the_format_is_still_refused_through_search():
    with pytest.raises(NotImplementedError):
        with pyrisearch_tauso.search(
            queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE, extra_args=["-p3"]
        ):
            pass


def test_a_matrix_reaches_risearch_as_its_value():
    """A member is a string, so it goes on the command line as one."""
    captured = []
    real = subprocess.Popen

    def spy(args, *a, **kw):
        captured.append([str(x) for x in args])
        return real(args, *a, **kw)

    subprocess.Popen = spy
    try:
        with pyrisearch_tauso.search(
            queries=QUERIES, targets=TARGETS, min_score=MIN_SCORE,
            matrix=pyrisearch_tauso.Matrix.T04,
        ) as batches:
            list(batches)
    finally:
        subprocess.Popen = real

    assert "-m" in captured[0]
    assert captured[0][captured[0].index("-m") + 1] == "t04"


def test_every_matrix_risearch_names_is_a_member():
    """RIsearch lists these in its own -m message."""
    assert {m.value for m in pyrisearch_tauso.Matrix} == {
        "t99", "t04", "su95", "su95_noGU", "slh04_noGU"
    }
