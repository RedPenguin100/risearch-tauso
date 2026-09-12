import subprocess

import pytest

import pyarrow as pa

import pyrisearch_tauso

QUERY_NAME = "q1"
TARGET_NAME = "t1"
QUERY = f">{QUERY_NAME}\nGCTAGCTAGCTAGCTAGCTA\n"
TARGET = f">{TARGET_NAME}\nTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\n"


@pytest.fixture
def pair(tmp_path):
    q = tmp_path / "q.fa"
    t = tmp_path / "t.fa"
    q.write_text(QUERY)
    t.write_text(TARGET)
    return q, t


def args_for(query, target, score="500"):
    return ["-q", str(query), "-t", str(target), "-s", score]


def test_batches_carry_the_hit_columns(pair):
    query, target = pair
    with pyrisearch_tauso.stream(args_for(query, target)) as batches:
        table = pa.Table.from_batches(list(batches))

    assert table.num_rows > 0
    assert table.column_names == list(pyrisearch_tauso.HIT_COLUMNS)


def test_streamed_rows_match_the_buffered_run(pair):
    """stream() and run() read the same hits out of the same invocation."""
    query, target = pair
    with pyrisearch_tauso.stream(args_for(query, target)) as batches:
        streamed = pa.Table.from_batches(list(batches)).num_rows

    result = pyrisearch_tauso.run(args_for(query, target) + ["-p2"])

    assert streamed == len(result.stdout.splitlines())


def test_a_column_subset_skips_the_rest(pair):
    query, target = pair
    with pyrisearch_tauso.stream(args_for(query, target), columns=("query", "energy")) as batches:
        table = pa.Table.from_batches(list(batches))

    assert table.column_names == ["query", "energy"]
    assert table.column("query")[0].as_py() == QUERY_NAME


def test_no_hits_yields_no_batches(pair):
    """Empty output is not an error: RIsearch searched and found nothing."""
    query, target = pair
    with pyrisearch_tauso.stream(args_for(query, target, score="100000")) as batches:
        assert list(batches) == []


def test_a_failed_run_raises_after_the_batches(tmp_path):
    missing = tmp_path / "nope.fa"
    with pytest.raises(pyrisearch_tauso.RIsearchError) as excinfo:
        with pyrisearch_tauso.stream(args_for(missing, missing)) as batches:
            list(batches)

    assert "not readable" in str(excinfo.value)


def test_the_failure_is_a_called_process_error(tmp_path):
    missing = tmp_path / "nope.fa"
    with pytest.raises(subprocess.CalledProcessError):
        with pyrisearch_tauso.stream(args_for(missing, missing)) as batches:
            list(batches)


def test_the_readable_format_may_be_passed_outright(pair):
    """-p2 is what stream() reads, so saying so explicitly is not an error."""
    query, target = pair
    with pyrisearch_tauso.stream(args_for(query, target) + ["-p2"]) as batches:
        assert pa.Table.from_batches(list(batches)).num_rows > 0


def test_the_readable_format_may_be_written_separately(pair):
    """RIsearch takes -p2 and -p 2 alike, so both have to be recognised."""
    query, target = pair
    with pyrisearch_tauso.stream(args_for(query, target) + ["-p", "2"]) as batches:
        assert pa.Table.from_batches(list(batches)).num_rows > 0


@pytest.mark.parametrize("fmt", [["-p1"], ["-p3"], ["-p", "1"], ["-p", "3"]])
def test_a_format_nothing_parses_yet_is_refused(pair, fmt):
    query, target = pair
    with pytest.raises(NotImplementedError, match="nothing parses yet"):
        with pyrisearch_tauso.stream(args_for(query, target) + fmt):
            pass


def test_a_value_that_is_not_a_format_is_refused(pair):
    query, target = pair
    with pytest.raises(ValueError, match="not a RIsearch output format"):
        with pyrisearch_tauso.stream(args_for(query, target) + ["-p9"]):
            pass


def test_an_unknown_column_is_refused(pair):
    query, target = pair
    with pytest.raises(ValueError, match="not RIsearch hit columns"):
        with pyrisearch_tauso.stream(args_for(query, target), columns=("query", "nonsense")):
            pass


def test_the_callers_error_is_the_one_raised(pair):
    """Leaving the block early must not bury the reason in a teardown failure."""
    query, target = pair
    with pytest.raises(ZeroDivisionError):
        with pyrisearch_tauso.stream(args_for(query, target)) as batches:
            next(iter(batches))
            1 / 0


def test_warnings_do_not_stall_the_stream(tmp_path):
    """A nonstandard base warns once per occurrence. Enough of them fill a pipe,
    so the warnings must not be going to one."""
    query = tmp_path / "q.fa"
    target = tmp_path / "t.fa"
    query.write_text(QUERY)
    target.write_text(">t1\n" + "Z" * 200_000 + "\n")

    with pyrisearch_tauso.stream(args_for(query, target)) as batches:
        list(batches)

