"""What stream() does when the output is not a clean run of hits."""

import io
import subprocess

import pyarrow as pa
import pyarrow.csv as pacsv
import pytest

import pyrisearch_tauso
from pyrisearch_tauso import HIT_COLUMNS, HIT_TYPES, _batches

QUERY_NAME = "q1"
TARGET_NAME = "t1"
QUERIES = {QUERY_NAME: "GCTAGCTAGCTAGCTAGCTA"}
TARGETS = {TARGET_NAME: "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"}
MIN_SCORE = 500

A_WELL_FORMED_HIT = b"q1\t1\t20\tt1\t13\t32\t4345\t-37.86\n"
# RIsearch prints eight columns; three is not a hit table.
A_MALFORMED_HIT = b"q1\t1\t20\n"


def drain(data):
    read_options = pacsv.ReadOptions(
        column_names=list(HIT_COLUMNS), use_threads=False, block_size=1 << 20
    )
    parse_options = pacsv.ParseOptions(delimiter="\t")
    convert_options = pacsv.ConvertOptions(
        include_columns=["query", "energy"],
        column_types={"query": HIT_TYPES["query"], "energy": HIT_TYPES["energy"]},
    )
    stdout = io.BufferedReader(io.BytesIO(data))
    reached_the_end = []
    rows = sum(
        b.num_rows
        for b in _batches(stdout, read_options, parse_options, convert_options, reached_the_end)
    )
    return rows, reached_the_end


def test_output_that_cannot_be_parsed_raises():
    """Nothing written and something unreadable both raise ArrowInvalid from the
    reader, and only the first of them means no hits."""
    with pytest.raises(pa.ArrowInvalid):
        drain(A_MALFORMED_HIT)


def test_nothing_written_is_no_hits():
    rows, reached_the_end = drain(b"")

    assert rows == 0
    assert reached_the_end


def test_well_formed_output_still_reads():
    rows, reached_the_end = drain(A_WELL_FORMED_HIT)

    assert rows == 1
    assert reached_the_end


def test_leaving_early_does_not_read_the_rest_out(tmp_path):
    """A caller that stops gets RIsearch stopped too, rather than the remaining
    output being read into memory to be tidy about it."""
    query = tmp_path / "q.fa"
    target = tmp_path / "t.fa"
    query.write_text("".join(f">q{i}\nGCTAGCTAGCTAGCTAGCTA\n" for i in range(16)))
    target.write_text("".join(f">t{i}\n{'TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC' * 40}\n" for i in range(400)))

    started = []
    real = subprocess.Popen

    def spy(*a, **kw):
        proc = real(*a, **kw)
        started.append(proc)
        return proc

    subprocess.Popen = spy
    try:
        with pyrisearch_tauso.stream(
            ["-q", str(query), "-t", str(target), "-s", str(MIN_SCORE)], block_size=4096
        ) as batches:
            next(iter(batches))  # take one and stop
    finally:
        subprocess.Popen = real

    assert started[0].returncode is not None, "RIsearch should have been stopped"
    assert started[0].returncode < 0, "stopped by a signal, not left to finish writing"


def test_reading_every_batch_still_checks_the_exit_code(tmp_path):
    missing = tmp_path / "nope.fa"
    with pytest.raises(pyrisearch_tauso.RIsearchError):
        with pyrisearch_tauso.stream(
            ["-q", str(missing), "-t", str(missing), "-s", str(MIN_SCORE)]
        ) as batches:
            list(batches)
