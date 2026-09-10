"""Tests for `risearch_tauso.run`, which drives the bundled binary."""

from __future__ import annotations

import subprocess

import pytest

import risearch_tauso

QUERY = ">q1\nGCTAGCTAGCTAGCTAGCTA\n"
TARGET = ">t1\nTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\n"


@pytest.fixture
def pair(tmp_path):
    """A query and a target file that produce hits at -s 500."""
    q = tmp_path / "q.fa"
    t = tmp_path / "t.fa"
    q.write_text(QUERY)
    t.write_text(TARGET)
    return q, t


def test_hits_are_eight_tab_separated_columns(pair):
    q, t = pair
    result = risearch_tauso.run(["-q", str(q), "-t", str(t), "-s", "500", "-p2"])

    assert result.returncode == 0
    lines = result.stdout.splitlines()
    assert lines, "expected at least one hit"
    for line in lines:
        assert len(line.split("\t")) == 8


def test_query_and_target_names_come_back(pair):
    q, t = pair
    result = risearch_tauso.run(["-q", str(q), "-t", str(t), "-s", "500", "-p2"])

    qname, _qbeg, _qend, tname = result.stdout.splitlines()[0].split("\t")[:4]
    assert (qname, tname) == ("q1", "t1")


def test_stderr_is_kept_out_of_stdout(pair):
    """A run that fails puts nothing on stdout, so a caller parsing it sees no rows."""
    q, _t = pair
    result = risearch_tauso.run(["-q", str(q)], check=False)

    assert result.returncode != 0
    assert result.stdout == ""


def test_check_raises_on_failure(pair):
    q, _t = pair
    with pytest.raises(subprocess.CalledProcessError):
        risearch_tauso.run(["-q", str(q)])


def test_check_false_returns_the_failure(pair):
    q, _t = pair
    result = risearch_tauso.run(["-q", str(q)], check=False)

    assert result.returncode != 0


def test_executable_path_is_the_bundled_binary():
    assert risearch_tauso.executable_path().endswith("/bin/RIsearch")
