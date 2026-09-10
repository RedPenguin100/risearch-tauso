import subprocess

import pytest

import risearch_tauso

QUERY_NAME = "q1"
TARGET_NAME = "t1"
QUERY = f">{QUERY_NAME}\nGCTAGCTAGCTAGCTAGCTA\n"
TARGET = f">{TARGET_NAME}\nTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\n"


@pytest.fixture
def pair(tmp_path):
    """A query and a target file that produce hits at -s 500."""
    q = tmp_path / "q.fa"
    t = tmp_path / "t.fa"
    q.write_text(QUERY)
    t.write_text(TARGET)
    return q, t


def test_hit_parsing(pair):
    query, target = pair
    result = risearch_tauso.run(["-q", str(query), "-t", str(target), "-s", "500", "-p2"])

    assert result.returncode == 0
    lines = result.stdout.splitlines()
    assert lines, "expected at least one hit"
    for line in lines:
        assert len(line.split("\t")) == 8


def test_query_target_name_parsing(pair):
    query, target = pair
    result = risearch_tauso.run(["-q", str(query), "-t", str(target), "-s", "500", "-p2"])

    qname, _qbeg, _qend, tname = result.stdout.splitlines()[0].split("\t")[:4]
    assert (qname, tname) == (QUERY_NAME, TARGET_NAME)


def test_no_hits_is_not_a_failure(pair):
    """A threshold nothing reaches: RIsearch searched and found none."""
    query, target = pair
    result = risearch_tauso.run(["-q", str(query), "-t", str(target), "-s", "100000", "-p2"])

    assert result.returncode == 0
    assert result.stdout == ""


def test_check_raises_on_failure(pair):
    query, _ = pair
    with pytest.raises(risearch_tauso.RIsearchError):
        risearch_tauso.run(["-q", str(query)])


def test_failure_is_a_called_process_error(pair):
    """Callers already catching subprocess.CalledProcessError keep working."""
    query, _ = pair
    with pytest.raises(subprocess.CalledProcessError):
        risearch_tauso.run(["-q", str(query)])


def test_failure_message_names_the_reason(tmp_path):
    missing = tmp_path / "nope.fa"
    with pytest.raises(risearch_tauso.RIsearchError) as excinfo:
        risearch_tauso.run(["-q", str(missing), "-t", str(missing), "-s", "500", "-p2"])

    assert "not readable" in str(excinfo.value)


def test_check_false_returns_the_failure(pair):
    query, _ = pair
    result = risearch_tauso.run(["-q", str(query)], check=False)

    assert result.returncode != 0


def test_warning_stays_out_of_stdout(tmp_path):
    """A nonstandard base warns on stderr; stdout stays a parseable hit table."""
    query = tmp_path / "q.fa"
    target = tmp_path / "t.fa"
    query.write_text(QUERY)
    target.write_text(f">{TARGET_NAME}\nZAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\n")

    result = risearch_tauso.run(["-q", str(query), "-t", str(target), "-s", "500", "-p2"])

    assert result.returncode == 0
    assert "Nonstandard" in result.stderr
    lines = result.stdout.splitlines()
    assert lines, "expected hits despite the replaced base"
    for line in lines:
        assert len(line.split("\t")) == 8


def test_executable_path_is_the_bundled_binary():
    assert risearch_tauso.executable_path().endswith("/bin/RIsearch")
