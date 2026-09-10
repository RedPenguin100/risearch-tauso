"""A Python API over the RIsearch binary that risearch-tauso ships.

`run()` drives one search and hands back the finished process. `stream()` reads
the hits as pyarrow record batches while the search is still going, which is what
a large one needs: RIsearch prints a line per hit, and that count grows with the
product of the inputs.
"""

from __future__ import annotations

import subprocess
import tempfile
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _installed_version
from pathlib import Path

import pyarrow as pa
import pyarrow.csv as pacsv

from risearch_tauso import executable_path

__all__ = [
    "DEFAULT_BLOCK_SIZE",
    "HIT_COLUMNS",
    "RIsearchError",
    "executable_path",
    "run",
    "stream",
    "__version__",
]

try:
    __version__ = _installed_version("pyrisearch-tauso")
except PackageNotFoundError:  # running from a source tree, not installed
    __version__ = "0.0.0+unknown"

# The table `-p2` prints, in order.
HIT_COLUMNS = ("qname", "qbeg", "qend", "tname", "tbeg", "tend", "score", "energy")

# Coordinates and scores are C ints; energy keeps the width it is printed at.
HIT_TYPES = {
    "qname": pa.string(),
    "qbeg": pa.int32(),
    "qend": pa.int32(),
    "tname": pa.string(),
    "tbeg": pa.int32(),
    "tend": pa.int32(),
    "score": pa.int32(),
    "energy": pa.float64(),
}

# Peak memory follows this rather than the size of the run, and parsing does not
# get cheaper as it grows: 16 MB blocks held a 670 MB run to 282 MB, 64 MB blocks
# to 738 MB, for the same wall time.
DEFAULT_BLOCK_SIZE = 16 << 20


class RIsearchError(subprocess.CalledProcessError):
    """A RIsearch run that exited non-zero.

    Subclasses CalledProcessError, so `except subprocess.CalledProcessError`
    catches it. What it adds is RIsearch's own stderr in the message, which is
    where the reason for the failure is written.
    """

    def __str__(self) -> str:
        reason = (self.stderr or "").strip()
        return f"{super().__str__()} {reason}" if reason else super().__str__()


def run(
    args: Sequence[str],
    *,
    cwd: str | Path | None = None,
    timeout: float | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run RIsearch and return the finished process.

    `args` carries RIsearch's own options -- the binary path is prepended here.
    `cwd` is what RIsearch resolves relative paths against, the `-m` energy
    matrix among them.

    A run that fails raises RIsearchError, so an empty result means RIsearch
    found no hits rather than that it never searched. `check=False` hands the
    failure back instead of raising.

    Output is buffered in memory, and RIsearch prints a line per hit, which a
    large run reaches tens of millions of. Use stream() for those.
    """
    completed = subprocess.run(
        [executable_path(), *args],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        cwd=cwd,
        timeout=timeout,
    )
    if check and completed.returncode != 0:
        raise RIsearchError(
            completed.returncode, completed.args, completed.stdout, completed.stderr
        )
    return completed


def _batches(stdout, read_options, parse_options, convert_options):
    try:
        reader = pacsv.open_csv(
            stdout,
            read_options=read_options,
            parse_options=parse_options,
            convert_options=convert_options,
        )
    except pa.ArrowInvalid:
        # RIsearch wrote nothing, which is what no hits looks like. Raised here
        # rather than during the walk, so a malformed row still surfaces.
        return
    for batch in reader:
        if batch.num_rows:
            yield batch


@contextmanager
def stream(
    args: Sequence[str],
    *,
    columns: Sequence[str] = HIT_COLUMNS,
    block_size: int = DEFAULT_BLOCK_SIZE,
    cwd: str | Path | None = None,
) -> Iterator[Iterator[pa.RecordBatch]]:
    """Run RIsearch and read its hits as pyarrow record batches.

    `args` carries RIsearch's own options without an output format -- this reads
    the eight-column table and passes `-p2` itself. `columns` picks the subset to
    parse, and the rest are never converted.

    Batches arrive while RIsearch is still running, so what is held at once
    follows `block_size` and not the length of the run. Reducing each batch as it
    arrives keeps it that way; collecting them into one table does not.

    A run that fails raises RIsearchError once the batches are done. RIsearch's
    stderr goes to a temporary file, not a pipe: it writes a warning per
    nonstandard base, and a pipe nobody drains fills up and stops the process.
    """
    if any(a.startswith("-p") for a in args):
        raise ValueError("stream() sets the output format itself; leave -p out of args")
    unknown = [c for c in columns if c not in HIT_COLUMNS]
    if unknown:
        raise ValueError(f"not RIsearch hit columns: {unknown}")

    read_options = pacsv.ReadOptions(
        column_names=list(HIT_COLUMNS), use_threads=False, block_size=block_size
    )
    parse_options = pacsv.ParseOptions(delimiter="\t")
    convert_options = pacsv.ConvertOptions(
        include_columns=list(columns), column_types={c: HIT_TYPES[c] for c in columns}
    )

    with tempfile.TemporaryFile() as errors:
        proc = subprocess.Popen(
            [executable_path(), *args, "-p2"],
            stdout=subprocess.PIPE,
            stderr=errors,
            stdin=subprocess.DEVNULL,
            cwd=cwd,
            bufsize=-1,
        )
        try:
            yield _batches(proc.stdout, read_options, parse_options, convert_options)
        except BaseException:
            # The caller stopped early or failed; theirs is the error worth
            # seeing, so end the run without reading it out.
            proc.kill()
            proc.wait()
            raise
        else:
            proc.stdout.read()
            if proc.wait() != 0:
                errors.seek(0)
                raise RIsearchError(
                    proc.returncode,
                    proc.args,
                    None,
                    errors.read().decode("utf-8", "replace"),
                )
        finally:
            proc.stdout.close()
