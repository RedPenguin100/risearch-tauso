"""A Python API over the RIsearch binary that risearch-tauso ships.

`search()` takes sequences and gives back the hits as pyarrow record batches,
reading them while RIsearch is still going -- which is what a large search needs,
since RIsearch prints a line per hit and that count grows with the product of the
inputs. `stream()` is the same thing for inputs that are already FASTA files, and
`run()` drives one search and hands back the finished process.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from enum import Enum
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _installed_version
from pathlib import Path

import pyarrow as pa
import pyarrow.csv as pacsv

from risearch_tauso import executable_path

__all__ = [
    "DEFAULT_BLOCK_SIZE",
    "EnergyStats",
    "energy_stats",
    "fasta_targets",
    "HIT_COLUMNS",
    "Matrix",
    "RIsearchError",
    "Reduction",
    "executable_path",
    "hits_table",
    "run",
    "search",
    "search_reduced",
    "stream",
    "__version__",
]

try:
    __version__ = _installed_version("pyrisearch-tauso")
except PackageNotFoundError:  # running from a source tree, not installed
    __version__ = "0.0.0+unknown"

class Matrix(str, Enum):
    """The scoring matrices RIsearch takes for `-m`.

    Each is a string, so it can be passed to `matrix=` as it is or written out
    in full. The interaction each one scores is RIsearch's own, from the table
    in its src/dsm.h.
    """

    T99 = "t99"  # RNA-RNA
    T04 = "t04"  # RNA-RNA, and what RIsearch uses when not told otherwise
    SU95 = "su95"  # RNA-DNA
    SU95_NO_GU = "su95_noGU"  # RNA-DNA, wobble pairs left out
    SLH04_NO_GU = "slh04_noGU"  # DNA-DNA, wobble pairs left out

    # Without this a member reads as "Matrix.T04" wherever it is printed, and
    # what RIsearch wants is the value.
    __str__ = str.__str__


# The table `-p2` prints, in order.
HIT_COLUMNS = ("query", "query_start", "query_end", "target", "target_start", "target_end", "score", "energy")

# Coordinates and scores are C ints; energy keeps the width it is printed at.
HIT_TYPES = {
    "query": pa.string(),
    "query_start": pa.int32(),
    "query_end": pa.int32(),
    "target": pa.string(),
    "target_start": pa.int32(),
    "target_end": pa.int32(),
    "score": pa.int32(),
    "energy": pa.float64(),
}

# What RIsearch's -p values print. Only the eight-column table is parsed here;
# the others have a shape of their own and nothing reads them yet.
KNOWN_FORMATS = ("1", "2", "3")
READABLE_FORMAT = "2"

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


def _output_format(args):
    """The -p value in args, or None. RIsearch takes it attached or separate."""
    args = list(args)
    for i, arg in enumerate(args):
        if arg == "-p":
            return args[i + 1] if i + 1 < len(args) else ""
        if arg.startswith("-p"):
            return arg[2:]
    return None


def _with_format(args):
    """args as RIsearch should get them, once the output format is settled."""
    fmt = _output_format(args)
    if fmt is None:
        return [*args, f"-p{READABLE_FORMAT}"]
    if fmt == READABLE_FORMAT:
        return list(args)
    if fmt in KNOWN_FORMATS:
        raise NotImplementedError(
            f"stream() reads the -p{READABLE_FORMAT} table; -p{fmt} prints a different "
            "shape, which nothing parses yet"
        )
    raise ValueError(f"not a RIsearch output format: -p{fmt}")


def _batches(stdout, read_options, parse_options, convert_options, reached_the_end):
    """The hits, a batch at a time, appending to `reached_the_end` when there are
    no more.

    A search that found nothing writes nothing, and the reader cannot tell that
    from a first block it could not parse -- both raise ArrowInvalid. Asking the
    pipe whether anything arrived separates them, so unparseable output raises
    rather than passing for no hits.
    """
    if not stdout.peek(1):
        reached_the_end.append(True)
        return
    reader = pacsv.open_csv(
        stdout,
        read_options=read_options,
        parse_options=parse_options,
        convert_options=convert_options,
    )
    for batch in reader:
        if batch.num_rows:
            yield batch
    reached_the_end.append(True)


@contextmanager
def stream(
    args: Sequence[str],
    *,
    columns: Sequence[str] = HIT_COLUMNS,
    block_size: int = DEFAULT_BLOCK_SIZE,
    cwd: str | Path | None = None,
) -> Iterator[Iterator[pa.RecordBatch]]:
    """Run RIsearch and read its hits as pyarrow record batches.

    `args` carries RIsearch's own options. The output format may be left out, in
    which case `-p2` is added; `-p2` may also be passed outright. The other
    formats raise NotImplementedError -- RIsearch prints them, this does not read
    them yet. `columns` picks the subset to parse, and the rest are never
    converted.

    Batches arrive while RIsearch is still running, so what is held at once
    follows `block_size` and not the length of the run. Reducing each batch as it
    arrives keeps it that way; collecting them into one table does not.

    A run that fails raises RIsearchError once the batches are done. RIsearch's
    stderr goes to a temporary file, not a pipe: it writes a warning per
    nonstandard base, and a pipe nobody drains fills up and stops the process.
    """
    args = _with_format(args)
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
            [executable_path(), *args],
            stdout=subprocess.PIPE,
            stderr=errors,
            stdin=subprocess.DEVNULL,
            cwd=cwd,
            bufsize=-1,
        )
        reached_the_end = []
        try:
            yield _batches(
                proc.stdout, read_options, parse_options, convert_options, reached_the_end
            )
        except BaseException:
            # The caller stopped early or failed; theirs is the error worth
            # seeing, so end the run without reading it out.
            proc.kill()
            proc.wait()
            raise
        else:
            if not reached_the_end:
                # The caller has what it wanted and RIsearch has more to write.
                # Reading the rest out to be polite would hold all of it.
                proc.kill()
                proc.wait()
            elif proc.wait() != 0:
                errors.seek(0)
                raise RIsearchError(
                    proc.returncode,
                    proc.args,
                    None,
                    errors.read().decode("utf-8", "replace"),
                )
        finally:
            proc.stdout.close()


def _scratch_dir():
    """Where the FASTA files a search writes go.

    tempfile.gettempdir() reads TMPDIR, which is how a machine says where its
    scratch space is. What is written stays in the page cache, so RIsearch reads
    it back out of memory, and the kernel is free to drop it when memory is
    wanted elsewhere -- which a tmpfs like /dev/shm cannot do, since what is
    written there is held until it is deleted.
    """
    return Path(tempfile.gettempdir()) / "pyrisearch_tauso"


# RIsearch reads a header line into a 512-byte buffer: the '>', the name, the
# newline, and the NUL fgets adds. A header that does not fit is read as two
# lines, the second as sequence, which moves every coordinate it reports.
MAX_NAME_BYTES = 512 - 3


def _fasta_name(name):
    """`name` as it goes after the '>', refused when RIsearch would misread it.

    RIsearch takes the name to be the first whitespace-separated token after the
    '>', so whitespace cuts it short and the hits come back under a different
    name; a header with nothing after the '>' crashes it.
    """
    name = str(name)
    if not name or any(c.isspace() for c in name):
        raise ValueError(f"a FASTA name must have no whitespace and not be empty: {name!r}")
    if len(name.encode()) > MAX_NAME_BYTES:
        raise ValueError(f"a FASTA name must be at most {MAX_NAME_BYTES} bytes: {name[:32]!r}...")
    return name


def _write_fasta(sequences, path):
    with open(path, "w") as out:
        if isinstance(sequences, Mapping):
            sequences = sequences.items()
        for name, sequence in sequences:
            out.write(f">{_fasta_name(name)}\n{sequence}\n")


@contextmanager
def _as_fasta(sequences, role):
    """A path to `sequences`, writing one only when they are not already a file.

    The file is named uniquely and removed on the way out, so searches running
    beside each other in threads cannot land on the same name.
    """
    if isinstance(sequences, (str, Path)):
        yield str(sequences)
        return

    scratch = _scratch_dir()
    scratch.mkdir(parents=True, exist_ok=True)
    fd, path = tempfile.mkstemp(prefix=f"{role}-", suffix=".fa", dir=str(scratch))
    os.close(fd)
    try:
        _write_fasta(sequences, path)
        yield path
    finally:
        try:
            os.remove(path)
        except OSError:
            pass


def _search_args(query_path, target_path, min_score, matrix, extension_penalty,
                 neighborhood, transpose, extra_args):
    """Only the options that were asked for; the rest are RIsearch's own defaults."""
    args = ["-q", query_path, "-t", target_path, "-s", str(min_score)]
    if matrix is not None:
        args += ["-m", matrix]
    if extension_penalty is not None:
        args += ["-d", str(extension_penalty)]
    if neighborhood is not None:
        args += ["-n", str(neighborhood)]
    if transpose:
        args.append("-R")
    return args + list(extra_args)


@contextmanager
def search(
    queries,
    targets,
    *,
    min_score: int,
    matrix: str | None = None,
    extension_penalty: int | None = None,
    neighborhood: int | None = None,
    transpose: bool = False,
    columns: Sequence[str] = HIT_COLUMNS,
    block_size: int = DEFAULT_BLOCK_SIZE,
    extra_args: Sequence[str] = (),
) -> Iterator[Iterator[pa.RecordBatch]]:
    """Search `queries` against `targets` and read the hits as record batches.

    Either side is a mapping of name to sequence, a sequence of (name, sequence)
    pairs, or a path to a FASTA file. Sequences are written to a file of their
    own under the machine's scratch directory and removed afterwards; a path is
    used where it lies.

    The sequences go to RIsearch as given. A query meant to bind its target has
    to be handed over already reverse complemented -- that belongs to whoever
    knows what the sequences are for.

    `matrix`, `extension_penalty` and `neighborhood` are left out of the command
    line when they are None, so RIsearch's own defaults govern rather than any
    chosen here. `extra_args` carries options this signature does not name.

    What is held at once follows `block_size`, so the number of queries is
    limited by time rather than memory: 2000 queries producing 385 million hits
    ran in 99 MB.
    """
    with _as_fasta(queries, "query") as query_path, _as_fasta(targets, "target") as target_path:
        args = _search_args(query_path, target_path, min_score, matrix,
                            extension_penalty, neighborhood, transpose, extra_args)
        with stream(args, columns=columns, block_size=block_size) as batches:
            yield batches


def _empty_hits(columns):
    """A table with the hit schema and no rows, for a search that found none."""
    return pa.table({column: pa.array([], HIT_TYPES[column]) for column in columns})


def hits_table(queries, targets, **search_options):
    """Every hit of one search, as a single pyarrow table.

    Takes what search() takes. The whole result is held at once, so this is for
    a search whose hits fit in memory -- checking one against another, or a run
    small enough to look at. search_reduced() is what a large one wants.

    A search that found nothing gives back a table with the hit schema and no
    rows, rather than nothing at all, so a caller can read its columns either
    way.
    """
    columns = search_options.get("columns", HIT_COLUMNS)
    with search(queries, targets, **search_options) as batches:
        tables = [pa.Table.from_batches([batch]) for batch in batches]
    return pa.concat_tables(tables) if tables else _empty_hits(columns)


@dataclass(frozen=True)
class Reduction:
    """How to turn the hits of a search into one answer.

    `columns` names what the reduction reads, and is what the search parses.
    `combine` takes one record batch and gives back a partial table; `finalize`
    takes those partials concatenated and gives back the answer. `empty` is the
    answer when the search found no hits at all.

    Two stages rather than one because a batch holds part of the hits: the
    smallest energy within a batch is the smallest of that batch, and the
    smallest overall is only known once the partials are put together. The
    second pass is what does that, and leaving it out gives one row per batch
    where there should be one.

    `min_score_at_most` is the highest `min_score` the search may run with: a
    reduction that counts every hit above some score has to be handed all of
    them, and a search run above that score leaves some out with nothing to say
    so. None puts no limit on it.
    """

    columns: Sequence[str]
    combine: object
    finalize: object
    empty: object = None
    min_score_at_most: int | None = None


def search_reduced(queries, targets, *, reduction: Reduction, **search_options):
    """Search, and reduce the hits as they arrive.

    Takes what search() takes, apart from `columns`, which the reduction names.
    Each batch is reduced as it comes and only the partials are kept, so a
    search whose hits would not fit in memory still has an answer that does.
    """
    if "columns" in search_options:
        raise TypeError("search_reduced reads the columns off the reduction")
    limit = reduction.min_score_at_most
    if limit is not None and search_options.get("min_score", limit) > limit:
        raise ValueError(
            f"min_score={search_options['min_score']} leaves out hits the reduction "
            f"counts: it must be at most {limit}"
        )

    parts = []
    with search(queries, targets, columns=reduction.columns, **search_options) as batches:
        for batch in batches:
            parts.append(reduction.combine(batch))

    if not parts:
        return reduction.empty
    return reduction.finalize(pa.concat_tables(parts))


@contextmanager
def fasta_targets(sequences):
    """Prepare a target FASTA once for reuse across searches, cleaning it up on exit.

    Accepts the same targets as search(). An existing path is borrowed and never
    deleted. Keep this context open until all searches using the path have finished.
    """
    with _as_fasta(sequences, "target") as path:
        yield path


# Last, not at the top: energy.py builds on Reduction, so it can only be
# imported once this module has defined it.
from .energy import EnergyStats, energy_stats  # noqa: E402
