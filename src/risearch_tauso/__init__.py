"""RIsearch1 (tauso fork) packaged for Python.

The C binary is bundled as package data. `run()` invokes it with RIsearch's own
command-line options; `executable_path()` hands back the path for callers that
drive the process themselves.
"""

from __future__ import annotations

import subprocess
from collections.abc import Sequence
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _installed_version
from importlib.resources import files
from pathlib import Path

__all__ = ["executable_path", "run", "__version__"]

# Read back from the installed distribution, whose version comes from
# pyproject.toml -- the same place CMake reads it for the binary's banner. A
# literal here is a second copy that nothing checks against the first.
try:
    __version__ = _installed_version("risearch-tauso")
except PackageNotFoundError:  # running from a source tree, not installed
    __version__ = "0.0.0+unknown"

_BINARY_NAME = "RIsearch"


def executable_path() -> str:
    """Absolute path to the bundled RIsearch binary, as a string.

    Raises FileNotFoundError if the binary is missing — that means the wheel
    was built incorrectly or an editable install hasn't compiled it yet.
    """
    p = files(__name__) / "bin" / _BINARY_NAME
    if not Path(str(p)).is_file():
        raise FileNotFoundError(
            f"Bundled RIsearch binary not found at {p}. "
            "If installing in editable mode, ensure `make` succeeded during install."
        )
    return str(p)


def run(
    args: Sequence[str],
    *,
    cwd: str | Path | None = None,
    timeout: float | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run the bundled RIsearch and return the finished process.

    `args` carries RIsearch's own options; the binary path is prepended here, so
    a caller writes `run(["-q", queries, "-t", targets, "-s", "900", "-p2"])`.

    stdout and stderr are captured separately as text, so a warning on stderr
    cannot land in the middle of the hit table. `cwd` sets the directory
    RIsearch resolves relative paths against, the `-m` energy matrix among them.
    `check` raises subprocess.CalledProcessError on a non-zero exit.

    Everything RIsearch writes is held in memory. It prints one line per hit, and
    that count grows with the product of the inputs: sixteen 20-nt queries against
    4000 records of 600 nt at `-s 900` come to 16.4 million lines. Read the
    process's stdout incrementally for runs of that size.
    """
    return subprocess.run(
        [executable_path(), *args],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        cwd=cwd,
        timeout=timeout,
        check=check,
    )
