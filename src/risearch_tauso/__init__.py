"""RIsearch1 (tauso fork) packaged for Python.

The C binary is bundled as package data. `run()` invokes it with RIsearch's own
command-line options; `executable_path()` hands back the path for callers that
drive the process themselves.
"""

import subprocess
from collections.abc import Sequence
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _installed_version
from importlib.resources import files
from pathlib import Path
from typing import Optional, Union

__all__ = ["RIsearchError", "executable_path", "run", "__version__"]


class RIsearchError(subprocess.CalledProcessError):
    """A RIsearch run that exited non-zero.

    Subclasses CalledProcessError, so `except subprocess.CalledProcessError`
    catches it. What it adds is RIsearch's own stderr in the message, which is
    where the reason for the failure is written.
    """

    def __str__(self) -> str:
        reason = (self.stderr or "").strip()
        return f"{super().__str__()} {reason}" if reason else super().__str__()

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
    cwd: Optional[Union[str, Path]] = None,
    timeout: Optional[float] = None,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run the bundled RIsearch and return the finished process.

    `args` carries RIsearch's own options -- the binary path is prepended here.
    `cwd` is what RIsearch resolves relative paths against, the `-m` energy
    matrix among them.

    A run that fails raises RIsearchError, so an empty result means RIsearch
    found no hits rather than that it never searched. `check=False` hands the
    failure back instead of raising.

    Output is buffered in memory, and RIsearch prints a line per hit, which a
    large run reaches tens of millions of. Read stdout incrementally for those.
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
