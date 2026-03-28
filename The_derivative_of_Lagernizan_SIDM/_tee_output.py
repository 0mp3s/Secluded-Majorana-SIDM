"""Tee stdout to both terminal and output/{name}.txt."""
from __future__ import annotations

import contextlib
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_OUTPUT_DIR = _HERE / "output"


class _Tee:
    """Write to multiple streams simultaneously."""
    __slots__ = ("_streams", "encoding", "newline")

    def __init__(self, *streams):
        self._streams = streams
        # Mimic sys.stdout attributes that libraries may read
        self.encoding = getattr(streams[0], "encoding", "utf-8")
        self.newline = getattr(streams[0], "newline", None)

    def write(self, data):
        for s in self._streams:
            s.write(data)

    def flush(self):
        for s in self._streams:
            s.flush()

    def fileno(self):
        return self._streams[0].fileno()

    def isatty(self):
        return False


@contextlib.contextmanager
def tee_to_output(name: str):
    """Context manager: tees stdout to output/{name}.txt while it's active."""
    _OUTPUT_DIR.mkdir(exist_ok=True)
    fpath = _OUTPUT_DIR / f"{name}.txt"
    old_stdout = sys.stdout
    with open(fpath, "w", encoding="utf-8") as fh:
        sys.stdout = _Tee(old_stdout, fh)
        try:
            yield fpath
        finally:
            sys.stdout = old_stdout
