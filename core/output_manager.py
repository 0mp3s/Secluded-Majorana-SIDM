#!/usr/bin/env python3
"""
core/output_manager.py
======================
Pipeline output management: timestamped writes + latest-file reads.

Every pipeline script that writes a CSV or NPY file consumed by a
downstream script should use:

    from output_manager import timestamped_path, get_latest

    # WRITE — always creates a new file, never overwrites previous runs
    out = timestamped_path("all_viable_raw_v8")
    # → data/archive/all_viable_raw_v8_2026_03_24.csv

    # READ — always loads the most recent run
    inp = get_latest("all_viable_raw_v8")
    # → data/archive/all_viable_raw_v8_2026_03_24.csv  (latest)

Rules
-----
- Date format: yyyy_mm_dd  (lexicographic sort == chronological sort)
- Default archive directory: <repo_root>/data/archive/
- MCMC-specific files (npy/csv):  pass  archive=_MCMC_ARCHIVE
- Downstream scripts never hardcode filenames — they always call get_latest()
"""
from __future__ import annotations

import glob
import pathlib
from datetime import date
from typing import Optional

# ---------------------------------------------------------------------------
#  Repo-relative archive locations
# ---------------------------------------------------------------------------
_REPO_ROOT    = pathlib.Path(__file__).resolve().parent.parent
_ARCHIVE      = _REPO_ROOT / "data" / "archive"          # shared pipeline CSVs
_MCMC_ARCHIVE = _REPO_ROOT / "stats_mcmc" / "output" / "archive"  # mcmc npy/csv


def _today() -> str:
    """Return today's date as yyyy_mm_dd."""
    return date.today().strftime("%Y_%m_%d")


def timestamped_path(
    stem: str,
    ext: str = ".csv",
    archive: Optional[pathlib.Path] = None,
) -> pathlib.Path:
    """
    Return a new, timestamped path inside the archive directory.

    Parameters
    ----------
    stem    : base filename without date/extension, e.g. "all_viable_raw_v8"
    ext     : file extension including dot, default ".csv"
    archive : archive directory (default: data/archive/)

    Returns
    -------
    pathlib.Path  e.g.  data/archive/all_viable_raw_v8_2026_03_24.csv

    Notes
    -----
    The archive directory is created automatically if it does not exist.
    Two runs on the same calendar day produce the same path — which is
    intentional: re-running on the same day replaces that day's output
    while keeping all previous days intact.
    """
    d = archive if archive is not None else _ARCHIVE
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{stem}_{_today()}{ext}"


def get_latest(
    stem: str,
    ext: str = ".csv",
    archive: Optional[pathlib.Path] = None,
) -> pathlib.Path:
    """
    Return the path of the most-recent archive file matching
    ``stem_yyyy_mm_dd{ext}``.

    Parameters
    ----------
    stem    : base filename without date/extension, e.g. "all_viable_raw_v8"
    ext     : file extension including dot, default ".csv"
    archive : archive directory (default: data/archive/)

    Returns
    -------
    pathlib.Path  e.g.  data/archive/all_viable_raw_v8_2026_03_24.csv

    Raises
    ------
    FileNotFoundError
        If no matching file exists in the archive.  This means the upstream
        pipeline script has not been run yet.
    """
    d = archive if archive is not None else _ARCHIVE
    pattern = str(d / f"{stem}_*{ext}")
    matches = sorted(glob.glob(pattern))   # yyyy_mm_dd → lexicographic == chronological
    if not matches:
        raise FileNotFoundError(
            f"No archive file found matching:\n  {pattern}\n"
            "Run the upstream pipeline script first to generate this data."
        )
    return pathlib.Path(matches[-1])
