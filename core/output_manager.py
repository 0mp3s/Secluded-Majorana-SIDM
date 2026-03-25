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
    # → data/archive/all_viable_raw_v8_2026_03_24_r12.csv

    # READ — always loads the most recent run
    inp = get_latest("all_viable_raw_v8")
    # → data/archive/all_viable_raw_v8_2026_03_24_r12.csv  (latest)

Rules
-----
- Filename format: stem_yyyy_mm_dd_rNNN.ext  (run_id from RunLogger)
- Fallback (no active RunLogger):  stem_yyyy_mm_dd.ext  (backward compat)
- Default archive directory: <repo_root>/data/archive/
- MCMC-specific files (npy/csv):  pass  archive=_MCMC_ARCHIVE
- Downstream scripts never hardcode filenames — they always call get_latest()
"""
from __future__ import annotations

import glob
import pathlib
import subprocess
import warnings
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


def _git_hash_short() -> str:
    """Return first 7 chars of HEAD, or 'no-git'."""
    try:
        r = subprocess.run(
            ["git", "rev-parse", "--short=7", "HEAD"],
            capture_output=True, text=True, cwd=str(_REPO_ROOT), timeout=3,
        )
        h = r.stdout.strip()
        return h if h else "no-git"
    except Exception:
        return "no-git"


def timestamped_path(
    stem: str,
    ext: str = ".csv",
    archive: Optional[pathlib.Path] = None,
) -> pathlib.Path:
    """
    Return a new, timestamped path inside the archive directory.

    If called inside an active RunLogger context, the run_id is appended
    (e.g. ``_r12``), guaranteeing uniqueness even for multiple runs on
    the same calendar day.  Without an active RunLogger the old format
    ``stem_yyyy_mm_dd.ext`` is used for backward compatibility.

    Parameters
    ----------
    stem    : base filename without date/extension, e.g. "all_viable_raw_v8"
    ext     : file extension including dot, default ".csv"
    archive : archive directory (default: data/archive/)
    """
    from core.run_logger import get_active_logger  # local import to avoid circular

    d = archive if archive is not None else _ARCHIVE
    d.mkdir(parents=True, exist_ok=True)

    rl = get_active_logger()
    if rl is not None and rl.run_id is not None:
        name = f"{stem}_{_today()}_r{rl.run_id}{ext}"
    else:
        name = f"{stem}_{_today()}{ext}"
    return d / name


def get_latest(
    stem: str,
    ext: str = ".csv",
    archive: Optional[pathlib.Path] = None,
) -> pathlib.Path:
    """
    Return the path of the most-recent archive file matching *stem*.

    Matches both old (``stem_yyyy_mm_dd.ext``) and new
    (``stem_yyyy_mm_dd_rNNN.ext``) filenames, sorted lexicographically
    so the newest is always last.

    Side-effects
    -------------
    * If an active RunLogger exists, the returned path is automatically
      registered as a data source (``add_data_source``).
    * If the current git hash differs from the hash recorded in the run
      that *produced* the file, a warning is emitted.  (Requires the
      run log; silently skipped if lookup fails.)
    """
    from core.run_logger import get_active_logger  # local import

    d = archive if archive is not None else _ARCHIVE
    pattern = str(d / f"{stem}_*{ext}")
    matches = sorted(glob.glob(pattern))   # lexicographic == chronological
    if not matches:
        raise FileNotFoundError(
            f"No archive file found matching:\n  {pattern}\n"
            "Run the upstream pipeline script first to generate this data."
        )

    latest = pathlib.Path(matches[-1])

    # ── auto-register as data source ────────────────────────────────────
    rl = get_active_logger()
    if rl is not None:
        rl.add_data_source(str(latest))

    # ── git hash consistency warning ────────────────────────────────────
    _warn_if_hash_mismatch(latest)

    return latest


def _warn_if_hash_mismatch(path: pathlib.Path) -> None:
    """Emit a warning if *path* was produced by a different git commit."""
    try:
        from core.run_logger import _LOG_PATH
        import csv

        if not _LOG_PATH.exists():
            return

        fname = path.name
        producing_hash = None
        with open(_LOG_PATH, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                if fname in row.get("output_files", ""):
                    producing_hash = row.get("git_hash", "")

        if not producing_hash or producing_hash == "no-git":
            return

        current = _git_hash_short()
        if current != "no-git" and current != producing_hash:
            warnings.warn(
                f"Git hash mismatch for {fname}:\n"
                f"  produced at {producing_hash}, current HEAD is {current}.\n"
                "  The code may have changed since this data was generated.",
                stacklevel=3,
            )
    except Exception:
        pass  # non-critical — never break the pipeline for a warning
