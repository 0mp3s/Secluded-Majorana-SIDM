#!/usr/bin/env python3
"""
core/output_manager.py
======================
Timestamped output paths for dark-energy-T-breaking tests.

Every test that writes a file (txt/csv/png) should use:

    from core.output_manager import timestamped_path, get_latest

    # WRITE — new timestamped file, never overwrites
    out = timestamped_path("test19_qcd_coincidence", ext=".txt")
    # → data/output/test19_qcd_coincidence_2026_03_27_r05.txt

    # READ — load the most recent run of a given stem
    inp = get_latest("test06_sidm_consistency")
    # → data/output/test06_sidm_consistency_2026_03_26_r03.txt

Rules
-----
- Format: stem_yyyy_mm_dd_rNNN.ext  (run_id from active RunLogger)
- Fallback (no active RunLogger): stem_yyyy_mm_dd.ext
- All outputs go to data/output/ by default
"""
from __future__ import annotations

import pathlib
from datetime import date
from typing import Optional

from run_logger import get_active_logger

_REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
_OUTPUT    = _REPO_ROOT / "data" / "output"


def _today() -> str:
    return date.today().strftime("%Y_%m_%d")


def timestamped_path(
    stem: str,
    ext: str = ".txt",
    output_dir: Optional[pathlib.Path] = None,
) -> pathlib.Path:
    """
    Return a new timestamped path for writing.

    If called inside an active RunLogger context, appends run_id:
        stem_2026_03_27_r05.txt
    Otherwise falls back to:
        stem_2026_03_27.txt
    """
    d = output_dir or _OUTPUT
    d.mkdir(parents=True, exist_ok=True)

    rl = get_active_logger()
    if rl is not None and rl.run_id is not None:
        name = f"{stem}_{_today()}_r{rl.run_id:02d}{ext}"
    else:
        name = f"{stem}_{_today()}{ext}"

    return d / name


def get_latest(
    stem: str,
    ext: str = ".txt",
    output_dir: Optional[pathlib.Path] = None,
) -> Optional[pathlib.Path]:
    """
    Return the most recent file matching stem_*.ext in output_dir.
    Returns None if no match found.
    """
    d = output_dir or _OUTPUT
    matches = sorted(d.glob(f"{stem}_*{ext}"))
    if not matches:
        return None
    return matches[-1]
