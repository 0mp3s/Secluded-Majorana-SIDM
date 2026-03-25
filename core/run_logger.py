#!/usr/bin/env python3
"""
core/run_logger.py
==================
Append-only run log — records every pipeline script execution to
docs/runs_log.csv.

Each row = one script run. The file is created on first use and
never overwritten (only appended to).

────────────────────────────────────────────────────────────────
Usage — context manager (recommended)
────────────────────────────────────────────────────────────────
    from run_logger import RunLogger

    with RunLogger("core/v22_raw_scan.py",
                   stage="0 - VPM Scan",
                   params={"m_chi_min_GeV": 5.0, "n_grid": 500},
                   data_source="grid scan from scratch") as rl:

        # ... your computation ...

        rl.add_output("data/archive/all_viable_raw_v8_2026_03_25.csv")
        rl.set_n_viable(8743)
        rl.set_notes("First clean run — Harvey 0.47 bound")

    # On __exit__: status set to OK (or FAILED if exception),
    # duration captured, row appended to runs_log.csv.

────────────────────────────────────────────────────────────────
Usage — one-shot call (for simple/one-liner scripts)
────────────────────────────────────────────────────────────────
    from run_logger import log_run
    log_run(
        script="vpm_scan/sidm_cuts_sensitivity.py",
        stage="Sensitivity",
        params={"sd_lo": 0.5, "sd_hi": 10.0, "sc_hi": 0.47},
        output_files=["docs/sidm_cuts_sensitivity_2026_03_25_0225.txt"],
        n_viable=131,
        status="OK",
        duration_sec=4.3,
        notes="Fine sweep 0.01→1.01",
    )

────────────────────────────────────────────────────────────────
CSV columns
────────────────────────────────────────────────────────────────
  run_id          — auto-incremented (max existing + 1)
  timestamp_start — YYYY-MM-DD HH:MM:SS
  script          — relative path from repo root
  stage           — pipeline stage label (from execution_pipeline.csv)
  params_json     — JSON object with every relevant parameter
  data_source     — input file(s) read, comma-separated
  output_files    — output files written (full timestamped names)
  n_viable        — # viable BPs found (empty if not applicable)
  status          — OK | FAILED | PARTIAL
  duration_sec    — wall-clock seconds (2 dp)
  git_hash        — 7-char HEAD hash (or "no-git")
  notes           — free text
"""
from __future__ import annotations

import csv
import json
import pathlib
import subprocess
import threading
import time
from datetime import datetime
from typing import Any, Dict, List, Optional

# ── thread-local active RunLogger instance ────────────────────────────────────
_active_lock = threading.Lock()
_active_logger: Optional["RunLogger"] = None


def get_active_logger() -> Optional["RunLogger"]:
    """Return the currently active RunLogger (inside a `with` block), or None."""
    return _active_logger

# ── paths ────────────────────────────────────────────────────────────────────
_REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
_LOG_PATH  = _REPO_ROOT / "docs" / "runs_log.csv"

_COLUMNS = [
    "run_id",
    "timestamp_start",
    "script",
    "stage",
    "params_json",
    "data_source",
    "output_files",
    "n_viable",
    "status",
    "duration_sec",
    "git_hash",
    "notes",
]


# ── helpers ───────────────────────────────────────────────────────────────────

def _git_hash() -> str:
    """Return first 7 chars of HEAD commit hash, or 'no-git'."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short=7", "HEAD"],
            capture_output=True, text=True, cwd=str(_REPO_ROOT), timeout=3,
        )
        h = result.stdout.strip()
        return h if h else "no-git"
    except Exception:
        return "no-git"


def _next_run_id() -> int:
    """Return max existing run_id + 1, or 1 if log is empty/missing."""
    if not _LOG_PATH.exists():
        return 1
    max_id = 0
    with open(_LOG_PATH, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                max_id = max(max_id, int(row["run_id"]))
            except (ValueError, KeyError):
                pass
    return max_id + 1


def _ensure_header():
    """Create docs/runs_log.csv with header row if it doesn't exist yet."""
    _LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    if not _LOG_PATH.exists():
        with open(_LOG_PATH, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(_COLUMNS)


def _append_row(row: Dict[str, Any]):
    """Append one row dict to runs_log.csv (creates file+header if needed)."""
    _ensure_header()
    with open(_LOG_PATH, "a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=_COLUMNS)
        writer.writerow(row)


# ── context manager ───────────────────────────────────────────────────────────

class RunLogger:
    """
    Context manager that auto-records a pipeline run to runs_log.csv.

    Parameters
    ----------
    script      : relative path from repo root, e.g. "core/v22_raw_scan.py"
    stage       : pipeline stage label, e.g. "0 - VPM Scan"
    params      : dict of key parameter values (will be JSON-serialised)
    data_source : string describing input data (file paths or "grid scan…")
    """

    def __init__(
        self,
        script: str,
        stage: str = "",
        params: Optional[Dict[str, Any]] = None,
        data_source: str = "",
    ):
        self.script      = script
        self.stage       = stage
        self.params      = params or {}
        self.data_source = data_source

        self._outputs: List[str] = []
        self._data_sources: List[str] = []
        self._n_viable: Optional[int] = None
        self._notes: str = ""
        self._status: str = "OK"

        self._start_ts: Optional[str] = None
        self._start_t:  Optional[float] = None
        self._run_id: Optional[int] = None

    # ── mutators called inside `with` block ──────────────────────────────────

    @property
    def run_id(self) -> Optional[int]:
        """Return the run_id assigned at __enter__, or None if not yet entered."""
        return self._run_id

    def add_output(self, path: str):
        """Register an output file (call once per file written)."""
        self._outputs.append(str(path))

    def add_data_source(self, path: str):
        """Register an input file (auto-called by get_latest)."""
        s = str(path)
        if s not in self._data_sources:
            self._data_sources.append(s)

    def set_n_viable(self, n: int):
        """Record the number of viable benchmark points found."""
        self._n_viable = int(n)

    def set_notes(self, text: str):
        """Set free-text notes for this run."""
        self._notes = text

    def set_status(self, status: str):
        """Manually override status (OK | FAILED | PARTIAL)."""
        self._status = status

    # ── context protocol ─────────────────────────────────────────────────────

    def __enter__(self) -> "RunLogger":
        global _active_logger
        self._start_ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self._start_t  = time.monotonic()
        self._run_id   = _next_run_id()
        with _active_lock:
            _active_logger = self
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> bool:
        global _active_logger
        with _active_lock:
            _active_logger = None

        duration = round(time.monotonic() - self._start_t, 2)

        if exc_type is not None:
            self._status = "FAILED"
            if not self._notes:
                self._notes = f"{exc_type.__name__}: {exc_val}"

        # Merge auto-collected data sources with manually provided one
        all_sources = []
        if self.data_source:
            all_sources.append(self.data_source)
        all_sources.extend(self._data_sources)

        row = {
            "run_id":          self._run_id,
            "timestamp_start": self._start_ts,
            "script":          self.script,
            "stage":           self.stage,
            "params_json":     json.dumps(self.params, separators=(",", ":")),
            "data_source":     " | ".join(all_sources),
            "output_files":    " | ".join(self._outputs),
            "n_viable":        "" if self._n_viable is None else self._n_viable,
            "status":          self._status,
            "duration_sec":    duration,
            "git_hash":        _git_hash(),
            "notes":           self._notes,
        }
        _append_row(row)

        # Do not suppress exceptions
        return False


# ── one-shot convenience function ─────────────────────────────────────────────

def log_run(
    script: str,
    stage: str = "",
    params: Optional[Dict[str, Any]] = None,
    data_source: str = "",
    output_files: Optional[List[str]] = None,
    n_viable: Optional[int] = None,
    status: str = "OK",
    duration_sec: float = 0.0,
    notes: str = "",
):
    """
    Append a single log row without a context manager.

    Useful for scripts that cannot be wrapped in a `with` block,
    or for retroactively logging a known run.
    """
    row = {
        "run_id":          _next_run_id(),
        "timestamp_start": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "script":          script,
        "stage":           stage,
        "params_json":     json.dumps(params or {}, separators=(",", ":")),
        "data_source":     data_source,
        "output_files":    " | ".join(output_files or []),
        "n_viable":        "" if n_viable is None else n_viable,
        "status":          status,
        "duration_sec":    round(duration_sec, 2),
        "git_hash":        _git_hash(),
        "notes":           notes,
    }
    _append_row(row)
