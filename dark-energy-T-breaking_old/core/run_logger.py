#!/usr/bin/env python3
"""
core/run_logger.py
==================
Append-only run log — records every test execution to docs/runs_log.csv.

Each row = one test run. The file is created on first use and never
overwritten (only appended to).

────────────────────────────────────────────────────────────────
Usage — context manager (recommended)
────────────────────────────────────────────────────────────────
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core'))
    from run_logger import RunLogger

    with RunLogger("sigma_radiative_stability.py",
                   test="Test 1 - Radiative Stability",
                   params={"m_chi_GeV": 94.07, "m_phi_MeV": 11.10}) as rl:
        # ... your computation ...
        rl.add_output("data/output/test01_sigma_radiative_2026_03_27.txt")
        rl.set_notes("Option D (dark axion) is only viable coupling")

    # On __exit__: status set to OK (or FAILED if exception),
    # duration captured, row appended to docs/runs_log.csv.

────────────────────────────────────────────────────────────────
Usage — one-shot call
────────────────────────────────────────────────────────────────
    from run_logger import log_run
    log_run(
        script="qcd_scale_coincidence.py",
        test="Test 19 - QCD Coincidence",
        params={"m_chi_MeV": 94.07, "Lambda_QCD_MeV": 200.0},
        output_files=["data/output/test19_qcd_2026_03_27.txt"],
        status="OK",
        duration_sec=1.2,
        notes="delta_N_eff=0.153 at T_D=200 MeV",
    )

────────────────────────────────────────────────────────────────
CSV columns
────────────────────────────────────────────────────────────────
  run_id          — auto-incremented (max existing + 1)
  timestamp_start — YYYY-MM-DD HH:MM:SS
  script          — filename (relative to repo root)
  test            — test label, e.g. "Test 19 - QCD Coincidence"
  params_json     — JSON object with key parameters
  output_files    — output files written (pipe-separated)
  status          — OK | FAILED | PARTIAL
  duration_sec    — wall-clock seconds (2 dp)
  git_hash        — 7-char HEAD hash (or "no-git")
  key_result      — one-line result summary (set via set_key_result)
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

# ── thread-local active RunLogger instance ───────────────────────────────────
_active_lock = threading.Lock()
_active_logger: Optional["RunLogger"] = None


def get_active_logger() -> Optional["RunLogger"]:
    """Return the currently active RunLogger (inside a `with` block), or None."""
    return _active_logger


# ── paths ─────────────────────────────────────────────────────────────────────
_REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
_LOG_PATH  = _REPO_ROOT / "docs" / "runs_log.csv"

_COLUMNS = [
    "run_id",
    "timestamp_start",
    "script",
    "test",
    "params_json",
    "output_files",
    "status",
    "duration_sec",
    "git_hash",
    "key_result",
    "notes",
]


# ── helpers ───────────────────────────────────────────────────────────────────

def _git_hash() -> str:
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
    if not _LOG_PATH.exists():
        return 1
    max_id = 0
    try:
        with open(_LOG_PATH, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                try:
                    max_id = max(max_id, int(row.get("run_id", 0)))
                except ValueError:
                    pass
    except Exception:
        pass
    return max_id + 1


def _ensure_log() -> None:
    _LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    if not _LOG_PATH.exists():
        with open(_LOG_PATH, "w", newline="", encoding="utf-8") as f:
            csv.DictWriter(f, fieldnames=_COLUMNS).writeheader()


def _append_row(row: dict) -> None:
    _ensure_log()
    with open(_LOG_PATH, "a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=_COLUMNS, extrasaction="ignore")
        w.writerow(row)


# ── RunLogger class ───────────────────────────────────────────────────────────

class RunLogger:
    """Context manager that records a single test execution to runs_log.csv."""

    def __init__(
        self,
        script: str,
        test: str = "",
        params: Optional[Dict[str, Any]] = None,
    ):
        self.script = script
        self.test   = test
        self.params = params or {}

        self._outputs: List[str]    = []
        self._key_result: str       = ""
        self._notes: str            = ""
        self._status: str           = "OK"
        self._start_ts: Optional[str]   = None
        self._start_t:  Optional[float] = None
        self._run_id:   Optional[int]   = None

    # ── mutators ─────────────────────────────────────────────────────────────

    @property
    def run_id(self) -> Optional[int]:
        return self._run_id

    def add_output(self, path: str):
        self._outputs.append(str(path))

    def set_key_result(self, text: str):
        """One-line summary of the main finding (shown in runs_log)."""
        self._key_result = text

    def set_notes(self, text: str):
        self._notes = text

    def set_status(self, status: str):
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

        row = {
            "run_id":          self._run_id,
            "timestamp_start": self._start_ts,
            "script":          self.script,
            "test":            self.test,
            "params_json":     json.dumps(self.params, separators=(",", ":")),
            "output_files":    " | ".join(self._outputs),
            "status":          self._status,
            "duration_sec":    duration,
            "git_hash":        _git_hash(),
            "key_result":      self._key_result,
            "notes":           self._notes,
        }
        _append_row(row)
        return False   # don't suppress exceptions


# ── one-shot helper ───────────────────────────────────────────────────────────

def log_run(
    script: str,
    test: str = "",
    params: Optional[Dict[str, Any]] = None,
    output_files: Optional[List[str]] = None,
    status: str = "OK",
    duration_sec: float = 0.0,
    key_result: str = "",
    notes: str = "",
) -> None:
    row = {
        "run_id":          _next_run_id(),
        "timestamp_start": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "script":          script,
        "test":            test,
        "params_json":     json.dumps(params or {}, separators=(",", ":")),
        "output_files":    " | ".join(output_files or []),
        "status":          status,
        "duration_sec":    round(duration_sec, 2),
        "git_hash":        _git_hash(),
        "key_result":      key_result,
        "notes":           notes,
    }
    _append_row(row)
