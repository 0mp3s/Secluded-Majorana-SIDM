#!/usr/bin/env python3
r"""
run_pipeline.py — The_derivative_of_Lagernizan_SIDM/
=====================================================
Pipeline runner for the Lagrangian derivative investigation.

Reads docs/execution_pipeline.csv and runs scripts in order,
checking dependencies and I/O between steps.

Mirrors the SIDM pipeline runner pattern but simplified:
  - No FastAPI/WebSocket (overkill for 6 scripts)
  - Sequential execution with dependency checks
  - Timestamped outputs → data/archive/
  - Full run log via RunLogger

Usage:
  python run_pipeline.py                 # run all steps
  python run_pipeline.py --step 4        # run step 4 only
  python run_pipeline.py --from 4        # run steps 4–6
  python run_pipeline.py --list          # show pipeline status
  python run_pipeline.py --step 4 5 6    # run specific steps
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

# ── path setup ────────────────────────────────────────────────────────
_HERE      = Path(__file__).resolve().parent
_SIDM_ROOT = _HERE.parent
_CORE      = _SIDM_ROOT / "core"
sys.path.insert(0, str(_SIDM_ROOT))
sys.path.insert(0, str(_CORE))

from core.run_logger import RunLogger

# ── constants ─────────────────────────────────────────────────────────
CSV_PATH     = _HERE / "docs" / "execution_pipeline.csv"
ARCHIVE_DIR  = _HERE / "data" / "archive"


# ══════════════════════════════════════════════════════════════════════
# PIPELINE STEP DATACLASS
# ══════════════════════════════════════════════════════════════════════

@dataclass
class PipelineStep:
    order:        int
    stage:        str
    dependencies: str          # "2;4" means depends on steps 2 and 4
    script:       str          # filename relative to _HERE
    reads:        str
    writes:       str          # stem for data/archive/ lookup
    description:  str

    # resolved at runtime
    missing_inputs: list[str] = field(default_factory=list)
    last_run_at:    Optional[str] = None
    status:         str = "pending"   # pending | running | passed | failed | skipped

    @property
    def label(self) -> str:
        return Path(self.script).name

    @property
    def dep_orders(self) -> list[int]:
        """Parse dependency string → list of step order numbers."""
        if not self.dependencies or self.dependencies.lower() == "none":
            return []
        parts = self.dependencies.replace(",", ";").split(";")
        return [int(p.strip()) for p in parts if p.strip().isdigit()]


# ══════════════════════════════════════════════════════════════════════
# LOAD PIPELINE FROM CSV
# ══════════════════════════════════════════════════════════════════════

def load_pipeline() -> list[PipelineStep]:
    """Read execution_pipeline.csv → list of PipelineStep."""
    if not CSV_PATH.exists():
        print(f"ERROR: Pipeline CSV not found: {CSV_PATH}")
        sys.exit(1)

    steps = []
    with open(CSV_PATH, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            step = PipelineStep(
                order=int(row["order"]),
                stage=row["stage"].strip(),
                dependencies=row.get("dependencies", "").strip(),
                script=row["script"].strip(),
                reads=row.get("reads", "").strip(),
                writes=row.get("writes", "").strip(),
                description=row.get("description", "").strip(),
            )
            steps.append(step)

    # Sort by order
    steps.sort(key=lambda s: s.order)

    # Resolve output timestamps
    for step in steps:
        stem = Path(step.writes).stem if step.writes else ""
        if stem:
            matches = sorted(ARCHIVE_DIR.glob(f"{stem}_*.csv"))
            if matches:
                # Extract timestamp from filename
                step.last_run_at = matches[-1].stat().st_mtime
                step.status = "passed"  # has output → assume passed

    # Check missing inputs
    for step in steps:
        _check_inputs(step)

    return steps


def _check_inputs(step: PipelineStep) -> None:
    """Check if all declared inputs exist."""
    if not step.reads or step.reads.lower() == "none":
        return

    missing = []
    for token in step.reads.split(","):
        token = token.strip()
        if not token:
            continue

        # Check as-is relative to _HERE
        if (_HERE / token).exists():
            continue

        # Check in archive (timestamped)
        stem = Path(token).stem
        ext = Path(token).suffix or ".csv"
        if list(ARCHIVE_DIR.glob(f"{stem}_*{ext}")):
            continue

        # Check relative to SIDM root
        if (_SIDM_ROOT / token).exists():
            continue

        missing.append(token)

    step.missing_inputs = missing


# ══════════════════════════════════════════════════════════════════════
# DISPLAY PIPELINE STATUS
# ══════════════════════════════════════════════════════════════════════

def show_pipeline(steps: list[PipelineStep]) -> None:
    """Print pipeline status table."""
    print()
    print("=" * 85)
    print("  DERIVATIVE INVESTIGATION — PIPELINE STATUS")
    print("=" * 85)
    print()

    header = f"  {'#':>2}  {'Stage':<28} {'Script':<35} {'Status':<10} {'Last Run'}"
    print(header)
    print("  " + "─" * 80)

    for step in steps:
        # Status icon
        if step.status == "passed":
            icon = "✓"
            ts = datetime.fromtimestamp(step.last_run_at).strftime("%Y-%m-%d %H:%M") \
                 if step.last_run_at else "—"
        elif step.status == "failed":
            icon = "✗"
            ts = "FAILED"
        elif step.status == "running":
            icon = "►"
            ts = "running..."
        else:
            icon = "○"
            ts = "—"

        # Missing inputs warning
        warn = ""
        if step.missing_inputs:
            warn = f"  ⚠ missing: {', '.join(step.missing_inputs)}"

        print(f"  {step.order:>2}  {step.stage:<28} {step.label:<35} {icon:<10} {ts}{warn}")

    print()
    n_done = sum(1 for s in steps if s.status == "passed")
    print(f"  Progress: {n_done}/{len(steps)} steps completed")
    print("=" * 85)
    print()


# ══════════════════════════════════════════════════════════════════════
# RUN A SINGLE STEP
# ══════════════════════════════════════════════════════════════════════

def run_step(step: PipelineStep, steps: list[PipelineStep],
             force: bool = False) -> bool:
    """Run one pipeline step. Returns True if passed."""

    script_path = _HERE / step.script
    if not script_path.exists():
        print(f"  ✗ Script not found: {step.script}")
        step.status = "failed"
        return False

    # Check dependencies
    if not force:
        for dep_order in step.dep_orders:
            dep_step = next((s for s in steps if s.order == dep_order), None)
            if dep_step and dep_step.status not in ("passed",):
                print(f"  ⚠ Dependency step {dep_order} ({dep_step.label}) "
                      f"not yet passed. Use --force to override.")
                step.status = "skipped"
                return False

    # Run
    step.status = "running"
    print(f"\n{'━' * 70}")
    print(f"  ► Step {step.order}: {step.description}")
    print(f"    Script: {step.script}")
    print(f"{'━' * 70}\n")

    t0 = time.time()

    result = subprocess.run(
        [sys.executable, str(script_path)],
        cwd=str(_HERE),
        timeout=1800,   # 30 min max
    )

    elapsed = time.time() - t0

    if result.returncode == 0:
        step.status = "passed"
        print(f"\n  ✓ Step {step.order} PASSED ({elapsed:.1f}s)")
        return True
    else:
        step.status = "failed"
        print(f"\n  ✗ Step {step.order} FAILED (exit code {result.returncode}, {elapsed:.1f}s)")
        return False


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Derivative investigation pipeline runner"
    )
    parser.add_argument("--list", action="store_true",
                        help="Show pipeline status and exit")
    parser.add_argument("--step", type=int, nargs="+",
                        help="Run specific step(s) by order number")
    parser.add_argument("--from", type=int, dest="from_step",
                        help="Run from step N to the end")
    parser.add_argument("--force", action="store_true",
                        help="Ignore dependency checks")
    args = parser.parse_args()

    steps = load_pipeline()

    # ── List mode ─────────────────────────────────────────────────────
    if args.list:
        show_pipeline(steps)
        return

    # ── Determine which steps to run ──────────────────────────────────
    if args.step:
        to_run = [s for s in steps if s.order in args.step]
        if not to_run:
            print(f"ERROR: No steps found with order {args.step}")
            sys.exit(1)
    elif args.from_step:
        to_run = [s for s in steps if s.order >= args.from_step]
    else:
        # Run all
        to_run = steps

    # ── Banner ────────────────────────────────────────────────────────
    print()
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  LAGRANGIAN DERIVATIVE INVESTIGATION — PIPELINE RUNNER         ║")
    print("║  Secluded Majorana SIDM: L = ½χ̄(i∂̸-m)χ + ½(∂φ)²-½m²φ²-(y/2)χ̄χφ  ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    print()
    print(f"  Steps to run: {[s.order for s in to_run]}")
    print(f"  Working dir:  {_HERE}")
    print()

    # ── Run with RunLogger ────────────────────────────────────────────
    with RunLogger(
        script="The_derivative_of_Lagernizan_SIDM/run_pipeline.py",
        stage="Pipeline Run",
        params={"steps": [s.order for s in to_run],
                "force": args.force},
        data_source="docs/execution_pipeline.csv",
    ) as rl:
        t_start = time.time()
        n_pass = 0
        n_fail = 0
        n_skip = 0

        for step in to_run:
            ok = run_step(step, steps, force=args.force)
            if ok:
                n_pass += 1
            elif step.status == "skipped":
                n_skip += 1
            else:
                n_fail += 1
                if not args.force:
                    print(f"\n  ⚠ Pipeline stopped at step {step.order}. "
                          f"Use --force to continue past failures.")
                    break

        elapsed = time.time() - t_start

        # ── Final summary ─────────────────────────────────────────────
        print()
        print("=" * 70)
        print("  PIPELINE SUMMARY")
        print("=" * 70)
        print(f"  Passed:  {n_pass}")
        print(f"  Failed:  {n_fail}")
        print(f"  Skipped: {n_skip}")
        print(f"  Total:   {elapsed:.1f}s")
        print()

        if n_fail == 0 and n_skip == 0:
            print("  ✓ ALL STEPS PASSED")
        elif n_fail > 0:
            print("  ✗ PIPELINE INCOMPLETE — see failures above")
        print("=" * 70)

        # Show updated status
        steps_reloaded = load_pipeline()
        show_pipeline(steps_reloaded)

        rl.set_notes(f"pass={n_pass}, fail={n_fail}, skip={n_skip}, "
                     f"elapsed={elapsed:.1f}s")


if __name__ == "__main__":
    main()
