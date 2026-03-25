"""
Auto-patcher: adds Telegram notify block to all pipeline scripts with __main__.
Run once from project root.
"""
import os

PROJECT_ROOT = r"C:\Users\omerp\source\omer_mind\Secluded-Majorana-SIDM"
CORE_DIR = os.path.join(PROJECT_ROOT, "core")

# Files to skip (modules/utilities, not runnable pipeline scripts)
SKIP_FILES = {
    "__init__.py",
    "config_loader.py",
    "global_config.py",
    "output_manager.py",
    "run_logger.py",
    "tg_notify.py",
    "setup.py",
}

SKIP_DIRS = {".venv", "__pycache__", ".git"}


def get_core_path_expr(file_path):
    """Return the os.path.join(...) expression that points to core/ at runtime."""
    file_dir = os.path.dirname(os.path.abspath(file_path))
    rel = os.path.relpath(CORE_DIR, file_dir)  # e.g. '..\\core' or '..\\..\\core'
    parts = rel.split(os.sep)                   # ['..', 'core'] or ['..', '..', 'core']
    parts_str = ", ".join(repr(p) for p in parts)
    return f"_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), {parts_str})"


patched = []
skipped_reason = {}

for root, dirs, files in os.walk(PROJECT_ROOT):
    # Prune unwanted dirs in-place
    dirs[:] = [d for d in dirs if d not in SKIP_DIRS]

    for fname in sorted(files):
        if not fname.endswith(".py"):
            continue
        if fname in SKIP_FILES:
            continue

        fpath = os.path.join(root, fname)

        with open(fpath, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()

        # Already patched
        if "tg_notify" in content:
            skipped_reason[fpath] = "already has tg_notify"
            continue

        # No __main__ block → not a runnable script
        has_main = (
            "if __name__ == '__main__':" in content
            or 'if __name__ == "__main__":' in content
        )
        if not has_main:
            skipped_reason[fpath] = "no __main__ block"
            continue

        script_name = fname[:-3]  # strip .py

        file_dir = os.path.dirname(os.path.abspath(fpath))
        is_in_core = os.path.normpath(file_dir) == os.path.normpath(CORE_DIR)

        if is_in_core:
            block = (
                "\n\nif __name__ == '__main__':\n"
                "    try:\n"
                "        from tg_notify import notify\n"
                f'        notify("\\u2705 {script_name} done!")\n'
                "    except Exception:\n"
                "        pass\n"
            )
        else:
            core_expr = get_core_path_expr(fpath)
            block = (
                "\n\nif __name__ == '__main__':\n"
                "    try:\n"
                "        import sys as _sys, os as _os\n"
                f"        _sys.path.insert(0, {core_expr})\n"
                "        from tg_notify import notify\n"
                f'        notify("\\u2705 {script_name} done!")\n'
                "    except Exception:\n"
                "        pass\n"
            )

        with open(fpath, "a", encoding="utf-8") as f:
            f.write(block)

        patched.append(fpath)

# Report
print(f"\n{'='*60}")
print(f"PATCHED  {len(patched)} files:")
for p in patched:
    rel = os.path.relpath(p, PROJECT_ROOT)
    print(f"  + {rel}")

print(f"\nSKIPPED  {len(skipped_reason)} files:")
for p, reason in sorted(skipped_reason.items()):
    rel = os.path.relpath(p, PROJECT_ROOT)
    print(f"  - {rel}  ({reason})")
