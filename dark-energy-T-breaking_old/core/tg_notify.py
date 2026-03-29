"""
tg_notify.py — Telegram notifications for dark-energy-T-breaking tests.

Reads credentials from core/tg_config.json (gitignored).
If the file is missing, notifications are silently skipped.

Usage:
    from core.tg_notify import notify
    notify("✅ Test 19 done! ΔN_eff=0.153")

tg_config.json format:
    {"token": "YOUR_BOT_TOKEN", "chat_id": YOUR_CHAT_ID}
"""
import json
import pathlib
import requests

_CONFIG_PATH = pathlib.Path(__file__).parent / "tg_config.json"


def notify(msg: str, silent: bool = False) -> bool:
    """Send a Telegram message. Returns True on success, False on failure."""
    try:
        cfg = json.loads(_CONFIG_PATH.read_text())
        api = f"https://api.telegram.org/bot{cfg['token']}/sendMessage"
        r = requests.post(api, json={
            "chat_id": cfg["chat_id"],
            "text": msg,
            "disable_notification": silent,
        }, timeout=10)
        return r.ok
    except Exception:
        return False


if __name__ == "__main__":
    ok = notify("✅ tg_notify.py עובד! הודעת בדיקה מ-dark-energy-T-breaking")
    print("sent" if ok else "FAILED (no tg_config.json or network error)")
