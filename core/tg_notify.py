"""
tg_notify.py — Telegram notifications for long-running scans.

Reads credentials from core/tg_config.json (gitignored).
If the file is missing, notifications are silently skipped.

Usage:
    from tg_notify import notify
    notify("Scan done! raw=157k rep=580 elapsed=9000s")

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
    ok = notify("✅ tg_notify.py עובד! הודעת בדיקה מ-Secluded-Majorana-SIDM")
    print("sent" if ok else "FAILED")
