"""
tg_notify.py — Telegram notifications for long-running scans.

Usage:
    from tg_notify import notify
    notify("Scan done! raw=157k rep=580 elapsed=9000s")
"""
import requests

_TOKEN   = "8637696133:AAFmJgrPJzBZZk_Wp91qz1psYq3rp1nXZ3M"
_CHAT_ID = 5810172011
_API     = f"https://api.telegram.org/bot{_TOKEN}/sendMessage"


def notify(msg: str, silent: bool = False) -> bool:
    """Send a Telegram message. Returns True on success, False on failure."""
    try:
        r = requests.post(_API, json={
            "chat_id": _CHAT_ID,
            "text": msg,
            "disable_notification": silent,
        }, timeout=10)
        return r.ok
    except Exception:
        return False


if __name__ == "__main__":
    ok = notify("✅ tg_notify.py עובד! הודעת בדיקה מ-Secluded-Majorana-SIDM")
    print("sent" if ok else "FAILED")
