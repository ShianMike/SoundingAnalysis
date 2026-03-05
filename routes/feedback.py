"""
Feedback routes: submit and retrieve user feedback.
"""
import html
import json
import os
import re
from datetime import datetime, timezone

from flask import Blueprint, jsonify, request

bp = Blueprint("feedback", __name__)

# Use /tmp on Cloud Run (ephemeral but always writable)
_fb_dir = "/tmp" if os.environ.get("K_SERVICE") else os.path.dirname(os.path.dirname(__file__))
FEEDBACK_FILE = os.path.join(_fb_dir, "feedback.json")

# Limits
_MAX_TEXT_LEN = 2000
_MAX_UA_LEN = 500
_ALLOWED_TYPES = {"suggestion", "bug", "other"}


def _sanitize(text: str, max_len: int = _MAX_TEXT_LEN) -> str:
    """Strip control chars, HTML entities, and enforce length."""
    text = text[:max_len]
    text = html.escape(text)
    text = re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f]", "", text)
    return text.strip()


def _load_feedback():
    if os.path.exists(FEEDBACK_FILE):
        try:
            with open(FEEDBACK_FILE, "r", encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            return []
    return []


def _save_feedback(data):
    with open(FEEDBACK_FILE, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


@bp.route("/api/feedback", methods=["POST"])
def submit_feedback():
    """Store user feedback/suggestions."""
    body = request.get_json(force=True)
    text = (body.get("text") or "").strip()
    if not text:
        return jsonify({"error": "Feedback text is required"}), 400

    fb_type = body.get("type", "suggestion")
    if fb_type not in _ALLOWED_TYPES:
        fb_type = "other"

    entry = {
        "id": int(datetime.now(timezone.utc).timestamp() * 1000),
        "type": fb_type,
        "text": _sanitize(text),
        "date": datetime.now(timezone.utc).isoformat(),
        "userAgent": _sanitize(body.get("userAgent", ""), _MAX_UA_LEN),
    }

    print(f"[FEEDBACK] type={entry['type']} | {entry['text'][:200]}")

    feedback = _load_feedback()
    feedback.append(entry)
    _save_feedback(feedback)

    return jsonify({"ok": True, "id": entry["id"]})


@bp.route("/api/feedback", methods=["GET"])
def get_feedback():
    """Retrieve all feedback (for admin/developer review)."""
    return jsonify(_load_feedback())
