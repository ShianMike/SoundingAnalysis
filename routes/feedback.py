"""
Feedback routes: submit and retrieve user feedback.
"""
import json
import os
from datetime import datetime, timezone

from flask import Blueprint, jsonify, request

bp = Blueprint("feedback", __name__)

# Use /tmp on Cloud Run (ephemeral but always writable)
_fb_dir = "/tmp" if os.environ.get("K_SERVICE") else os.path.dirname(os.path.dirname(__file__))
FEEDBACK_FILE = os.path.join(_fb_dir, "feedback.json")


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

    entry = {
        "id": int(datetime.now(timezone.utc).timestamp() * 1000),
        "type": body.get("type", "suggestion"),
        "text": text,
        "date": datetime.now(timezone.utc).isoformat(),
        "userAgent": body.get("userAgent", ""),
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
