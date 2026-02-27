import { useState } from "react";
import { createPortal } from "react-dom";
import { Github, MessageSquarePlus, X, Send, Loader2 } from "lucide-react";
import "./Header.css";

export default function Header({ showFeedback, onCloseFeedback }) {
  const [feedbackText, setFeedbackText] = useState("");
  const [feedbackType, setFeedbackType] = useState("suggestion");
  const [submitted, setSubmitted] = useState(false);
  const [sending, setSending] = useState(false);

  const handleSubmitFeedback = async (e) => {
    e.preventDefault();
    if (!feedbackText.trim()) return;
    setSending(true);
    try {
      const API_BASE = import.meta.env.VITE_API_URL || "";
      const res = await fetch(`${API_BASE}/api/feedback`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          type: feedbackType,
          text: feedbackText.trim(),
          userAgent: navigator.userAgent,
        }),
      });
      if (!res.ok) throw new Error("Failed to send");
      setSubmitted(true);
      setFeedbackText("");
      setTimeout(() => {
        setSubmitted(false);
        onCloseFeedback();
      }, 2000);
    } catch {
      // Fallback to localStorage if backend is unreachable
      const feedbackStore = JSON.parse(localStorage.getItem("sounding_feedback") || "[]");
      feedbackStore.push({
        id: Date.now(),
        type: feedbackType,
        text: feedbackText.trim(),
        date: new Date().toISOString(),
      });
      localStorage.setItem("sounding_feedback", JSON.stringify(feedbackStore));
      setSubmitted(true);
      setFeedbackText("");
      setTimeout(() => {
        setSubmitted(false);
        onCloseFeedback();
      }, 2000);
    } finally {
      setSending(false);
    }
  };

  return (
    <>
      {showFeedback && createPortal(
        <div className="fb-overlay" onClick={() => onCloseFeedback()}>
          <div className="fb-modal" onClick={(e) => e.stopPropagation()}>
            <div className="fb-header">
              <h3>Feedback & Suggestions</h3>
              <button className="fb-close" onClick={() => onCloseFeedback()}>
                <X size={16} />
              </button>
            </div>

            {submitted ? (
              <div className="fb-success">
                <span className="fb-check">✓</span>
                <p>Thanks for your feedback!</p>
              </div>
            ) : (
              <form onSubmit={handleSubmitFeedback} className="fb-form">
                <div className="fb-type-row">
                  {[
                    { id: "suggestion", label: "💡 Suggestion" },
                    { id: "bug", label: "🐛 Bug Report" },
                    { id: "feature", label: "✨ Feature" },
                  ].map((t) => (
                    <button
                      key={t.id}
                      type="button"
                      className={`fb-type-btn ${feedbackType === t.id ? "active" : ""}`}
                      onClick={() => setFeedbackType(t.id)}
                    >
                      {t.label}
                    </button>
                  ))}
                </div>
                <textarea
                  className="fb-textarea"
                  rows={4}
                  placeholder="Tell us what you think, report a bug, or suggest a feature..."
                  value={feedbackText}
                  onChange={(e) => setFeedbackText(e.target.value)}
                  autoFocus
                />
                <button
                  type="submit"
                  className="fb-submit"
                  disabled={!feedbackText.trim() || sending}
                >
                  {sending ? (
                    <><Loader2 size={14} className="spin" /> Sending...</>
                  ) : (
                    <><Send size={14} /> Send Feedback</>
                  )}
                </button>
              </form>
            )}
          </div>
        </div>,
        document.body
      )}
    </>
  );
}
