import { useState } from "react";
import { Github, MessageSquarePlus, X, Send, Loader2 } from "lucide-react";
import "./Header.css";

export default function Header() {
  const [showFeedback, setShowFeedback] = useState(false);
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
        setShowFeedback(false);
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
        setShowFeedback(false);
      }, 2000);
    } finally {
      setSending(false);
    }
  };

  return (
    <header className="header">
      <div className="header-inner">
        <div className="header-brand">
          <div className="header-logo">
            <svg width="28" height="28" viewBox="0 0 28 28" fill="none" xmlns="http://www.w3.org/2000/svg">
              <rect x="2" y="2" width="24" height="24" rx="4" fill="rgba(59,130,246,0.12)" stroke="#3b82f6" strokeWidth="1.5"/>
              <path d="M8 22 L10 18 L11 16 L12 13 L14 10 L16 8 L19 5" stroke="#ef4444" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" fill="none"/>
              <path d="M6 22 L8 19 L9 17 L9.5 15 L10 13 L10 11 L10.5 9 L11 6" stroke="#22c55e" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" fill="none"/>
              <line x1="21" y1="10" x2="21" y2="22" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
              <line x1="21" y1="10" x2="24" y2="8" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
              <line x1="21" y1="13" x2="24" y2="11" stroke="#60a5fa" strokeWidth="1.2" strokeLinecap="round"/>
            </svg>
          </div>
          <div>
            <h1 className="header-title">Sounding Analysis</h1>
            <p className="header-subtitle">Atmospheric Vertical Profile Tool</p>
          </div>
        </div>

        <div className="header-right">
          <div className="header-tags">
            <span className="header-tag">IEM</span>
            <span className="header-tag">UWyo</span>
            <span className="header-tag">RAP</span>
            <span className="header-tag">BUFKIT</span>
            <span className="header-tag">ERA5</span>
            <span className="header-tag">ACARS</span>
          </div>

          <div className="header-actions">
            <button
              className={`header-action-btn ${showFeedback ? "active" : ""}`}
              onClick={() => setShowFeedback((v) => !v)}
              title="Send feedback or suggestion"
            >
              <MessageSquarePlus size={16} />
            </button>
            <a
              href="https://github.com/ShianMike/SoundingAnalysis"
              target="_blank"
              rel="noopener noreferrer"
              className="header-action-btn"
              title="View on GitHub"
            >
              <Github size={16} />
            </a>
          </div>
        </div>
      </div>

      {showFeedback && (
        <div className="fb-overlay" onClick={() => setShowFeedback(false)}>
          <div className="fb-modal" onClick={(e) => e.stopPropagation()}>
            <div className="fb-header">
              <h3>Feedback & Suggestions</h3>
              <button className="fb-close" onClick={() => setShowFeedback(false)}>
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
        </div>
      )}
    </header>
  );
}
