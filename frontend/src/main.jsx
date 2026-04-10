import React from 'react'
import ReactDOM from 'react-dom/client'
import App from './App.jsx'
import './index.css'

ReactDOM.createRoot(document.getElementById('root')).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>,
)

// Service worker: register in production, unregister in development
if ('serviceWorker' in navigator) {
  if (import.meta.env.PROD) {
    window.addEventListener('load', () => {
      const base = import.meta.env.BASE_URL || '/';
      navigator.serviceWorker.register(`${base}sw.js`).catch(() => {});
    });
  } else {
    // Unregister any stale SW from prior sessions so it doesn't intercept Vite dev requests
    navigator.serviceWorker.getRegistrations().then((regs) => {
      for (const r of regs) r.unregister();
    });
  }
}
