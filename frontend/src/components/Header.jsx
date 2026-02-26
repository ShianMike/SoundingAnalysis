import "./Header.css";

export default function Header() {
  return (
    <header className="header">
      <div className="header-inner">
        <div className="header-brand">
          <div>
            <h1 className="header-title">Sounding Analysis</h1>
            <p className="header-subtitle">Vertical Profile Analysis Tool</p>
          </div>
        </div>
        <div className="header-meta">
          <span className="header-tag">IEM</span>
          <span className="header-tag">UWyo</span>
          <span className="header-tag">RAP</span>
          <span className="header-tag">BUFKIT</span>
          <span className="header-tag">ERA5</span>
          <span className="header-tag">ACARS</span>
        </div>
      </div>
    </header>
  );
}
