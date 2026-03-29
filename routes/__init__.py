"""
Routes package — registers all Flask blueprints.
"""
from .meta import bp as meta_bp
from .sounding_routes import bp as sounding_bp
from .analysis import bp as analysis_bp
from .wind import bp as wind_bp
from .risk import bp as risk_bp
from .spc import bp as spc_bp
from .feedback import bp as feedback_bp
from .overlays import bp as overlays_bp

ALL_BLUEPRINTS = [
    meta_bp,
    sounding_bp,
    analysis_bp,
    wind_bp,
    risk_bp,
    spc_bp,
    feedback_bp,
    overlays_bp,
]


def register_all(app):
    """Register every blueprint on the Flask app."""
    for bp in ALL_BLUEPRINTS:
        app.register_blueprint(bp)
