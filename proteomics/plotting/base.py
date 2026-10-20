"""Shared matplotlib helpers: thread safety, figure creation, PNG rendering, labels."""
import io
import threading
from functools import wraps

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg


# ------------------------------------------------------------------
# Thread safety
# Streamlit runs script in threads. Every time a widget changes the whole script reruns. 
#matplotlib can't handle overlapping threads. All figure building / rendering therefore goes
# through one re-entrant lock. Figures are created with the object-oriented
# Figure class instead of pyplot, so there is no shared "current figure" and
# nothing to close afterwards (a figure is freed once it is no longer referenced). 
# ------------------------------------------------------------------
MPL_LOCK = threading.RLock()


def locked(func):
    """Run a plotting function while holding MPL_LOCK."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        with MPL_LOCK:
            return func(*args, **kwargs)
    return wrapper


def new_figure(**kwargs):
    """A standalone Figure with an Agg canvas (pyplot-free, so also usable from threads)."""
    fig = Figure(**kwargs)
    FigureCanvasAgg(fig)
    return fig


def fig_to_png(fig, dpi=150):
    """Render a figure to PNG bytes (under the lock). Used for both display and download."""
    with MPL_LOCK:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
        return buf.getvalue()


def fmt_fraction(f):
    """1.0 -> '1', 2 -> '2', 'x' -> 'x' (clean fraction labels for titles / axes / file names)."""
    try:
        if float(f).is_integer():
            return str(int(f))
    except (TypeError, ValueError):
        pass
    return str(f)


def with_method(title, method):
    """Prefix a plot title with the fractionation method, e.g. 'RP high pH | 3D Waterfall'."""
    return f"{method} | {title}" if method else title
