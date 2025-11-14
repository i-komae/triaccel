from __future__ import annotations

from typing import Optional, Any, Dict, Sequence

from . import viz as viz
from . import _core as _core


def simulate(
    M: int,
    sites: Sequence[dict],
    d_deg: float,
    seed: Optional[int] = None,
    bins: Optional[int | Sequence[int]] = None,
    return_histograms: bool = False,
    threads: Optional[int] = None,
    progress: bool = True,
    debug: bool = False,
    cluster_size: int = 3,
) -> Dict[str, Any]:
    """Public API wrapper around _core.simulate with optional debug rendering.

    When debug=True, automatically renders hammer PDFs from the generated logs.
    """
    # Validate and adapt M to C++ int range. pybind11 cannot pass Python ints
    # larger than 32-bit 'int' to our C++ signature. When debug=True, keep
    # behavior that the backend clamps to 1000 by ensuring we pass a valid
    # 32-bit integer (INT_MAX) so the C++ side can emit the clamp message.
    INT_MAX = 2_147_483_647
    try:
        M_in = int(M)
    except Exception as e:
        raise ValueError(f"M must be an integer-like value: {e}")
    if M_in < 1:
        raise ValueError("M must be >= 1")
    if M_in > INT_MAX:
        if debug:
            M_in = INT_MAX  # allow backend to clamp and log
        else:
            raise ValueError(f"M is too large for 32-bit int (got {M_in}, max={INT_MAX})")

    res = _core.simulate(
        M_in,
        sites,
        d_deg,
        seed,
        bins,
        return_histograms,
        threads,
        progress,
        debug,
        cluster_size,
    )
    if debug:
        try:
            from ._debug_plot import render_all
            render_all(log_dir=res.get("log_dir"), base_log_dir="log", circle_radius_deg=2.5, fmt="pdf", overwrite=True)
        except Exception:
            # Rendering failure should not break simulation results
            pass
    return res


__all__ = ["simulate", "viz"]
