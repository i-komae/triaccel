from __future__ import annotations

import os
import re
from typing import Optional


def _latest_log_dir(base_log_dir: str) -> Optional[str]:
    if not os.path.isdir(base_log_dir):
        return None
    names = [n for n in os.listdir(base_log_dir) if os.path.isdir(os.path.join(base_log_dir, n))]
    if not names:
        return None
    names.sort(reverse=True)
    return os.path.join(base_log_dir, names[0])


def _wrap_pi(x: float) -> float:
    import math
    y = (x + math.pi) % (2 * math.pi)
    if y < 0:
        y += 2 * math.pi
    return y - math.pi


def _read_events(path: str):
    """Read sim_*_events.txt -> (ra_deg, dec_deg, in_cluster)
    Returns tuples of floats and ints for drawing.
    """
    ra, dec, inc = [], [], []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            # idx, site, ra_deg, dec_deg, in_cluster
            try:
                ra.append(float(parts[2]))
                dec.append(float(parts[3]))
                inc.append(int(parts[4]))
            except Exception:
                continue
    return ra, dec, inc


def _read_clique_members(path: str):
    """Read clique file and return a set of vertex indices that appear in any clique."""
    members = set()
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                try:
                    for tok in parts:
                        if tok:
                            members.add(int(tok))
                except Exception:
                    continue
    except FileNotFoundError:
        pass
    return members


def _plot_wrapped(ax, lam, lat, **kwargs):
    """Plot a polyline on a wrapped longitude domain by splitting at seams.
    Expects lam, lat as numpy arrays (radians), lam in [-pi, pi).
    """
    import numpy as np
    if len(lam) == 0:
        return
    dif = np.abs(np.diff(lam))
    # seam where jump is large (crossing +-pi). threshold slightly below pi
    seam = np.where(dif > (np.pi * 0.9))[0]
    start = 0
    for idx in seam:
        end = idx + 1
        if end > start:
            ax.plot(lam[start:end], lat[start:end], **kwargs)
        start = end
    ax.plot(lam[start:], lat[start:], **kwargs)


def _draw_small_circle(ax, ra_deg: float, dec_deg: float, radius_deg: float, *, color: str, lw: float = 0.5, segments: int = 180):
    import numpy as np
    import math

    alpha = math.radians(ra_deg)
    delta = math.radians(dec_deg)
    # Adjust radius so that the OUTER edge lies at 'radius_deg' including line width.
    # Approximate angular half-width (deg) from line width in points using axis width in pixels.
    try:
        renderer = ax.figure.canvas.get_renderer()
        width_px = ax.get_window_extent(renderer=renderer).width
        # Hammer x-range spans [-pi, pi] in data coords; rough angular per pixel at equator
        deg_per_px = (360.0) / float(width_px) if width_px > 0 else 0.0
        lw_ang_deg = lw * (deg_per_px / 72.0)  # points -> inches -> pixels already accounted; treat points as 1/72 inch
    except Exception:
        lw_ang_deg = 0.0
    rho = math.radians(max(1e-6, radius_deg - 0.5 * lw_ang_deg))
    phi = np.linspace(0.0, 2 * math.pi, segments, endpoint=True)
    sin_delta = math.sin(delta)
    cos_delta = math.cos(delta)
    sin_rho = math.sin(rho)
    cos_rho = math.cos(rho)
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)

    lat = np.arcsin(sin_delta * cos_rho + cos_delta * sin_rho * cos_phi)
    # central longitude alpha + atan2(...)
    lon = alpha + np.arctan2(sin_phi * sin_rho * cos_delta, cos_rho - sin_delta * np.sin(lat))
    # wrap lon to [-pi, pi) and apply RA inversion (lambda = -lon)
    lon_wrapped = np.vectorize(_wrap_pi)(lon)
    lam = -lon_wrapped
    _plot_wrapped(ax, lam, lat, color=color, linewidth=lw, alpha=1.0)


def render_all(log_dir: Optional[str] = None, base_log_dir: str = "log", circle_radius_deg: float = 2.5, fmt: str = "pdf", overwrite: bool = True) -> None:
    """Render hammer figures for all events logs in a debug directory.

    - Picks latest directory under base_log_dir when log_dir is None.
    - Reads text/sim_*_events.txt and writes fig/sim_*_hammer.pdf.
    - Uses matplotlib Agg backend lazily.
    """
    if log_dir is None:
        log_dir = _latest_log_dir(base_log_dir)
    if log_dir is None:
        return

    text_dir = os.path.join(log_dir, "text")
    fig_dir = os.path.join(log_dir, "fig")
    os.makedirs(fig_dir, exist_ok=True)

    os.environ.setdefault("MPLBACKEND", "Agg")
    # Lazy import matplotlib with Agg
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'cm'
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.linestyle'] = "--"
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.major.width'] = 1.0
    plt.rcParams['ytick.major.width'] = 1.0
    plt.rcParams['axes.linewidth'] = 1.2

    # Print selected directory
    print(f"[debug-plot] using log_dir: {log_dir}")

    # Collect event files
    files = [f for f in os.listdir(text_dir) if f.startswith("sim_") and f.endswith("_events.txt")]
    files.sort()

    for fn in files:
        path = os.path.join(text_dir, fn)
        m = re.match(r"sim_(\d{5})_events\.txt\Z", fn)
        tlabel = m.group(1) if m else "00000"
        out_path = os.path.join(fig_dir, f"sim_{tlabel}_hammer.{fmt}")
        if (not overwrite) and os.path.exists(out_path):
            continue

        ra_deg, dec_deg, inc = _read_events(path)
        if not ra_deg:
            continue

        # Try to read clique membership and reconcile
        clique_file = None
        # pick any k file for this t
        for name in os.listdir(text_dir):
            if name.startswith(f"sim_{tlabel}_cliques_k") and name.endswith(".txt"):
                clique_file = os.path.join(text_dir, name)
                break
        members = _read_clique_members(clique_file) if clique_file else set()
        inc_from_events = set(i for i, v in enumerate(inc) if v)
        if members and inc_from_events != members:
            print(f"[debug-plot] warn: in_cluster mismatch: events={len(inc_from_events)} vs cliques={len(members)}; using union for coloring")
        use_members = members if members else inc_from_events

        fig = plt.figure(figsize=(8, 4))
        ax = fig.add_subplot(111, projection="hammer")

        # Draw in two passes to ensure red (members) stays on top
        n_all = len(ra_deg)
        n_mem = len(use_members)
        n_black = 0
        for i in range(n_all):
            if i not in use_members:
                _draw_small_circle(ax, ra_deg[i], dec_deg[i], circle_radius_deg, color="black", lw=0.5, segments=180)
                n_black += 1
        n_red = 0
        for i in sorted(use_members):
            if 0 <= i < n_all:
                _draw_small_circle(ax, ra_deg[i], dec_deg[i], circle_radius_deg, color="red", lw=0.5, segments=180)
                n_red += 1

        # Grid and ticks
        ax.grid(True)
        ra_ticks_deg = np.array([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150])
        lam_ticks = np.deg2rad(-ra_ticks_deg)
        ax.set_xticks(lam_ticks)
        ax.set_xticklabels([str(int(x)) for x in ra_ticks_deg])
        yticks_deg = np.arange(-60, 61, 30)
        ax.set_yticks(np.deg2rad(yticks_deg))
        ax.set_yticklabels([str(int(x)) for x in yticks_deg])
        ax.set_xlabel("RA [deg]")
        ax.set_ylabel("Dec [deg]")

        fig.tight_layout()
        fig.savefig(out_path, format=fmt, bbox_inches="tight")
        print(f"[debug-plot] wrote {out_path} (radius={circle_radius_deg} deg, total={n_all}, red={n_red}, black={n_black})")
        plt.close(fig)

    return None
