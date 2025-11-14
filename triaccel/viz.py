import numpy as np
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


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


def _ensure_outdir(outdir: str | None) -> str:
    """Return a directory path to save figures. Create default if None."""
    if outdir is None:
        outdir = os.path.join(os.getcwd(), "triaccel_plots")
    os.makedirs(outdir, exist_ok=True)
    return outdir

def _save_fig(fig, *, outdir: str | None, filename: str):
    outdir = _ensure_outdir(outdir)
    if not filename.lower().endswith(".pdf"):
        filename = filename + ".pdf"
    path = os.path.join(outdir, filename)
    fig.tight_layout()
    print(f"Saving figure to {path}")
    fig.savefig(path, format="pdf", bbox_inches="tight")
    return path


def count_hist(rec, *, log=True, ax=None, title=None, outdir=None, filename="count_hist"):
    counts = rec["counts"]
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    bins = np.arange(counts.min(), counts.max()+2) - 0.5
    ax.hist(counts, bins=bins, edgecolor="black", log=log, alpha=0.7)
    ax.set_xlabel("Triplet count per trial")
    ax.set_ylabel("Frequency" + (" (log)" if log else ""))
    if title:
        ax.set_title(title)
    fig = ax.figure
    _save_fig(fig, outdir=outdir, filename=filename)
    return ax


def _require_keys(res, keys):
    for k in keys:
        if k not in res:
            raise ValueError(
                f"simulate(..., bins=..., return_histograms=True) で '{k}' を返す設定にしてください。")


def hist1d(res, what="events", *, ax_ra=None, ax_dec=None, title_prefix=None, outdir=None, filename=None):
    """1次元ヒスト（等間隔bin）。what: 'events' or 'triplets'

    2D ヒスト（Dec×RA）の周辺和を取り、RA/Dec ごとの 1D 分布を描画します。
    """
    key_2d = f"{what}_hist2d"
    _require_keys(res, [key_2d, "ra_edges", "dec_edges"])

    H2d = np.asarray(res[key_2d])
    Era = np.asarray(res["ra_edges"])
    Ede = np.asarray(res["dec_edges"])
    Hra = H2d.sum(axis=0)
    Hde = H2d.sum(axis=1)

    if ax_ra is None or ax_dec is None:
        fig, (ax_ra, ax_dec) = plt.subplots(1, 2, figsize=(10, 3))

    ax_ra.stairs(Hra, Era)
    ax_ra.set_xlim(float(Era.min()), float(Era.max()))
    ax_ra.set_xlabel("RA [deg]")
    ax_ra.set_ylabel("Counts")
    ax_ra.set_title(
        f"{title_prefix+' ' if title_prefix else ''}{what.capitalize()} RA (1D)")

    ax_dec.stairs(Hde, Ede)
    ax_dec.set_xlim(float(Ede.min()), float(Ede.max()))
    ax_dec.set_xlabel("Dec [deg]")
    ax_dec.set_ylabel("Counts")
    ax_dec.set_title(
        f"{title_prefix+' ' if title_prefix else ''}{what.capitalize()} Dec (1D)")

    # decide filename
    if filename is None:
        filename = f"{what}_hist1d"
    fig = ax_ra.figure
    _save_fig(fig, outdir=outdir, filename=filename)

    return ax_ra, ax_dec


def hammer(res, what="events", *, ax=None, title=None, outdir=None, filename=None,
           colorbar=True, cbar_kwargs=None, **pcolormesh_kwargs):
    """Hammer 図（等間隔binを pcolormesh で描画）
    - RAエッジは [-180, 180) を想定し、Hammer の経度 λ = -RA [rad] に変換
    - pcolormesh の単調増加要件に合わせ、必要ならエッジと列を反転するだけ（対応を厳密維持）
    """
    key = f"{what}_hist2d"
    _require_keys(res, [key, "ra_edges", "dec_edges"])

    H = np.asarray(res[key])
    ra_edges_deg = np.asarray(res["ra_edges"])
    dec_edges_deg = np.asarray(res["dec_edges"])

    if H.shape[1] + 1 != ra_edges_deg.size:
        raise ValueError(f"shape mismatch: H.shape[1]+1={H.shape[1]+1} vs ra_edges.size={ra_edges_deg.size}")
    if H.shape[0] + 1 != dec_edges_deg.size:
        raise ValueError(f"shape mismatch: H.shape[0]+1={H.shape[0]+1} vs dec_edges.size={dec_edges_deg.size}")

    lon_edges_deg = -ra_edges_deg
    d = np.diff(lon_edges_deg)
    if np.all(d > 0):
        pass
    elif np.all(d < 0):
        lon_edges_deg = lon_edges_deg[::-1]
        H = H[:, ::-1]
    else:
        raise ValueError("RA edge array is not monotonic; cannot map H to cells safely.")

    lat_edges = np.deg2rad(dec_edges_deg)
    lon_edges = np.deg2rad(lon_edges_deg)

    if ax is None:
        fig = plt.figure(figsize=(8, 4))
        ax = fig.add_subplot(111, projection="hammer")
    else:
        fig = ax.figure

    pkw = dict(shading="auto", antialiased=False, linewidth=0)
    pkw.update(pcolormesh_kwargs)
    im = ax.pcolormesh(lon_edges, lat_edges, H, **pkw)

    ax.grid(True)
    ax.set_title(title if title else f"{what.capitalize()} (Hammer)")

    ra_ticks_deg = np.arange(-150, 151, 30)
    lam_ticks = np.deg2rad(-ra_ticks_deg)
    yticks_deg = np.arange(-60,  61, 30)
    ax.set_xticks(lam_ticks)
    ax.set_xticklabels([rf"${d}^\circ$" for d in ra_ticks_deg])
    ax.set_yticks(np.deg2rad(yticks_deg))
    ax.set_yticklabels([rf"${d}^\circ$" for d in yticks_deg])
    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("Dec [deg]")

    if colorbar:
        if cbar_kwargs is None:
            cbar_kwargs = dict(orientation="vertical", pad=0.02, shrink=0.9)
        cbar = fig.colorbar(im, ax=ax, **cbar_kwargs)
        cbar.set_label("Counts")

    if filename is None:
        filename = f"{what}_hammer"
    _save_fig(fig, outdir=outdir, filename=filename)
    return ax


def _norm_ppf(p: float) -> float:
    import numpy as _np
    a = [-3.969683028665376e+01,  2.209460984245205e+02,
         -2.759285104469687e+02,  1.383577518672690e+02,
         -3.066479806614716e+01,  2.506628277459239e+00]
    b = [-5.447609879822406e+01,  1.615858368580409e+02,
         -1.556989798598866e+02,  6.680131188771972e+01,
         -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01,
         -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00,  2.938163982698783e+00]
    d = [7.784695709041462e-03,  3.224671290700398e-01,
         2.445134137142996e+00,  3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p <= 0.0:
        return float("-inf")
    if p >= 1.0:
        return float("inf")
    if p < plow:
        q = (-2.0 * _np.log(p))**0.5
        num = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5])
        den = ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0)
        return num / den
    if p > phigh:
        q = (-2.0 * _np.log(1.0 - p))**0.5
        num = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5])
        den = ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0)
        return -(num / den)
    q = p - 0.5
    r = q*q
    num = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q
    den = (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0)
    return num / den


def write_triplet_prob_sigma(
    res: dict,
    n: int = 10,
    *,
    outdir: str | None = None,
    filename: str = "counts.npy",
) -> str:
    counts = np.asarray(res.get("counts"))
    if counts is None:
        raise ValueError("res does not contain 'counts'")
    if n < 1:
        raise ValueError("n must be >= 1")
    n_min, n_max = 1, int(n)

    outdir = _ensure_outdir(outdir)
    counts_path = os.path.join(outdir, filename)
    try:
        if counts_path.lower().endswith(".npy"):
            np.save(counts_path, counts)
        else:
            np.savetxt(counts_path, counts, fmt="%d")
    except Exception as e:
        print(f"Failed to save counts to {counts_path}: {e}")

    lines = []
    header = "\nProbability P(triplets >= n) and one-sided Gaussian sigma (n=%d..%d):" % (n_min, n_max)
    lines.append(header)

    counts_int = counts.astype(np.int64, copy=False)
    total = counts_int.size
    maxv = int(counts_int.max()) if total > 0 else 0
    L = max(n_max, maxv) + 1
    # freq[v] = 値 v の出現回数
    freq = np.bincount(counts_int, minlength=L)
    # cum_ge[n] = 値が n 以上のサンプル数（後ろ向き累積）
    cum_ge = np.cumsum(freq[::-1])[::-1]

    for n in range(n_min, n_max + 1):
        ge = int(cum_ge[n]) if n < cum_ge.size else 0  # 安全策: 範囲外は 0
        p = ge / total if total > 0 else 0.0           # P(counts >= n)
        if p <= 0.0:
            z = float('inf')                           # 観測上ゼロ → 片側∞σ相当
        elif p >= 1.0:
            z = 0.0                                    # 常に達成 → 0σ
        else:
            z = float(_norm_ppf(1.0 - p))              # 片側: Z such that P(Z>=z)=p
        lines.append(f"n={n:2d}  p={p:<12.6g}  sigma={z:<10.6f}")
    lines.append("")

    for ln in lines:
        print(ln)

    base, ext = os.path.splitext(os.path.basename(filename))
    summary_name = f"{base}_prob_sigma.txt"
    summary_path = os.path.join(outdir, summary_name)
    try:
        with open(summary_path, "w", encoding="utf-8") as f:
            for ln in lines:
                f.write(ln + "\n")
    except Exception as e:
        print(f"Failed to write summary to {summary_path}: {e}")

    return summary_path
