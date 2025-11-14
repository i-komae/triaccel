from __future__ import annotations

from typing import Any, Mapping, Optional, Tuple, Literal

from ._core import SimResult

import matplotlib.axes as _mpl_axes
import matplotlib.figure as _mpl_figure

Axes = _mpl_axes.Axes
Figure = _mpl_figure.Figure


def _ensure_outdir(outdir: Optional[str]) -> str: ...


def _save_fig(fig: Figure, *, outdir: Optional[str], filename: str) -> str: ...


def count_hist(
    rec: Mapping[str, Any] | SimResult,
    *,
    log: bool = ...,
    ax: Optional[Axes] = ...,
    title: Optional[str] = ...,
    outdir: Optional[str] = ...,
    filename: str = ...,
) -> Axes: ...


def hist1d(
    res: Mapping[str, Any] | SimResult,
    what: Literal["events", "triplets"] = ...,
    *,
    ax_ra: Optional[Axes] = ...,
    ax_dec: Optional[Axes] = ...,
    title_prefix: Optional[str] = ...,
    outdir: Optional[str] = ...,
    filename: Optional[str] = ...,
) -> tuple[Axes, Axes]: ...


def hammer(
    res: Mapping[str, Any] | SimResult,
    what: Literal["events", "triplets"] = ...,
    *,
    ax: Optional[Axes] = ...,
    title: Optional[str] = ...,
    outdir: Optional[str] = ...,
    filename: Optional[str] = ...,
    colorbar: bool = ...,
    cbar_kwargs: Optional[Mapping[str, Any]] = ...,
    **pcolormesh_kwargs: Any,
) -> Axes: ...


def write_triplet_prob_sigma(
    res: Mapping[str, Any] | SimResult,
    n: int = ...,
    *,
    outdir: Optional[str] = ...,
    filename: str = ...,
) -> str: ...
