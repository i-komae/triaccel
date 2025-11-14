from __future__ import annotations

from typing import Optional, Sequence, TypedDict, NotRequired

import numpy as np
import numpy.typing as npt


class Site(TypedDict):
    lat: float
    lon: float
    zmax: float
    n: int


class SimResult(TypedDict):
    # Always present
    counts: npt.NDArray[np.int32]

    # Present when returning histograms/edges
    ra_edges: NotRequired[npt.NDArray[np.float64]]
    dec_edges: NotRequired[npt.NDArray[np.float64]]

    # 2D histograms (shape: (dec_bins, ra_bins))
    events_hist2d: NotRequired[npt.NDArray[np.int64]]
    triplets_hist2d: NotRequired[npt.NDArray[np.int64]]


def simulate(
    M: int,
    sites: Sequence[Site],
    d_deg: float,
    seed: Optional[int] = ...,
    bins: Optional[int | Sequence[int]] = ...,
    return_histograms: bool = ...,
    threads: Optional[int] = ...,
    progress: bool = ...,
    debug: bool = ...,
    cluster_size: int = ...,
) -> SimResult: ...
