"""Lp-type scan statistic (centred empirical process), ported from the R
reference implementation. Loop over window starts with the window-length
axis vectorised via cumulative sums; identical arithmetic per window."""

from __future__ import annotations

import numpy as np

from ._weights import holder_weights

__all__ = ["lp_statistic", "lp_scan"]


def _ecdf_transform(x: np.ndarray) -> np.ndarray:
    """R's ecdf(x)(x): the fraction of the sample <= x_i."""
    return np.searchsorted(np.sort(x), x, side="right") / x.size


def _prepare(x, gamma: float, p: float, cdf, t_grid: int):
    x = np.asarray(x, dtype=float)
    n = x.size
    Fx = _ecdf_transform(x) if cdf is None else np.asarray(cdf(x), dtype=float)
    w = holder_weights(n, gamma)
    ts = np.linspace(0.0, 1.0, t_grid + 1)
    ind = (Fx[:, None] <= ts[None, :]).astype(float) - ts[None, :]
    tot = ind.sum(axis=0)
    return n, w, ind, tot


def lp_statistic(x, gamma: float, p: float = 2, cdf=None,
                 t_grid: int = 100) -> float:
    """Weighted Lp scan of the F_n-centred empirical process.

    ``cdf=None`` (default) uses the empirical CDF of the sample itself -
    the calibration-safe choice for real data; pass the true null CDF when
    known (e.g. ``scipy.stats.norm.cdf`` for a standard normal null).
    """
    n, w, ind, tot = _prepare(x, gamma, p, cdf, t_grid)
    stat = 0.0
    for k in range(n):                       # window start (0-based)
        cs = np.cumsum(ind[k:], axis=0)      # rows: l = 1..n-k
        ls = np.arange(1, n - k + 1)
        cc = cs - (ls[:, None] / n) * tot[None, :]
        sums = (np.abs(cc) ** p).sum(axis=1) / t_grid
        vals = sums ** (1.0 / p) * w[ls - 1]
        m = float(vals.max())
        if m > stat:
            stat = m
    return stat / np.sqrt(n)


def lp_scan(x, gamma: float, p: float = 2, cdf=None, t_grid: int = 100):
    """Lp scan with localization (R convention: window {k+1..m})."""
    n, w, ind, tot = _prepare(x, gamma, p, cdf, t_grid)
    stat, bk, bm = 0.0, None, None
    for k in range(n):
        cs = np.cumsum(ind[k:], axis=0)
        ls = np.arange(1, n - k + 1)
        cc = cs - (ls[:, None] / n) * tot[None, :]
        sums = (np.abs(cc) ** p).sum(axis=1) / t_grid
        vals = sums ** (1.0 / p) * w[ls - 1]
        i = int(np.argmax(vals))             # first max within the block,
        if vals[i] > stat:                   # strict >: keeps the earliest
            stat = float(vals[i])            # (k, l) exactly like R
            bk, bm = k, k + i + 1
    return {"stat": stat / np.sqrt(n), "k": bk, "m": bm,
            "interval": (None, None) if bk is None else (bk / n, bm / n)}
