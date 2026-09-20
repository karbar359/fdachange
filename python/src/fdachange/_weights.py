"""Hoelder-type scan weights and Monte-Carlo critical values.

Mirrors the fdachange R package (the reference implementation): the scan
statistics MULTIPLY window increments by ``w(t) = t^-gamma * log(0.3/t)^-beta``,
with ``w(0) = 0`` by convention.
"""

from __future__ import annotations

import numpy as np

__all__ = ["holder_weight", "holder_weights", "critical_value"]


def holder_weight(t, gamma: float, beta: float = 0.0):
    """Hoelder-type weight ``w(t) = t^-gamma * log(0.3/t)^-beta`` (``w(0)=0``)."""
    t = np.asarray(t, dtype=float)
    out = np.zeros_like(t)
    nz = t != 0
    denom = t[nz] ** gamma * np.log(0.3 / t[nz]) ** beta
    out[nz] = 1.0 / denom
    if out.ndim == 0:
        return float(out)
    return out


def holder_weights(n: int, gamma: float, beta: float = 0.0,
                   symmetric: bool = False) -> np.ndarray:
    """Weight per window length ``l = 1..n-1``; entry ``n`` is 0.

    ``w[l-1] = w(l/n)`` (RKHS/Lp), or ``w((l/n)(1-l/n))`` when
    ``symmetric`` (Wilcoxon). Returned array has length ``n`` and is
    0-indexed by ``l - 1`` (the R vector is 1-indexed by ``l``).
    """
    l = np.arange(1, n, dtype=float)
    t = l / n
    if symmetric:
        t = t * (1.0 - t)
    w = np.zeros(n)
    w[: n - 1] = holder_weight(t, gamma, beta)
    return w


def critical_value(draws, alpha: float = 0.05, splits: int = 1) -> float:
    """Upper quantile of null draws with the split-aggregation adjustment.

    Linear-interpolation quantile (numpy ``method="linear"``), which is R's
    type 7 - published thresholds depend on it, so it is not configurable.
    """
    draws = np.asarray(draws, dtype=float)
    return float(np.quantile(draws, (1.0 - alpha) ** (1.0 / splits),
                             method="linear"))
