"""Kuiper's law: the exact gamma = 0 null of the Wilcoxon scan."""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

__all__ = ["kuiper_cdf", "kuiper_quantile"]

_K = np.arange(1, 1001, dtype=float)


def _cdf_scalar(v: float) -> float:
    if v <= 0.3:
        return 0.0
    kv2 = (_K * v) ** 2
    s = 1.0 - float(np.sum(2.0 * (4.0 * kv2 - 1.0) * np.exp(-2.0 * kv2)))
    return min(1.0, max(0.0, s))


def kuiper_cdf(v):
    """P(V <= v) for Kuiper's V (the Brownian-bridge range)."""
    arr = np.asarray(v, dtype=float)
    out = np.vectorize(_cdf_scalar, otypes=[float])(arr)
    return float(out) if arr.ndim == 0 else out


def kuiper_quantile(p):
    """Inverse of :func:`kuiper_cdf` by root finding on [0.3, 6]."""
    def q(pp: float) -> float:
        pp = min(max(pp, 1e-9), 1.0 - 1e-12)
        return brentq(lambda x: _cdf_scalar(x) - pp, 0.3, 6.0,
                      xtol=np.finfo(float).eps ** 0.5)
    arr = np.asarray(p, dtype=float)
    out = np.vectorize(q, otypes=[float])(arr)
    return float(out) if arr.ndim == 0 else out


_GRID_CACHE: dict[int, np.ndarray] = {}


def _kuiper_grid(n: int = 200_000) -> np.ndarray:
    """Deterministic quantile grid qK((i - 0.5)/N) standing in for MC draws.

    Same construction as the R package: dense tabulated CDF on
    ``linspace(0.3, 6, 20001)``, clamped and made monotone, inverted by
    linear interpolation; error ~1e-5, far below the 1e-3 the 3-decimal
    critical values need.
    """
    if n not in _GRID_CACHE:
        v = np.linspace(0.3, 6.0, 20001)
        vk = np.outer(v ** 2, _K[:100] ** 2)
        cdf = 1.0 - np.sum(2.0 * (4.0 * vk - 1.0) * np.exp(-2.0 * vk), axis=1)
        cdf = np.maximum.accumulate(np.clip(cdf, 0.0, 1.0))
        keep = np.concatenate(([True], np.diff(cdf) > 0))
        probs = (np.arange(1, n + 1) - 0.5) / n
        _GRID_CACHE[n] = np.interp(probs, cdf[keep], v[keep])
    return _GRID_CACHE[n]
