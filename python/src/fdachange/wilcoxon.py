"""Wilcoxon-type (antisymmetric kernel) scan statistic.

Faithful port of the R reference implementation's corrected family: every
window I = {a..b}, 1 <= a <= b <= n, weight MULTIPLIED and indexed by the
true window length, sigma-hat normalisation per the printed theorem. The
window enumeration order matches R (b-major), so argmax tie-breaking - and
therefore the reported (k, m) - is identical.
"""

from __future__ import annotations

import numpy as np

from ._weights import holder_weights

__all__ = ["wilcoxon_statistic", "wilcoxon_scan"]


def _feature(x: np.ndarray, kernel: str) -> np.ndarray:
    if kernel == "h_w":
        return x
    if kernel == "h_0":
        return np.abs(x)
    raise ValueError(f"unknown wilcoxon kernel: {kernel!r}")


def _windows_b_major(n: int):
    """(a_vec, b_vec, l_vec) over all 1 <= a <= b <= n in R's ordering."""
    counts = np.arange(1, n + 1)
    b_vec = np.repeat(counts, counts)
    starts = np.cumsum(counts) - counts
    a_vec = np.arange(b_vec.size) - np.repeat(starts, counts) + 1
    return a_vec, b_vec, b_vec - a_vec + 1


def _core(x, kernel: str):
    x = _feature(np.asarray(x, dtype=float), kernel)
    n = x.size
    H = np.sign(x[:, None] - x[None, :])
    sigma = float(np.sqrt(np.mean(H.mean(axis=0) ** 2)))
    if not np.isfinite(sigma) or sigma == 0.0:
        return None
    P = np.zeros((n + 1, n + 1))
    P[1:, 1:] = H.cumsum(axis=0).cumsum(axis=1)

    a_vec, b_vec, l_vec = _windows_b_major(n)
    rows = P[b_vec, n] - P[a_vec - 1, n]
    block = (P[b_vec, b_vec] - P[b_vec, a_vec - 1]
             - P[a_vec - 1, b_vec] + P[a_vec - 1, a_vec - 1])
    return {
        "diff_abs": np.abs(rows - block), "l_vec": l_vec,
        "a_vec": a_vec, "b_vec": b_vec, "sigma": sigma, "n": n,
    }


def wilcoxon_statistic(x, gamma, beta: float = 0.0, kernel: str = "h_w"):
    """Weighted Wilcoxon scan statistic; ``gamma`` may be an array."""
    gam = np.atleast_1d(np.asarray(gamma, dtype=float))
    co = _core(x, kernel)
    if co is None:
        out = np.zeros(gam.size)
        return float(out[0]) if np.isscalar(gamma) or np.ndim(gamma) == 0 else out
    norm = co["n"] ** 1.5 * co["sigma"]
    li = co["l_vec"] - 1
    vals = np.array([
        np.max(co["diff_abs"] * holder_weights(co["n"], g, beta,
                                               symmetric=True)[li]) / norm
        for g in gam
    ])
    if np.isscalar(gamma) or np.ndim(gamma) == 0:
        return float(vals[0])
    return vals


def wilcoxon_scan(x, gamma: float, beta: float = 0.0, kernel: str = "h_w"):
    """Scan with localization: dict with stat, k, m, interval (R convention
    I(k, m) = {k+1..m}, so k = a - 1, m = b)."""
    co = _core(x, kernel)
    if co is None:
        return {"stat": 0.0, "k": None, "m": None, "interval": (None, None)}
    w = holder_weights(co["n"], float(gamma), beta, symmetric=True)
    vals = co["diff_abs"] * w[co["l_vec"] - 1]
    i = int(np.argmax(vals))
    n = co["n"]
    k = int(co["a_vec"][i] - 1)
    m = int(co["b_vec"][i])
    return {"stat": float(vals[i] / (n ** 1.5 * co["sigma"])),
            "k": k, "m": m, "interval": (k / n, m / n)}
