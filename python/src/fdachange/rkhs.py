"""RKHS (kernel mean embedding) scan statistic, ported from the R
reference implementation: windows of length l <= n/2, all starts, full
cross term; l-major enumeration matching R's argmax tie-breaking."""

from __future__ import annotations

import numpy as np

from ._kernels import gram
from ._weights import holder_weights

__all__ = ["rkhs_statistic", "rkhs_scan"]


def _windows_l_major(n: int):
    l_max = n // 2
    ls = np.arange(1, l_max + 1)
    starts_per_l = n - ls + 1
    l_rep = np.repeat(ls, starts_per_l)
    offsets = np.cumsum(starts_per_l) - starts_per_l
    a_vec = np.arange(l_rep.size) - np.repeat(offsets, starts_per_l) + 1
    b_vec = a_vec + l_rep - 1
    return a_vec, b_vec, l_rep


def _core(x, gamma: float, kernel, beta: float = 0.0):
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 2:
        return None
    K = gram(x, kernel)
    P = np.zeros((n + 1, n + 1))
    P[1:, 1:] = K.cumsum(axis=0).cumsum(axis=1)
    K_sum = P[n, n]
    w = holder_weights(n, gamma, beta)

    a_vec, b_vec, l_rep = _windows_l_major(n)
    rows = P[b_vec, n] - P[a_vec - 1, n]
    block = (P[b_vec, b_vec] - P[b_vec, a_vec - 1]
             - P[a_vec - 1, b_vec] + P[a_vec - 1, a_vec - 1])
    l_n = l_rep / n
    s = block + l_n ** 2 * K_sum - 2.0 * l_n * rows
    np.maximum(s, 0.0, out=s)                # clamp floating-point negatives
    return {"vals": np.sqrt(s) * w[l_rep - 1], "a_vec": a_vec,
            "b_vec": b_vec, "n": n}


def rkhs_statistic(x, gamma: float, kernel="gaussian",
                   beta: float = 0.0) -> float:
    """Weighted RKHS scan statistic; ``kernel`` is a name or a callable."""
    co = _core(x, gamma, kernel, beta)
    if co is None:
        return 0.0
    return float(co["vals"].max() / np.sqrt(co["n"]))


def rkhs_scan(x, gamma: float, kernel="gaussian", beta: float = 0.0):
    """RKHS scan with localization (R convention: window {k+1..m})."""
    co = _core(x, gamma, kernel, beta)
    if co is None:
        return {"stat": 0.0, "k": None, "m": None, "interval": (None, None)}
    i = int(np.argmax(co["vals"]))
    n = co["n"]
    k = int(co["a_vec"][i] - 1)
    m = int(co["b_vec"][i])
    return {"stat": float(co["vals"][i] / np.sqrt(n)),
            "k": k, "m": m, "interval": (k / n, m / n)}
