"""Birbilas's functional epidemic-change statistic E_{n,beta,gamma},
ported from the R reference implementation (the corrected printed
denominator (m - k)^gamma, all window lengths including 1)."""

from __future__ import annotations

import numpy as np

from ._interop import _as_matrix

__all__ = ["far1_statistic", "far1_scan", "far1_betas"]

far1_betas = {
    "beta1": lambda t: 0.7 * np.sqrt(t),
    "beta2": lambda t: -0.7 * np.sqrt(t),
    "beta3": lambda t: 0.7 * np.sqrt(1.0 - t),
    "beta4": lambda t: -0.7 * np.sqrt(1.0 - t),
}


def _prefix(X, beta_t):
    X = _as_matrix(X)
    L, n = X.shape
    s = np.linspace(0.0, 1.0, L)
    W = (1.0 - beta_t(s))[:, None] * (X - X.mean(axis=1, keepdims=True))
    CS = np.zeros((L, n + 1))
    CS[:, 1:] = np.cumsum(W, axis=1)
    return CS, L, n


def far1_statistic(X, gamma: float = 0.4,
                       beta_t=far1_betas["beta1"]) -> float:
    """E_{n,beta,gamma} on the test scale n^(-1/2+gamma) E.

    ``X``: (grid, n) matrix or skfda FDataGrid, one curve per observation.
    """
    CS, L, n = _prefix(X, beta_t)
    if n < 1:
        return 0.0
    E = 0.0
    for d in range(n):                        # window of d + 1 columns
        block = np.abs(CS[:, d + 1:n + 1] - CS[:, 0:n - d])
        v = block.max() / (d + 1) ** gamma
        if v > E:
            E = float(v)
    return E * n ** (-0.5 + gamma)


def far1_scan(X, gamma: float = 0.4, beta_t=far1_betas["beta1"]):
    """E statistic with the maximising window (R convention {k+1..m})."""
    CS, L, n = _prefix(X, beta_t)
    if n < 1:
        return {"stat": 0.0, "k": None, "m": None, "interval": (None, None)}
    E, b_start, b_end = 0.0, None, None
    for d in range(n):
        block = np.abs(CS[:, d + 1:n + 1] - CS[:, 0:n - d])
        i = int(np.argmax(block))             # row-major first max
        v = block.flat[i] / (d + 1) ** gamma
        if v > E:
            E = float(v)
            col = i % block.shape[1]          # start column (0-based)
            b_start, b_end = col + 1, col + 1 + d
    k, m = b_start - 1, b_end
    return {"stat": E * n ** (-0.5 + gamma), "k": k, "m": m,
            "interval": (k / n, m / n)}
