"""Monte-Carlo null simulators for cells the packaged caches do not cover.

Serial, with an explicit ``rng`` argument for reproducibility. For heavy
calibration runs prefer the R package's parallel simulators; store results
with :func:`fdachange.save_null_draws` either way.
"""

from __future__ import annotations

import numpy as np

from ._interop import _as_rng
from ._kernels import _resolve_kernel
from ._weights import holder_weight
from .rkhs import rkhs_statistic

__all__ = [
    "simulate_wilcoxon_null", "simulate_lp_null",
    "simulate_rkhs_null", "simulate_birbilas_null",
]


def simulate_wilcoxon_null(n_draws: int, gamma: float, beta: float = 0.0,
                           grid: int = 2000, rng=None) -> np.ndarray:
    """Weighted Brownian-bridge range draws (the Wilcoxon scan's null)."""
    rng = _as_rng(rng)
    N = int(grid)
    ls = np.arange(1, N + 1)
    w = holder_weight(ls / N * (1.0 - ls / N), gamma, beta)
    out = np.empty(n_draws)
    for i in range(n_draws):
        Wp = np.concatenate(([0.0],
                             np.cumsum(rng.normal(0.0, np.sqrt(1.0 / N), N))))
        B = Wp - np.arange(N + 1) / N * Wp[-1]
        vals = np.fromiter(
            (np.max(np.abs(B[l:] - B[:N - l + 1])) * w[l - 1]
             for l in range(1, N + 1)), dtype=float, count=N)
        out[i] = vals.max()
    return out


def simulate_lp_null(n_draws: int, gamma: float, p: float = 2,
                     grid: int = 200, rng=None) -> np.ndarray:
    """Weighted Lp functional of the centred Wiener sheet."""
    rng = _as_rng(rng)
    N = int(grid)
    w = holder_weight(np.arange(1, N) / N, gamma, 0.0)
    us = np.arange(1, N + 1) / N
    out = np.empty(n_draws)
    for i in range(n_draws):
        Z = rng.normal(0.0, 1.0, (N, N)) / N
        W = Z.cumsum(axis=0).cumsum(axis=1)
        end_row = W[N - 1, :]
        end_col = W[:, N - 1]
        end_end = W[N - 1, N - 1]
        best = 0.0
        for k in range(N - 1):
            for l in range(k + 1, N):
                delta = ((W[:, l] - W[:, k]) - us * end_row[l]
                         - (l - k) / N * end_col + us * end_row[k]
                         + (l - k) / N * us * end_end)
                val = (np.mean(np.abs(delta) ** p)) ** (1.0 / p) * w[l - k - 1]
                if val > best:
                    best = float(val)
        out[i] = best
    return out


def simulate_rkhs_null(n_draws: int, gamma: float, kernel="gaussian",
                       rdist=None, n: int = 2000, beta: float = 0.0,
                       rng=None) -> np.ndarray:
    """Draws of the RKHS statistic under an i.i.d. H0 sampler ``rdist``
    (default standard normal)."""
    rng = _as_rng(rng)
    kfn = _resolve_kernel(kernel)
    if rdist is None:
        rdist = lambda size: rng.normal(0.0, 1.0, size)
    return np.array([rkhs_statistic(rdist(n), gamma, kernel=kfn, beta=beta)
                     for _ in range(n_draws)])


def simulate_birbilas_null(n_draws: int, gamma: float, grid: int = 1200,
                           rng=None) -> np.ndarray:
    """Draws of the asymptotic limit E_{q,gamma} on a Brownian lattice."""
    rng = _as_rng(rng)
    N = int(grid)
    out = np.empty(n_draws)
    for i in range(n_draws):
        x1 = np.cumsum(rng.normal(0.0, 1.0, (N, N)), axis=1)  # rows: BM paths
        CS = np.zeros((N + 1, N))
        CS[1:] = np.cumsum(x1, axis=0)
        E = float(np.max(np.abs(x1)))          # single-row boundary terms
        for d in range(1, N):
            block = np.abs(CS[d + 1:N + 1] - CS[0:N - d])
            v = block.max() / d ** gamma
            if v > E:
                E = float(v)
        out[i] = E * N ** (-1.0 + gamma)
    return out
