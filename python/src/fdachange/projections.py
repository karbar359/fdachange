"""Functional-to-scalar projections, mirroring the R reference pipeline.

The basis-backed projections deliberately hand-roll the R pipeline's exact
semantics (least-squares fit onto order-2 B-splines on the dyadic J = 4
mesh - the Schauder basis - then variance-critical-point evaluation /
FPCA with trapezoid inner products) rather than delegating to scikit-fda:
faithfulness to the thesis pipeline is their purpose. Inputs may be
(grid, n) arrays or skfda FDataGrid objects.
"""

from __future__ import annotations

import numpy as np

from ._interop import _as_matrix, _as_rng

__all__ = [
    "project_l2", "project_ito", "project_wiener", "project_c_max_var",
    "project_pc1", "project_pc_d", "projection_fn",
    "projection_multiplicity",
]


def project_l2(X, **_ignored) -> np.ndarray:
    """Trapezoid L2 norm per curve (deterministic)."""
    X = _as_matrix(X)
    L = X.shape[0]
    w = np.ones(L)
    w[0] = w[-1] = 0.5
    return np.sqrt((w @ X ** 2) / (L - 1))


def project_ito(X, rng=None, **_ignored) -> np.ndarray:
    """Midpoint integral of every curve against one shared random Wiener
    path per call (fresh random direction each call)."""
    X = _as_matrix(X)
    rng = _as_rng(rng)
    L = X.shape[0]
    l = L - 1
    X_shift = np.vstack([np.zeros((1, X.shape[1])), X[:-1]])
    B = rng.normal(0.0, 1.0 / l, L)
    return B @ ((X + X_shift) / 2.0)


# ---- Schauder (order-2 B-spline) basis machinery ---------------------------

def _dyadic_grid(K: int = 4) -> np.ndarray:
    return np.arange(0.0, 1.0 + 0.5 ** K, 0.5 ** K)


def _hat_design(argvals: np.ndarray, knots: np.ndarray) -> np.ndarray:
    """Design matrix of the piecewise-linear hat basis on ``knots``:
    B[i, j] = hat_j(argvals[i]). Equivalent to order-2 B-splines."""
    t = np.asarray(argvals, dtype=float)
    B = np.zeros((t.size, knots.size))
    for j, kj in enumerate(knots):
        left = knots[j - 1] if j > 0 else kj
        right = knots[j + 1] if j < knots.size - 1 else kj
        up = (t >= left) & (t <= kj)
        if j > 0:
            B[up, j] = (t[up] - left) / (kj - left)
        else:
            B[t <= kj, j] = 1.0
        down = (t > kj) & (t <= right)
        if j < knots.size - 1:
            B[down, j] = (right - t[down]) / (right - kj)
    return B


def _fit_vertices(X: np.ndarray, K: int = 4,
                  argvals: np.ndarray | None = None):
    """Least-squares fit of each curve onto the level-K polygonal space;
    returns (vertex grid, fitted vertex values (2^K+1, n)). The hat basis
    evaluates to the identity at the vertices, so coefficients ARE the
    fitted vertex values (the Data2fd map)."""
    L = X.shape[0]
    t = np.linspace(0.0, 1.0, L) if argvals is None else np.asarray(argvals)
    knots = _dyadic_grid(K)
    B = _hat_design(t, knots)
    coefs, *_ = np.linalg.lstsq(B, X, rcond=None)
    return knots, coefs


def project_c_max_var(X, rng=None, K: int = 4, **_ignored) -> np.ndarray:
    """Evaluate every fitted curve at a randomly chosen critical point of
    the (fitted) variance function; falls back to a uniform point when the
    variance derivative never changes sign."""
    X = _as_matrix(X)
    rng = _as_rng(rng)
    grid, V = _fit_vertices(X, K)
    u = rng.uniform(0.0, 1.0)
    point = u
    var = V.var(axis=1, ddof=1)
    d = np.diff(var)
    idx = np.where(d[:-1] * d[1:] < 0)[0]
    if idx.size:
        point = grid[idx[int(np.ceil(u * idx.size)) - 1]]
    # piecewise-linear evaluation of the fitted curves at `point`
    return np.array([np.interp(point, grid, V[:, j])
                     for j in range(V.shape[1])])


def project_wiener(X, rng=None, K: int = 4, **_ignored) -> np.ndarray:
    """Trapezoid inner product of each fitted curve with one random
    Brownian path on the vertex grid (fresh path per call)."""
    X = _as_matrix(X)
    rng = _as_rng(rng)
    grid, V = _fit_vertices(X, K)
    Lg = grid.size
    W = np.concatenate(([0.0],
                        np.cumsum(np.sqrt(np.diff(grid))
                                  * rng.normal(0.0, 1.0, Lg - 1))))
    w = np.ones(Lg)
    w[0] = w[-1] = 0.5
    dx = grid[1] - grid[0]
    return ((V * W[:, None]) * w[:, None]).sum(axis=0) * dx


def _fpca(V: np.ndarray, grid: np.ndarray, ncomp: int):
    """FPCA of fitted curves via SVD under the trapezoid inner product."""
    Lg = grid.size
    w = np.ones(Lg)
    w[0] = w[-1] = 0.5
    dx = grid[1] - grid[0]
    sw = np.sqrt(w * dx)
    C = (V - V.mean(axis=1, keepdims=True)) * sw[:, None]
    U, s, _ = np.linalg.svd(C, full_matrices=False)
    comps = (U[:, :ncomp] / sw[:, None])      # orthonormal in L2(trapezoid)
    return comps


def project_pc1(X, K: int = 4, **_ignored) -> np.ndarray:
    """Inner product with the leading functional principal component."""
    return project_pc_d(X, d=1, K=K)[0]


def project_pc_d(X, d: int = 3, K: int = 4, **_ignored) -> np.ndarray:
    """Inner products with the leading ``d`` components; (d, n) matrix
    (one scalar series per row) for the outer-max aggregation of
    ``make_run``; adjust critical values by ``projection_multiplicity``."""
    X = _as_matrix(X)
    grid, V = _fit_vertices(X, K)
    comps = _fpca(V, grid, d)
    Lg = grid.size
    w = np.ones(Lg)
    w[0] = w[-1] = 0.5
    dx = grid[1] - grid[0]
    return (comps * (w * dx)[:, None]).T @ V


_PROJECTIONS = {
    "none": lambda X, **kw: X,
    "l2": project_l2,
    "ito": project_ito,
    "wiener": project_wiener,
    "c_max_var": project_c_max_var,
    "pc1": project_pc1,
    "pc_d": project_pc_d,
}


def projection_fn(name: str):
    """Resolve a projection id to its function."""
    try:
        return _PROJECTIONS[name]
    except KeyError:
        raise ValueError(f"unknown projection: {name!r}") from None


def projection_multiplicity(name: str, pc_d_n: int = 3) -> int:
    """Number of scalar series a projection emits (CV adjustment factor)."""
    return pc_d_n if name == "pc_d" else 1
