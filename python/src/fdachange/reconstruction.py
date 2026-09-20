"""Curve reconstruction (smoothing) as an explicit pipeline step.

Three maps into related spaces, preserving the R reference package's
distinction. The level-K Faber-Schauder truncation and order-2 B-splines
on the dyadic mesh 2^-K span the SAME space (continuous piecewise-linear
functions, dimension 2^K + 1) but are DIFFERENT maps into it:

- ``fs_interp``: the Faber-Schauder truncation proper - interpolation at
  the dyadic points (the map the thesis text describes);
- ``fs_ls``: least squares over the full sampling grid onto that space
  (the ``fda::Data2fd`` map, which produced the printed application
  tables; NOT the FS truncation - the two differ materially on rough
  curves);
- ``bspline_pen``: order-4 B-splines with an L2 second-derivative
  penalty - the only genuinely smooth arm, outside the polygonal family;
  provided via scikit-fda (the ``fdachange[fda]`` extra).
"""

from __future__ import annotations

import numpy as np

from ._interop import _as_matrix
from .projections import _dyadic_grid, _fit_vertices

__all__ = ["reconstruct"]

_METHODS = ("none", "fs_interp", "fs_ls", "bspline_pen")


def reconstruct(X, method: str = "none", K: int = 4, grid: str = "vertices",
                lam: float = 1e-4, argvals=None) -> dict:
    """Reconstruct a sample of curves (one per column).

    Returns ``{"X": matrix on the evaluation grid, "argvals": that grid,
    "method": ..., "K": ..., "grid": ...}``. ``grid="vertices"`` evaluates
    on the 2^K + 1 dyadic points (reproduces the printed application
    tables); ``"input"`` evaluates back on the sampling grid.
    """
    if method not in _METHODS:
        raise ValueError(f"method must be one of {_METHODS}")
    X = _as_matrix(X)
    L = X.shape[0]
    t_in = np.linspace(0.0, 1.0, L) if argvals is None else np.asarray(argvals)

    if method == "none":
        return {"X": X, "argvals": t_in, "method": method, "K": None,
                "grid": "input"}

    vertices = _dyadic_grid(K)
    out_grid = t_in if grid == "input" else vertices

    if method == "fs_interp":
        V = np.column_stack([np.interp(vertices, t_in, X[:, j])
                             for j in range(X.shape[1])])
        Xr = V if grid == "vertices" else np.column_stack(
            [np.interp(t_in, vertices, V[:, j]) for j in range(V.shape[1])])
    elif method == "fs_ls":
        _, V = _fit_vertices(X, K, argvals=t_in)
        Xr = V if grid == "vertices" else np.column_stack(
            [np.interp(t_in, vertices, V[:, j]) for j in range(V.shape[1])])
    else:  # bspline_pen
        try:
            import skfda
            from skfda.misc.operators import LinearDifferentialOperator
            from skfda.misc.regularization import L2Regularization
            from skfda.preprocessing.smoothing import BasisSmoother
            from skfda.representation.basis import BSplineBasis
        except ImportError as e:
            raise ImportError(
                "reconstruct(method='bspline_pen') needs scikit-fda; "
                "install the extra: pip install 'fdachange[fda]'") from e
        basis = BSplineBasis(domain_range=(0.0, 1.0), order=4,
                             knots=vertices)
        fd = skfda.FDataGrid(data_matrix=X.T, grid_points=t_in)
        smoother = BasisSmoother(
            basis=basis, smoothing_parameter=lam,
            regularization=L2Regularization(LinearDifferentialOperator(2)),
            return_basis=True)
        fdb = smoother.fit_transform(fd)
        Xr = fdb(out_grid)[:, :, 0].T

    return {"X": np.asarray(Xr, dtype=float), "argvals": out_grid,
            "method": method, "K": int(K), "grid": grid}
