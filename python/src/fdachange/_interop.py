"""Input interop: accept plain arrays or scikit-fda FDataGrid objects.

The package convention (mirroring the R reference implementation) is a
matrix with the evaluation grid on the ROWS and one observation per
COLUMN. scikit-fda's ``FDataGrid`` stores ``data_matrix`` as
``(n_samples, n_points[, 1])`` - duck-typed here and transposed, so any
matrix-consuming entry point takes either representation.
"""

from __future__ import annotations

import numpy as np

__all__ = ["_as_matrix", "_as_rng"]


def _as_matrix(X) -> np.ndarray:
    """Coerce to the (grid, n_observations) float matrix convention."""
    dm = getattr(X, "data_matrix", None)
    if dm is not None and hasattr(X, "grid_points"):
        arr = np.asarray(dm, dtype=float)
        if arr.ndim == 3:
            if arr.shape[2] != 1:
                raise ValueError("only univariate FDataGrid objects are supported")
            arr = arr[:, :, 0]
        return np.ascontiguousarray(arr.T)
    arr = np.asarray(X, dtype=float)
    if arr.ndim != 2:
        raise ValueError("expected a (grid, n) matrix or an skfda FDataGrid")
    return arr


def _as_rng(rng) -> np.random.Generator:
    """None | int seed | Generator -> Generator."""
    if isinstance(rng, np.random.Generator):
        return rng
    return np.random.default_rng(rng)
