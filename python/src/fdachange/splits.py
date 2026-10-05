"""Split-sample aggregation (paper 1's T^(d) statistic)."""

from __future__ import annotations

import numpy as np

from ._interop import _as_matrix

__all__ = ["scan_over_splits", "make_run"]


def scan_over_splits(X, m: int, run) -> float:
    """Maximum of ``run(chunk)`` over ``m`` equal consecutive chunks.

    Calibrate with ``critical_value(draws, alpha, splits=m)``. ``run`` must
    encapsulate projection + statistic so stochastic projections re-draw
    per chunk.
    """
    if m == 1:
        return float(run(X))
    if hasattr(X, "data_matrix") or (np.ndim(X) == 2):
        Xm = _as_matrix(X)
        n = Xm.shape[1]
        if n % m:
            raise ValueError(f"n = {n} not divisible by m = {m}")
        N = n // m
        return float(max(run(Xm[:, i * N:(i + 1) * N]) for i in range(m)))
    x = np.asarray(X, dtype=float)
    n = x.size
    if n % m:
        raise ValueError(f"n = {n} not divisible by m = {m}")
    N = n // m
    return float(max(run(x[i * N:(i + 1) * N]) for i in range(m)))


def make_run(project, stat_fn, gamma, stat_kwargs=None, matrix_in=False):
    """Build ``run(chunk) -> float`` composing projection and statistic.

    Matrix-valued projections (one scalar series per row, e.g.
    ``project_pc_d``) are aggregated by the inner maximum; adjust the
    critical value by ``projection_multiplicity``. ``matrix_in=True``
    forwards the projection output wholesale (Birbilas's statistic).
    """
    stat_kwargs = dict(stat_kwargs or {})

    def run(chunk) -> float:
        proj = project(chunk)
        if matrix_in:
            return float(stat_fn(proj, gamma, **stat_kwargs))
        proj = np.asarray(proj, dtype=float)
        if proj.ndim == 2:
            return float(max(stat_fn(row, gamma, **stat_kwargs)
                             for row in proj))
        return float(stat_fn(proj, gamma, **stat_kwargs))

    return run
