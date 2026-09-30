"""scikit-fda interop: FDataGrid inputs equal ndarray inputs bit for bit,
and the penalized reconstruction matches the R fda::smooth.basis fixture."""

from __future__ import annotations

import numpy as np
import pytest

import fdachange as fc

from conftest import load_mat

skfda = pytest.importorskip("skfda")


def _as_fdatagrid(X):
    return skfda.FDataGrid(data_matrix=X.T,
                           grid_points=np.linspace(0, 1, X.shape[0]))


def test_fdatagrid_inputs_match_ndarray(F5):
    fd = _as_fdatagrid(F5)
    assert fc.far1_statistic(fd, 0.4) == fc.far1_statistic(F5, 0.4)
    np.testing.assert_array_equal(fc.project_l2(fd), fc.project_l2(F5))
    a = fc.project_ito(fd, rng=11)
    b = fc.project_ito(F5, rng=11)
    np.testing.assert_array_equal(a, b)
    r_fd = fc.reconstruct(fd, "fs_ls", K=4)["X"]
    r_np = fc.reconstruct(F5, "fs_ls", K=4)["X"]
    np.testing.assert_array_equal(r_fd, r_np)


def test_bspline_pen_matches_r_fixture(F5):
    # Same penalized LS problem as fda::smooth.basis (order-4 dyadic basis,
    # second-derivative L2 penalty, lambda = 1e-4); different solvers, so a
    # loose numerical tolerance.
    r = fc.reconstruct(F5, "bspline_pen", K=4, lam=1e-4)
    expected = load_mat("expected_bsplinepen.csv")
    np.testing.assert_allclose(r["X"], expected, rtol=0, atol=1e-6)
