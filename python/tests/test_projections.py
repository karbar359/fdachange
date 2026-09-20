from __future__ import annotations

import numpy as np
import pytest

import fdachange as fc


def test_l2_deterministic(F5):
    a = fc.project_l2(F5)
    assert a.shape == (F5.shape[1],)
    assert np.all(a >= 0)
    np.testing.assert_array_equal(a, fc.project_l2(F5))


def test_ito_stochastic(F5):
    a = fc.project_ito(F5, rng=2)
    b = fc.project_ito(F5, rng=2)
    np.testing.assert_array_equal(a, b)      # same path, same projection
    c = fc.project_ito(F5, rng=3)
    assert not np.array_equal(a, c)          # fresh random direction


def test_basis_projections_shapes(F5):
    n = F5.shape[1]
    assert fc.project_c_max_var(F5, rng=4).shape == (n,)
    assert fc.project_wiener(F5, rng=4).shape == (n,)
    assert fc.project_pc1(F5).shape == (n,)
    assert fc.project_pc_d(F5, d=2).shape == (2, n)


def test_c_max_var_selects_the_variance_peak():
    # Deterministic mechanism check: low-noise curves where half the sample
    # carries the F5 spike (a = 4). The variance function peaks at t = 1/2,
    # the selected critical point lands within one dyadic vertex of it, and
    # point evaluation there separates the groups by ~a * f_alpha(point).
    rng = np.random.default_rng(5)
    L, n, a = 101, 40, 4.0
    s = np.arange(L) / (L - 1)
    spike = a * (1 - np.abs(s - 0.5) ** 0.45) ** (1 / 0.45)
    X = rng.normal(0.0, 0.05, (L, n))
    X[:, : n // 2] += spike[:, None]
    p = fc.project_c_max_var(X, rng=6)
    sep = p[: n // 2].mean() - p[n // 2:].mean()
    # f_alpha at the vertex one step off the peak (0.4375) is ~0.47
    assert sep > a * 0.4


def test_projection_fn_registry():
    assert fc.projection_fn("l2") is fc.project_l2
    x = np.arange(5.0)
    np.testing.assert_array_equal(fc.projection_fn("none")(x), x)
    with pytest.raises(ValueError):
        fc.projection_fn("nope")
    assert fc.projection_multiplicity("l2") == 1
    assert fc.projection_multiplicity("pc_d", 4) == 4


def test_make_run_and_splits(F5):
    run = fc.make_run(fc.project_l2, fc.rkhs_statistic, 0.25,
                      stat_kwargs={"kernel": "gaussian"})
    assert fc.scan_over_splits(F5, 1, run) == pytest.approx(run(F5))
    v2 = fc.scan_over_splits(F5, 2, run)
    assert v2 == pytest.approx(max(run(F5[:, :25]), run(F5[:, 25:])))
    with pytest.raises(ValueError, match="divisible"):
        fc.scan_over_splits(F5, 7, run)
    run_b = fc.make_run(lambda X: X, fc.birbilas_statistic, 0.4,
                        matrix_in=True)
    assert run_b(F5) == pytest.approx(fc.birbilas_statistic(F5, 0.4))


def test_reconstruct_contracts(F5):
    r = fc.reconstruct(F5, "none")
    np.testing.assert_array_equal(r["X"], F5)
    r = fc.reconstruct(F5, "fs_interp", K=4)
    assert r["X"].shape == (17, F5.shape[1])
    # polygonal input on the dyadic mesh is a fixed point
    grid = np.arange(0, 1.0001, 2 ** -4)
    V = np.random.default_rng(2).normal(size=(17, 8))
    out = fc.reconstruct(V, "fs_interp", K=4, argvals=grid)["X"]
    np.testing.assert_allclose(out, V, atol=1e-12)
    ri = fc.reconstruct(F5, "fs_interp", K=4, grid="input")
    assert ri["X"].shape == F5.shape
