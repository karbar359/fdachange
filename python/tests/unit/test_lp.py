from __future__ import annotations

import numpy as np
import pytest
from hamcrest import assert_that, close_to, is_
from scipy.stats import norm

import fdachange as fc

# Printed-formula brute force vs the vectorised implementation.
TOL = 1e-10


@pytest.fixture
def x() -> np.ndarray:
    return np.random.default_rng(42).normal(size=60)


def _brute_lp(x, gamma, p=2, t_grid=100):
    n = x.size
    Fx = norm.cdf(x)
    w = fc.holder_weights(n, gamma)
    ts = np.linspace(0, 1, t_grid + 1)
    ind = (Fx[:, None] <= ts[None, :]).astype(float) - ts[None, :]
    tot = ind.sum(axis=0)
    best = 0.0
    for k in range(n):
        for l in range(1, n - k + 1):
            cs = ind[k:k + l].sum(axis=0) - (l / n) * tot
            v = (np.sum(np.abs(cs) ** p) / t_grid) ** (1 / p) * w[l - 1]
            best = max(best, v)
    return best / np.sqrt(n)


def describe_lp_statistic():
    def describe_given_the_normal_cdf_transform():
        params = (0.0, 0.25, 0.45)

        @pytest.mark.parametrize("gamma", params)
        def test_it_matches_the_printed_formula_brute_force(gamma, x):
            # Setup.
            exp_stat = _brute_lp(x, gamma)

            # Test.
            stat = fc.lp_statistic(x, gamma, cdf=norm.cdf)

            # Assert.
            assert_that(stat, is_(close_to(exp_stat, TOL)))
