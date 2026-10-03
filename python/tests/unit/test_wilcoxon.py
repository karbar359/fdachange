from __future__ import annotations

import numpy as np
import pytest
from hamcrest import assert_that, close_to, is_

import fdachange as fc

# Printed-formula brute force vs the vectorised implementation: an
# R-independent proof on arbitrary data, so the sample is generated here.
TOL = 1e-10


@pytest.fixture
def x() -> np.ndarray:
    return np.random.default_rng(42).normal(size=60)


def _brute_wilcoxon(x, gamma, kernel="h_w"):
    xv = np.abs(x) if kernel == "h_0" else x
    n = xv.size
    H = np.sign(xv[:, None] - xv[None, :])
    sigma = np.sqrt(np.mean(H.mean(axis=0) ** 2))
    w = fc.holder_weights(n, gamma, 0.0, symmetric=True)
    best = 0.0
    for k in range(n):                       # window {k+1..m}, 1-based R
        for m in range(k + 1, n + 1):
            I = np.arange(k, m)
            comp = np.setdiff1d(np.arange(n), I)
            v = abs(H[np.ix_(I, comp)].sum()) * w[m - k - 1]
            best = max(best, v)
    return best / (n ** 1.5 * sigma)


def describe_wilcoxon_statistic():
    def describe_given_the_location_kernel_h_w():
        params = (0.0, 0.25, 0.45)

        @pytest.mark.parametrize("gamma", params)
        def test_it_matches_the_printed_formula_brute_force(gamma, x):
            # Setup.
            exp_stat = _brute_wilcoxon(x, gamma)

            # Test.
            stat = fc.wilcoxon_statistic(x, gamma)

            # Assert.
            assert_that(stat, is_(close_to(exp_stat, TOL)))

    def describe_given_the_scale_kernel_h_0():
        def test_it_matches_the_printed_formula_brute_force(x):
            # Setup.
            exp_stat = _brute_wilcoxon(x, 0.25, "h_0")

            # Test.
            stat = fc.wilcoxon_statistic(x, 0.25, kernel="h_0")

            # Assert.
            assert_that(stat, is_(close_to(exp_stat, TOL)))
