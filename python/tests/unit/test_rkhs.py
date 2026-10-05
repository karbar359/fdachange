from __future__ import annotations

import numpy as np
import pytest
from hamcrest import assert_that, close_to, is_

import fdachange as fc

# Printed-formula brute force vs the vectorised implementation.
TOL = 1e-10


@pytest.fixture
def x() -> np.ndarray:
    return np.random.default_rng(42).normal(size=60)


def _brute_rkhs(x, gamma, kernel):
    n = x.size
    K = fc.gram(x, kernel)
    w = fc.holder_weights(n, gamma)
    best = 0.0
    for l in range(1, n // 2 + 1):
        for a in range(0, n - l + 1):
            I = np.arange(a, a + l)
            ln = l / n
            d2 = K[np.ix_(I, I)].sum() + ln ** 2 * K.sum() \
                - 2 * ln * K[I, :].sum()
            best = max(best, np.sqrt(max(d2, 0.0)) * w[l - 1])
    return best / np.sqrt(n)


def describe_rkhs_statistic():
    def describe_given_the_gaussian_kernel():
        params = (0.0, 0.25, 0.45)

        @pytest.mark.parametrize("gamma", params)
        def test_it_matches_the_printed_formula_brute_force(gamma, x):
            # Setup.
            exp_stat = _brute_rkhs(x, gamma, fc.kernel_gaussian)

            # Test.
            stat = fc.rkhs_statistic(x, gamma, kernel="gaussian")

            # Assert.
            assert_that(stat, is_(close_to(exp_stat, TOL)))
