from __future__ import annotations

import numpy as np
import pytest
from hamcrest import assert_that, close_to, is_

import fdachange as fc

# Printed-formula brute force vs the vectorised implementation.
TOL = 1e-10


@pytest.fixture
def X() -> np.ndarray:
    return fc.generate_f4(40, a=0, ec_interval=(0.4, 0.6), rng=7)


def _brute_far1(X, gamma, beta_t):
    L, n = X.shape
    s = np.linspace(0, 1, L)
    W = (1 - beta_t(s))[:, None] * (X - X.mean(axis=1, keepdims=True))
    best = 0.0
    for k in range(n):
        for m in range(k + 1, n + 1):
            v = np.abs(W[:, k:m].sum(axis=1)).max() / (m - k) ** gamma
            best = max(best, v)
    return best * n ** (-0.5 + gamma)


def describe_far1_statistic():
    def describe_given_the_default_beta1_dependence():
        params = (0.25, 0.4)

        @pytest.mark.parametrize("gamma", params)
        def test_it_matches_the_printed_formula_brute_force(gamma, X):
            # Setup.
            exp_stat = _brute_far1(X, gamma, fc.far1_betas["beta1"])

            # Test.
            stat = fc.far1_statistic(X, gamma)

            # Assert.
            assert_that(stat, is_(close_to(exp_stat, TOL)))
