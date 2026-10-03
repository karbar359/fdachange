from __future__ import annotations

import numpy as np
import pytest
from hamcrest import assert_that, close_to, equal_to, is_

import fdachange as fc

# Brute-force agreement tolerance: identical float formulas, no
# cancellation, so anything above machine noise is a real defect.
TOL = 1e-12


def describe_holder_weight():
    def describe_given_a_nonzero_t():
        params = (
            (0.1, 0.25, 0.0),
            (0.25, 0.45, 0.0),
            (0.2, 0.5, 1.0),
        )

        @pytest.mark.parametrize("t, gamma, beta", params)
        def test_it_computes_the_printed_formula(t, gamma, beta):
            # Setup.
            exp_weight = 1.0 / (t ** gamma * np.log(0.3 / t) ** beta)

            # Test.
            weight = fc.holder_weight(t, gamma, beta)

            # Assert.
            assert_that(weight, is_(close_to(exp_weight, TOL)))

    def describe_when_t_is_zero():
        def test_it_returns_zero_by_convention():
            # Setup.

            # Test.
            weight = fc.holder_weight(0.0, 0.25)

            # Assert.
            assert_that(weight, is_(equal_to(0.0)))


def describe_holder_weights():
    def test_it_evaluates_the_weight_on_the_window_length_grid():
        # Setup.
        n, gamma = 40, 0.25
        l = np.arange(1, n, dtype=float)

        # Test.
        w = fc.holder_weights(n, gamma)
        w_symmetric = fc.holder_weights(n, gamma, symmetric=True)

        # Assert.
        np.testing.assert_array_equal(w[: n - 1], fc.holder_weight(l / n, gamma))
        np.testing.assert_array_equal(
            w_symmetric[: n - 1], fc.holder_weight(l / n * (1 - l / n), gamma))
        assert_that(w[n - 1], is_(equal_to(0.0)))


def describe_critical_value():
    def test_it_is_the_type_7_quantile_with_the_split_adjustment():
        # Setup. Any draws work: the definition is what is under test.
        draws = np.random.default_rng(42).normal(size=501)
        exp_cv = float(np.quantile(draws, 0.95, method="linear"))
        exp_cv_split = float(np.quantile(draws, 0.95 ** 0.5, method="linear"))

        # Test.
        cv = fc.critical_value(draws, alpha=0.05)
        cv_split = fc.critical_value(draws, alpha=0.05, splits=2)

        # Assert.
        assert_that(cv, is_(equal_to(exp_cv)))
        assert_that(cv_split, is_(equal_to(exp_cv_split)))
