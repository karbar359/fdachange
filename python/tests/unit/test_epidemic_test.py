from __future__ import annotations

import numpy as np
import pytest
from hamcrest import (assert_that, close_to, contains_string, equal_to,
                      greater_than_or_equal_to, instance_of, is_, less_than,
                      less_than_or_equal_to, none)

import fdachange as fc

# Relative tolerance for identity-style checks (same arithmetic path).
RTOL = 1e-12


def _close_rel(exp, rtol=RTOL):
    return close_to(exp, abs(exp) * rtol)


def describe_epidemic_test():
    def describe_given_a_moderate_epidemic_signal():
        params = ("wilcoxon", "lp", "rkhs")

        @pytest.mark.parametrize("method", params)
        def test_it_returns_a_coherent_result(method):
            # Setup.
            x = fc.generate_s1(200, a=1, ec_interval=(0.4, 0.7), rng=3)

            # Test.
            r = fc.epidemic_test(x, method, gamma=0.25)

            # Assert.
            assert_that(r.p_value, is_(greater_than_or_equal_to(0.0)))
            assert_that(r.p_value, is_(less_than_or_equal_to(1.0)))
            assert_that(r.p_value < r.alpha,
                        is_(equal_to(r.statistic > r.critical_value)))
            assert_that(r.reject,
                        is_(equal_to(r.statistic > r.critical_value)))
            assert_that(r.k, is_(greater_than_or_equal_to(0)))
            assert_that(r.k, is_(less_than(r.m)))
            assert_that(r.m, is_(less_than_or_equal_to(200)))

    def describe_when_the_signal_is_strong():
        def test_it_rejects_and_localises_inside_the_window():
            # Setup.
            x = fc.generate_s1(200, a=3, ec_interval=(0.4, 0.7), rng=4)

            # Test.
            r = fc.epidemic_test(x, "wilcoxon", gamma=0.25)

            # Assert.
            assert_that(r.p_value, is_(less_than(0.01)))
            assert_that(r.m, is_(greater_than_or_equal_to(101)))
            assert_that(r.k, is_(less_than(140)))

    def describe_when_gamma_is_zero_for_the_wilcoxon_scan():
        def test_it_calibrates_from_the_shipped_cell():
            # Setup.
            x = fc.generate_s1(150, rng=3)
            draws = fc.wilcoxon_null_draws(0.0)

            # Test.
            r = fc.epidemic_test(x, "wilcoxon", gamma=0.0)

            # Assert.
            assert_that(r.n_null, is_(equal_to(30_000.0)))
            assert_that(r.p_value, is_(_close_rel(
                1.0 - float(np.mean(draws <= r.statistic)))))
            assert_that(r.critical_value,
                        is_(_close_rel(fc.critical_value(draws, 0.05, 1))))

    def describe_when_null_draws_are_injected():
        def test_it_calibrates_against_them_with_the_type_7_quantile():
            # Setup.
            x = fc.generate_s1(100, rng=5)
            draws = np.linspace(0.0, 10.0, 1001)
            exp_cv = float(np.quantile(draws, 0.95, method="linear"))

            # Test.
            r = fc.epidemic_test(x, "rkhs", gamma=0.25, null=draws)

            # Assert.
            assert_that(r.critical_value, is_(_close_rel(exp_cv)))
            assert_that(r.n_null, is_(equal_to(1001)))
            exp_p = 1.0 - float(np.mean(draws <= r.statistic))
            assert_that(r.p_value, is_(close_to(exp_p, 1e-12)))

    def describe_when_splits_aggregation_is_requested():
        def test_it_takes_the_chunk_maximum_and_adjusts_the_calibration():
            # Setup.
            x = fc.generate_s1(200, rng=6)
            exp_stat = max(fc.wilcoxon_statistic(x[:100], 0.25),
                           fc.wilcoxon_statistic(x[100:], 0.25))

            # Test.
            r = fc.epidemic_test(x, "wilcoxon", gamma=0.25, splits=2)

            # Assert.
            assert_that(r.statistic, is_(_close_rel(exp_stat)))
            assert_that(r.critical_value, is_(_close_rel(
                fc.wilcoxon_critical_value(0.25, splits=2))))
            assert_that(r.k, is_(none()))
            d = fc.wilcoxon_null_draws(0.25)
            exp_p = 1.0 - float(np.mean(d <= exp_stat)) ** 2
            assert_that(r.p_value, is_(_close_rel(exp_p)))

    def describe_given_a_functional_matrix_and_the_far1_method():
        def test_it_consumes_the_curves_directly():
            # Setup.
            exp_n = 60
            X = fc.generate_f4(exp_n, a=2, l=50, rng=7)

            # Test.
            r = fc.epidemic_test(X, "far1", gamma=0.4)

            # Assert.
            assert_that(r.statistic, is_(_close_rel(fc.far1_statistic(X, 0.4))))
            assert_that(r.k, is_(greater_than_or_equal_to(0)))
            assert_that(r.k, is_(less_than(r.m)))
            assert_that(r.m, is_(less_than_or_equal_to(exp_n)))

    def describe_when_a_custom_rkhs_kernel_has_no_cached_null():
        def test_it_demands_injected_null_draws():
            # Setup.
            x = fc.generate_s1(50, rng=8)
            my_kernel = lambda a, b: np.exp(-np.abs(a - b))  # noqa: E731

            # Test.
            with pytest.raises(ValueError, match="null"):
                fc.epidemic_test(x, "rkhs", gamma=0.25, kernel=my_kernel)
            r = fc.epidemic_test(x, "rkhs", gamma=0.25, kernel=my_kernel,
                                 null=np.linspace(0, 5, 200))

            # Assert.
            assert_that(r, is_(instance_of(fc.EpidemicTestResult)))

    def describe_the_string_report():
        def test_it_names_the_critical_value_and_the_detected_window():
            # Setup.
            x = fc.generate_s1(100, a=2, rng=9)

            # Test.
            out = str(fc.epidemic_test(x, "wilcoxon", gamma=0.25))

            # Assert.
            assert_that(out, contains_string("critical value"))
            assert_that(out, contains_string("detected window"))
