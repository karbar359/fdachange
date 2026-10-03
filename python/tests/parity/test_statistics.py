"""Statistics and scan windows vs the R reference implementation on
identical inputs."""

from __future__ import annotations

import numpy as np
import pytest
from hamcrest import assert_that, equal_to, is_
from scipy.stats import norm

import fdachange as fc

from parity.helpers import close_rel

GAMMAS = (0.0, 0.25, 0.45)
# Same formulas, different vectorisation order: float-noise agreement.
RTOL = 1e-10


def describe_wilcoxon_statistic():
    @pytest.mark.parametrize("gamma", GAMMAS)
    def test_it_matches_r_per_gamma(gamma, expected, x60):
        # Setup.
        exp_stat = expected[f"wilcoxon|g={gamma:g}"]

        # Test.
        stat = fc.wilcoxon_statistic(x60, gamma)

        # Assert.
        assert_that(stat, is_(close_rel(exp_stat, RTOL)))

    def test_it_matches_r_for_the_h_0_kernel(expected, x60):
        # Setup.
        exp_stat = expected["wilcoxon_h0|g=0.25"]

        # Test.
        stat = fc.wilcoxon_statistic(x60, 0.25, kernel="h_0")

        # Assert.
        assert_that(stat, is_(close_rel(exp_stat, RTOL)))

    def describe_given_a_vector_of_gammas():
        def test_it_equals_the_per_gamma_calls(expected, x60):
            # Setup.

            # Test.
            v = fc.wilcoxon_statistic(x60, np.array(GAMMAS))

            # Assert.
            for gamma, value in zip(GAMMAS, v):
                assert_that(value, is_(close_rel(
                    expected[f"wilcoxon|g={gamma:g}"], RTOL)))


def describe_lp_statistic():
    @pytest.mark.parametrize("gamma", GAMMAS)
    def test_it_matches_r_for_both_cdf_paths(gamma, expected, x60):
        # Setup.
        exp_pnorm = expected[f"lp_pnorm|g={gamma:g}"]
        exp_ecdf = expected[f"lp_ecdf|g={gamma:g}"]

        # Test.
        stat_pnorm = fc.lp_statistic(x60, gamma, p=2, cdf=norm.cdf)
        stat_ecdf = fc.lp_statistic(x60, gamma, p=2)

        # Assert.
        assert_that(stat_pnorm, is_(close_rel(exp_pnorm, RTOL)))
        assert_that(stat_ecdf, is_(close_rel(exp_ecdf, RTOL)))

    def test_it_matches_r_at_p_4(expected, x60):
        # Setup.
        exp_stat = expected["lp_p4|g=0.25"]

        # Test.
        stat = fc.lp_statistic(x60, 0.25, p=4, cdf=norm.cdf)

        # Assert.
        assert_that(stat, is_(close_rel(exp_stat, RTOL)))


def describe_rkhs_statistic():
    @pytest.mark.parametrize("gamma", GAMMAS)
    def test_it_matches_r_for_the_gaussian_kernel(gamma, expected, x60):
        # Setup.
        exp_stat = expected[f"rkhs_gaussian|g={gamma:g}"]

        # Test.
        stat = fc.rkhs_statistic(x60, gamma, kernel="gaussian")

        # Assert.
        assert_that(stat, is_(close_rel(exp_stat, RTOL)))

    def test_it_matches_r_for_the_min_and_2max_kernels(expected, x60):
        # Setup.
        exp_min = expected["rkhs_min|g=0.25"]
        exp_2max = expected["rkhs_2max|g=0.25"]

        # Test.
        stat_min = fc.rkhs_statistic(x60, 0.25, kernel="min")
        stat_2max = fc.rkhs_statistic(x60, 0.25, kernel="2max")

        # Assert.
        assert_that(stat_min, is_(close_rel(exp_min, RTOL)))
        assert_that(stat_2max, is_(close_rel(exp_2max, RTOL)))


def describe_far1_statistic():
    @pytest.mark.parametrize("gamma", (0.25, 0.4))
    def test_it_matches_r_per_gamma(gamma, expected, F4):
        # Setup.
        exp_stat = expected[f"far1|g={gamma:g}"]

        # Test.
        stat = fc.far1_statistic(F4, gamma)

        # Assert.
        assert_that(stat, is_(close_rel(exp_stat, RTOL)))

    def test_it_matches_r_for_the_beta3_dependence(expected, F4):
        # Setup.
        exp_stat = expected["far1_beta3|g=0.4"]

        # Test.
        stat = fc.far1_statistic(F4, 0.4, beta_t=fc.far1_betas["beta3"])

        # Assert.
        assert_that(stat, is_(close_rel(exp_stat, RTOL)))


def describe_the_scan_windows():
    def test_it_reports_r_identical_windows_and_statistics(expected, x60, F4):
        # Setup. Window indices must be EXACT: the implementations replicate
        # R's window enumeration order so argmax ties break identically.
        scans = (
            ("scan_wilcoxon", fc.wilcoxon_scan(x60, 0.25)),
            ("scan_lp", fc.lp_scan(x60, 0.25, cdf=norm.cdf)),
            ("scan_rkhs", fc.rkhs_scan(x60, 0.25, kernel="gaussian")),
            ("scan_far1", fc.far1_scan(F4, 0.4)),
        )

        # Test.

        # Assert.
        for key, sc in scans:
            assert_that(sc["stat"],
                        is_(close_rel(expected[f"{key}|stat"], RTOL)))
            assert_that(sc["k"], is_(equal_to(int(expected[f"{key}|k"]))))
            assert_that(sc["m"], is_(equal_to(int(expected[f"{key}|m"]))))
