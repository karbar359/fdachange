"""Critical values vs R: same frozen draws and the same quantile
definition, so agreement is at bit-parity level."""

from __future__ import annotations

import pytest
from hamcrest import assert_that, is_

import fdachange as fc

from parity.helpers import close_rel

GAMMAS = (0.0, 0.25, 0.45)
# Same frozen draws + same type-7 quantile: only representation noise left.
RTOL = 1e-12


def describe_wilcoxon_critical_value():
    @pytest.mark.parametrize(
        "gamma", (0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.495))
    def test_it_matches_r_on_every_cached_cell(gamma, expected):
        # Setup.
        exp_cv = expected[f"cv_wilcoxon|g={gamma:g}|b=0"]

        # Test.
        cv = fc.wilcoxon_critical_value(gamma)

        # Assert.
        assert_that(cv, is_(close_rel(exp_cv, RTOL)))

    def test_it_matches_r_for_the_beta_and_splits_variants(expected):
        # Setup.
        exp_beta = expected["cv_wilcoxon|g=0.5|b=1"]
        exp_splits = expected["cv_wilcoxon_s2|g=0.25"]

        # Test.
        cv_beta = fc.wilcoxon_critical_value(0.5, beta=1)
        cv_splits = fc.wilcoxon_critical_value(0.25, splits=2)

        # Assert.
        assert_that(cv_beta, is_(close_rel(exp_beta, RTOL)))
        assert_that(cv_splits, is_(close_rel(exp_splits, RTOL)))


def describe_lp_critical_value():
    params = tuple((g, p) for g in GAMMAS for p in (2, 4, 10))

    @pytest.mark.parametrize("gamma, p", params)
    def test_it_matches_r_on_every_cached_cell(gamma, p, expected):
        # Setup.
        exp_cv = expected[f"cv_lp|g={gamma:g}|p={p}"]

        # Test.
        cv = fc.lp_critical_value(gamma, p)

        # Assert.
        assert_that(cv, is_(close_rel(exp_cv, RTOL)))


def describe_rkhs_critical_value():
    cells = (("gaussian", ""), ("min", "rnorm.sq"), ("min", "runif"),
             ("2max", "runif"), ("2max", "rbeta.1.9"))
    params = tuple((g, k, b) for g in GAMMAS for k, b in cells)

    @pytest.mark.parametrize("gamma, kernel, baseline", params)
    def test_it_matches_r_on_every_cached_cell(gamma, kernel, baseline,
                                               expected):
        # Setup.
        exp_cv = expected[f"cv_rkhs|g={gamma:g}|k={kernel}|b={baseline}"]

        # Test.
        cv = fc.rkhs_critical_value(gamma, kernel, baseline)

        # Assert.
        assert_that(cv, is_(close_rel(exp_cv, RTOL)))


def describe_far1_critical_value():
    @pytest.mark.parametrize("gamma", (0, 0.25, 0.4, 0.45))
    def test_it_matches_r_on_every_cached_cell(gamma, expected):
        # Setup.
        exp_cv = expected[f"cv_far1|g={gamma:g}"]

        # Test.
        cv = fc.far1_critical_value(gamma)

        # Assert.
        assert_that(cv, is_(close_rel(exp_cv, RTOL)))


def describe_the_kuiper_closed_form():
    def test_it_matches_r_in_both_directions(expected):
        # Setup.
        exp_quantile = expected["kuiper_q|p=0.95"]
        exp_cdf = expected["kuiper_cdf|v=1.7473"]

        # Test.
        quantile = fc.kuiper_quantile(0.95)
        cdf = fc.kuiper_cdf(1.7473)

        # Assert.
        assert_that(quantile, is_(close_rel(exp_quantile, 1e-9)))
        assert_that(cdf, is_(close_rel(exp_cdf, RTOL)))
