from __future__ import annotations

import numpy as np
import pytest
from hamcrest import assert_that, close_to, equal_to, is_

import fdachange as fc
from fdachange import nulls

# Identity-style checks on the same arithmetic path.
RTOL = 1e-15


def describe_the_shipped_null_caches():
    def test_it_serves_the_published_cells_with_their_draw_counts():
        # Setup.
        params = (
            (fc.wilcoxon_null_draws(0.25), 30_000),
            (fc.wilcoxon_null_draws(0.5, beta=1), 30_000),
            (fc.lp_null_draws(0.45, p=10), 300),
            (fc.rkhs_null_draws(0.45, "2max", "rbeta.1.9"), 1_000),
            (fc.far1_null_draws(0.4), 1_000),
        )

        # Test.

        # Assert.
        for draws, exp_size in params:
            assert_that(draws.size, is_(equal_to(exp_size)))

    def describe_when_gamma_is_zero_for_the_wilcoxon_scan():
        def test_it_serves_a_cell_just_below_the_bridge_range_law():
            # Setup.
            exact_q95 = 1.7473  # 95% point of the Brownian-bridge range

            # Test.
            d = fc.wilcoxon_null_draws(0.0)

            # Assert.
            assert_that(d.size, is_(equal_to(30_000)))
            assert_that(float(np.quantile(d, 0.95)),
                        is_(close_to(exact_q95 - 0.01, 0.01)))

    def describe_when_a_cell_is_not_cached():
        def test_it_names_the_simulator_in_the_error():
            # Setup.

            # Test.
            with pytest.raises(KeyError, match="simulate_wilcoxon_null"):
                fc.wilcoxon_null_draws(0.33)
            with pytest.raises(KeyError, match="simulate_lp_null"):
                fc.lp_null_draws(0.25, p=3)
            with pytest.raises(KeyError, match="simulate_rkhs_null"):
                fc.rkhs_null_draws(0.25, "min")

            # Assert.


def describe_the_user_cache():
    def test_it_round_trips_saved_draws_and_overrides_shipped_cells(
            tmp_path, monkeypatch):
        # Setup.
        monkeypatch.setenv("FDACHANGE_CACHE_DIR", str(tmp_path))
        nulls._invalidate()
        draws = np.linspace(0.0, 5.0, 500)

        try:
            # Test.
            key = fc.save_null_draws(draws, "wilcoxon", 0.33)
            fc.save_null_draws(np.arange(1, 101) / 100, "lp", 0.25, p=2)
            fc.save_null_draws(draws, "rkhs", 0.1, kernel="2max",
                               baseline="runif")

            # Assert.
            assert_that(fc.null_cache_dir(), is_(equal_to(tmp_path)))
            assert_that(key, is_(equal_to("g0.33_b0")))
            np.testing.assert_array_equal(fc.wilcoxon_null_draws(0.33), draws)
            # a user entry overrides a shipped cell
            assert_that(fc.lp_null_draws(0.25, p=2).size, is_(equal_to(100)))
            # rkhs keying round-trips kernel + baseline
            np.testing.assert_array_equal(
                fc.rkhs_null_draws(0.1, "2max", "runif"), draws)
        finally:
            nulls._invalidate()


def describe_the_critical_value_accessors():
    def test_it_applies_the_type_7_quantile_and_split_adjustment():
        # Setup.
        d = fc.lp_null_draws(0.25, p=2)
        exp_cv = float(np.quantile(d, 0.95, method="linear"))
        exp_cv_split = float(np.quantile(d, np.sqrt(0.95), method="linear"))

        # Test.
        cv = fc.lp_critical_value(0.25, p=2)
        cv_split = fc.lp_critical_value(0.25, p=2, splits=2)

        # Assert.
        assert_that(cv, is_(close_to(exp_cv, abs(exp_cv) * RTOL)))
        assert_that(cv_split,
                    is_(close_to(exp_cv_split, abs(exp_cv_split) * RTOL)))
