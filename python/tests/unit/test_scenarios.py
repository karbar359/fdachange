from __future__ import annotations

import numpy as np
import pytest
from hamcrest import assert_that, close_to, equal_to, greater_than, is_
from scipy.stats import norm

import fdachange as fc


def describe_the_scenario_generators():
    def test_it_returns_the_documented_shapes():
        # Setup.
        exp_scalar_shape = (50,)
        exp_curve_shape = (31, 20)            # l + 1 grid points, one curve per column

        # Test.
        x = fc.generate_s1(50, rng=1)

        # Assert.
        assert_that(x.shape, is_(equal_to(exp_scalar_shape)))
        for gen in (fc.generate_f1, fc.generate_f2, fc.generate_f3,
                    fc.generate_f4, fc.generate_f5):
            assert_that(gen(20, a=2, l=30, rng=1).shape,
                        is_(equal_to(exp_curve_shape)))

    def describe_the_rng_argument():
        def test_it_reproduces_streams_for_equal_seeds_only():
            # Setup.

            # Test.
            a = fc.generate_f1(10, a=3, rng=5)
            b = fc.generate_f1(10, a=3, rng=5)
            c = fc.generate_f1(10, a=3, rng=6)

            # Assert.
            np.testing.assert_array_equal(a, b)
            assert_that(bool(np.array_equal(a, c)), is_(False))


def describe_generate_s1():
    def describe_given_the_short_window_regime():
        params = (
            (True, 6),    # exact_window floors to l = 6 at n = 400
            (False, 7),   # the legacy rounding keeps a 7th point
        )

        @pytest.mark.parametrize("exact, exp_changed", params)
        def test_it_puts_the_documented_count_of_points_in_the_window(
                exact, exp_changed):
            # Setup.
            x = fc.generate_s1(400, a=50, ec_interval=(0.4925, 0.5075),
                               exact_window=exact, rng=1)

            # Test.
            changed = int(np.sum(x > 25))

            # Assert.
            assert_that(changed, is_(equal_to(exp_changed)))


def describe_generate_s2():
    def test_it_hits_the_noncentral_chi_square_mean_inside_the_window():
        # Setup.
        exp_mean = 1 + 2 ** 2

        # Test.
        x = fc.generate_s2(200_000, a=2, ec_interval=(0.1, 0.9), rng=7)

        # Assert.
        assert_that(float(np.mean(x[20_000:180_000])),
                    is_(close_to(exp_mean, 0.06)))


def describe_generate_f5():
    def test_it_adds_exactly_the_f_alpha_spike_on_the_window():
        # Setup. Same stream, different spike_alpha: outside the window the
        # curves are identical, inside they differ by the spike difference.
        s = np.arange(101) / 100

        def f_alpha(t, al):
            return (1 - np.abs(t - 0.5) ** al) ** (1 / al)

        exp_diff = 3 * (f_alpha(s, 0.45) - f_alpha(s, 0.2))
        win = slice(20, 35)                   # exact window at n = 50
        out = np.r_[0:20, 35:50]

        # Test.
        Xa = fc.generate_f5(50, a=3, spike_alpha=0.45, rng=9)
        Xb = fc.generate_f5(50, a=3, spike_alpha=0.2, rng=9)

        # Assert.
        np.testing.assert_allclose(
            Xa[:, win] - Xb[:, win],
            np.broadcast_to(exp_diff[:, None], Xa[:, win].shape),
            rtol=0, atol=1e-12)
        np.testing.assert_array_equal(Xa[:, out], Xb[:, out])


def describe_scenario_cdf():
    def test_it_returns_the_documented_transforms():
        # Setup.
        x = np.array([-1.0, 0.0, 1.5])
        pareto_x = np.array([1.5, 3.0])

        # Test.
        s1 = fc.scenario_cdf("S1")(x)
        s3 = fc.scenario_cdf("S3")(pareto_x)

        # Assert.
        np.testing.assert_allclose(s1, norm.cdf(x))
        np.testing.assert_allclose(s3, 1 - pareto_x ** -3.0)
        assert_that(fc.scenario_cdf("S3")(np.array([0.5]))[0],
                    is_(equal_to(0.0)))

    def describe_when_the_scenario_uses_the_ecdf():
        def test_it_maps_data_into_the_unit_interval():
            # Setup.
            y = fc.generate_s1(30, rng=2)

            # Test.
            e = fc.scenario_cdf("F5")(y)

            # Assert.
            assert_that(float(e.min()), is_(greater_than(0.0)))
            assert_that(float(e.max()), is_(equal_to(1.0)))

    def describe_when_the_scenario_id_is_unknown():
        def test_it_raises_a_value_error():
            # Setup.

            # Test.
            with pytest.raises(ValueError):
                fc.scenario_cdf("X9")

            # Assert.
