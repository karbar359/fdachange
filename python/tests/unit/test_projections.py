from __future__ import annotations

import numpy as np
import pytest
from hamcrest import assert_that, equal_to, greater_than, is_, same_instance

import fdachange as fc


@pytest.fixture
def F5() -> np.ndarray:
    return fc.generate_f5(50, a=3, rng=9)


def describe_project_l2():
    def test_it_is_deterministic_and_nonnegative(F5):
        # Setup.
        exp_shape = (F5.shape[1],)

        # Test.
        a = fc.project_l2(F5)

        # Assert.
        assert_that(a.shape, is_(equal_to(exp_shape)))
        assert_that(bool(np.all(a >= 0)), is_(True))
        np.testing.assert_array_equal(a, fc.project_l2(F5))


def describe_project_ito():
    def test_it_redraws_the_direction_per_seed(F5):
        # Setup.

        # Test.
        a = fc.project_ito(F5, rng=2)
        b = fc.project_ito(F5, rng=2)
        c = fc.project_ito(F5, rng=3)

        # Assert.
        np.testing.assert_array_equal(a, b)   # same path, same projection
        assert_that(bool(np.array_equal(a, c)), is_(False))


def describe_the_basis_projections():
    def test_it_returns_one_value_per_curve(F5):
        # Setup.
        n = F5.shape[1]

        # Test.

        # Assert.
        assert_that(fc.project_c_max_var(F5, rng=4).shape, is_(equal_to((n,))))
        assert_that(fc.project_wiener(F5, rng=4).shape, is_(equal_to((n,))))
        assert_that(fc.project_pc1(F5).shape, is_(equal_to((n,))))
        assert_that(fc.project_pc_d(F5, d=2).shape, is_(equal_to((2, n))))


def describe_project_c_max_var():
    def describe_given_low_noise_curves_with_a_spike_in_half_the_sample():
        def test_it_evaluates_at_the_variance_peak_and_separates_the_groups():
            # Setup. The variance function peaks at t = 1/2; the selected
            # critical point lands within one dyadic vertex of it, and point
            # evaluation there separates the groups by ~a * f_alpha(point)
            # (f_alpha one vertex off the peak, at 0.4375, is ~0.47).
            rng = np.random.default_rng(5)
            L, n, a = 101, 40, 4.0
            s = np.arange(L) / (L - 1)
            spike = a * (1 - np.abs(s - 0.5) ** 0.45) ** (1 / 0.45)
            X = rng.normal(0.0, 0.05, (L, n))
            X[:, : n // 2] += spike[:, None]

            # Test.
            p = fc.project_c_max_var(X, rng=6)
            separation = p[: n // 2].mean() - p[n // 2:].mean()

            # Assert.
            assert_that(float(separation), is_(greater_than(a * 0.4)))


def describe_projection_fn():
    def test_it_resolves_names_to_the_projection_functions():
        # Setup.
        x = np.arange(5.0)

        # Test.

        # Assert.
        assert_that(fc.projection_fn("l2"), is_(same_instance(fc.project_l2)))
        np.testing.assert_array_equal(fc.projection_fn("none")(x), x)
        with pytest.raises(ValueError):
            fc.projection_fn("nope")


def describe_projection_multiplicity():
    def test_it_counts_the_series_a_projection_produces():
        # Setup.

        # Test.

        # Assert.
        assert_that(fc.projection_multiplicity("l2"), is_(equal_to(1)))
        assert_that(fc.projection_multiplicity("pc_d", 4), is_(equal_to(4)))
