from __future__ import annotations

import numpy as np
from hamcrest import assert_that, equal_to, is_

import fdachange as fc

# Small draw counts: these tests check contracts, not distributions.
N_DRAWS = 5


def describe_the_null_simulators():
    def test_it_returns_the_requested_number_of_finite_positive_draws():
        # Setup.
        params = (
            fc.simulate_wilcoxon_null(N_DRAWS, 0.1, grid=200, rng=1),
            fc.simulate_lp_null(N_DRAWS, 0.1, grid=50, rng=1),
            fc.simulate_rkhs_null(N_DRAWS, 0.1, n=100, rng=1),
            fc.simulate_far1_null(N_DRAWS, 0.1, grid=100, rng=1),
        )

        # Test.

        # Assert.
        for draws in params:
            assert_that(draws.shape, is_(equal_to((N_DRAWS,))))
            assert_that(bool(np.all(np.isfinite(draws))), is_(True))
            assert_that(bool(np.all(draws > 0)), is_(True))

    def describe_given_equal_seeds():
        def test_it_reproduces_the_draw_stream():
            # Setup.

            # Test.
            a = fc.simulate_wilcoxon_null(N_DRAWS, 0.25, grid=200, rng=3)
            b = fc.simulate_wilcoxon_null(N_DRAWS, 0.25, grid=200, rng=3)
            c = fc.simulate_wilcoxon_null(N_DRAWS, 0.25, grid=200, rng=4)

            # Assert.
            np.testing.assert_array_equal(a, b)
            assert_that(bool(np.array_equal(a, c)), is_(False))
