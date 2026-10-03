from __future__ import annotations

import numpy as np
import pytest
from hamcrest import assert_that, close_to, is_

import fdachange as fc

# Identity-style checks on the same arithmetic path.
RTOL = 1e-12


def _close_rel(exp):
    return close_to(exp, abs(exp) * RTOL)


@pytest.fixture
def F5() -> np.ndarray:
    return fc.generate_f5(50, a=3, rng=9)


def describe_scan_over_splits():
    def describe_given_a_projection_statistic_run():
        def test_it_reduces_to_the_plain_run_at_one_split(F5):
            # Setup.
            run = fc.make_run(fc.project_l2, fc.rkhs_statistic, 0.25,
                              stat_kwargs={"kernel": "gaussian"})

            # Test.
            value = fc.scan_over_splits(F5, 1, run)

            # Assert.
            assert_that(value, is_(_close_rel(run(F5))))

        def test_it_takes_the_outer_maximum_over_the_chunks(F5):
            # Setup.
            run = fc.make_run(fc.project_l2, fc.rkhs_statistic, 0.25,
                              stat_kwargs={"kernel": "gaussian"})
            exp_value = max(run(F5[:, :25]), run(F5[:, 25:]))

            # Test.
            value = fc.scan_over_splits(F5, 2, run)

            # Assert.
            assert_that(value, is_(_close_rel(exp_value)))

    def describe_when_the_sample_size_is_not_divisible():
        def test_it_refuses_the_split(F5):
            # Setup.
            run = fc.make_run(fc.project_l2, fc.rkhs_statistic, 0.25,
                              stat_kwargs={"kernel": "gaussian"})

            # Test.
            with pytest.raises(ValueError, match="divisible"):
                fc.scan_over_splits(F5, 7, run)

            # Assert.


def describe_make_run():
    def describe_given_a_matrix_consuming_statistic():
        def test_it_skips_the_projection(F5):
            # Setup.
            run = fc.make_run(lambda X: X, fc.far1_statistic, 0.4,
                              matrix_in=True)

            # Test.
            value = run(F5)

            # Assert.
            assert_that(value, is_(_close_rel(fc.far1_statistic(F5, 0.4))))
