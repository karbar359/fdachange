from __future__ import annotations

import numpy as np
import pytest
from hamcrest import assert_that, equal_to, is_, less_than

import fdachange as fc

CLUSTERS = tuple(f"cl{i}" for i in range(1, 12))


@pytest.fixture(scope="module")
def telco():
    return fc.load_telco_far1()


def describe_load_telco_far1():
    def test_it_returns_the_eleven_documented_clusters(telco):
        # Setup.
        exp_shape = (50, 28)                  # months x days, as in R

        # Test.

        # Assert.
        assert_that(tuple(telco), is_(equal_to(CLUSTERS)))
        for cl, c in telco.items():
            assert_that(c.y.shape, is_(equal_to(exp_shape)))
            assert_that(c.y.dtype, is_(equal_to(np.float64)))
            assert_that(c.multiplier.shape, is_(equal_to((50,))))
            assert_that(c.jump.shape, is_(equal_to((50,))))
            assert_that(c.outliers.ndim, is_(equal_to(2)))
            assert_that(c.outliers.shape[1], is_(equal_to(2)))
            assert_that(bool(np.all(c.y >= 0)), is_(True))

    def describe_the_designed_change_windows():
        def test_it_marks_exactly_the_generator_script_windows(telco):
            # Setup. The design constants live in data-raw/telco_far1.R.

            # Test.

            # Assert.
            assert_that(telco["cl1"].change_months.size, is_(equal_to(0)))
            np.testing.assert_array_equal(telco["cl2"].change_months,
                                          np.arange(37, 47))
            np.testing.assert_array_equal(telco["cl3"].change_months,
                                          np.arange(37, 50))
            np.testing.assert_array_equal(telco["cl4"].change_months,
                                          np.arange(37, 51))
            for cl in CLUSTERS[4:]:
                assert_that(telco[cl].change_months.size, is_(equal_to(0)))

    def describe_the_jump_and_level_designs():
        def test_it_carries_the_generator_script_parameters(telco):
            # Setup.
            cl3 = telco["cl3"]
            jump_months = np.setdiff1d(cl3.change_months,
                                       cl3.outlier_months) - 1

            # Test.

            # Assert.
            assert_that(bool(np.all(telco["cl2"].jump[36:46] == 0.25)),
                        is_(True))
            assert_that(bool(np.all(cl3.jump[jump_months] == -0.25)),
                        is_(True))
            assert_that(telco["cl4"].jump[36], is_(equal_to(0.0)))
            assert_that(float(telco["cl4"].jump[49]), is_(less_than(-0.3)))
            assert_that(bool(np.all(telco["cl7"].multiplier == 1.5)),
                        is_(True))
            assert_that(bool(np.all(telco["cl8"].multiplier == 0.6)),
                        is_(True))

    def describe_the_outlier_designs():
        def test_it_places_the_documented_outlier_days(telco):
            # Setup.
            exp_cl9_outliers = round(
                0.05 * (50 - telco["cl9"].outlier_months.size) * 28)

            # Test.

            # Assert.
            assert_that(telco["cl3"].outliers.shape, is_(equal_to((3, 2))))
            assert_that(telco["cl9"].outliers.shape[0],
                        is_(equal_to(exp_cl9_outliers)))
