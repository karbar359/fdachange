from __future__ import annotations

import pytest
from hamcrest import assert_that, close_to, greater_than, is_

import fdachange as fc


def describe_kuiper_quantile():
    def describe_given_an_upper_tail_probability():
        params = (0.9, 0.95, 0.99)

        @pytest.mark.parametrize("prob", params)
        def test_it_inverts_the_closed_form_cdf(prob):
            # Setup.

            # Test.
            quantile = fc.kuiper_quantile(prob)

            # Assert.
            assert_that(fc.kuiper_cdf(quantile), is_(close_to(prob, 1e-9)))

    def test_it_is_increasing_in_the_probability():
        # Setup.

        # Test.
        q_low = fc.kuiper_quantile(0.9)
        q_high = fc.kuiper_quantile(0.99)

        # Assert.
        assert_that(q_high, is_(greater_than(q_low)))
