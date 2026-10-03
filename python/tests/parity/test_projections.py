"""Deterministic projections vs R on the fixed F5 input."""

from __future__ import annotations

import numpy as np

import fdachange as fc

from parity.helpers import load_vec


def describe_project_l2():
    def test_it_matches_r_on_the_fixed_curves(F5):
        # Setup.
        exp_values = load_vec("expected_l2.csv")

        # Test.
        values = fc.project_l2(F5)

        # Assert.
        np.testing.assert_allclose(values, exp_values, rtol=1e-12)
