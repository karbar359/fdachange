"""Curve reconstruction vs R on the fixed F5 input."""

from __future__ import annotations

import os

import numpy as np
import pytest

import fdachange as fc

from parity.helpers import load_mat


def describe_reconstruct():
    def describe_given_the_faber_schauder_methods():
        params = (
            ("fs_interp", "expected_fsinterp.csv", 1e-10),
            ("fs_ls", "expected_fsls.csv", 1e-8),
        )

        @pytest.mark.parametrize("method, fixture, exp_rtol", params)
        def test_it_matches_r_within_the_solver_tolerance(method, fixture,
                                                          exp_rtol, F5):
            # Setup.
            exp_X = load_mat(fixture)

            # Test.
            r = fc.reconstruct(F5, method, K=4)

            # Assert.
            np.testing.assert_allclose(r["X"], exp_X, rtol=exp_rtol)

    def describe_given_the_penalized_bspline_method():
        def test_it_matches_the_r_fda_smooth_basis_fixture(F5):
            # Setup. Same penalized LS problem as fda::smooth.basis (order-4
            # dyadic basis, second-derivative L2 penalty, lambda = 1e-4);
            # different solvers, so a loose absolute tolerance.
            try:
                import skfda  # noqa: F401 - bspline_pen needs the [fda] extra
            except ImportError:
                if os.environ.get("CI"):
                    pytest.fail("scikit-fda is in the [test] extra and must "
                                "be installed on CI", pytrace=False)
                pytest.skip("scikit-fda not installed")
            exp_X = load_mat("expected_bsplinepen.csv")

            # Test.
            r = fc.reconstruct(F5, "bspline_pen", K=4, lam=1e-4)

            # Assert.
            np.testing.assert_allclose(r["X"], exp_X, rtol=0, atol=1e-6)
