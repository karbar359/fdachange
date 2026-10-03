from __future__ import annotations

import numpy as np
import pytest
from hamcrest import assert_that, equal_to, is_

import fdachange as fc


@pytest.fixture
def F5() -> np.ndarray:
    return fc.generate_f5(50, a=3, rng=9)


def describe_reconstruct():
    def describe_given_method_none():
        def test_it_passes_the_curves_through(F5):
            # Setup.

            # Test.
            r = fc.reconstruct(F5, "none")

            # Assert.
            np.testing.assert_array_equal(r["X"], F5)

    def describe_given_the_faber_schauder_truncation():
        def test_it_evaluates_on_the_dyadic_vertices(F5):
            # Setup.
            exp_shape = (17, F5.shape[1])     # 2^K + 1 vertices at K = 4

            # Test.
            r = fc.reconstruct(F5, "fs_interp", K=4)

            # Assert.
            assert_that(r["X"].shape, is_(equal_to(exp_shape)))

        def test_it_fixes_polygonal_input_on_the_dyadic_mesh():
            # Setup.
            grid = np.arange(0, 1.0001, 2 ** -4)
            V = np.random.default_rng(2).normal(size=(17, 8))

            # Test.
            out = fc.reconstruct(V, "fs_interp", K=4, argvals=grid)["X"]

            # Assert.
            np.testing.assert_allclose(out, V, atol=1e-12)

        def test_it_can_evaluate_back_on_the_input_grid(F5):
            # Setup.

            # Test.
            r = fc.reconstruct(F5, "fs_interp", K=4, grid="input")

            # Assert.
            assert_that(r["X"].shape, is_(equal_to(F5.shape)))
