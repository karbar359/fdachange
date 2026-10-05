from __future__ import annotations

import os

import numpy as np
import pytest
from hamcrest import assert_that, equal_to, is_

import fdachange as fc

# scikit-fda ships in the [test] extra, so a CI runner without it is
# misconfigured and must go red rather than silently skip.
try:
    import skfda
except ImportError:  # pragma: no cover - environment-dependent
    if os.environ.get("CI"):
        pytest.fail("scikit-fda is in the [test] extra and must be "
                    "installed on CI", pytrace=False)
    pytest.skip("scikit-fda not installed", allow_module_level=True)


@pytest.fixture
def F5() -> np.ndarray:
    return fc.generate_f5(50, a=3, rng=9)


def _as_fdatagrid(X):
    return skfda.FDataGrid(data_matrix=X.T,
                           grid_points=np.linspace(0, 1, X.shape[0]))


def describe_fdatagrid_interop():
    def describe_given_the_same_curves_as_ndarray_and_fdatagrid():
        def test_it_computes_bit_identical_results(F5):
            # Setup.
            fd = _as_fdatagrid(F5)

            # Test.
            stat_fd = fc.far1_statistic(fd, 0.4)
            proj_fd = fc.project_l2(fd)
            ito_fd = fc.project_ito(fd, rng=11)
            rec_fd = fc.reconstruct(fd, "fs_ls", K=4)["X"]

            # Assert.
            assert_that(stat_fd, is_(equal_to(fc.far1_statistic(F5, 0.4))))
            np.testing.assert_array_equal(proj_fd, fc.project_l2(F5))
            np.testing.assert_array_equal(ito_fd, fc.project_ito(F5, rng=11))
            np.testing.assert_array_equal(rec_fd,
                                          fc.reconstruct(F5, "fs_ls", K=4)["X"])
