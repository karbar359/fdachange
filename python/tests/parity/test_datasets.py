"""telco_far1 conversion parity vs the R dataset, plus the end-to-end
scan on it."""

from __future__ import annotations

import pytest
from hamcrest import assert_that, equal_to, is_

import fdachange as fc

from parity.helpers import close_rel

CLUSTERS = tuple(f"cl{i}" for i in range(1, 12))


@pytest.fixture(scope="module")
def telco():
    return fc.load_telco_far1()


def describe_load_telco_far1():
    def test_it_matches_the_r_values_per_cluster(telco, expected):
        # Setup. Sums at 1e-12: R and numpy accumulate in different orders.

        # Test.

        # Assert.
        for cl, c in telco.items():
            assert_that(float(c.y.sum()),
                        is_(close_rel(expected[f"telco|{cl}|sum"], 1e-12)))

    def test_it_matches_r_spot_cells_exactly(telco, expected):
        # Setup. Full-precision conversion: spot cells are bit-identical.
        params = (
            ("cl2", 0, 0, "telco|cl2|y[1,1]"),
            ("cl2", 36, 27, "telco|cl2|y[37,28]"),
            ("cl9", 49, 13, "telco|cl9|y[50,14]"),
        )

        # Test.

        # Assert.
        for cl, i, j, key in params:
            assert_that(telco[cl].y[i, j], is_(equal_to(expected[key])))


def describe_the_end_to_end_scan_on_cl2():
    def test_it_matches_the_r_statistic_and_window(telco, expected):
        # Setup.
        exp_stat = expected["telco|scan_cl2|stat"]
        exp_k = int(expected["telco|scan_cl2|k"])
        exp_m = int(expected["telco|scan_cl2|m"])

        # Test.
        sc = fc.wilcoxon_scan(fc.project_l2(telco["cl2"].y.T), 0.25)

        # Assert.
        assert_that(sc["stat"], is_(close_rel(exp_stat, 1e-10)))
        assert_that(sc["k"], is_(equal_to(exp_k)))
        assert_that(sc["m"], is_(equal_to(exp_m)))
