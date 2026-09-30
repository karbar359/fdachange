"""telco_far1 dataset: structure, R-conversion parity, end-to-end scan."""

from __future__ import annotations

import numpy as np
import pytest

import fdachange as fc

CLUSTERS = [f"cl{i}" for i in range(1, 12)]


@pytest.fixture(scope="module")
def telco():
    return fc.load_telco_far1()


def test_structure(telco):
    assert list(telco) == CLUSTERS
    for cl, c in telco.items():
        assert c.y.shape == (50, 28), cl
        assert c.y.dtype == np.float64
        assert c.multiplier.shape == (50,)
        assert c.jump.shape == (50,)
        assert c.outliers.ndim == 2 and c.outliers.shape[1] == 2
        assert np.all(c.y >= 0), cl


def test_design(telco):
    # designed change windows (data-raw/telco_far1.R)
    assert telco["cl1"].change_months.size == 0
    assert np.array_equal(telco["cl2"].change_months, np.arange(37, 47))
    assert np.array_equal(telco["cl3"].change_months, np.arange(37, 50))
    assert np.array_equal(telco["cl4"].change_months, np.arange(37, 51))
    for cl in CLUSTERS[4:]:
        assert telco[cl].change_months.size == 0, cl
    # jumps only where designed
    assert np.all(telco["cl2"].jump[36:46] == 0.25)
    assert np.all(telco["cl3"].jump[np.setdiff1d(
        telco["cl3"].change_months, telco["cl3"].outlier_months) - 1] == -0.25)
    assert telco["cl4"].jump[36] == 0 and telco["cl4"].jump[49] < -0.3
    # outliers: cl3 has 3 outlier days, cl9 has 5% of its non-outlier months
    assert telco["cl3"].outliers.shape == (3, 2)
    assert telco["cl9"].outliers.shape[0] == round(
        0.05 * (50 - telco["cl9"].outlier_months.size) * 28)
    # level designs
    assert np.all(telco["cl7"].multiplier == 1.5)
    assert np.all(telco["cl8"].multiplier == 0.6)


def test_r_parity(telco, expected):
    for cl, c in telco.items():
        assert c.y.sum() == pytest.approx(
            expected[f"telco|{cl}|sum"], rel=1e-12)
    assert telco["cl2"].y[0, 0] == expected["telco|cl2|y[1,1]"]
    assert telco["cl2"].y[36, 27] == expected["telco|cl2|y[37,28]"]
    assert telco["cl9"].y[49, 13] == expected["telco|cl9|y[50,14]"]


def test_end_to_end_scan(telco, expected):
    sc = fc.wilcoxon_scan(fc.project_l2(telco["cl2"].y.T), 0.25)
    assert sc["stat"] == pytest.approx(expected["telco|scan_cl2|stat"],
                                       rel=1e-10)
    assert sc["k"] == expected["telco|scan_cl2|k"]
    assert sc["m"] == expected["telco|scan_cl2|m"]
