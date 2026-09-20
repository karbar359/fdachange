from __future__ import annotations

import numpy as np
import pytest

import fdachange as fc
from fdachange import nulls


def test_shipped_cells():
    assert fc.wilcoxon_null_draws(0.25).size == 30_000
    assert fc.wilcoxon_null_draws(0.5, beta=1).size == 30_000
    assert fc.lp_null_draws(0.45, p=10).size == 300
    assert fc.rkhs_null_draws(0.45, "2max", "rbeta.1.9").size == 1_000
    assert fc.birbilas_null_draws(0.4).size == 1_000


def test_kuiper_grid_at_gamma_zero():
    d = fc.wilcoxon_null_draws(0.0)
    assert d.size == 200_000
    assert np.all(np.diff(d) >= 0)            # deterministic quantile grid
    assert np.quantile(d, 0.95) == pytest.approx(1.7473, abs=5e-4)


def test_missing_cell_names_the_simulator():
    with pytest.raises(KeyError, match="simulate_wilcoxon_null"):
        fc.wilcoxon_null_draws(0.33)
    with pytest.raises(KeyError, match="simulate_lp_null"):
        fc.lp_null_draws(0.25, p=3)
    with pytest.raises(KeyError, match="simulate_rkhs_null"):
        fc.rkhs_null_draws(0.25, "min")


def test_user_cache_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setenv("FDACHANGE_CACHE_DIR", str(tmp_path))
    nulls._invalidate()
    try:
        assert fc.null_cache_dir() == tmp_path
        draws = np.linspace(0.0, 5.0, 500)
        key = fc.save_null_draws(draws, "wilcoxon", 0.33)
        assert key == "g0.33_b0"
        np.testing.assert_array_equal(fc.wilcoxon_null_draws(0.33), draws)
        # user entry overrides a shipped cell
        fc.save_null_draws(np.arange(1, 101) / 100, "lp", 0.25, p=2)
        assert fc.lp_null_draws(0.25, p=2).size == 100
        # rkhs keying round-trips kernel + baseline
        fc.save_null_draws(draws, "rkhs", 0.1, kernel="2max",
                           baseline="runif")
        np.testing.assert_array_equal(
            fc.rkhs_null_draws(0.1, "2max", "runif"), draws)
    finally:
        nulls._invalidate()


def test_critical_value_definition():
    d = fc.lp_null_draws(0.25, p=2)
    assert fc.lp_critical_value(0.25, p=2) == pytest.approx(
        float(np.quantile(d, 0.95, method="linear")), rel=1e-15)
    assert fc.lp_critical_value(0.25, p=2, splits=2) == pytest.approx(
        float(np.quantile(d, np.sqrt(0.95), method="linear")), rel=1e-15)
