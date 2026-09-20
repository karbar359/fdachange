from __future__ import annotations

import numpy as np
import pytest

import fdachange as fc


def test_result_is_coherent():
    x = fc.generate_s1(200, a=1, ec_interval=(0.4, 0.7), rng=3)
    for m in ("wilcoxon", "lp", "rkhs"):
        r = fc.epidemic_test(x, m, gamma=0.25)
        assert 0.0 <= r.p_value <= 1.0
        assert (r.p_value < r.alpha) == (r.statistic > r.critical_value)
        assert 0 <= r.k < r.m <= 200
        assert r.reject == (r.statistic > r.critical_value)


def test_strong_signal_detected_and_localised():
    x = fc.generate_s1(200, a=3, ec_interval=(0.4, 0.7), rng=4)
    r = fc.epidemic_test(x, "wilcoxon", gamma=0.25)
    assert r.p_value < 0.01
    assert r.m > 100 and r.k < 140


def test_kuiper_exact_path():
    x = fc.generate_s1(150, rng=3)
    r = fc.epidemic_test(x, "wilcoxon", gamma=0.0)
    assert np.isinf(r.n_null)
    assert 1.0 - r.p_value == pytest.approx(fc.kuiper_cdf(r.statistic),
                                            rel=1e-12)
    assert r.critical_value == pytest.approx(fc.kuiper_quantile(0.95),
                                             rel=1e-9)


def test_null_injection():
    x = fc.generate_s1(100, rng=5)
    draws = np.linspace(0.0, 10.0, 1001)
    r = fc.epidemic_test(x, "rkhs", gamma=0.25, null=draws)
    assert r.critical_value == pytest.approx(
        float(np.quantile(draws, 0.95, method="linear")))
    assert r.n_null == 1001
    assert 1.0 - r.p_value == pytest.approx(float(np.mean(draws <= r.statistic)))


def test_splits_aggregation():
    x = fc.generate_s1(200, rng=6)
    r = fc.epidemic_test(x, "wilcoxon", gamma=0.25, splits=2)
    s1 = fc.wilcoxon_statistic(x[:100], 0.25)
    s2 = fc.wilcoxon_statistic(x[100:], 0.25)
    assert r.statistic == pytest.approx(max(s1, s2), rel=1e-12)
    assert r.critical_value == pytest.approx(
        fc.wilcoxon_critical_value(0.25, splits=2), rel=1e-12)
    assert r.k is None
    d = fc.wilcoxon_null_draws(0.25)
    assert r.p_value == pytest.approx(
        1.0 - float(np.mean(d <= max(s1, s2))) ** 2, rel=1e-12)


def test_birbilas_on_matrix():
    X = fc.generate_f4(60, a=2, l=50, rng=7)
    r = fc.epidemic_test(X, "birbilas", gamma=0.4)
    assert r.statistic == pytest.approx(fc.birbilas_statistic(X, 0.4),
                                        rel=1e-12)
    assert 0 <= r.k < r.m <= 60


def test_custom_rkhs_kernel_requires_null():
    x = fc.generate_s1(50, rng=8)
    my_kernel = lambda a, b: np.exp(-np.abs(a - b))
    with pytest.raises(ValueError, match="null"):
        fc.epidemic_test(x, "rkhs", gamma=0.25, kernel=my_kernel)
    r = fc.epidemic_test(x, "rkhs", gamma=0.25, kernel=my_kernel,
                         null=np.linspace(0, 5, 200))
    assert isinstance(r, fc.EpidemicTestResult)


def test_str_report():
    x = fc.generate_s1(100, a=2, rng=9)
    out = str(fc.epidemic_test(x, "wilcoxon", gamma=0.25))
    assert "critical value" in out and "detected window" in out
