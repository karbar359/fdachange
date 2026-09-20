from __future__ import annotations

import numpy as np
import pytest

import fdachange as fc


def test_shapes():
    assert fc.generate_s1(50, rng=1).shape == (50,)
    for gen in (fc.generate_f1, fc.generate_f2, fc.generate_f3,
                fc.generate_f4, fc.generate_f5):
        assert gen(20, a=2, l=30, rng=1).shape == (31, 20)


def test_exact_vs_legacy_window_counts():
    # short regime: l = 6 at n = 400 with interval (0.4925, 0.5075)
    def changed(exact):
        x = fc.generate_s1(400, a=50, ec_interval=(0.4925, 0.5075),
                           exact_window=exact, rng=1)
        return int(np.sum(x > 25))
    assert changed(True) == 6
    assert changed(False) == 7


def test_s2_noncentrality():
    x = fc.generate_s2(200_000, a=2, ec_interval=(0.1, 0.9), rng=7)
    assert np.mean(x[20_000:180_000]) == pytest.approx(1 + 2 ** 2, abs=0.06)


def test_f5_spike_is_deterministic_given_stream():
    s = np.arange(101) / 100

    def f_alpha(t, al):
        return (1 - np.abs(t - 0.5) ** al) ** (1 / al)

    Xa = fc.generate_f5(50, a=3, spike_alpha=0.45, rng=9)
    Xb = fc.generate_f5(50, a=3, spike_alpha=0.2, rng=9)
    diff_expected = 3 * (f_alpha(s, 0.45) - f_alpha(s, 0.2))
    win = slice(20, 35)                       # exact window at n = 50
    np.testing.assert_allclose(
        Xa[:, win] - Xb[:, win],
        np.broadcast_to(diff_expected[:, None], Xa[:, win].shape),
        rtol=0, atol=1e-12)
    out = np.r_[0:20, 35:50]
    assert np.array_equal(Xa[:, out], Xb[:, out])


def test_scenario_cdf():
    from scipy.stats import norm
    x = np.array([-1.0, 0.0, 1.5])
    np.testing.assert_allclose(fc.scenario_cdf("S1")(x), norm.cdf(x))
    np.testing.assert_allclose(fc.scenario_cdf("S3")(np.array([1.5, 3.0])),
                               1 - np.array([1.5, 3.0]) ** -3.0)
    assert fc.scenario_cdf("S3")(np.array([0.5]))[0] == 0.0
    y = fc.generate_s1(30, rng=2)
    e = fc.scenario_cdf("F5")(y)
    assert e.min() > 0 and e.max() == 1.0
    with pytest.raises(ValueError):
        fc.scenario_cdf("X9")


def test_rng_reproducibility():
    a = fc.generate_f1(10, a=3, rng=5)
    b = fc.generate_f1(10, a=3, rng=5)
    assert np.array_equal(a, b)
    c = fc.generate_f1(10, a=3, rng=6)
    assert not np.array_equal(a, c)
