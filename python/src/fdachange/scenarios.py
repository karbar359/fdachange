"""Benchmark scenarios S1-S5 (scalar) and F1-F5 (functional).

Same data-generating processes as the R reference implementation. RNG is
numpy's ``Generator`` (pass a ``Generator``, a seed, or ``None``); draw
streams are NOT reproducible across languages - cross-implementation
parity is established on fixed inputs by the test suite instead.
"""

from __future__ import annotations

import numpy as np

from ._interop import _as_rng

__all__ = [
    "generate_s1", "generate_s2", "generate_s3", "generate_s4", "generate_s5",
    "generate_f1", "generate_f2", "generate_f3", "generate_f4", "generate_f5",
    "scenario_cdf",
]


def _rpareto(rng, n, location, shape):
    return location * (1.0 - rng.random(n)) ** (-1.0 / shape)


def _ppareto(q, location, shape):
    q = np.asarray(q, dtype=float)
    return np.where(q >= location, 1.0 - (location / q) ** shape, 0.0)


def _window_slice(n, ec_interval, exact_window):
    lo = round(n * ec_interval[0])
    hi = round(n * ec_interval[1])
    # R (1-based): exact {lo+1..hi} -> 0-based [lo, hi); legacy {lo..hi} ->
    # 0-based [lo-1, hi)
    return slice(lo, hi) if exact_window else slice(lo - 1, hi)


def generate_s1(n, a=0.0, ec_interval=(0.4, 0.7), exact_window=True,
                rng=None):
    """N(0, 1) baseline, epidemic mean shift to N(a, 1)."""
    rng = _as_rng(rng)
    x = rng.normal(0.0, 1.0, n)
    if a != 0:
        w = _window_slice(n, ec_interval, exact_window)
        x[w] = rng.normal(a, 1.0, len(range(*w.indices(n))))
    return x


def generate_s2(n, a=0.0, ec_interval=(0.4, 0.7), exact_window=True,
                rng=None):
    """Squared-normal baseline; in-window N(a,1)^2 (ncp = a^2)."""
    rng = _as_rng(rng)
    x = rng.normal(0.0, 1.0, n) ** 2
    if a != 0:
        w = _window_slice(n, ec_interval, exact_window)
        x[w] = rng.normal(a, 1.0, len(range(*w.indices(n)))) ** 2
    return x


def generate_s3(n, a=0.0, ec_interval=(0.4, 0.7), exact_window=True,
                rng=None):
    """Pareto(1, 3) baseline; location shifts to 1 + a."""
    rng = _as_rng(rng)
    x = _rpareto(rng, n, 1.0, 3.0)
    if a != 0:
        w = _window_slice(n, ec_interval, exact_window)
        x[w] = _rpareto(rng, len(range(*w.indices(n))), 1.0 + a, 3.0)
    return x


def generate_s4(n, a=0.0, ec_interval=(0.4, 0.7), exact_window=True,
                rng=None):
    """Beta(1, 9) baseline; shape shifts to Beta(1 + a, 9 - a)."""
    rng = _as_rng(rng)
    x = rng.beta(1.0, 9.0, n)
    if a != 0:
        w = _window_slice(n, ec_interval, exact_window)
        x[w] = rng.beta(1.0 + a, 9.0 - a, len(range(*w.indices(n))))
    return x


def generate_s5(n, a=0.0, ec_interval=(0.4, 0.7), exact_window=True,
                rng=None):
    """Pareto(1, 3) baseline; shape tightens to 3 + 10a."""
    rng = _as_rng(rng)
    x = _rpareto(rng, n, 1.0, 3.0)
    if a != 0:
        w = _window_slice(n, ec_interval, exact_window)
        x[w] = _rpareto(rng, len(range(*w.indices(n))), 1.0, 3.0 + 10.0 * a)
    return x


def _cosine_curve(s, coefs):
    k = np.arange(1, coefs.size + 1)
    return np.sqrt(2.0) * (np.cos(np.outer(s, k) * np.pi) @ (coefs / k))


def _window_flags(n, ec_interval, exact_window, a):
    lo = round(n * ec_interval[0])
    hi = round(n * ec_interval[1])
    j = np.arange(1, n + 1)                    # 1-based observation index
    j_lo = lo + 1 if exact_window else lo
    return (a != 0) & (j >= j_lo) & (j <= hi)


def generate_f1(n, a=0.0, ec_interval=(0.4, 0.7), l=100, exact_window=True,
                rng=None):
    """Cosine basis; N(0,1) coefficients, centred-Pareto in the window."""
    rng = _as_rng(rng)
    s = np.arange(l + 1) / l
    inw = _window_flags(n, ec_interval, exact_window, a)
    X = np.empty((l + 1, n))
    for j in range(n):
        if inw[j]:
            coefs = _rpareto(rng, 10, 1.0, a) - _rpareto(rng, 10, 1.0, a)
        else:
            coefs = rng.normal(0.0, 1.0, 10)
        X[:, j] = _cosine_curve(s, coefs)
    return X


def generate_f2(n, a=0.0, ec_interval=(0.4, 0.7), l=100, exact_window=True,
                rng=None):
    """Cosine basis; Pareto-difference coefficients, in-window shape a."""
    rng = _as_rng(rng)
    s = np.arange(l + 1) / l
    inw = _window_flags(n, ec_interval, exact_window, a)
    X = np.empty((l + 1, n))
    for j in range(n):
        theta = a if inw[j] else 1.5
        coefs = _rpareto(rng, 10, 1.0, theta) - _rpareto(rng, 10, 1.0, theta)
        X[:, j] = _cosine_curve(s, coefs)
    return X


def generate_f3(n, a=0.0, ec_interval=(0.4, 0.7), l=100, exact_window=True,
                rng=None):
    """Integrated partial-sum (random polygonal) curves, F1's change."""
    rng = _as_rng(rng)
    s = np.arange(l + 1) / l
    d = 10
    inw = _window_flags(n, ec_interval, exact_window, a)
    dt = d * s
    k_floor = np.minimum(np.floor(dt).astype(int), d)
    frac = dt - k_floor
    extra_idx = np.minimum(k_floor, d - 1)     # 0-based xi index
    X = np.empty((l + 1, n))
    for j in range(n):
        if inw[j]:
            xi = _rpareto(rng, d, 1.0, a) - _rpareto(rng, d, 1.0, a)
        else:
            xi = rng.normal(0.0, 1.0, d)
        cs = np.concatenate(([0.0], np.cumsum(xi)))
        X[:, j] = cs[k_floor] + frac * xi[extra_idx]
    return X


def generate_f4(n, a=0.0, ec_interval=(0.25, 0.45), l=100, beta_t=None,
                rng=None):
    """Birbilas's FAR(1) series with epidemic mean jump a*sqrt(t)."""
    from .birbilas import birbilas_betas
    rng = _as_rng(rng)
    if beta_t is None:
        beta_t = birbilas_betas["beta1"]
    s = np.linspace(0.0, 1.0, l + 1)
    beta_grid = beta_t(s)
    k_lo = round(n * ec_interval[0])
    k_hi = round(n * ec_interval[1])
    jump = a * np.sqrt(s)
    X = np.empty((l + 1, n))
    prev = np.zeros(l + 1)
    for i in range(1, n + 1):
        Y = np.concatenate(([0.0],
                            np.cumsum(np.sqrt(1.0 / l)
                                      * rng.normal(0.0, 1.0, l))))
        cur = beta_grid * prev + Y
        if a != 0 and k_lo < i <= k_hi:
            cur = cur + jump
        X[:, i - 1] = cur
        prev = cur
    return X


def generate_f5(n, a=0.0, ec_interval=(0.4, 0.7), l=100, spike_alpha=0.45,
                exact_window=True, rng=None):
    """I.i.d. Wiener curves; in-window spike a*f_alpha(t) with the
    dissertation's spike function (1 - |t - 1/2|^alpha)^(1/alpha)."""
    rng = _as_rng(rng)
    s = np.arange(l + 1) / l
    spike = a * (1.0 - np.abs(s - 0.5) ** spike_alpha) ** (1.0 / spike_alpha)
    inw = _window_flags(n, ec_interval, exact_window, a)
    X = np.empty((l + 1, n))
    for j in range(n):
        W = np.concatenate(([0.0],
                            np.cumsum(np.sqrt(1.0 / l)
                                      * rng.normal(0.0, 1.0, l))))
        X[:, j] = W + spike if inw[j] else W
    return X


def scenario_cdf(scenario: str):
    """Null CDF matching a scenario, for the ``cdf`` argument of the Lp
    scan; the functional scenarios use the empirical CDF (also ``cdf=None``)."""
    from scipy import stats

    s = scenario.upper()
    if s == "S1":
        return lambda x: stats.norm.cdf(x)
    if s == "S2":
        return lambda x: stats.chi2.cdf(x, df=1)
    if s in ("S3", "S5"):
        return lambda x: _ppareto(x, 1.0, 3.0)
    if s == "S4":
        return lambda x: stats.beta.cdf(x, 1.0, 9.0)
    if s in ("F1", "F2", "F3", "F4", "F5"):
        def ecdf(x):
            x = np.asarray(x, dtype=float)
            return np.searchsorted(np.sort(x), x, side="right") / x.size
        return ecdf
    raise ValueError(f"unknown scenario: {scenario!r}")
