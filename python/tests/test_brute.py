"""Printed-formula brute forces in pure Python: an R-independent proof
that the vectorised implementations compute the dissertation's formulas."""

from __future__ import annotations

import numpy as np
import pytest

import fdachange as fc

TOL = 1e-10


def brute_wilcoxon(x, gamma, kernel="h_w"):
    xv = np.abs(x) if kernel == "h_0" else x
    n = xv.size
    H = np.sign(xv[:, None] - xv[None, :])
    sigma = np.sqrt(np.mean(H.mean(axis=0) ** 2))
    w = fc.holder_weights(n, gamma, 0.0, symmetric=True)
    best = 0.0
    for k in range(n):                       # window {k+1..m}, 1-based R
        for m in range(k + 1, n + 1):
            I = np.arange(k, m)
            comp = np.setdiff1d(np.arange(n), I)
            v = abs(H[np.ix_(I, comp)].sum()) * w[m - k - 1]
            best = max(best, v)
    return best / (n ** 1.5 * sigma)


def brute_lp(x, gamma, p=2, t_grid=100):
    from scipy.stats import norm
    n = x.size
    Fx = norm.cdf(x)
    w = fc.holder_weights(n, gamma)
    ts = np.linspace(0, 1, t_grid + 1)
    ind = (Fx[:, None] <= ts[None, :]).astype(float) - ts[None, :]
    tot = ind.sum(axis=0)
    best = 0.0
    for k in range(n):
        for l in range(1, n - k + 1):
            cs = ind[k:k + l].sum(axis=0) - (l / n) * tot
            v = (np.sum(np.abs(cs) ** p) / t_grid) ** (1 / p) * w[l - 1]
            best = max(best, v)
    return best / np.sqrt(n)


def brute_rkhs(x, gamma, kernel):
    n = x.size
    K = fc.gram(x, kernel)
    w = fc.holder_weights(n, gamma)
    best = 0.0
    for l in range(1, n // 2 + 1):
        for a in range(0, n - l + 1):
            I = np.arange(a, a + l)
            ln = l / n
            d2 = K[np.ix_(I, I)].sum() + ln ** 2 * K.sum() \
                - 2 * ln * K[I, :].sum()
            best = max(best, np.sqrt(max(d2, 0.0)) * w[l - 1])
    return best / np.sqrt(n)


def brute_far1(X, gamma, beta_t):
    L, n = X.shape
    s = np.linspace(0, 1, L)
    W = (1 - beta_t(s))[:, None] * (X - X.mean(axis=1, keepdims=True))
    best = 0.0
    for k in range(n):
        for m in range(k + 1, n + 1):
            v = np.abs(W[:, k:m].sum(axis=1)).max() / (m - k) ** gamma
            best = max(best, v)
    return best * n ** (-0.5 + gamma)


def test_wilcoxon_brute(x60):
    for g in (0.0, 0.25, 0.45):
        assert fc.wilcoxon_statistic(x60, g) == pytest.approx(
            brute_wilcoxon(x60, g), abs=TOL)
    assert fc.wilcoxon_statistic(x60, 0.25, kernel="h_0") == pytest.approx(
        brute_wilcoxon(x60, 0.25, "h_0"), abs=TOL)


def test_lp_brute(x60):
    from scipy.stats import norm
    for g in (0.0, 0.25, 0.45):
        assert fc.lp_statistic(x60, g, cdf=norm.cdf) == pytest.approx(
            brute_lp(x60, g), abs=TOL)


def test_rkhs_brute(x60):
    for g in (0.0, 0.25, 0.45):
        assert fc.rkhs_statistic(x60, g, kernel="gaussian") == pytest.approx(
            brute_rkhs(x60, g, fc.kernel_gaussian), abs=TOL)


def test_far1_brute(F4):
    for g in (0.25, 0.4):
        assert fc.far1_statistic(F4, g) == pytest.approx(
            brute_far1(F4, g, fc.far1_betas["beta1"]), abs=TOL)
