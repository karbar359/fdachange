"""Calibrated epidemic change-point test, mirroring the R epidemic_test()."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._kuiper import kuiper_cdf, kuiper_quantile
from ._weights import critical_value
from .wilcoxon import wilcoxon_scan
from .lp import lp_scan
from .rkhs import rkhs_scan
from .far1 import far1_scan, far1_betas
from .splits import scan_over_splits
from . import nulls

__all__ = ["epidemic_test", "EpidemicTestResult"]


@dataclass
class EpidemicTestResult:
    statistic: float
    p_value: float
    critical_value: float
    alpha: float
    k: int | None
    m: int | None
    interval: tuple
    method: str
    gamma: float
    splits: int
    n_null: float
    parameters: dict = field(default_factory=dict)

    @property
    def reject(self) -> bool:
        return self.statistic > self.critical_value

    def __str__(self) -> str:  # htest-flavoured report
        lines = [
            f"\t{self.method}",
            "",
            f"T = {self.statistic:.4f}, gamma = {self.gamma:g}, "
            f"p-value = {self.p_value:.4g}",
            "alternative hypothesis: an epidemic (changed then restored) segment",
            f"critical value (alpha = {self.alpha:g}): "
            f"{self.critical_value:.4f}  "
            f"[{'closed-form' if np.isinf(self.n_null) else int(self.n_null)}"
            f"{'' if np.isinf(self.n_null) else ' null draws'}]",
        ]
        if self.k is not None:
            lines.append(
                f"detected window: observations {self.k + 1}..{self.m} "
                f"(fractions {self.interval[0]:.3f}..{self.interval[1]:.3f})")
        return "\n".join(lines)


_LABELS = {
    "wilcoxon": "Epidemic change scan test (Wilcoxon-type antisymmetric kernel)",
    "lp": "Epidemic change scan test (Lp empirical process)",
    "rkhs": "Epidemic change scan test (RKHS kernel mean embedding)",
    "far1": "Epidemic change test for FAR(1) (Birbilas E_n,beta,gamma)",
}


def epidemic_test(x, method: str, gamma: float, *, beta: float = 0.0,
                  kernel=None, p: float = 2, cdf=None, t_grid: int = 100,
                  beta_t=None, baseline: str = "", alpha: float = 0.05,
                  splits: int = 1, null=None) -> EpidemicTestResult:
    """Run one scan statistic against its null calibration.

    Calibration resolution mirrors the R package: the closed-form Kuiper
    law for Wilcoxon at gamma = 0 (beta = 0), otherwise the packaged
    Monte-Carlo caches, unless explicit draws are passed via ``null=``
    (required for custom RKHS kernels or uncached baselines).
    """
    if method not in _LABELS:
        raise ValueError(f"method must be one of {tuple(_LABELS)}")
    if beta_t is None:
        beta_t = far1_betas["beta1"]

    if method == "wilcoxon":
        scan = lambda z: wilcoxon_scan(z, gamma, beta,
                                       kernel=kernel or "h_w")
    elif method == "lp":
        scan = lambda z: lp_scan(z, gamma, p=p, cdf=cdf, t_grid=t_grid)
    elif method == "rkhs":
        scan = lambda z: rkhs_scan(z, gamma, kernel=kernel or "gaussian",
                                   beta=beta)
    else:
        scan = lambda z: far1_scan(z, gamma, beta_t=beta_t)

    if splits == 1:
        sc = scan(x)
        stat, k, m, interval = sc["stat"], sc["k"], sc["m"], sc["interval"]
    else:
        stat = scan_over_splits(x, splits, lambda chunk: scan(chunk)["stat"])
        k = m = None
        interval = (None, None)

    exact_kuiper = (null is None and method == "wilcoxon"
                    and float(gamma) == 0.0 and float(beta) == 0.0)
    if exact_kuiper:
        p_value = 1.0 - float(kuiper_cdf(stat)) ** splits
        cv = float(kuiper_quantile((1.0 - alpha) ** (1.0 / splits)))
        n_null = float("inf")
    else:
        if null is not None:
            draws = np.asarray(null, dtype=float)
            if draws.size < 2:
                raise ValueError("null= needs a vector of draws")
        elif method == "wilcoxon":
            draws = nulls.wilcoxon_null_draws(gamma, beta)
        elif method == "lp":
            draws = nulls.lp_null_draws(gamma, p)
        elif method == "rkhs":
            if callable(kernel):
                raise ValueError(
                    "rkhs with a custom kernel function needs explicit "
                    "null draws: pass null=simulate_rkhs_null(...)")
            draws = nulls.rkhs_null_draws(gamma, kernel or "gaussian",
                                          baseline)
        else:
            draws = nulls.far1_null_draws(gamma)
        p_value = 1.0 - float(np.mean(draws <= stat)) ** splits
        cv = critical_value(draws, alpha, splits)
        n_null = float(draws.size)

    params = {"gamma": gamma}
    if beta:
        params["beta"] = beta
    if method == "lp":
        params["p"] = p

    return EpidemicTestResult(
        statistic=float(stat), p_value=float(p_value), critical_value=cv,
        alpha=alpha, k=k, m=m, interval=interval, method=_LABELS[method],
        gamma=float(gamma), splits=int(splits), n_null=n_null,
        parameters=params)
