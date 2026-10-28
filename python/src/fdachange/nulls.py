"""Null-draw storage: packaged caches plus a user-side cache directory.

The packaged ``.npz`` files are converted from the same frozen CSV sources
as the R package's extdata (``data-raw/cache-src`` in the repository), so
critical values are bit-identical across the two implementations. Cell
keys match the R package exactly (``g0.25_b0``, ``g0.25_p2``,
``g0.45_2max_rbeta.1.9``, ``g0.4``, ...).
"""

from __future__ import annotations

import os
from importlib import resources
from pathlib import Path

import numpy as np

from ._kernels import _kernel_name
from ._weights import critical_value

__all__ = [
    "null_cache_dir", "save_null_draws",
    "wilcoxon_null_draws", "wilcoxon_critical_value",
    "lp_null_draws", "lp_critical_value",
    "rkhs_null_draws", "rkhs_critical_value",
    "far1_null_draws", "far1_critical_value",
]

_METHODS = ("wilcoxon", "lp", "rkhs", "far1")
_TABLES: dict[str, dict[str, np.ndarray]] = {}


def _fmt(x) -> str:
    s = f"{float(x):g}"
    return s


def _null_key(method: str, gamma, beta=0, p=2, kernel="gaussian",
              baseline="") -> str:
    if method == "wilcoxon":
        return f"g{_fmt(gamma)}_b{_fmt(beta)}"
    if method == "lp":
        return f"g{_fmt(gamma)}_p{_fmt(p)}"
    if method == "rkhs":
        suffix = f"_{baseline}" if baseline else ""
        return f"g{_fmt(gamma)}_{_kernel_name(kernel)}{suffix}"
    if method == "far1":
        return f"g{_fmt(gamma)}"
    raise ValueError(f"unknown method: {method}")


def null_cache_dir() -> Path:
    """Directory for user-generated null draws.

    ``FDACHANGE_CACHE_DIR`` overrides; otherwise
    ``$XDG_CACHE_HOME/fdachange`` or ``~/.cache/fdachange``. User entries
    take precedence over the packaged caches.
    """
    env = os.environ.get("FDACHANGE_CACHE_DIR")
    if env:
        return Path(env)
    base = os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache")
    return Path(base) / "fdachange"


def _load_table(method: str) -> dict[str, np.ndarray]:
    if method in _TABLES:
        return _TABLES[method]
    table: dict[str, np.ndarray] = {}
    pkg_file = resources.files("fdachange").joinpath(f"_data/{method}_null.npz")
    if pkg_file.is_file():
        with pkg_file.open("rb") as fh:
            with np.load(fh) as z:
                table.update({k: np.array(z[k]) for k in z.files})
    user_file = null_cache_dir() / f"{method}_null.npz"
    if user_file.is_file():
        with np.load(user_file) as z:
            table.update({k: np.array(z[k]) for k in z.files})
    _TABLES[method] = table
    return table


def _invalidate():
    _TABLES.clear()


def _null_get(method: str, key: str, hint: str) -> np.ndarray:
    table = _load_table(method)
    if key not in table:
        cells = ", ".join(sorted(table)) or "(none)"
        raise KeyError(
            f"No {method} null draws for cell {key!r}. Available cells: "
            f"{cells}. Generate the missing cell with {hint} and store it "
            f"with save_null_draws().")
    return table[key]


def save_null_draws(draws, method: str, gamma, *, beta=0, p=2,
                    kernel="gaussian", baseline="") -> str:
    """Store draws in :func:`null_cache_dir` under the standard cell key."""
    if method not in _METHODS:
        raise ValueError(f"method must be one of {_METHODS}")
    draws = np.asarray(draws, dtype=float)
    if draws.size == 0:
        raise ValueError("draws must be non-empty")
    key = _null_key(method, gamma, beta=beta, p=p, kernel=kernel,
                    baseline=baseline)
    path = null_cache_dir() / f"{method}_null.npz"
    existing: dict[str, np.ndarray] = {}
    if path.is_file():
        with np.load(path) as z:
            existing = {k: np.array(z[k]) for k in z.files}
    existing[key] = draws
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **existing)
    _invalidate()
    return key


def wilcoxon_null_draws(gamma, beta=0) -> np.ndarray:
    """Null draws for the Wilcoxon scan: the cached cell for (gamma, beta)."""
    return _null_get("wilcoxon", _null_key("wilcoxon", gamma, beta=beta),
                     f"simulate_wilcoxon_null(n_draws, gamma={gamma}, beta={beta})")


def lp_null_draws(gamma, p=2) -> np.ndarray:
    return _null_get("lp", _null_key("lp", gamma, p=p),
                     f"simulate_lp_null(n_draws, gamma={gamma}, p={p})")


def rkhs_null_draws(gamma, kernel="gaussian", baseline="") -> np.ndarray:
    return _null_get("rkhs",
                     _null_key("rkhs", gamma, kernel=kernel, baseline=baseline),
                     f"simulate_rkhs_null(n_draws, gamma={gamma}, "
                     f"kernel={kernel!r}, rdist=<H0 sampler>)")


def far1_null_draws(gamma) -> np.ndarray:
    return _null_get("far1", _null_key("far1", gamma),
                     f"simulate_far1_null(n_draws, gamma={gamma})")


def wilcoxon_critical_value(gamma, beta=0, alpha=0.05, splits=1) -> float:
    return critical_value(wilcoxon_null_draws(gamma, beta), alpha, splits)


def lp_critical_value(gamma, p=2, alpha=0.05, splits=1) -> float:
    return critical_value(lp_null_draws(gamma, p), alpha, splits)


def rkhs_critical_value(gamma, kernel="gaussian", baseline="",
                        alpha=0.05, splits=1) -> float:
    return critical_value(rkhs_null_draws(gamma, kernel, baseline),
                          alpha, splits)


def far1_critical_value(gamma, alpha=0.05, splits=1) -> float:
    return critical_value(far1_null_draws(gamma), alpha, splits)
