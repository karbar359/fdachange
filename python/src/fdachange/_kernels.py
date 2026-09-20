"""Built-in RKHS kernels and the Gram matrix."""

from __future__ import annotations

import numpy as np

__all__ = ["kernel_gaussian", "kernel_min", "kernel_2max", "gram"]


def kernel_gaussian(x, y):
    """Gaussian kernel ``exp(-(x - y)^2 / 2)``."""
    return np.exp(-(np.asarray(x, float) - np.asarray(y, float)) ** 2 / 2.0)


def kernel_min(x, y):
    """Min (Wiener) kernel ``min(x, y)``."""
    return np.minimum(np.asarray(x, float), np.asarray(y, float))


def kernel_2max(x, y):
    """2-max kernel ``2 - max(x, y)``."""
    return 2.0 - np.maximum(np.asarray(x, float), np.asarray(y, float))


_KERNELS = {
    "gaussian": kernel_gaussian,
    "min": kernel_min,
    "2max": kernel_2max,
    "2_max": kernel_2max,
    "2.max": kernel_2max,
}


def _resolve_kernel(kernel):
    if callable(kernel):
        return kernel
    try:
        return _KERNELS[str(kernel)]
    except KeyError:
        raise ValueError(f"unknown kernel: {kernel!r}") from None


def _kernel_name(kernel) -> str:
    if isinstance(kernel, str):
        k = "2max" if kernel in ("2_max", "2.max") else kernel
        if k in ("gaussian", "min", "2max"):
            return k
        raise ValueError(f"unknown kernel name: {kernel!r}")
    for name, fn in (("gaussian", kernel_gaussian), ("min", kernel_min),
                     ("2max", kernel_2max)):
        if kernel is fn:
            return name
    raise ValueError("custom kernel functions have no cache name; "
                     "pass explicit null draws instead")


def gram(x, kernel) -> np.ndarray:
    """Gram matrix ``K[i, j] = kernel(x_i, x_j)`` of a numeric sample."""
    x = np.asarray(x, dtype=float)
    fn = _resolve_kernel(kernel)
    return np.asarray(fn(x[:, None], x[None, :]), dtype=float)
