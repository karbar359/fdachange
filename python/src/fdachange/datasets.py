"""Shipped datasets (converted from the R package's data/)."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources

import numpy as np

__all__ = ["TelcoCluster", "load_telco_far1"]


@dataclass(frozen=True)
class TelcoCluster:
    """One cluster of the ``telco_far1`` dataset.

    ``y`` is 50 months x 28 days, as in the R object; transpose (``y.T``)
    to get the package's ``(grid, n)`` one-curve-per-column layout.
    """

    y: np.ndarray               #: (50, 28) average daily voice minutes
    change_months: np.ndarray   #: designed epidemic window (empty if none)
    multiplier: np.ndarray      #: (50,) level multiplier f_j
    jump: np.ndarray            #: (50,) sqrt(t)-jump size c_j (level share)
    outliers: np.ndarray        #: (k, 2) int (month, day) outlier days
    outlier_months: np.ndarray  #: months with a raised whole-month level


def load_telco_far1() -> dict[str, TelcoCluster]:
    """Simulated FAR(1) telco usage curves (the R dataset ``telco_far1``).

    Daily voice-usage curves for 11 simulated customer clusters, 50 months
    each, following a functional autoregressive FAR(1) process with
    ``beta(t) = 0.7 * sqrt(t)`` (``far1_betas["beta1"]``). Clusters
    ``cl2``-``cl4`` contain a designed epidemic change (sqrt(t)-shaped
    jumps as in ``generate_f4``); ``cl5``-``cl11`` have no change but
    trends, level shifts or outliers, for checking behaviour when the null
    is violated in other ways; ``cl1`` is a clean null. The full design
    table and the generating model are documented in the R package
    (``?telco_far1``) and built by ``data-raw/telco_far1.R``.

    Values are bit-identical to the R dataset (full-precision conversion,
    asserted by ``tools/build_telco_data.py``).

    Calibration caveat: months are serially dependent, so the shipped null
    caches (which assume independent observations) are not exact for this
    data. Simulate the null from the same FAR(1) model and pass it via
    ``null=`` in :func:`fdachange.epidemic_test` for exact calibration.

    Returns a dict ``{"cl1": TelcoCluster, ..., "cl11": TelcoCluster}``.
    """
    pkg_file = resources.files("fdachange").joinpath("_data/telco_far1.npz")
    out: dict[str, TelcoCluster] = {}
    with pkg_file.open("rb") as fh:
        with np.load(fh) as z:
            for i in range(1, 12):
                cl = f"cl{i}"
                out[cl] = TelcoCluster(
                    y=np.array(z[f"{cl}_y"]),
                    change_months=np.array(z[f"{cl}_change_months"]),
                    multiplier=np.array(z[f"{cl}_multiplier"]),
                    jump=np.array(z[f"{cl}_jump"]),
                    outliers=np.array(z[f"{cl}_outliers"]),
                    outlier_months=np.array(z[f"{cl}_outlier_months"]),
                )
    return out
