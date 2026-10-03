"""Loaders and matchers shared by the parity suite: the fixtures are
R-generated (tools/make_fixtures.R) full-precision CSVs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from hamcrest import close_to

FIXTURES = Path(__file__).parent / "fixtures"


def load_vec(name: str) -> np.ndarray:
    return np.loadtxt(FIXTURES / name, dtype=float)


def load_mat(name: str) -> np.ndarray:
    return np.loadtxt(FIXTURES / name, dtype=float, delimiter=",")


def close_rel(exp: float, rtol: float):
    """Hamcrest matcher for a RELATIVE tolerance (close_to takes a delta)."""
    return close_to(exp, abs(exp) * rtol)
