from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

FIXTURES = Path(__file__).parent / "fixtures"


def load_vec(name: str) -> np.ndarray:
    return np.loadtxt(FIXTURES / name, dtype=float)


def load_mat(name: str) -> np.ndarray:
    return np.loadtxt(FIXTURES / name, dtype=float, delimiter=",")


@pytest.fixture(scope="session")
def expected() -> dict[str, float]:
    out = {}
    lines = (FIXTURES / "expected_scalars.csv").read_text().splitlines()[1:]
    for line in lines:
        key, value = line.rsplit(",", 1)
        out[key.strip('"')] = float(value.strip('"'))
    return out


@pytest.fixture(scope="session")
def x60() -> np.ndarray:
    return load_vec("input_x60.csv")


@pytest.fixture(scope="session")
def F4() -> np.ndarray:
    return load_mat("input_F4.csv")


@pytest.fixture(scope="session")
def F5() -> np.ndarray:
    return load_mat("input_F5.csv")
