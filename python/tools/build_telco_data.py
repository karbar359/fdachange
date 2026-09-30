#!/usr/bin/env python3
"""Pack the telco_far1 CSVs exported by export_telco.R (tools/telco-src)
into the packaged src/fdachange/_data/telco_far1.npz, with exact float64
round-trip asserts on every cell.

Run:  Rscript python/tools/export_telco.R
      python python/tools/build_telco_data.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
SRC = HERE.parent / "telco-src"
OUT = HERE.parents[1] / "src" / "fdachange" / "_data"

CLUSTERS = [f"cl{i}" for i in range(1, 12)]


def read_float_mat(path: Path) -> np.ndarray:
    with path.open(newline="") as fh:
        rows = [[float(c) for c in row] for row in csv.reader(fh)]
    return np.array(rows, dtype=np.float64)


def read_float_vec(path: Path) -> np.ndarray:
    return np.array([float(l) for l in path.read_text().split()],
                    dtype=np.float64)


def read_int_vec(path: Path) -> np.ndarray:
    text = path.read_text().split()
    return np.array([int(l) for l in text], dtype=np.int64)


def read_outliers(path: Path) -> np.ndarray:
    with path.open(newline="") as fh:
        rows = list(csv.reader(fh))
    assert rows[0] == ["month", "day"], rows[0]
    out = np.array([[int(c) for c in row] for row in rows[1:]],
                   dtype=np.int64)
    return out.reshape(-1, 2)


def main() -> None:
    table: dict[str, np.ndarray] = {}
    for cl in CLUSTERS:
        y = read_float_mat(SRC / f"{cl}_y.csv")
        assert y.shape == (50, 28), (cl, y.shape)
        table[f"{cl}_y"] = y
        table[f"{cl}_change_months"] = read_int_vec(
            SRC / f"{cl}_change_months.csv")
        table[f"{cl}_multiplier"] = read_float_vec(
            SRC / f"{cl}_multiplier.csv")
        table[f"{cl}_jump"] = read_float_vec(SRC / f"{cl}_jump.csv")
        table[f"{cl}_outliers"] = read_outliers(SRC / f"{cl}_outliers.csv")
        table[f"{cl}_outlier_months"] = read_int_vec(
            SRC / f"{cl}_outlier_months.csv")
        assert table[f"{cl}_multiplier"].shape == (50,), cl
        assert table[f"{cl}_jump"].shape == (50,), cl

    out_file = OUT / "telco_far1.npz"
    np.savez_compressed(out_file, **table)

    # exact round-trip on every array
    with np.load(out_file) as z:
        assert sorted(z.files) == sorted(table), "key mismatch"
        for k in table:
            back = z[k]
            assert back.dtype == table[k].dtype, k
            assert np.array_equal(back, table[k]), k
    n = sum(v.size for v in table.values())
    print(f"telco_far1: {len(CLUSTERS)} clusters, {n} values -> {out_file} "
          f"(exact round-trip OK)")


if __name__ == "__main__":
    main()
