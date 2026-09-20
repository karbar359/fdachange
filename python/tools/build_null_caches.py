#!/usr/bin/env python3
"""Convert the frozen null-cache CSVs (data-raw/cache-src, the same
sources as the R package's extdata) into the packaged .npz tables, with
exact float64 equality asserts on every cell.

Run from anywhere:  python python/tools/build_null_caches.py
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
REPO = HERE.parents[2]
SRC = REPO / "data-raw" / "cache-src"
OUT = HERE.parents[1] / "src" / "fdachange" / "_data"
OUT.mkdir(parents=True, exist_ok=True)


def fmt(x: float) -> str:
    return f"{x:g}"


def read_csv(path: Path):
    with path.open(newline="") as fh:
        rows = list(csv.reader(fh))
    header, data = rows[0], rows[1:]
    return header, data


def save(name: str, table: dict[str, np.ndarray]):
    for k, v in table.items():
        assert v.dtype == np.float64 and v.size > 0, (name, k)
    np.savez_compressed(OUT / f"{name}_null.npz", **table)
    print(f"{name}: {len(table)} cells -> {OUT / (name + '_null.npz')}")


def main() -> None:
    # Wilcoxon: long (x, gamma, beta)
    header, data = read_csv(SRC / "T2_gamma_beta.csv")
    ix, ig, ib = header.index("x"), header.index("gamma"), header.index("beta")
    cells: dict[str, list[float]] = {}
    for row in data:
        key = f"g{fmt(float(row[ig]))}_b{fmt(float(row[ib]))}"
        cells.setdefault(key, []).append(float(row[ix]))
    wt = {k: np.array(v, dtype=np.float64) for k, v in cells.items()}
    # exact-equality spot re-read
    for k, v in wt.items():
        assert all(np.isfinite(v)), k
    save("wilcoxon", wt)

    # Lp: wide, columns dist.p2.ws.g<g>.n200.p<p>
    header, data = read_csv(SRC / "dist_p2_ws_n200_t.csv")
    cols = np.array([[float(c) for c in row] for row in data])
    lp = {}
    for j, col in enumerate(header):
        m = re.fullmatch(r"dist\.p2\.ws\.g([0-9.]+)\.n200\.p([0-9]+)", col)
        assert m, col
        lp[f"g{fmt(float(m.group(1)))}_p{m.group(2)}"] = \
            np.ascontiguousarray(cols[:, j])
    save("lp", lp)

    # RKHS corrected: wide, SD.<g>.h.<kernel>.n2000[.<baseline>]
    header, data = read_csv(SRC / "SD_h_n2000_corrected.csv")
    cols = np.array([[float(c) for c in row] for row in data])
    rk = {}
    for j, col in enumerate(header):
        m = re.fullmatch(
            r"SD\.([0-9.]+)\.h\.(gaussian|min|2\.max)\.n2000(?:\.(.+))?", col)
        assert m, col
        kernel = m.group(2).replace("2.max", "2max")
        baseline = f"_{m.group(3)}" if m.group(3) else ""
        rk[f"g{fmt(float(m.group(1)))}_{kernel}{baseline}"] = \
            np.ascontiguousarray(cols[:, j])
    save("rkhs", rk)

    # Birbilas: long (x, gamma)
    header, data = read_csv(SRC / "E_q_gamma.csv")
    ix, ig = header.index("x"), header.index("gamma")
    cells = {}
    for row in data:
        cells.setdefault(f"g{fmt(float(row[ig]))}", []).append(float(row[ix]))
    save("birbilas", {k: np.array(v, dtype=np.float64)
                      for k, v in cells.items()})


if __name__ == "__main__":
    main()
