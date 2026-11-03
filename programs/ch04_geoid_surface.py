"""ch04_geoid_surface.py - Heights and vertical datums

Companion program to "Datums, Coordinate Systems and Map Projections: A Practical
Guide for Surveyors and Engineers, with SVY21 Worked Examples" by Yam Khoon Tor.

Copyright (c) 2026 Yam Khoon Tor. All rights reserved. Readers who have bought the
book may run, copy and modify this program for personal study and their own
professional work; it may not be republished, distributed or sold without the
author's written permission. Provided as a teaching aid, without warranty: verify
results independently before use in real work.

Contents: Listings 4.1
Run with:  python ch04_geoid_surface.py
Requires:  Python 3.9 or later and NumPy
"""


# ============================================================================
# Listing 4.1  Geometric geoid surface by least squares
# ============================================================================
import numpy as np

# Benchmarks: name, northing (m), easting (m), ellipsoidal height h (m),
# levelled reduced level H (m). Illustrative values.
BENCHMARKS = [
    ("BM1", 30500.0, 22400.0, 16.918,  7.132),
    ("BM2", 36200.0, 24800.0, 25.451, 15.640),
    ("BM3", 41800.0, 23100.0, 33.777, 23.862),
    ("BM4", 43300.0, 29700.0, 21.342, 11.562),
    ("BM5", 38100.0, 31200.0, 48.809, 39.115),
    ("BM6", 32600.0, 30100.0, 15.160,  5.523),
    ("BM7", 35400.0, 36900.0, 12.371,  2.871),
    ("BM8", 40600.0, 37600.0, 30.033, 20.469),
]
N0, E0 = 38744.572, 28001.642        # reference point of the surface


def fit_plane(benchmarks):
    """Least squares fit of N = c0 + c1*dn + c2*de, with dn, de in km."""
    A, L = [], []
    for _, n, e, h, H in benchmarks:
        A.append([1.0, (n - N0) / 1000.0, (e - E0) / 1000.0])
        L.append(h - H)                        # observed undulation
    A, L = np.array(A), np.array(L)
    x = np.linalg.solve(A.T @ A, A.T @ L)      # normal equations
    v = A @ x - L                              # residuals
    dof = len(L) - len(x)                      # degrees of freedom
    s0 = float(np.sqrt(v @ v / dof))           # standard deviation of fit
    return x, v, s0


def undulation(x, n, e):
    """Interpolated undulation at grid position (n, e)."""
    return x[0] + x[1] * (n - N0) / 1000.0 + x[2] * (e - E0) / 1000.0


if __name__ == "__main__":
    x, v, s0 = fit_plane(BENCHMARKS)
    print(f"c0 = {x[0]:.4f} m  c1 = {x[1]:.5f} m/km  c2 = {x[2]:.5f} m/km")
    for (name, n, e, h, H), vi in zip(BENCHMARKS, v):
        print(f"{name}  observed N = {h - H:.3f}  residual = {vi * 1000:+.1f} mm")
    print(f"standard deviation of fit = {s0 * 1000:.1f} mm")

    n, e, h = 37000.0, 28500.0, 27.818         # new point
    Nu = undulation(x, n, e)
    print(f"new point: N = {Nu:.3f} m  H = {h - Nu:.3f} m")
