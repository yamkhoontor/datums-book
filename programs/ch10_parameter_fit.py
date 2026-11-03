"""ch10_parameter_fit.py - Deriving transformation parameters

Companion program to "Datums, Coordinate Systems and Map Projections: A Practical
Guide for Surveyors and Engineers, with SVY21 Worked Examples" by Yam Khoon Tor.

Copyright (c) 2026 Yam Khoon Tor. All rights reserved. Readers who have bought the
book may run, copy and modify this program for personal study and their own
professional work; it may not be republished, distributed or sold without the
author's written permission. Provided as a teaching aid, without warranty: verify
results independently before use in real work.

Contents: Listings 10.1, 10.2
Run with:  python ch10_parameter_fit.py
Requires:  Python 3.9 or later and NumPy
"""


# ============================================================================
# Listing 10.1  Least squares core and two-dimensional conformal fit
# ============================================================================
import math
import numpy as np

ARCSEC = math.radians(1.0 / 3600.0)

# Common points: name, site e, site n, SVY21 E, SVY21 N (metres). Illustrative.
SITE_POINTS = [
    ("C1", 1000.000, 1000.000, 29170.138, 32443.589),
    ("C2", 1420.350, 1065.200, 29566.420, 32598.231),
    ("C3", 1510.800, 1390.450, 29584.336, 32935.352),
    ("C4", 1180.250, 1525.900, 29232.288, 32996.048),
    ("C5",  980.600, 1310.750, 29083.934, 32742.784),
]
CHECK_2D = ("K1", 1265.400, 1240.300, 29377.246, 32735.645)


def least_squares(A, L):
    """Solve A x = L by least squares.

    Returns the parameters, the residuals, s0 and the covariance matrix.
    """
    N = A.T @ A
    x = np.linalg.solve(N, A.T @ L)
    v = A @ x - L
    dof = len(L) - len(x)
    s0 = math.sqrt(float(v @ v) / dof)
    Cx = s0 ** 2 * np.linalg.inv(N)
    return x, v, s0, Cx


def fit_conformal_2d(points):
    """Four-parameter fit: E = a e' - b n' + c,  N = b e' + a n' + d,
    where e', n' are site coordinates reduced to their centroid."""
    e = np.array([p[1] for p in points])
    n = np.array([p[2] for p in points])
    e0, n0 = e.mean(), n.mean()
    A, L = [], []
    for p, ei, ni in zip(points, e - e0, n - n0):
        A.append([ei, -ni, 1.0, 0.0]); L.append(p[3])
        A.append([ni, ei, 0.0, 1.0]); L.append(p[4])
    x, v, s0, Cx = least_squares(np.array(A), np.array(L))
    return x, v, s0, Cx, (e0, n0)


def apply_conformal_2d(x, centroid, e, n):
    """Transform one site point to the target grid."""
    a, b, c, d = x
    de, dn = e - centroid[0], n - centroid[1]
    return a * de - b * dn + c, b * de + a * dn + d


if __name__ == "__main__":
    x, v, s0, Cx, cen = fit_conformal_2d(SITE_POINTS)
    a, b, c, d = x
    scale = math.hypot(a, b)
    theta = math.degrees(math.atan2(b, a))
    print(f"centroid e0 = {cen[0]:.3f}  n0 = {cen[1]:.3f}")
    print(f"a = {a:.8f}  b = {b:.8f}  c = {c:.4f}  d = {d:.4f}")
    print(f"scale = {scale:.8f} ({(scale - 1) * 1e6:+.1f} ppm)  "
          f"rotation = {theta:.5f} deg")
    for i, p in enumerate(SITE_POINTS):
        print(f"{p[0]}  vE = {v[2 * i] * 1000:+.1f} mm  "
              f"vN = {v[2 * i + 1] * 1000:+.1f} mm")
    print(f"s0 = {s0 * 1000:.1f} mm")
    E, N = apply_conformal_2d(x, cen, CHECK_2D[1], CHECK_2D[2])
    print(f"check {CHECK_2D[0]}: dE = {(E - CHECK_2D[3]) * 1000:+.1f} mm  "
          f"dN = {(N - CHECK_2D[4]) * 1000:+.1f} mm")


# ============================================================================
# Listing 10.2  Seven-parameter fit in either form
# ============================================================================
# Common points: name, source X Y Z, target X Y Z (metres). Illustrative.
POINTS_3D = [
    ("A", -1510204.669, 6195108.523, 143734.927,
     -1510190.226, 6194246.224, 143728.183),
    ("B", -1515519.310, 6193416.229, 160316.815,
     -1515504.883, 6192553.902, 160310.032),
    ("C", -1542589.194, 6186987.328, 149262.091,
     -1542574.317, 6186123.253, 149255.102),
    ("D", -1526445.287, 6191243.328, 138207.319,
     -1526430.522, 6190380.020, 138200.428),
    ("E", -1504606.825, 6195637.499, 176897.706,
     -1504592.720, 6194775.995, 176890.984),
    ("F", -1553254.927, 6183767.066, 171370.725,
     -1553240.034, 6182902.602, 171363.664),
    ("G", -1531511.918, 6188552.533, 193477.334,
     -1531497.553, 6187689.623, 193470.401),
]
CHECK_3D = ("H", -1529585.343, 6190032.629, 157000.340,
            -1529570.682, 6189169.401, 156993.429)


def fit_helmert7(points, pivot="centroid"):
    """Seven-parameter fit, position vector convention.

    pivot="centroid" gives the Molodensky-Badekas form;
    pivot="origin" gives the Bursa-Wolf form.
    Unknowns: Tx, Ty, Tz (m), Rx, Ry, Rz (rad), s (pure number).
    """
    S = np.array([p[1:4] for p in points])
    T = np.array([p[4:7] for p in points])
    X0 = S.mean(axis=0) if pivot == "centroid" else np.zeros(3)
    A, L = [], []
    for (x, y, z), d in zip(S - X0, T - S):
        A.append([1, 0, 0, 0.0, z, -y, x]); L.append(d[0])
        A.append([0, 1, 0, -z, 0.0, x, y]); L.append(d[1])
        A.append([0, 0, 1, y, -x, 0.0, z]); L.append(d[2])
    x, v, s0, Cx = least_squares(np.array(A, dtype=float), np.array(L))
    return x, v, s0, Cx, X0


def apply_helmert7_fit(x, X0, X):
    """Apply fitted parameters to one source point."""
    Tx, Ty, Tz, rx, ry, rz, s = x
    dx, dy, dz = np.array(X) - X0
    return (X[0] + Tx + s * dx - rz * dy + ry * dz,
            X[1] + Ty + s * dy + rz * dx - rx * dz,
            X[2] + Tz + s * dz - ry * dx + rx * dy)


if __name__ == "__main__":
    np.set_printoptions(suppress=True)
    for pivot in ("centroid", "origin"):
        x, v, s0, Cx, X0 = fit_helmert7(POINTS_3D, pivot)
        sd = np.sqrt(np.diag(Cx))
        corr = Cx / np.outer(sd, sd)
        print(f"--- pivot = {pivot}   s0 = {s0 * 1000:.1f} mm")
        print("T (m)      ", np.round(x[:3], 3), " sd ", np.round(sd[:3], 3))
        print("R (arcsec) ", np.round(x[3:6] / ARCSEC, 3),
              " sd ", np.round(sd[3:6] / ARCSEC, 3))
        print(f"s (ppm)     {x[6] * 1e6:.2f}  sd {sd[6] * 1e6:.2f}")
        print(f"correlation Tx-Rz = {corr[0, 5]:+.2f}")
        H = apply_helmert7_fit(x, X0, CHECK_3D[1:4])
        diff = (np.array(H) - np.array(CHECK_3D[4:7])) * 1000
        print("check H (mm):", np.round(diff, 1))
