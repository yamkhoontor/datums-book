"""appendix_e_tests.py - test points with known answers (Appendix E)

Companion program to "Datums, Coordinate Systems and Map Projections: A Practical
Guide for Surveyors and Engineers, with SVY21 Worked Examples" by Yam Khoon Tor.

Copyright (c) 2026 Yam Khoon Tor. All rights reserved. Readers who have bought the
book may run, copy and modify this program for personal study and their own
professional work; it may not be republished, distributed or sold without the
author's written permission. Provided as a teaching aid, without warranty: verify
results independently before use in real work.

Contents: tests of the book's functions against the published values of Appendix E
Run with:  python appendix_e_tests.py
Requires:  Python 3.9 or later and NumPy
"""


# Functions from earlier chapters, exactly as printed in the book.


# ---------------- from Listing 2.1: Ellipsoid constants and radii of curvature ----------------
import math

# Defining parameters: semi-major axis a (m) and reciprocal flattening 1/f
ELLIPSOIDS = {
    "WGS84": (6378137.0, 298.257223563),
    "GRS80": (6378137.0, 298.257222101),
    "Everest 1830 (Modified)": (6377304.063, 300.8017),
}


def ellipsoid_constants(a, inv_f):
    """Derived constants of an ellipsoid from a and 1/f."""
    f = 1.0 / inv_f
    b = a * (1.0 - f)
    e2 = 2.0 * f - f * f          # first eccentricity squared
    ep2 = e2 / (1.0 - e2)         # second eccentricity squared
    n = f / (2.0 - f)             # third flattening
    return {"a": a, "f": f, "b": b, "e2": e2, "ep2": ep2, "n": n}


def dms_to_deg(d, m, s):
    """Degrees, minutes, seconds to decimal degrees."""
    sign = -1.0 if d < 0 else 1.0
    return sign * (abs(d) + m / 60.0 + s / 3600.0)


def radii_of_curvature(lat_deg, a, inv_f):
    """Meridian radius M, prime vertical radius nu, mean radius R."""
    c = ellipsoid_constants(a, inv_f)
    phi = math.radians(lat_deg)
    w2 = 1.0 - c["e2"] * math.sin(phi) ** 2
    nu = a / math.sqrt(w2)
    M = a * (1.0 - c["e2"]) / w2 ** 1.5
    R = math.sqrt(M * nu)
    return M, nu, R


# ---------------- from Listing F.1: Utility functions ----------------
import math


def deg_to_dms(x, places=4):
    """Format decimal degrees as degrees, minutes and seconds."""
    sign = "-" if x < 0 else ""
    x = abs(x)
    t = round(x * 3600, places)              # total seconds, rounded once
    d = int(t // 3600)
    m = int((t - d * 3600) // 60)
    s = t - d * 3600 - m * 60
    return f"{sign}{d}° {m:02d}' {s:0{3 + places}.{places}f}\""


def vincenty_inverse(lat1, lon1, lat2, lon2, a, inv_f):
    """Geodesic distance (m) and forward azimuths at both ends (degrees)."""
    c = ellipsoid_constants(a, inv_f)
    f, b = c["f"], c["b"]
    U1 = math.atan((1 - f) * math.tan(math.radians(lat1)))
    U2 = math.atan((1 - f) * math.tan(math.radians(lat2)))
    L = math.radians(lon2 - lon1)
    lam = L
    for _ in range(100):
        sl, cl = math.sin(lam), math.cos(lam)
        ss = math.hypot(math.cos(U2) * sl,
                        math.cos(U1) * math.sin(U2)
                        - math.sin(U1) * math.cos(U2) * cl)
        cs = math.sin(U1) * math.sin(U2) + math.cos(U1) * math.cos(U2) * cl
        sig = math.atan2(ss, cs)
        sa = math.cos(U1) * math.cos(U2) * sl / ss
        c2a = 1 - sa * sa
        c2sm = cs - 2 * math.sin(U1) * math.sin(U2) / c2a if c2a else 0.0
        C = f / 16 * c2a * (4 + f * (4 - 3 * c2a))
        lam_new = L + (1 - C) * f * sa * (
            sig + C * ss * (c2sm + C * cs * (-1 + 2 * c2sm ** 2)))
        converged = abs(lam_new - lam) < 1e-13
        lam = lam_new
        if converged:
            break
    u2 = c2a * (a * a - b * b) / (b * b)
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    ds = B * ss * (c2sm + B / 4 * (cs * (-1 + 2 * c2sm ** 2)
                                    - B / 6 * c2sm * (-3 + 4 * ss ** 2)
                                    * (-3 + 4 * c2sm ** 2)))
    s = b * A * (sig - ds)
    sl, cl = math.sin(lam), math.cos(lam)
    az1 = math.degrees(math.atan2(
        math.cos(U2) * sl,
        math.cos(U1) * math.sin(U2) - math.sin(U1) * math.cos(U2) * cl)) % 360
    az2 = math.degrees(math.atan2(
        math.cos(U1) * sl,
        -math.sin(U1) * math.cos(U2) + math.cos(U1) * math.sin(U2) * cl)) % 360
    return s, az1, az2


# ---------------- from Listing 3.1: Conversion between geodetic and Cartesian coordinates ----------------
import math


def geodetic_to_ecef(lat_deg, lon_deg, h, a, inv_f):
    """Latitude, longitude (degrees) and height (m) to X, Y, Z (m)."""
    e2 = ellipsoid_constants(a, inv_f)["e2"]
    phi, lam = math.radians(lat_deg), math.radians(lon_deg)
    nu = a / math.sqrt(1.0 - e2 * math.sin(phi) ** 2)
    X = (nu + h) * math.cos(phi) * math.cos(lam)
    Y = (nu + h) * math.cos(phi) * math.sin(lam)
    Z = (nu * (1.0 - e2) + h) * math.sin(phi)
    return X, Y, Z


def ecef_to_geodetic(X, Y, Z, a, inv_f, tol=1e-12):
    """X, Y, Z (m) to latitude, longitude (degrees) and height (m)."""
    e2 = ellipsoid_constants(a, inv_f)["e2"]
    lam = math.atan2(Y, X)
    p = math.hypot(X, Y)
    phi = math.atan2(Z, p * (1.0 - e2))          # first approximation
    for _ in range(20):
        nu = a / math.sqrt(1.0 - e2 * math.sin(phi) ** 2)
        phi_new = math.atan2(Z + e2 * nu * math.sin(phi), p)
        converged = abs(phi_new - phi) < tol
        phi = phi_new
        if converged:
            break
    nu = a / math.sqrt(1.0 - e2 * math.sin(phi) ** 2)
    if abs(math.cos(phi)) > 1e-3:
        h = p / math.cos(phi) - nu
    else:                                        # close to a pole
        h = Z / math.sin(phi) - nu * (1.0 - e2)
    return math.degrees(phi), math.degrees(lam), h


# ---------------- from Listing 9.1: Seven-parameter transformation, rigorous inverse and Molodensky-Badekas translations ----------------
import math
import numpy as np

ARCSEC = math.radians(1.0 / 3600.0)          # one second of arc in radians


def rotation_matrix(R_arcsec, convention="position_vector"):
    """Linearised rotation matrix from rotations in arc-seconds."""
    rx, ry, rz = (r * ARCSEC for r in R_arcsec)
    if convention == "coordinate_frame":
        rx, ry, rz = -rx, -ry, -rz
    elif convention != "position_vector":
        raise ValueError("convention must be position_vector or coordinate_frame")
    return np.array([[1.0, -rz, ry],
                     [rz, 1.0, -rx],
                     [-ry, rx, 1.0]])


def helmert7(X, T, R_arcsec, s_ppm, convention="position_vector"):
    """Seven-parameter transformation of Cartesian coordinates (Bursa-Wolf)."""
    m = 1.0 + s_ppm * 1.0e-6
    return np.array(T) + m * rotation_matrix(R_arcsec, convention) @ np.array(X)


def helmert7_inverse(X, T, R_arcsec, s_ppm, convention="position_vector"):
    """Rigorous inverse: recover source coordinates from target coordinates."""
    m = 1.0 + s_ppm * 1.0e-6
    A = m * rotation_matrix(R_arcsec, convention)
    return np.linalg.solve(A, np.array(X) - np.array(T))


def bw_to_mb(T, R_arcsec, s_ppm, X0, convention="position_vector"):
    """Translations of the Molodensky-Badekas form for pivot point X0."""
    m = 1.0 + s_ppm * 1.0e-6
    X0 = np.array(X0)
    return np.array(T) - X0 + m * rotation_matrix(R_arcsec, convention) @ X0


def transform_geodetic(lat, lon, h, src, tgt, T, R_arcsec, s_ppm,
                       convention="position_vector"):
    """Geodetic on the source datum to geodetic on the target datum."""
    X = geodetic_to_ecef(lat, lon, h, *ELLIPSOIDS[src])
    X2 = helmert7(X, T, R_arcsec, s_ppm, convention)
    return ecef_to_geodetic(*X2, *ELLIPSOIDS[tgt])


# WGS84 (PMGSN94) to Kertau (MRT68), after JUPEM (2009a)
T_MRT = (379.776, -775.384, 86.609)        # metres
R_MRT = (-2.59674, -2.10213, 12.11377)     # arc-seconds, position vector
S_MRT = 0.0                                # ppm


# ---------------- from Listing 13.1: Transverse Mercator forward, inverse, scale factor and convergence (Krüger series) ----------------
import math


def _kruger(a, inv_f):
    """Constants of the Kruger series for an ellipsoid."""
    f = 1.0 / inv_f
    n = f / (2.0 - f)
    B = a / (1.0 + n) * (1.0 + n ** 2 / 4.0 + n ** 4 / 64.0 + n ** 6 / 256.0)
    alpha = (n / 2 - 2 * n ** 2 / 3 + 5 * n ** 3 / 16 + 41 * n ** 4 / 180,
             13 * n ** 2 / 48 - 3 * n ** 3 / 5 + 557 * n ** 4 / 1440,
             61 * n ** 3 / 240 - 103 * n ** 4 / 140,
             49561 * n ** 4 / 161280)
    beta = (n / 2 - 2 * n ** 2 / 3 + 37 * n ** 3 / 96 - n ** 4 / 360,
            n ** 2 / 48 + n ** 3 / 15 - 437 * n ** 4 / 1440,
            17 * n ** 3 / 480 - 37 * n ** 4 / 840,
            4397 * n ** 4 / 161280)
    e = math.sqrt(2.0 * f - f * f)
    return B, alpha, beta, e


def _xi_eta(lat_deg, dlon_rad, e):
    """Coordinates xi', eta' on the conformal sphere, and the auxiliary t."""
    phi = math.radians(lat_deg)
    t = math.sinh(math.atanh(math.sin(phi)) - e * math.atanh(e * math.sin(phi)))
    xi = math.atan2(t, math.cos(dlon_rad))
    eta = math.atanh(math.sin(dlon_rad) / math.sqrt(1.0 + t * t))
    return xi, eta, t


def tm_forward(lat_deg, lon_deg, lat0, lon0, k0, FN, FE, a, inv_f):
    """Latitude, longitude (degrees) to northing, easting (m)."""
    B, alpha, beta, e = _kruger(a, inv_f)
    xi, eta, t = _xi_eta(lat_deg, math.radians(lon_deg - lon0), e)
    X = xi + sum(alpha[j] * math.sin(2 * (j + 1) * xi)
                 * math.cosh(2 * (j + 1) * eta) for j in range(4))
    Y = eta + sum(alpha[j] * math.cos(2 * (j + 1) * xi)
                  * math.sinh(2 * (j + 1) * eta) for j in range(4))
    xi0, _, _ = _xi_eta(lat0, 0.0, e)
    X0 = xi0 + sum(alpha[j] * math.sin(2 * (j + 1) * xi0) for j in range(4))
    return FN + k0 * B * (X - X0), FE + k0 * B * Y


def tm_inverse(N, E, lat0, lon0, k0, FN, FE, a, inv_f):
    """Northing, easting (m) to latitude, longitude (degrees)."""
    B, alpha, beta, e = _kruger(a, inv_f)
    xi0, _, _ = _xi_eta(lat0, 0.0, e)
    X0 = xi0 + sum(alpha[j] * math.sin(2 * (j + 1) * xi0) for j in range(4))
    X = (N - FN) / (k0 * B) + X0
    Y = (E - FE) / (k0 * B)
    xi = X - sum(beta[j] * math.sin(2 * (j + 1) * X)
                 * math.cosh(2 * (j + 1) * Y) for j in range(4))
    eta = Y - sum(beta[j] * math.cos(2 * (j + 1) * X)
                  * math.sinh(2 * (j + 1) * Y) for j in range(4))
    chi = math.asin(math.sin(xi) / math.cosh(eta))        # conformal latitude
    dlon = math.atan2(math.sinh(eta), math.cos(xi))
    phi = chi                                     # iterate to geodetic latitude
    for _ in range(10):
        phi_new = math.asin(math.tanh(math.atanh(math.sin(chi))
                                      + e * math.atanh(e * math.sin(phi))))
        converged = abs(phi_new - phi) < 1e-14
        phi = phi_new
        if converged:
            break
    return math.degrees(phi), lon0 + math.degrees(dlon)


def tm_scale_convergence(lat_deg, lon_deg, lon0, k0, a, inv_f):
    """Point scale factor and grid convergence (degrees)."""
    B, alpha, beta, e = _kruger(a, inv_f)
    f = 1.0 / inv_f
    n = f / (2.0 - f)
    dl = math.radians(lon_deg - lon0)
    xi, eta, t = _xi_eta(lat_deg, dl, e)
    sigma = 1.0 + sum(2 * (j + 1) * alpha[j] * math.cos(2 * (j + 1) * xi)
                      * math.cosh(2 * (j + 1) * eta) for j in range(4))
    tau = sum(2 * (j + 1) * alpha[j] * math.sin(2 * (j + 1) * xi)
              * math.sinh(2 * (j + 1) * eta) for j in range(4))
    phi = math.radians(lat_deg)
    q = (1.0 - n) / (1.0 + n) * math.tan(phi)
    k = k0 * B / a * math.sqrt((1.0 + q ** 2) * (sigma ** 2 + tau ** 2)
                               / (t * t + math.cos(dl) ** 2))
    gamma = math.atan2(tau * math.sqrt(1.0 + t * t) + sigma * t * math.tan(dl),
                       sigma * math.sqrt(1.0 + t * t) - tau * t * math.tan(dl))
    return k, math.degrees(gamma)


# ---------------- from Listing 14.1: Cassini-Soldner projection ----------------
import math


def meridian_arc(lat_deg, a, inv_f):
    """Length of the meridian from the equator to the given latitude (m)."""
    e2 = ellipsoid_constants(a, inv_f)["e2"]
    p = math.radians(lat_deg)
    e4, e6 = e2 * e2, e2 ** 3
    return a * ((1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256) * p
                - (3 * e2 / 8 + 3 * e4 / 32 + 45 * e6 / 1024) * math.sin(2 * p)
                + (15 * e4 / 256 + 45 * e6 / 1024) * math.sin(4 * p)
                - (35 * e6 / 3072) * math.sin(6 * p))


def footpoint_latitude(m1, a, inv_f):
    """Latitude (radians) at which the meridian arc from the equator is m1."""
    e2 = ellipsoid_constants(a, inv_f)["e2"]
    e4, e6 = e2 * e2, e2 ** 3
    mu = m1 / (a * (1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256))
    e1 = (1 - math.sqrt(1 - e2)) / (1 + math.sqrt(1 - e2))
    return (mu + (3 * e1 / 2 - 27 * e1 ** 3 / 32) * math.sin(2 * mu)
            + (21 * e1 ** 2 / 16 - 55 * e1 ** 4 / 32) * math.sin(4 * mu)
            + (151 * e1 ** 3 / 96) * math.sin(6 * mu)
            + (1097 * e1 ** 4 / 512) * math.sin(8 * mu))


def cassini_forward(lat_deg, lon_deg, lat0, lon0, FN, FE, a, inv_f):
    """Cassini-Soldner forward: latitude, longitude to northing, easting."""
    c = ellipsoid_constants(a, inv_f)
    e2, ep2 = c["e2"], c["ep2"]
    p = math.radians(lat_deg)
    nu = a / math.sqrt(1 - e2 * math.sin(p) ** 2)
    A = math.radians(lon_deg - lon0) * math.cos(p)
    T = math.tan(p) ** 2
    C = ep2 * math.cos(p) ** 2
    E = FE + nu * (A - T * A ** 3 / 6 - (8 - T + 8 * C) * T * A ** 5 / 120)
    N = (FN + meridian_arc(lat_deg, a, inv_f) - meridian_arc(lat0, a, inv_f)
         + nu * math.tan(p) * (A ** 2 / 2 + (5 - T + 6 * C) * A ** 4 / 24))
    return N, E


def cassini_inverse(N, E, lat0, lon0, FN, FE, a, inv_f):
    """Cassini-Soldner inverse: northing, easting to latitude, longitude."""
    e2 = ellipsoid_constants(a, inv_f)["e2"]
    p1 = footpoint_latitude(meridian_arc(lat0, a, inv_f) + (N - FN), a, inv_f)
    nu1 = a / math.sqrt(1 - e2 * math.sin(p1) ** 2)
    M1 = a * (1 - e2) / (1 - e2 * math.sin(p1) ** 2) ** 1.5
    T1 = math.tan(p1) ** 2
    D = (E - FE) / nu1
    phi = p1 - (nu1 * math.tan(p1) / M1) * (D ** 2 / 2 - (1 + 3 * T1) * D ** 4 / 24)
    lam = math.radians(lon0) + (D - T1 * D ** 3 / 3
                                + (1 + 3 * T1) * T1 * D ** 5 / 15) / math.cos(p1)
    return math.degrees(phi), math.degrees(lam)


# ---------------- from Listing 14.2: Hotine oblique Mercator projection, variant A ----------------
def hom_setup(latc, lonc, azc, gammac, kc, FE, FN, a, inv_f):
    """Derived constants of the oblique Mercator from its defining parameters.

    latc, lonc: projection centre; azc: azimuth of the central line;
    gammac: angle from rectified to skew grid (all in degrees);
    kc: scale factor on the central line; FE, FN: false easting and northing.
    """
    e2 = ellipsoid_constants(a, inv_f)["e2"]
    e = math.sqrt(e2)
    pc, lc = math.radians(latc), math.radians(lonc)
    ac, gc = math.radians(azc), math.radians(gammac)
    B = math.sqrt(1 + e2 * math.cos(pc) ** 4 / (1 - e2))
    A = a * B * kc * math.sqrt(1 - e2) / (1 - e2 * math.sin(pc) ** 2)
    t0 = math.tan(math.pi / 4 - pc / 2) / (
        (1 - e * math.sin(pc)) / (1 + e * math.sin(pc))) ** (e / 2)
    D = B * math.sqrt(1 - e2) / (math.cos(pc)
                                 * math.sqrt(1 - e2 * math.sin(pc) ** 2))
    D2 = max(D * D, 1.0)
    F = D + math.sqrt(D2 - 1) * (1 if latc >= 0 else -1)
    H = F * t0 ** B
    G = (F - 1 / F) / 2
    g0 = math.asin(math.sin(ac) / D)
    l0 = lc - math.asin(G * math.tan(g0)) / B
    return dict(e=e, A=A, B=B, H=H, g0=g0, l0=l0, gc=gc, FE=FE, FN=FN)


def hom_forward(lat_deg, lon_deg, P):
    """Oblique Mercator forward: latitude, longitude to northing, easting."""
    e, A, B, H, g0, l0, gc = (P[k] for k in
                              ("e", "A", "B", "H", "g0", "l0", "gc"))
    p, l = math.radians(lat_deg), math.radians(lon_deg)
    t = math.tan(math.pi / 4 - p / 2) / (
        (1 - e * math.sin(p)) / (1 + e * math.sin(p))) ** (e / 2)
    Q = H / t ** B
    S, T = (Q - 1 / Q) / 2, (Q + 1 / Q) / 2
    V = math.sin(B * (l - l0))
    U = (-V * math.cos(g0) + S * math.sin(g0)) / T
    v = A * math.log((1 - U) / (1 + U)) / (2 * B)        # across the central line
    u = A * math.atan2(S * math.cos(g0) + V * math.sin(g0),
                       math.cos(B * (l - l0))) / B      # along the central line
    E = v * math.cos(gc) + u * math.sin(gc) + P["FE"]    # rectify
    N = u * math.cos(gc) - v * math.sin(gc) + P["FN"]
    return N, E


def hom_inverse(N, E, P):
    """Oblique Mercator inverse: northing, easting to latitude, longitude."""
    e, A, B, H, g0, l0, gc = (P[k] for k in
                              ("e", "A", "B", "H", "g0", "l0", "gc"))
    dE, dN = E - P["FE"], N - P["FN"]
    v = dE * math.cos(gc) - dN * math.sin(gc)
    u = dN * math.cos(gc) + dE * math.sin(gc)
    Q = math.exp(-B * v / A)
    S, T = (Q - 1 / Q) / 2, (Q + 1 / Q) / 2
    V = math.sin(B * u / A)
    U = (V * math.cos(g0) + S * math.sin(g0)) / T
    t = (H / math.sqrt((1 + U) / (1 - U))) ** (1 / B)
    chi = math.pi / 2 - 2 * math.atan(t)                 # conformal latitude
    e2 = e * e
    c2 = e2 / 2 + 5 * e2 ** 2 / 24 + e2 ** 3 / 12 + 13 * e2 ** 4 / 360
    c4 = 7 * e2 ** 2 / 48 + 29 * e2 ** 3 / 240 + 811 * e2 ** 4 / 11520
    c6 = 7 * e2 ** 3 / 120 + 81 * e2 ** 4 / 1120
    c8 = 4279 * e2 ** 4 / 161280
    phi = (chi + c2 * math.sin(2 * chi) + c4 * math.sin(4 * chi)
           + c6 * math.sin(6 * chi) + c8 * math.sin(8 * chi))
    lam = l0 - math.atan2(S * math.cos(g0) - V * math.sin(g0),
                          math.cos(B * u / A)) / B
    return math.degrees(phi), math.degrees(lam)




# ============================================================================
# E.1 SiReNT reference stations: grid coordinates (Table E.1)
# ============================================================================
SIRENT = {   # latitude (d, m, s), longitude (d, m, s), h, published N, E, published X, Y, Z
    "SNYU": ((1, 20, 44.84810), (103, 40, 47.91672), 76.246, 36436.799, 10934.735, -1508025.227, 6195576.969, 148798.263),
    "SRPT": ((1, 26, 37.41926), (103, 47, 3.88829), 59.411, 47265.636, 22557.618, -1519248.959, 6192545.559, 159624.074),
    "SNPT": ((1, 22, 44.82064), (103, 50, 55.51059), 55.486, 40121.255, 29717.651, -1526243.018, 6191002.686, 152481.753),
    "SLYG": ((1, 22, 21.44707), (103, 58, 18.00928), 51.158, 39403.770, 43396.751, -1539524.135, 6187726.731, 151763.925),
    "SNSC": ((1, 18, 49.12095), (103, 57, 42.01188), 15.253, 32882.042, 42284.293, -1538472.450, 6188108.833, 145243.192),
    "SNUS": ((1, 17, 32.68755), (103, 46, 31.83287), 63.417, 30534.079, 21566.307, -1518382.637, 6193173.459, 142897.210),
    "SSMK": ((1, 12, 37.48752), (103, 46, 31.49934), 25.419, 21466.939, 21555.796, -1518410.730, 6193331.346, 133831.345),
    "SSTS": ((1, 14, 23.49077), (103, 38, 52.29482), 13.570, 24723.506, 7359.574, -1504599.629, 6196617.251, 137086.270),
    "SMS1": ((1, 20, 25.60771), (103, 53, 16.40533), 37.860, 35845.357, 34073.254, -1530491.719, 6190040.099, 148206.547),
}
SVY21_TM = (dms_to_deg(1, 22, 0), dms_to_deg(103, 50, 0), 1.0, 38744.572, 28001.642)
RESULTS = []


def report(name, worst, tolerance, unit="mm"):
    """Print one line and remember whether the test passed."""
    ok = worst <= tolerance
    RESULTS.append((name, ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}: largest difference {worst:.3g} {unit} "
          f"(tolerance {tolerance:g} {unit})")


a, inv_f = ELLIPSOIDS["WGS84"]
worst = 0.0
for name, (lat_dms, lon_dms, h, N_pub, E_pub, X_pub, Y_pub, Z_pub) in SIRENT.items():
    N, E = tm_forward(dms_to_deg(*lat_dms), dms_to_deg(*lon_dms), *SVY21_TM, a, inv_f)
    worst = max(worst, abs(N - N_pub), abs(E - E_pub))
report("geodetic to SVY21 grid, 9 stations", worst * 1000, 1.0)

worst = 0.0
for name, (lat_dms, lon_dms, h, N_pub, E_pub, X_pub, Y_pub, Z_pub) in SIRENT.items():
    lat, lon = tm_inverse(N_pub, E_pub, *SVY21_TM, a, inv_f)
    worst = max(worst, abs(lat - dms_to_deg(*lat_dms)) * 3600, abs(lon - dms_to_deg(*lon_dms)) * 3600)
report("SVY21 grid to geodetic, 9 stations", worst, 0.00005, "arc-seconds")


# ============================================================================
# E.1 SiReNT reference stations: Cartesian coordinates (Table E.2)
# ============================================================================
worst = 0.0
for name, (lat_dms, lon_dms, h, N_pub, E_pub, X_pub, Y_pub, Z_pub) in SIRENT.items():
    X, Y, Z = geodetic_to_ecef(dms_to_deg(*lat_dms), dms_to_deg(*lon_dms), h, a, inv_f)
    worst = max(worst, abs(X - X_pub), abs(Y - Y_pub), abs(Z - Z_pub))
report("geodetic to Cartesian, 9 stations", worst * 1000, 1.0)

worst = 0.0
for name, (lat_dms, lon_dms, h, N_pub, E_pub, X_pub, Y_pub, Z_pub) in SIRENT.items():
    lat, lon, h2 = ecef_to_geodetic(X_pub, Y_pub, Z_pub, a, inv_f)
    worst = max(worst, abs(h2 - h))
report("Cartesian to geodetic (height), 9 stations", worst * 1000, 1.0)


# ============================================================================
# E.2 Worked examples published by JUPEM
# ============================================================================
# Table E.3: WGS84 (PMGSN94) to Kertau (MRT68) with the parameters of Table B.5
lat, lon, h = transform_geodetic(dms_to_deg(1, 51, 27.08807), dms_to_deg(102, 56, 31.72178), 7.118,
                                 "WGS84", "Everest 1830 (Modified)", T_MRT, R_MRT, S_MRT)
worst = max(abs(lat - dms_to_deg(1, 51, 27.38567)), abs(lon - dms_to_deg(102, 56, 37.52660))) * 3600
report("PMGSN94 to MRT68, latitude and longitude", worst, 0.00005, "arc-seconds")
report("PMGSN94 to MRT68, height", abs(h - 2.338) * 1000, 20.0)

# Table E.4: Kedah and Perlis Cassini grid, GDM2000
a_g, f_g = ELLIPSOIDS["GRS80"]
lat, lon = dms_to_deg(6, 8, 22.98892), dms_to_deg(100, 23, 6.57684)
origin = (dms_to_deg(5, 57, 52.82155), dms_to_deg(100, 38, 10.93860))
N, E = cassini_forward(lat, lon, *origin, 0.0, 0.0, a_g, f_g)
report("Cassini-Soldner, Kedah and Perlis", max(abs(N - 19364.3195), abs(E + 27805.4063)) * 1000, 0.5)

# Table E.4: Peninsula RSO grid, GDM2000, with the parameters of Table B.3
MRSO = hom_setup(4.0, 102.25, dms_to_deg(323, 1, 32.8458), dms_to_deg(323, 7, 48.3685),
                 0.99984, 804671.0, 0.0, a_g, f_g)
N, E = hom_forward(lat, lon, MRSO)
report("Peninsula RSO, easting", abs(E - 266843.898) * 1000, 2.0)
print(f"      Peninsula RSO, northing: differs from the published value by "
      f"{(N - 679690.914) * 1000:+.0f} mm. This is the known discrepancy "
      f"discussed in Example 14.3; it is reported, not tested.")

# Cartesian coordinates on GRS80
X, Y, Z = geodetic_to_ecef(dms_to_deg(6, 27, 0.56909), dms_to_deg(100, 16, 47.05076), 18.078, a_g, f_g)
report("geodetic to Cartesian on GRS80",
       max(abs(X + 1131051.654), abs(Y - 6236311.800), abs(Z - 711748.112)) * 1000, 1.0)


# ============================================================================
# E.3 Values computed in this book (Table E.6 and the list after it)
# ============================================================================
# Bridge marks S and J of Table E.6
lat, lon, h = dms_to_deg(1, 27, 10.00803), dms_to_deg(103, 46, 9.97011), 12.014
latJ, lonJ = dms_to_deg(1, 27, 55.0), dms_to_deg(103, 46, 25.0)
X, Y, Z = geodetic_to_ecef(lat, lon, h, a, inv_f)
report("S: Cartesian",
       max(abs(X + 1517612.8594), abs(Y - 6192871.8908), abs(Z - 160623.5351)) * 1000, 0.1)
N, E = tm_forward(lat, lon, *SVY21_TM, a, inv_f)
report("S: SVY21 grid", max(abs(N - 48266.651), abs(E - 20890.909)) * 1000, 0.5)
N, E = tm_forward(lat, lon, 0.0, 105.0, 0.9996, 0.0, 500000.0, a, inv_f)
report("S: UTM zone 48 North", max(abs(N - 160613.463), abs(E - 363102.191)) * 1000, 0.5)
N, E = tm_forward(latJ, lonJ, *SVY21_TM, a, inv_f)
report("J: SVY21 grid", max(abs(N - 49648.580), abs(E - 21355.553)) * 1000, 0.5)

# Geodesic from S to J
s, az1, az2 = vincenty_inverse(lat, lon, latJ, lonJ, a, inv_f)
report("geodesic S to J", abs(s - 1457.950) * 1000, 0.5)

# S on Kertau (MRT68) by the six parameters of Table B.5
lat_k, lon_k, h_k = transform_geodetic(lat, lon, h, "WGS84", "Everest 1830 (Modified)",
                                       T_MRT, R_MRT, S_MRT)
worst = max(abs(lat_k - dms_to_deg(1, 27, 10.17821)),
            abs(lon_k - dms_to_deg(103, 46, 16.15825))) * 3600
report("S on Kertau (MRT68)", worst, 0.00001, "arc-seconds")


# ============================================================================
# Summary
# ============================================================================
passed = sum(ok for name, ok in RESULTS)
print(f"{passed} of {len(RESULTS)} tests passed.")
if passed == len(RESULTS):
    print("The functions reproduce the published and the book values.")
else:
    print("Not passed:", ", ".join(name for name, ok in RESULTS if not ok))
