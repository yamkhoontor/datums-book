"""ch13_svy21_tm.py - Transverse Mercator and the SVY21 grid

Companion program to "Datums, Coordinate Systems and Map Projections: A Practical
Guide for Surveyors and Engineers, with SVY21 Worked Examples" by Yam Khoon Tor.

Copyright (c) 2026 Yam Khoon Tor. All rights reserved. Readers who have bought the
book may run, copy and modify this program for personal study and their own
professional work; it may not be republished, distributed or sold without the
author's written permission. Provided as a teaching aid, without warranty: verify
results independently before use in real work.

Contents: Listings 13.1, 13.2 (with the functions of Listings 2.1, F.1 that they use)
Run with:  python ch13_svy21_tm.py
Requires:  Python 3.9 or later
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




# ============================================================================
# Listing 13.1  Transverse Mercator forward, inverse, scale factor and convergence (Krüger series)
# ============================================================================
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


# ============================================================================
# Listing 13.2  Known-value and round-trip tests against the published SiReNT coordinates
# ============================================================================
# Published SiReNT coordinates (SLA): latitude d m s, longitude d m s,
# SVY21 northing, SVY21 easting
SIRENT = {
    "SNYU": ((1, 20, 44.84810), (103, 40, 47.91672), 36436.799, 10934.735),
    "SRPT": ((1, 26, 37.41926), (103, 47, 3.88829), 47265.636, 22557.618),
    "SNPT": ((1, 22, 44.82064), (103, 50, 55.51059), 40121.255, 29717.651),
    "SLYG": ((1, 22, 21.44707), (103, 58, 18.00928), 39403.770, 43396.751),
    "SNSC": ((1, 18, 49.12095), (103, 57, 42.01188), 32882.042, 42284.293),
    "SNUS": ((1, 17, 32.68755), (103, 46, 31.83287), 30534.079, 21566.307),
    "SSMK": ((1, 12, 37.48752), (103, 46, 31.49934), 21466.939, 21555.796),
    "SSTS": ((1, 14, 23.49077), (103, 38, 52.29482), 24723.506, 7359.574),
    "SMS1": ((1, 20, 25.60771), (103, 53, 16.40533), 35845.357, 34073.254),
}
# SVY21: origin latitude, central meridian, k0, false northing, false easting
SVY21_TM = (dms_to_deg(1, 22, 0), dms_to_deg(103, 50, 0), 1.0,
            38744.572, 28001.642)

if __name__ == "__main__":
    a, inv_f = ELLIPSOIDS["WGS84"]
    worst = 0.0
    for name, (lat_dms, lon_dms, N_pub, E_pub) in SIRENT.items():
        lat, lon = dms_to_deg(*lat_dms), dms_to_deg(*lon_dms)
        N, E = tm_forward(lat, lon, *SVY21_TM, a, inv_f)         # known-value test
        lat2, lon2 = tm_inverse(N, E, *SVY21_TM, a, inv_f)       # round-trip test
        k, gamma = tm_scale_convergence(lat, lon, SVY21_TM[1], 1.0, a, inv_f)
        worst = max(worst, abs(N - N_pub), abs(E - E_pub))
        dN, dE = (N - N_pub) * 1000, (E - E_pub) * 1000
        print(f"{name}  dN = {dN:+.1f} mm  dE = {dE:+.1f} mm  "
              f"k = {(k - 1) * 1e6:+.2f} ppm  gamma = {gamma * 3600:+.1f}\"")
        assert abs(lat2 - lat) < 1e-11 and abs(lon2 - lon) < 1e-11
    print(f"largest difference from published values: {worst * 1000:.1f} mm")
