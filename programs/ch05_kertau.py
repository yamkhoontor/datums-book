"""ch05_kertau.py - Classical datums and their realisation

Companion program to "Datums, Coordinate Systems and Map Projections: A Practical
Guide for Surveyors and Engineers, with SVY21 Worked Examples" by Yam Khoon Tor.

Copyright (c) 2026 Yam Khoon Tor. All rights reserved. Readers who have bought the
book may run, copy and modify this program for personal study and their own
professional work; it may not be republished, distributed or sold without the
author's written permission. Provided as a teaching aid, without warranty: verify
results independently before use in real work.

Contents: Listings 5.1 (with the functions of Listings 2.1, F.1, 3.1 that they use)
Run with:  python ch05_kertau.py
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




# ============================================================================
# Listing 5.1  Three-parameter transformation from a geocentric datum to Kertau 1968
# ============================================================================
# Kertau 1968 to WGS84 geocentric translations (m), DMA 1987
KERTAU_TO_WGS84 = (-11.0, 851.0, 5.0)


def wgs84_to_kertau(lat_deg, lon_deg, h, source="WGS84"):
    """Three-parameter transformation to Kertau 1968."""
    a_s, inv_f_s = ELLIPSOIDS[source]
    a_t, inv_f_t = ELLIPSOIDS["Everest 1830 (Modified)"]
    X, Y, Z = geodetic_to_ecef(lat_deg, lon_deg, h, a_s, inv_f_s)
    dX, dY, dZ = KERTAU_TO_WGS84
    X, Y, Z = X - dX, Y - dY, Z - dZ      # reverse direction: subtract
    return ecef_to_geodetic(X, Y, Z, a_t, inv_f_t)


if __name__ == "__main__":
    # mark S of the running case
    lat, lon, h = wgs84_to_kertau(dms_to_deg(1, 27, 10.00803),
                                  dms_to_deg(103, 46, 9.97011), 12.014)
    print(deg_to_dms(lat), deg_to_dms(lon), f"{h:.3f}")

    lat, lon, h = wgs84_to_kertau(dms_to_deg(1, 51, 27.04635),
                                  dms_to_deg(102, 56, 31.75616),
                                  7.724, source="GRS80")
    print(deg_to_dms(lat), deg_to_dms(lon), f"{h:.3f}")
