"""ch02_ellipsoid.py - The Earth's figure

Companion program to "Datums, Coordinate Systems and Map Projections: A Practical
Guide for Surveyors and Engineers, with SVY21 Worked Examples" by Yam Khoon Tor.

Copyright (c) 2026 Yam Khoon Tor. All rights reserved. Readers who have bought the
book may run, copy and modify this program for personal study and their own
professional work; it may not be republished, distributed or sold without the
author's written permission. Provided as a teaching aid, without warranty: verify
results independently before use in real work.

Contents: Listings 2.1
Run with:  python ch02_ellipsoid.py
Requires:  Python 3.9 or later
"""


# ============================================================================
# Listing 2.1  Ellipsoid constants and radii of curvature
# ============================================================================
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


if __name__ == "__main__":
    a, inv_f = ELLIPSOIDS["WGS84"]
    c = ellipsoid_constants(a, inv_f)
    print(f"b   = {c['b']:.4f} m")
    print(f"e2  = {c['e2']:.12f}")
    print(f"ep2 = {c['ep2']:.12f}")
    print(f"n   = {c['n']:.12f}")

    lat = dms_to_deg(1, 22, 0)
    M, nu, R = radii_of_curvature(lat, a, inv_f)
    print(f"M   = {M:.3f} m")
    print(f"nu  = {nu:.3f} m")
    print(f"R   = {R:.3f} m")
