"""ch07_itrf_plate_motion.py - ITRF and a moving Earth

Companion program to "Datums, Coordinate Systems and Map Projections: A Practical
Guide for Surveyors and Engineers, with SVY21 Worked Examples" by Yam Khoon Tor.

Copyright (c) 2026 Yam Khoon Tor. All rights reserved. Readers who have bought the
book may run, copy and modify this program for personal study and their own
professional work; it may not be republished, distributed or sold without the
author's written permission. Provided as a teaching aid, without warranty: verify
results independently before use in real work.

Contents: Listings 7.1 (with the functions of Listings 2.1, F.1, 3.1, 3.2 that they use)
Run with:  python ch07_itrf_plate_motion.py
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


# ---------------- from Listing 3.2: Local topocentric components of a baseline ----------------
def ecef_to_enu(dX, dY, dZ, lat0_deg, lon0_deg):
    """Rotate a Cartesian difference into east, north, up at the origin."""
    phi, lam = math.radians(lat0_deg), math.radians(lon0_deg)
    sp, cp = math.sin(phi), math.cos(phi)
    sl, cl = math.sin(lam), math.cos(lam)
    e = -sl * dX + cl * dY
    n = -sp * cl * dX - sp * sl * dY + cp * dZ
    u = cp * cl * dX + cp * sl * dY + sp * dZ
    return e, n, u


def enu_to_ecef(e, n, u, lat0_deg, lon0_deg):
    """Rotate east, north, up at the origin into a Cartesian difference."""
    phi, lam = math.radians(lat0_deg), math.radians(lon0_deg)
    sp, cp = math.sin(phi), math.cos(phi)
    sl, cl = math.sin(lam), math.cos(lam)
    dX = -sl * e - sp * cl * n + cp * cl * u
    dY = cl * e - sp * sl * n + cp * sl * u
    dZ = cp * n + sp * u
    return dX, dY, dZ




# ============================================================================
# Listing 7.1  Plate velocity, epoch propagation and the 14-parameter transformation
# ============================================================================
import math

MAS = math.radians(1.0 / 3600.0 / 1000.0)     # one milliarc-second in radians

# Sunda plate Euler pole in ITRF2000 (Simons et al., 2007):
# latitude, longitude (degrees), rate (degrees per million years)
SUNDA_POLE = (49.0, -94.2, 0.336)

# ITRF2020 to ITRF2000 (IERS): T1, T2, T3 (mm), D (ppb), R1, R2, R3 (mas),
# then their rates per year, then the reference epoch
ITRF2020_TO_ITRF2000 = ((-0.2, 0.8, -34.2, 2.25, 0.0, 0.0, 0.0),
                        (0.1, 0.0, -1.7, 0.11, 0.0, 0.0, 0.0),
                        2015.0)


def euler_pole_velocity(lat_deg, lon_deg, pole_lat, pole_lon, rate_deg_myr,
                        R=6371000.0):
    """Velocity of a point on a rigid plate, from the plate's Euler pole.

    Returns east and north components and the Cartesian vector, in m/yr.
    """
    w = math.radians(rate_deg_myr) * 1.0e-6          # radians per year
    pl, pp = math.radians(pole_lat), math.radians(pole_lon)
    wx = w * math.cos(pl) * math.cos(pp)
    wy = w * math.cos(pl) * math.sin(pp)
    wz = w * math.sin(pl)
    phi, lam = math.radians(lat_deg), math.radians(lon_deg)
    x = R * math.cos(phi) * math.cos(lam)            # position on a sphere
    y = R * math.cos(phi) * math.sin(lam)
    z = R * math.sin(phi)
    vx = wy * z - wz * y                             # V = Omega x r
    vy = wz * x - wx * z
    vz = wx * y - wy * x
    ve, vn, vu = ecef_to_enu(vx, vy, vz, lat_deg, lon_deg)
    return ve, vn, (vx, vy, vz)


def propagate(X, V, t_from, t_to):
    """Move Cartesian coordinates X from epoch t_from to epoch t_to."""
    dt = t_to - t_from
    return tuple(X[i] + V[i] * dt for i in range(3))


def helmert14(X, params, rates, t0, t):
    """14-parameter transformation at epoch t, position vector convention."""
    p = [params[i] + rates[i] * (t - t0) for i in range(7)]
    T1, T2, T3 = p[0] / 1000.0, p[1] / 1000.0, p[2] / 1000.0
    D = p[3] * 1.0e-9
    R1, R2, R3 = p[4] * MAS, p[5] * MAS, p[6] * MAS
    x, y, z = X
    X2 = x + T1 + D * x - R3 * y + R2 * z
    Y2 = y + T2 + R3 * x + D * y - R1 * z
    Z2 = z + T3 - R2 * x + R1 * y + D * z
    return X2, Y2, Z2


if __name__ == "__main__":
    a, inv_f = ELLIPSOIDS["GRS80"]
    # mark S of the running case
    lat, lon, h = dms_to_deg(1, 27, 10.00803), dms_to_deg(103, 46, 9.97011), 12.014

    # Example 7.1: velocity from the pole
    ve, vn, V = euler_pole_velocity(lat, lon, *SUNDA_POLE)
    print(f"ve = {ve * 1000:.1f} mm/yr   vn = {vn * 1000:.1f} mm/yr")

    # Example 7.4: ITRF2020 at 2026.2 to ITRF2000 at 2000.0
    X = geodetic_to_ecef(lat, lon, h, a, inv_f)
    params, rates, t0 = ITRF2020_TO_ITRF2000
    X1 = helmert14(X, params, rates, t0, 2026.2)      # change frame
    X2 = propagate(X1, V, 2026.2, 2000.0)             # change epoch
    print(f"ITRF2000 at 2026.2: {X1[0]:.4f} {X1[1]:.4f} {X1[2]:.4f}")
    print(f"ITRF2000 at 2000.0: {X2[0]:.4f} {X2[1]:.4f} {X2[2]:.4f}")

    lat2, lon2, h2 = ecef_to_geodetic(*X2, a, inv_f)
    print(deg_to_dms(lat2, 5), deg_to_dms(lon2, 5), f"{h2:.3f}")
