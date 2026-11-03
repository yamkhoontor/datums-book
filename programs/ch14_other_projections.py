"""ch14_other_projections.py - Other projections a surveyor meets

Companion program to "Datums, Coordinate Systems and Map Projections: A Practical
Guide for Surveyors and Engineers, with SVY21 Worked Examples" by Yam Khoon Tor.

Copyright (c) 2026 Yam Khoon Tor. All rights reserved. Readers who have bought the
book may run, copy and modify this program for personal study and their own
professional work; it may not be republished, distributed or sold without the
author's written permission. Provided as a teaching aid, without warranty: verify
results independently before use in real work.

Contents: Listings 14.1, 14.2, 14.3 (with the functions of Listings 2.1, F.1, 13.1 that they use)
Run with:  python ch14_other_projections.py
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




# ============================================================================
# Listing 14.1  Cassini-Soldner projection
# ============================================================================
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


# ============================================================================
# Listing 14.2  Hotine oblique Mercator projection, variant A
# ============================================================================
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
# Listing 14.3  The worked examples of this chapter
# ============================================================================
if __name__ == "__main__":
    a_g, f_g = ELLIPSOIDS["GRS80"]
    a_w, f_w = ELLIPSOIDS["WGS84"]

    # Example 14.1: bridge mark S on UTM zone 48 North
    lat, lon = dms_to_deg(1, 27, 10.00803), dms_to_deg(103, 46, 9.97011)
    zone = int((lon + 180.0) // 6) + 1
    cm = 6 * zone - 183
    N, E = tm_forward(lat, lon, 0.0, cm, 0.9996, 0.0, 500000.0, a_w, f_w)
    k, gamma = tm_scale_convergence(lat, lon, cm, 0.9996, a_w, f_w)
    print(f"UTM zone {zone}, central meridian {cm} E")
    print(f"N = {N:.3f}  E = {E:.3f}  k = {k:.7f}  gamma = {gamma * 3600:.1f}\"")

    # Example 14.2: JUPEM test point on the Kedah and Perlis Cassini grid
    lat = dms_to_deg(6, 8, 22.98892)
    lon = dms_to_deg(100, 23, 6.57684)
    origin = (dms_to_deg(5, 57, 52.82155), dms_to_deg(100, 38, 10.93860))
    N, E = cassini_forward(lat, lon, *origin, 0.0, 0.0, a_g, f_g)
    print(f"Cassini  N = {N:.4f}  E = {E:.4f}")
    print("         (published 19364.3195, -27805.4063)")
    lat2, lon2 = cassini_inverse(N, E, *origin, 0.0, 0.0, a_g, f_g)
    print("inverse:", deg_to_dms(lat2, 5), deg_to_dms(lon2, 5))

    # Example 14.3: the same point on the Peninsula RSO grid
    MRSO = hom_setup(4.0, 102.25, dms_to_deg(323, 1, 32.8458),
                     dms_to_deg(323, 7, 48.3685), 0.99984, 804671.0, 0.0,
                     a_g, f_g)
    N, E = hom_forward(lat, lon, MRSO)
    print(f"RSO      N = {N:.3f}  E = {E:.3f}")
    print("         (published 679690.914, 266843.898)")
    lat2, lon2 = hom_inverse(N, E, MRSO)
    print("inverse:", deg_to_dms(lat2, 5), deg_to_dms(lon2, 5))
