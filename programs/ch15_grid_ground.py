"""ch15_grid_ground.py - Grid to ground

Companion program to "Datums, Coordinate Systems and Map Projections: A Practical
Guide for Surveyors and Engineers, with SVY21 Worked Examples" by Yam Khoon Tor.

Copyright (c) 2026 Yam Khoon Tor. All rights reserved. Readers who have bought the
book may run, copy and modify this program for personal study and their own
professional work; it may not be republished, distributed or sold without the
author's written permission. Provided as a teaching aid, without warranty: verify
results independently before use in real work.

Contents: Listings 15.1
Run with:  python ch15_grid_ground.py
Requires:  Python 3.9 or later
"""


# ============================================================================
# Listing 15.1  Reduction between ground and grid distances
# ============================================================================
import math

R_SG = 6356777.0            # mean radius of curvature in Singapore (m)
E0_SVY21 = 28001.642        # false easting of the SVY21 grid (m)


def horizontal_distance(slope, h1, h2):
    """Horizontal distance at the mean height of the line, from the slope
    distance and the heights of its two ends."""
    return math.sqrt(slope ** 2 - (h2 - h1) ** 2)


def elevation_factor(h_ellipsoidal, R=R_SG):
    """Factor reducing a horizontal distance at ellipsoidal height h
    to the ellipsoid."""
    return R / (R + h_ellipsoidal)


def point_scale_factor(E, E0=E0_SVY21, k0=1.0, R=R_SG):
    """Point scale factor of a Transverse Mercator grid at easting E."""
    x = E - E0
    return k0 * (1.0 + x * x / (2.0 * k0 * k0 * R * R))


def line_scale_factor(E1, E2, E0=E0_SVY21, k0=1.0, R=R_SG):
    """Line scale factor between two eastings."""
    x1, x2 = E1 - E0, E2 - E0
    return k0 * (1.0 + (x1 * x1 + x1 * x2 + x2 * x2) / (6.0 * k0 * k0 * R * R))


def ground_to_grid(horizontal, h_mean, E1, E2, **grid):
    """Grid distance from a horizontal ground distance.
    Returns the distance, the elevation factor and the line scale factor."""
    ef = elevation_factor(h_mean)
    K = line_scale_factor(E1, E2, **grid)
    return horizontal * ef * K, ef, K


def grid_to_ground(grid_distance, h_mean, E1, E2, **grid):
    """Horizontal ground distance from a grid distance.
    Returns the distance, the elevation factor and the line scale factor."""
    ef = elevation_factor(h_mean)
    K = line_scale_factor(E1, E2, **grid)
    return grid_distance / (ef * K), ef, K


def ppm(factor):
    """Departure of a scale factor from unity, in parts per million."""
    return (factor - 1.0) * 1.0e6


if __name__ == "__main__":
    # Example 15.1: a measured slope distance reduced to the SVY21 grid
    slope = 1254.386
    h1 = 12.400 + 1.550 + 9.7     # reduced level + instrument height + undulation
    h2 = 51.814 + 1.400 + 9.7     # reduced level + target height + undulation
    D = horizontal_distance(slope, h1, h2)
    grid, ef, K = ground_to_grid(D, (h1 + h2) / 2.0, 8500.0, 9712.9)
    print(f"horizontal = {D:.4f} m")
    print(f"EF = {ppm(ef):+.2f} ppm   K = {ppm(K):+.2f} ppm   "
          f"combined = {ppm(ef * K):+.2f} ppm")
    print(f"grid distance = {grid:.4f} m")

    # Example 15.2: a grid distance converted for setting out
    ground, ef, K = grid_to_ground(350.000, 14.0 + 9.7, 4200.0, 4520.0)
    print(f"EF = {ppm(ef):+.2f} ppm   K = {ppm(K):+.2f} ppm   "
          f"ground distance = {ground:.4f} m")

    # Example 15.3: the line of Example 15.1 on UTM zone 48 North
    grid, ef, K = ground_to_grid(D, (h1 + h2) / 2.0, 350704.461, 351917.369,
                                 E0=500000.0, k0=0.9996)
    print(f"UTM: K = {ppm(K):+.1f} ppm   grid distance = {grid:.4f} m")
