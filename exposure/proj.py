"""Map projections needed for the EMA exposure analysis, numpy only.

* tm_forward : geographic (WGS84) -> UTM zone 35N (EPSG:32635), Krüger series
* albers_inverse : custom Albers equal-area (lon0 35, lat1 36.5, lat2 41.0,
  lat0 0, WGS84) -> geographic.  This is the CRS of the district and mask
  shapefiles used in the published analysis.
"""
import os
import numpy as np

A_WGS84 = 6378137.0
F_WGS84 = 1.0 / 298.257223563
E2 = 2 * F_WGS84 - F_WGS84 ** 2
E = np.sqrt(E2)


# ------------------------------------------------------------ transverse mercator
def tm_forward(lon_deg, lat_deg, lon0_deg=27.0, k0=0.9996, fe=500000.0, fn=0.0):
    """Krüger-series transverse Mercator (Karney 2011), accurate to < 1 mm in zone."""
    n = F_WGS84 / (2 - F_WGS84)
    A = A_WGS84 / (1 + n) * (1 + n ** 2 / 4 + n ** 4 / 64)
    a1 = n / 2 - 2 * n ** 2 / 3 + 5 * n ** 3 / 16 + 41 * n ** 4 / 180
    a2 = 13 * n ** 2 / 48 - 3 * n ** 3 / 5 + 557 * n ** 4 / 1440
    a3 = 61 * n ** 3 / 240 - 103 * n ** 4 / 140
    a4 = 49561 * n ** 4 / 161280
    phi = np.deg2rad(np.asarray(lat_deg, dtype=np.float64))
    lam = np.deg2rad(np.asarray(lon_deg, dtype=np.float64) - lon0_deg)
    sphi = np.sin(phi)
    c = 2 * np.sqrt(n) / (1 + n)
    t = np.sinh(np.arctanh(sphi) - c * np.arctanh(c * sphi))
    xi_ = np.arctan2(t, np.cos(lam))
    eta_ = np.arctanh(np.sin(lam) / np.sqrt(1 + t * t))
    xi = xi_.copy()
    eta = eta_.copy()
    for j, aj in enumerate((a1, a2, a3, a4), start=1):
        xi += aj * np.sin(2 * j * xi_) * np.cosh(2 * j * eta_)
        eta += aj * np.cos(2 * j * xi_) * np.sinh(2 * j * eta_)
    return fe + k0 * A * eta, fn + k0 * A * xi


# ------------------------------------------------------------ albers (ellipsoid)
def _q(phi):
    s = np.sin(phi)
    return (1 - E2) * (s / (1 - E2 * s * s)
                       - (1 / (2 * E)) * np.log((1 - E * s) / (1 + E * s)))


def _m(phi):
    return np.cos(phi) / np.sqrt(1 - E2 * np.sin(phi) ** 2)


def _albers_consts(lon0=35.0, lat1=36.5, lat2=41.0, lat0=0.0):
    p1, p2, p0 = np.deg2rad([lat1, lat2, lat0])
    m1, m2 = _m(p1), _m(p2)
    q1, q2, q0 = _q(p1), _q(p2), _q(p0)
    n = (m1 ** 2 - m2 ** 2) / (q2 - q1)
    C = m1 ** 2 + n * q1
    rho0 = A_WGS84 * np.sqrt(C - n * q0) / n
    return np.deg2rad(lon0), n, C, rho0


def albers_forward(lon_deg, lat_deg, **kw):
    lam0, n, C, rho0 = _albers_consts(**kw)
    phi = np.deg2rad(np.asarray(lat_deg, dtype=np.float64))
    lam = np.deg2rad(np.asarray(lon_deg, dtype=np.float64))
    rho = A_WGS84 * np.sqrt(C - n * _q(phi)) / n
    th = n * (lam - lam0)
    return rho * np.sin(th), rho0 - rho * np.cos(th)


def albers_inverse(x, y, **kw):
    lam0, n, C, rho0 = _albers_consts(**kw)
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    rho = np.hypot(x, rho0 - y)
    th = np.arctan2(x, rho0 - y)          # n > 0 in the northern hemisphere
    q = (C - rho * rho * n * n / A_WGS84 ** 2) / n
    lam = lam0 + th / n
    phi = np.arcsin(np.clip(q / 2, -1, 1))
    for _ in range(8):                   # Snyder (3-16)
        s = np.sin(phi)
        d = 1 - E2 * s * s
        phi = phi + (d * d / (2 * np.cos(phi))) * (
            q / (1 - E2) - s / d + (1 / (2 * E)) * np.log((1 - E * s) / (1 + E * s)))
    return np.rad2deg(lam), np.rad2deg(phi)


def albers_to_utm(x, y):
    lon, lat = albers_inverse(x, y)
    return tm_forward(lon, lat)


if __name__ == "__main__":
    # The published grid was built from bbox (26.0, 39.8, 30.5, 41.5) -> UTM 35N
    x_sw, y_sw = tm_forward(26.0, 39.8)
    x_ne, y_ne = tm_forward(30.5, 41.5)
    x_nw, y_nw = tm_forward(26.0, 41.5)
    x_se, y_se = tm_forward(30.5, 39.8)
    print("corner SW", x_sw, y_sw)
    print("corner NE", x_ne, y_ne)
    print("corner NW", x_nw, y_nw)
    print("corner SE", x_se, y_se)
    print("grid x0 = 414390.994, y0 = 4600180.695 (published)")
    # round trip for albers
    lon, lat = 28.9784, 41.0082
    ax, ay = albers_forward(lon, lat)
    print("albers fwd", ax, ay, "inv", albers_inverse(ax, ay))
