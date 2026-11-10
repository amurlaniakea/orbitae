# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Orbitales hidrogenoides en unidades atómicas (a0 = 1, energía en hartree).

Convención: armónicos esféricos con fase de Condon-Shortley (scipy.special.sph_harm_y),
y armónicos reales estándar, con lóbulo positivo hacia +x, +y, +z.
Solo válido para sistemas de un electrón (H, He+, Li2+, ...) con carga nuclear Z.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import eval_genlaguerre, sph_harm_y


def _check(n: int, l: int, m: int | None = None) -> None:
    if n < 1 or not (0 <= l <= n - 1):
        raise ValueError(f"Números cuánticos inválidos: n={n}, l={l}")
    if m is not None and abs(m) > l:
        raise ValueError(f"m={m} fuera de rango para l={l}")


def energy_hartree(n: int, Z: float = 1.0) -> float:
    """E_n = -Z^2 / (2 n^2) hartree (= -13.6057 Z^2/n^2 eV)."""
    _check(n, 0)
    return -(Z**2) / (2.0 * n**2)


def radial(r, n: int, l: int, Z: float = 1.0):
    """Parte radial normalizada: integral de r^2 R^2 dr = 1."""
    _check(n, l)
    r = np.asarray(r, dtype=float)
    rho = 2.0 * Z * r / n
    norm = math.sqrt(
        (2.0 * Z / n) ** 3 * math.factorial(n - l - 1) / (2.0 * n * math.factorial(n + l))
    )
    return norm * np.exp(-rho / 2.0) * rho**l * eval_genlaguerre(n - l - 1, 2 * l + 1, rho)


def real_harmonic(l: int, m: int, theta, phi):
    """Armónico esférico real. theta = ángulo polar, phi = azimutal.

    m=0: Y_l^0; m>0: sqrt2 (-1)^m Re Y_l^m; m<0: sqrt2 (-1)^m Im Y_l^|m|.
    """
    _check(l + 1, l, m)
    mm = abs(m)
    y = sph_harm_y(l, mm, theta, phi)
    if m == 0:
        return y.real
    sign = -1.0 if mm % 2 else 1.0
    return math.sqrt(2.0) * sign * (y.real if m > 0 else y.imag)


def psi_real(n: int, l: int, m: int, r, theta, phi, Z: float = 1.0):
    """Función de onda con armónico real (p_x, d_xy, ...)."""
    return radial(r, n, l, Z) * real_harmonic(l, m, theta, phi)


def psi_complex(n: int, l: int, m: int, r, theta, phi, Z: float = 1.0):
    """Función de onda con armónico complejo (m=±1 da densidad toroidal)."""
    _check(n, l, m)
    return radial(r, n, l, Z) * sph_harm_y(l, m, theta, phi)


def _radial_cdf(n: int, l: int, Z: float, npts: int = 40001):
    rmax = (6.0 * n * n + 30.0) / Z
    r = np.linspace(0.0, rmax, npts)
    p = r**2 * radial(r, n, l, Z) ** 2
    cdf = cumulative_trapezoid(p, r, initial=0.0)
    if cdf[-1] < 0.999:
        raise RuntimeError("Rejilla radial insuficiente para estos números cuánticos")
    return r, cdf / cdf[-1]


def radial_cutoff(n: int, l: int, p: float = 0.95, Z: float = 1.0) -> float:
    """Radio (a0) que encierra la fracción p de la distribución radial r^2 R^2."""
    _check(n, l)
    if not 0.0 < p < 1.0:
        raise ValueError("p debe estar en (0, 1)")
    r, cdf = _radial_cdf(n, l, Z)
    return float(np.interp(p, cdf, r))


def sample_orbital(
    n: int,
    l: int,
    m: int,
    size: int,
    Z: float = 1.0,
    real: bool = True,
    seed: int | None = None,
):
    """Muestrea puntos (x, y, z) distribuidos según |psi|^2.

    Como psi = R(r) Y(theta, phi), la densidad se separa: r por CDF inversa de
    r^2 R^2 y los ángulos por CDF discreta de |Y|^2 sin(theta) en rejilla fina.
    """
    _check(n, l, m)
    rng = np.random.default_rng(seed)

    r_grid, cdf = _radial_cdf(n, l, Z)
    r = np.interp(rng.random(size), cdf, r_grid)

    n_t, n_p = 361, 720
    et = np.linspace(0.0, np.pi, n_t + 1)
    ep = np.linspace(0.0, 2.0 * np.pi, n_p + 1)
    tc, pc = np.meshgrid(0.5 * (et[:-1] + et[1:]), 0.5 * (ep[:-1] + ep[1:]), indexing="ij")
    if real:
        dens = real_harmonic(l, m, tc, pc) ** 2
    else:
        dens = np.abs(sph_harm_y(l, m, tc, pc)) ** 2
    dens = (dens * np.sin(tc)).ravel()
    idx = rng.choice(dens.size, size=size, p=dens / dens.sum())
    it, ip = np.divmod(idx, n_p)
    theta = et[it] + rng.random(size) * (et[1] - et[0])
    phi = ep[ip] + rng.random(size) * (ep[1] - ep[0])

    st = np.sin(theta)
    return np.stack([r * st * np.cos(phi), r * st * np.sin(phi), r * np.cos(theta)], axis=1)
