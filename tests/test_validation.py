# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
import math

import numpy as np
import pytest
from scipy.integrate import quad

from orbitae import (
    energy_hartree,
    radial,
    radial_cutoff,
    real_harmonic,
    sample_orbital,
)

RADIAL_CASES = [(1, 0, 1.0), (2, 1, 1.0), (3, 2, 1.0), (4, 2, 1.0), (5, 3, 1.0), (2, 0, 2.0)]


@pytest.mark.parametrize("n,l,Z", RADIAL_CASES)
def test_radial_normalization(n, l, Z):
    val = quad(lambda r: r**2 * radial(r, n, l, Z) ** 2, 0, np.inf, limit=200)[0]
    assert val == pytest.approx(1.0, abs=1e-7)


@pytest.mark.parametrize("l,n1,n2", [(0, 1, 2), (0, 2, 3), (1, 2, 3), (2, 3, 4)])
def test_radial_orthogonality(l, n1, n2):
    val = quad(lambda r: r**2 * radial(r, n1, l) * radial(r, n2, l), 0, np.inf, limit=200)[0]
    assert abs(val) < 1e-7


def test_real_harmonics_orthonormal():
    x, w = np.polynomial.legendre.leggauss(40)
    theta = np.arccos(x)
    n_phi = 80
    phi = np.arange(n_phi) * (2 * np.pi / n_phi)
    T, P = np.meshgrid(theta, phi, indexing="ij")
    W = np.outer(w, np.full(n_phi, 2 * np.pi / n_phi))
    labels = [(l, m) for l in range(4) for m in range(-l, l + 1)]
    vals = [real_harmonic(l, m, T, P) for l, m in labels]
    gram = np.array([[np.sum(a * b * W) for b in vals] for a in vals])
    assert np.allclose(gram, np.eye(len(labels)), atol=1e-12)


def test_p_and_d_shapes_match_textbooks():
    th = np.linspace(0.1, 3.0, 7)
    ph = np.linspace(0.0, 6.0, 7)
    x, y, z = np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)
    c1 = math.sqrt(3 / (4 * math.pi))
    assert np.allclose(real_harmonic(1, 1, th, ph), c1 * x)    # p_x
    assert np.allclose(real_harmonic(1, -1, th, ph), c1 * y)   # p_y
    assert np.allclose(real_harmonic(1, 0, th, ph), c1 * z)    # p_z
    c2 = 0.5 * math.sqrt(15 / math.pi)
    assert np.allclose(real_harmonic(2, -2, th, ph), c2 * x * y)  # d_xy


@pytest.mark.parametrize("n,l,Z", [(1, 0, 1.0), (2, 1, 1.0), (3, 2, 2.0), (3, 1, 1.0)])
def test_radial_schrodinger_equation(n, l, Z):
    h = 1e-3
    r = np.arange(0.5, 12.0, h)
    u = r * radial(r, n, l, Z)
    upp = (u[2:] - 2 * u[1:-1] + u[:-2]) / h**2
    rc, uc = r[1:-1], u[1:-1]
    resid = -0.5 * upp + (l * (l + 1) / (2 * rc**2) - Z / rc) * uc - energy_hartree(n, Z) * uc
    assert np.max(np.abs(resid)) < 1e-4


def test_energy_values():
    assert energy_hartree(1) == pytest.approx(-0.5)
    assert energy_hartree(1) * 27.211386 == pytest.approx(-13.6057, abs=1e-3)
    assert energy_hartree(2, Z=2) == pytest.approx(-0.5)


def test_cutoff_1s_matches_closed_form():
    r95 = radial_cutoff(1, 0, 0.95)
    p_in = 1.0 - math.exp(-2 * r95) * (1 + 2 * r95 + 2 * r95**2)
    assert p_in == pytest.approx(0.95, abs=1e-4)


def test_sampling_1s_mean_radius():
    pts = sample_orbital(1, 0, 0, size=200_000, seed=1)
    assert np.linalg.norm(pts, axis=1).mean() == pytest.approx(1.5, abs=0.02)


def test_sampling_px_is_anisotropic():
    pts = sample_orbital(2, 1, 1, size=200_000, seed=2)  # 2p_x
    r2 = np.sum(pts**2, axis=1)
    assert np.mean(pts[:, 0] ** 2) / np.mean(r2) == pytest.approx(0.6, abs=0.01)
    assert np.mean(pts[:, 1] ** 2) / np.mean(r2) == pytest.approx(0.2, abs=0.01)
