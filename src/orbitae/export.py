# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Exporta nubes electrónicas y perfiles radiales precalculados para el visor web.

Formato de cada .bin: N posiciones int16 (x, y, z intercalados, little-endian,
escala = extent/32767) seguidas de N bytes de signo (1 si psi < 0).
index.json incluye ademas, por cada par (n, l), la distribucion radial r^2 R^2
normalizada a su maximo y la posicion de los nodos radiales (unidades a0).
Orbitales s, p, d y f (n <= 4).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .hydrogenic import energy_hartree, psi_real, radial, radial_cutoff, sample_orbital

A0_ANGSTROM = 0.529177210903
HARTREE_EV = 27.211386245988
_NAMES = {
    (1, -1): "p_y", (1, 0): "p_z", (1, 1): "p_x",
    (2, -2): "d_xy", (2, -1): "d_yz", (2, 0): "d_z²", (2, 1): "d_xz", (2, 2): "d_x²-y²",
    (3, -3): "f_y(3x²-y²)", (3, -2): "f_xyz", (3, -1): "f_yz²", (3, 0): "f_z³",
    (3, 1): "f_xz²", (3, 2): "f_z(x²-y²)", (3, 3): "f_x(x²-3y²)",
}


def orbital_label(n: int, l: int, m: int) -> str:
    if l == 0:
        return f"{n}s"
    if (l, m) not in _NAMES:
        raise ValueError(f"Orbital no soportado todavía: l={l}, m={m}")
    return f"{n}{_NAMES[(l, m)]}"


def orbital_file(n: int, l: int, m: int) -> str:
    return f"o{n}{l}{'m' if m < 0 else 'p'}{abs(m)}.bin"


def catalog(nmax: int):
    return [(n, l, m) for n in range(1, nmax + 1) for l in range(n) for m in range(-l, l + 1)]


def radial_nodes(n: int, l: int, Z: float = 1.0):
    """Radios (a0) donde R_nl cambia de signo; hay n - l - 1."""
    r = np.linspace(1e-6, (6.0 * n * n + 30.0) / Z, 60001)
    R = radial(r, n, l, Z)
    idx = np.where(np.sign(R[:-1]) * np.sign(R[1:]) < 0)[0]
    return [float(r[i] - R[i] * (r[i + 1] - r[i]) / (R[i + 1] - R[i])) for i in idx]


def radial_profile(n: int, l: int, extent: float, Z: float = 1.0, npts: int = 300):
    r = np.linspace(0.0, extent, npts)
    p = r**2 * radial(r, n, l, Z) ** 2
    p = p / p.max()
    nodes = [round(x, 4) for x in radial_nodes(n, l, Z) if x <= extent]
    return {"r": np.round(r, 4).tolist(), "p": np.round(p, 5).tolist(), "nodes": nodes}


def export_orbital(n, l, m, size=30000, coverage=0.995, Z=1.0, seed=0):
    extent = radial_cutoff(n, l, coverage, Z)
    pts = sample_orbital(n, l, m, size=int(size * 1.03) + 10, Z=Z, seed=seed)
    pts = pts[np.linalg.norm(pts, axis=1) <= extent][:size]
    if len(pts) < size:
        raise RuntimeError("Muestras insuficientes tras el truncado")
    r = np.maximum(np.linalg.norm(pts, axis=1), 1e-12)
    theta = np.arccos(np.clip(pts[:, 2] / r, -1.0, 1.0))
    phi = np.arctan2(pts[:, 1], pts[:, 0])
    sign = (psi_real(n, l, m, r, theta, phi, Z) < 0).astype(np.uint8)
    scale = extent / 32767.0
    payload = np.round(pts / scale).astype("<i2").tobytes() + sign.tobytes()
    r95 = radial_cutoff(n, l, 0.95, Z)
    meta = {
        "label": orbital_label(n, l, m), "n": n, "l": l, "m": m, "count": size,
        "file": orbital_file(n, l, m), "extent_a0": extent, "scale": scale,
        "r95_a0": r95, "r95_angstrom": r95 * A0_ANGSTROM / Z,
        "energy_ev": energy_hartree(n, Z) * HARTREE_EV,
        "radial_nodes": n - l - 1, "angular_nodes": l, "radial_key": f"{n},{l}",
    }
    return meta, payload


def main(argv=None):
    ap = argparse.ArgumentParser(description="Exporta orbitales para el visor web")
    ap.add_argument("--out", default="web/data")
    ap.add_argument("--nmax", type=int, default=4)
    ap.add_argument("--points", type=int, default=30000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    items, rad = [], {}
    for i, (n, l, m) in enumerate(catalog(a.nmax)):
        meta, payload = export_orbital(n, l, m, size=a.points, seed=a.seed + i)
        (out / meta["file"]).write_bytes(payload)
        items.append(meta)
        if meta["radial_key"] not in rad:
            rad[meta["radial_key"]] = radial_profile(n, l, meta["extent_a0"])
    (out / "index.json").write_text(
        json.dumps({"z": 1, "orbitals": items, "radial": rad}, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"{len(items)} orbitales exportados en {out}")


if __name__ == "__main__":
    main()
