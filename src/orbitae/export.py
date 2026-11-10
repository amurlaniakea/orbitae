# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Exporta nubes electrónicas precalculadas para el visor web.

Formato de cada .bin: N posiciones int16 (x, y, z intercalados, little-endian,
escala = extent/32767) seguidas de N bytes de signo (1 si psi < 0).
Solo orbitales s, p y d (n <= 3) por ahora.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .hydrogenic import energy_hartree, psi_real, radial_cutoff, sample_orbital

A0_ANGSTROM = 0.529177210903
HARTREE_EV = 27.211386245988
_NAMES = {
    (1, -1): "p_y", (1, 0): "p_z", (1, 1): "p_x",
    (2, -2): "d_xy", (2, -1): "d_yz", (2, 0): "d_z²", (2, 1): "d_xz", (2, 2): "d_x²-y²",
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
        "radial_nodes": n - l - 1, "angular_nodes": l,
    }
    return meta, payload


def main(argv=None):
    ap = argparse.ArgumentParser(description="Exporta orbitales para el visor web")
    ap.add_argument("--out", default="web/data")
    ap.add_argument("--nmax", type=int, default=3)
    ap.add_argument("--points", type=int, default=30000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    items = []
    for i, (n, l, m) in enumerate(catalog(a.nmax)):
        meta, payload = export_orbital(n, l, m, size=a.points, seed=a.seed + i)
        (out / meta["file"]).write_bytes(payload)
        items.append(meta)
    (out / "index.json").write_text(
        json.dumps({"z": 1, "orbitals": items}, ensure_ascii=False, indent=1), encoding="utf-8"
    )
    print(f"{len(items)} orbitales exportados en {out}")


if __name__ == "__main__":
    main()
