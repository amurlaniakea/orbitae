# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
import json

import numpy as np
import pytest

from orbitae import real_harmonic
from orbitae.export import catalog, export_orbital, main, orbital_label, radial_nodes


def decode(payload, count, scale):
    q = np.frombuffer(payload[: 6 * count], dtype="<i2").reshape(count, 3).astype(float) * scale
    return q, np.frombuffer(payload[6 * count:], dtype=np.uint8)


def test_labels_and_catalog():
    assert orbital_label(1, 0, 0) == "1s"
    assert orbital_label(2, 1, 1) == "2p_x"
    assert orbital_label(3, 2, -2) == "3d_xy"
    assert orbital_label(4, 3, 0) == "4f_z³"
    assert len(catalog(3)) == 14
    assert len(catalog(4)) == 30


def test_layout_and_extent():
    meta, payload = export_orbital(2, 1, 1, size=3000, seed=3)
    assert len(payload) == 7 * meta["count"]
    q, _ = decode(payload, meta["count"], meta["scale"])
    assert np.linalg.norm(q, axis=1).max() <= meta["extent_a0"] * (1 + 1e-3)


def test_sign_follows_px_lobes():
    meta, payload = export_orbital(2, 1, 1, size=3000, seed=3)
    q, sign = decode(payload, meta["count"], meta["scale"])
    mask = np.abs(q[:, 0]) > 2 * meta["scale"]
    assert np.all(sign[mask] == (q[mask, 0] < 0))


def test_1s_has_single_phase():
    meta, payload = export_orbital(1, 0, 0, size=2000, seed=1)
    _, sign = decode(payload, meta["count"], meta["scale"])
    assert sign.sum() == 0


F_FORMS = {
    -3: lambda x, y, z: y * (3 * x**2 - y**2),
    -2: lambda x, y, z: x * y * z,
    -1: lambda x, y, z: y * (5 * z**2 - 1),
    0: lambda x, y, z: z * (5 * z**2 - 3),
    1: lambda x, y, z: x * (5 * z**2 - 1),
    2: lambda x, y, z: z * (x**2 - y**2),
    3: lambda x, y, z: x * (x**2 - 3 * y**2),
}


@pytest.mark.parametrize("m", sorted(F_FORMS))
def test_f_harmonics_match_cartesian_forms(m):
    rng = np.random.default_rng(5)
    th = np.arccos(rng.uniform(-1, 1, 400))
    ph = rng.uniform(0, 2 * np.pi, 400)
    x, y, z = np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)
    poly = F_FORMS[m](x, y, z)
    keep = np.abs(poly) > 0.05
    ratio = real_harmonic(3, m, th, ph)[keep] / poly[keep]
    assert ratio.min() > 0
    assert np.allclose(ratio, ratio[0], rtol=1e-9)


def test_radial_nodes_known_values():
    assert radial_nodes(1, 0) == []
    assert radial_nodes(2, 0) == pytest.approx([2.0], abs=2e-3)
    assert radial_nodes(3, 0) == pytest.approx([1.9019, 7.0981], abs=2e-3)


def test_radial_nodes_count_is_n_minus_l_minus_1():
    for n in range(1, 5):
        for l in range(n):
            assert len(radial_nodes(n, l)) == n - l - 1


def test_main_writes_files_and_radial_profile(tmp_path):
    main(["--out", str(tmp_path), "--nmax", "1", "--points", "500"])
    idx = json.loads((tmp_path / "index.json").read_text(encoding="utf-8"))
    assert len(idx["orbitals"]) == 1
    assert (tmp_path / idx["orbitals"][0]["file"]).stat().st_size == 7 * 500
    prof = idx["radial"]["1,0"]
    assert max(prof["p"]) == 1.0 and prof["nodes"] == [] and len(prof["r"]) == len(prof["p"])
