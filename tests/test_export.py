# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
import json

import numpy as np

from orbitae.export import catalog, export_orbital, main, orbital_label


def decode(payload, count, scale):
    q = np.frombuffer(payload[: 6 * count], dtype="<i2").reshape(count, 3).astype(float) * scale
    return q, np.frombuffer(payload[6 * count:], dtype=np.uint8)


def test_labels_and_catalog():
    assert orbital_label(1, 0, 0) == "1s"
    assert orbital_label(2, 1, 1) == "2p_x"
    assert orbital_label(3, 2, -2) == "3d_xy"
    assert len(catalog(3)) == 14


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


def test_main_writes_files(tmp_path):
    main(["--out", str(tmp_path), "--nmax", "1", "--points", "500"])
    idx = json.loads((tmp_path / "index.json").read_text(encoding="utf-8"))
    assert len(idx["orbitals"]) == 1
    assert (tmp_path / idx["orbitals"][0]["file"]).stat().st_size == 7 * 500
