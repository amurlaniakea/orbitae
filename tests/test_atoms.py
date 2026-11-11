# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
import json

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import gamma

from orbitae.atoms import (
    NAMES, SYMBOLS, atom_record, configuration, export_angular_pool, main,
    sample_sto_radial, slater_zeff, sto_radial,
)


def decode(payload, n):
    q = np.frombuffer(payload[: 6 * n], dtype="<i2").reshape(n, 3) / 32767.0
    return q, np.frombuffer(payload[6 * n:], dtype=np.uint8)


def test_names_and_symbols():
    assert len(NAMES) == len(SYMBOLS) == 36


def test_configuration_totals():
    for Z in range(1, 37):
        assert sum(configuration(Z).values()) == Z


def test_configuration_exceptions():
    assert (configuration(24)[(4, 0)], configuration(24)[(3, 2)]) == (1, 5)
    assert (configuration(29)[(4, 0)], configuration(29)[(3, 2)]) == (1, 10)
    assert (configuration(26)[(4, 0)], configuration(26)[(3, 2)]) == (2, 6)


@pytest.mark.parametrize("Z,n,l,esperado", [
    (1, 1, 0, 1.0), (2, 1, 0, 1.70), (6, 2, 1, 3.25), (11, 3, 0, 2.20),
    (18, 3, 1, 6.75), (26, 3, 2, 6.25), (26, 4, 0, 3.75),
])
def test_slater_known_values(Z, n, l, esperado):
    assert slater_zeff(Z, configuration(Z), n, l) == pytest.approx(esperado, abs=1e-9)


def test_sto_normalization_and_gamma_form():
    ns, zeta = 3.7, 1.1
    val = quad(lambda r: r**2 * sto_radial(r, ns, zeta) ** 2, 0, np.inf)[0]
    assert val == pytest.approx(1.0, abs=1e-8)
    r = np.linspace(0.1, 15, 50)
    assert np.allclose(r**2 * sto_radial(r, ns, zeta) ** 2, gamma.pdf(r, a=2 * ns + 1, scale=1 / (2 * zeta)))


def test_sto_sampling_mean():
    x = sample_sto_radial(3.7, 1.1, 200_000, seed=4)
    assert x.mean() == pytest.approx((2 * 3.7 + 1) / (2 * 1.1), rel=0.01)


def test_atom_record_extent_encloses_electrons():
    for Z in (1, 6, 26, 36):
        a = atom_record(Z)
        inside = sum(
            s["occ"] * gamma.cdf(a["extent_a0"], a=2 * s["nstar"] + 1, scale=1 / (2 * s["zeta"]))
            for s in a["subshells"]
        )
        assert inside / Z >= 0.994
        assert max(a["radial"]["total"]) == 1.0


def test_angular_pool_px_signs_and_norm():
    q, s = decode(export_angular_pool(1, 1, 3000, 7), 3000)
    assert np.allclose(np.linalg.norm(q, axis=1), 1.0, atol=1e-3)
    mask = np.abs(q[:, 0]) > 0.01
    assert np.all(s[mask] == (q[mask, 0] < 0))


def test_angular_pool_s_is_uniform_and_single_phase():
    q, s = decode(export_angular_pool(0, 0, 20000, 2), 20000)
    assert s.sum() == 0
    assert np.abs(q.mean(axis=0)).max() < 0.02


def test_main_writes_everything(tmp_path):
    main(["--out", str(tmp_path), "--points", "300"])
    idx = json.loads((tmp_path / "atoms.json").read_text(encoding="utf-8"))
    assert len(idx["atoms"]) == 36 and idx["pools"]["n"] == 300
    assert len(idx["pools"]["ang"]) == 9 and len(idx["pools"]["rad"]) == 4
    for p in idx["pools"]["ang"]:
        assert (tmp_path / p["file"]).stat().st_size == 7 * 300
    for p in idx["pools"]["rad"]:
        assert (tmp_path / p["file"]).stat().st_size == 4 * 300
    assert all(sum(s["occ"] for s in a["subshells"]) == a["Z"] for a in idx["atoms"])
