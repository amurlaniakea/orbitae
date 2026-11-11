# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
import json

import pytest

from orbitae.atoms import IONS, atom_record, ion_configuration, ion_symbol, main, slater_zeff


def test_cations_lose_4s_before_3d():
    assert ion_configuration(26, 2) == {**{k: v for k, v in ion_configuration(26, 0).items() if k != (4, 0)}}
    assert ion_configuration(26, 2)[(3, 2)] == 6 and (4, 0) not in ion_configuration(26, 2)
    assert ion_configuration(26, 3)[(3, 2)] == 5
    assert ion_configuration(29, 1)[(3, 2)] == 10 and (4, 0) not in ion_configuration(29, 1)
    assert ion_configuration(29, 2)[(3, 2)] == 9
    assert ion_configuration(25, 2)[(3, 2)] == 5 and ion_configuration(30, 2)[(3, 2)] == 10


def test_isoelectronic_series_share_configuration():
    assert ion_configuration(19, 1) == ion_configuration(20, 2) == ion_configuration(17, -1)
    assert ion_configuration(11, 1) == ion_configuration(12, 2) == ion_configuration(9, -1)


def test_electron_counts_and_unique_keys():
    keys = []
    for Z, q in IONS:
        rec = atom_record(Z, q)
        assert sum(s["occ"] for s in rec["subshells"]) == Z - q == rec["electrons"]
        keys.append(rec["key"])
    assert len(set(keys)) == len(keys) == 14


def test_ion_symbols():
    assert ion_symbol("Na", 1) == "Na⁺" and ion_symbol("Fe", 3) == "Fe³⁺" and ion_symbol("Cl", -1) == "Cl⁻"


@pytest.mark.parametrize("Z,q,n,l,esperado", [(11, 1, 2, 1, 6.85), (26, 3, 3, 2, 6.60), (17, -1, 3, 1, 5.75)])
def test_slater_zeff_for_ions(Z, q, n, l, esperado):
    assert slater_zeff(Z, ion_configuration(Z, q), n, l) == pytest.approx(esperado, abs=1e-9)


def test_slater_extents_shrink_for_cations_and_swell_for_anions():
    assert atom_record(11, 1)["extent_a0"] < atom_record(11)["extent_a0"]
    assert atom_record(17, -1)["extent_a0"] > atom_record(17)["extent_a0"]
    assert atom_record(26, 3)["ref_extent_a0"] == atom_record(26)["extent_a0"]


def test_main_writes_ions(tmp_path):
    main(["--out", str(tmp_path), "--points", "300"])
    idx = json.loads((tmp_path / "atoms.json").read_text(encoding="utf-8"))
    assert len(idx["ions"]) == 14 and idx["atoms"][25]["key"] == "26"
