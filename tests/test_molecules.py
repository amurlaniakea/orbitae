# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
import json

import numpy as np
import pytest

pytest.importorskip("pyscf")

from orbitae.molecules import (  # noqa: E402
    ALL_KEYS, ENERGY_REF, KEYS, OPEN_KEYS, SPIN, _spin_density, export_cube, main, make_grid, orbital_coeff,
    orbital_values, read_cube, sample_orbital_points, select_orbitals, solve,
)


@pytest.fixture(scope="module")
def solved():
    return {k: solve(k) for k in KEYS}


@pytest.mark.parametrize("key", KEYS)
def test_hf_energy_near_reference(solved, key):
    assert solved[key][1].e_tot == pytest.approx(ENERGY_REF[key], abs=0.02)


def test_orbital_selection_counts_and_degeneracy(solved):
    sel = {k: select_orbitals(*solved[k]) for k in KEYS}
    assert {k: len(v) for k, v in sel.items()} == {"H2O": 5, "NH3": 5, "CO": 7, "N2": 7, "CO2": 9}
    en = lambda k, base: [o["energy_ev"] for o in sel[k] if o["base_label"] == base]
    assert len(en("CO", "HOMO")) == 1 and len(en("N2", "HOMO")) == 2 and len(en("CO", "LUMO")) == 2
    assert abs(en("N2", "HOMO")[0] - en("N2", "HOMO")[1]) < 1e-6


@pytest.mark.parametrize("key", ["H2O", "CO"])
def test_valence_orbitals_normalized_and_orthogonal_on_grid(solved, key):
    mol, mf = solved[key]
    occ = [o for o in select_orbitals(mol, mf) if o["occupied"]]
    X, dV = make_grid(mol)
    P = orbital_values(mol, mf.mo_coeff[:, [o["index"] for o in occ]], X)
    G = P.T @ P * dV
    assert np.allclose(np.diag(G), 1.0, atol=0.02)
    assert np.abs(G - np.diag(np.diag(G))).max() < 0.02


def test_water_homo_has_a_nodal_plane_and_opposite_signs(solved):
    mol, mf = solved["H2O"]
    h = next(o for o in select_orbitals(mol, mf) if o["base_label"] == "HOMO")
    assert h["irrep"] == "B1"
    X, dV = make_grid(mol)
    psi = orbital_values(mol, mf.mo_coeff[:, [h["index"]]], X)[:, 0]
    assert np.abs(psi[np.abs(X[:, 1]) < 1e-9]).max() < 1e-8 < np.abs(psi).max()
    pts, sign = sample_orbital_points(mol, mf.mo_coeff[:, [h["index"]]], X, psi, dV ** (1 / 3), 20000, np.random.default_rng(1))
    m = np.abs(pts[:, 1]) > 0.15
    frac = np.mean(sign[m] == (pts[m, 1] < 0))
    assert max(frac, 1 - frac) > 0.97


def test_co_homo_is_a_sigma_lone_pair_on_carbon(solved):
    mol, mf = solved["CO"]
    h = next(o for o in select_orbitals(mol, mf) if o["base_label"] == "HOMO")
    X, _ = make_grid(mol)
    psi = orbital_values(mol, mf.mo_coeff[:, [h["index"]]], X)[:, 0]
    mid = mol.atom_coords()[:, 2].mean()
    assert h["irrep"] == "A1"
    assert np.sum(psi[X[:, 2] < mid] ** 2) / np.sum(psi**2) > 0.8  # el carbono está en z = 0, el oxígeno en z > 0


def test_sampling_matches_grid_moments(solved):
    mol, mf = solved["H2O"]
    h = next(o for o in select_orbitals(mol, mf) if o["base_label"] == "HOMO")
    X, dV = make_grid(mol)
    psi = orbital_values(mol, mf.mo_coeff[:, [h["index"]]], X)[:, 0]
    pts, _ = sample_orbital_points(mol, mf.mo_coeff[:, [h["index"]]], X, psi, dV ** (1 / 3), 20000, np.random.default_rng(2))
    exact = np.sum(psi**2 * np.sum(X**2, axis=1)) / np.sum(psi**2)
    assert np.mean(np.sum(pts**2, axis=1)) == pytest.approx(exact, rel=0.04)


def test_cube_export_roundtrip(tmp_path):
    path = tmp_path / "h2o.cube"
    main(["--cube", "H2O", "--orbital", "HOMO", "--cube-out", str(path)])
    natm, dV, vals = read_cube(path)
    assert natm == 3
    assert np.sum(vals**2) * dV == pytest.approx(1.0, abs=0.02)
    with pytest.raises(ValueError):
        export_cube("H2O", "HOMO−9", tmp_path / "x.cube")


def test_main_writes_data_files(tmp_path):
    main(["--out", str(tmp_path), "--points", "300"])
    idx = json.loads((tmp_path / "molecules.json").read_text(encoding="utf-8"))
    mols = {m["key"]: m for m in idx["molecules"]}
    assert list(mols) == ALL_KEYS and sum(len(m["orbitals"]) for m in mols.values()) == 63
    assert len(mols["H2O"]["atoms"]) == 3 and len(mols["CO2"]["atoms"]) == 3
    for m in mols.values():
        assert m["name_en"] and m["bio_en"] and m["extent_a0"] > 0
        assert all(abs(o["norm_grid"] - 1) < 0.04 for o in m["orbitals"] if o["occupied"] and o["norm_grid"])  # error máximo medido: 3,2 % (CO₂)
        assert all(bool(o["note"]) == bool(o["note_en"]) for o in m["orbitals"])
        for o in m["orbitals"]:
            assert (tmp_path / o["file"]).stat().st_size == 7 * 300


@pytest.fixture(scope="module")
def solved_open():
    return {k: solve(k) for k in OPEN_KEYS}


def test_open_shell_scf_quality(solved_open):
    for key, (mol, mf) in solved_open.items():
        s2, _ = mf.spin_square()
        ideal = SPIN[key] / 2 * (SPIN[key] / 2 + 1)
        assert mf.converged and abs(s2 - ideal) < 0.06  # contaminación de espín moderada
        assert mol.nelec[0] - mol.nelec[1] == SPIN[key]
    assert -149.75 < solved_open["O2"][1].e_tot < -149.60
    assert -129.35 < solved_open["NO"][1].e_tot < -129.25


def test_open_shell_orbital_selection(solved_open):
    o2, no = (select_orbitals(*solved_open[k]) for k in ("O2", "NO"))
    assert len(o2) == 15 and len(no) == 13
    count = lambda sel, g: sum(o["group"] == g for o in sel)
    assert (count(o2, "occ_a"), count(o2, "occ_b")) == (7, 5) and (count(no, "occ_a"), count(no, "occ_b")) == (6, 5)
    ha = [o["energy_ev"] for o in o2 if o["base_label"] == "α HOMO"]
    assert len(ha) == 2 and abs(ha[0] - ha[1]) < 1e-6  # el par π* con los dos electrones desapareados
    assert len([o for o in o2 if o["base_label"] == "β LUMO"]) == 2


def test_spin_density_and_populations(solved_open):
    for key, (mol, mf) in solved_open.items():
        X, dV = make_grid(mol)
        assert _spin_density(mol, mf, X).sum() * dV == pytest.approx(SPIN[key], abs=0.02)
    pop_o2 = solved_open["O2"][1].mulliken_spin_pop(verbose=0)[1]
    pop_no = solved_open["NO"][1].mulliken_spin_pop(verbose=0)[1]
    assert np.allclose(pop_o2, [1.0, 1.0], atol=0.05)
    assert pop_no[0] > 0.6 and pop_no[0] > pop_no[1]  # el electrón desapareado del NO está sobre todo en el N


@pytest.mark.parametrize("key", ["O2", "NO"])
def test_open_shell_valence_orbitals_normalized(solved_open, key):
    mol, mf = solved_open[key]
    occ = [o for o in select_orbitals(mol, mf) if o["occupied"]]
    X, dV = make_grid(mol)
    P = orbital_values(mol, np.hstack([orbital_coeff(mf, o) for o in occ]), X)
    assert np.allclose(np.sum(P**2, axis=0) * dV, 1.0, atol=0.03)


def test_spin_density_cube_roundtrip(tmp_path):
    path = tmp_path / "o2_spin.cube"
    main(["--cube", "O2", "--orbital", "SPIN", "--cube-out", str(path)])
    natm, dV, vals = read_cube(path)
    assert natm == 2 and np.sum(vals) * dV == pytest.approx(2.0, abs=0.05)
