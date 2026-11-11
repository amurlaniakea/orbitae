# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
import json

import numpy as np
import pytest

pytest.importorskip("pyscf")

from orbitae.atoms import configuration  # noqa: E402
from orbitae.hf_compare import (  # noqa: E402
    compare_atom, hf_atom, main, radial_density, unpaired_electrons,
)


@pytest.mark.parametrize("Z,esperado", [(1, 1), (2, 0), (6, 2), (7, 3), (8, 2), (24, 6), (26, 4), (29, 1), (36, 0)])
def test_unpaired_electrons(Z, esperado):
    assert unpaired_electrons(configuration(Z)) == esperado


@pytest.mark.parametrize("Z,ref,tol", [(1, -0.5, 0.005), (2, -2.86168, 0.01), (10, -128.54710, 0.03), (18, -526.81751, 0.05)])
def test_hf_energy_close_to_hf_limit(Z, ref, tol):
    _, mf, _ = hf_atom(Z)
    assert mf.e_tot == pytest.approx(ref, abs=tol)


@pytest.mark.parametrize("Z", [1, 6, 26])
def test_hf_density_integrates_to_Z(Z):
    mol, _, dm = hf_atom(Z)
    r = 40.0 * np.linspace(0, 1, 1500) ** 3 + 1e-8
    D = radial_density(mol, dm, r)
    N = np.sum((D[1:] + D[:-1]) * (r[1:] - r[:-1])) / 2.0
    assert N == pytest.approx(Z, rel=2e-3)


def test_hydrogen_matches_exact_radius_and_slater():
    rec = compare_atom(1)
    assert rec["mean_r_hf_a0"] == pytest.approx(1.5, abs=0.01)
    assert abs(rec["mean_r_error_pct"]) < 1.0


def test_slater_overestimates_3d_metals():
    assert 15.0 < compare_atom(26)["mean_r_error_pct"] < 35.0


def test_main_writes_profile(tmp_path):
    main(["--out", str(tmp_path), "--atoms", "1", "2", "--sin-iones"])
    idx = json.loads((tmp_path / "hf_compare.json").read_text(encoding="utf-8"))
    assert set(idx["atoms"]) == {"1", "2"}
    he = idx["atoms"]["2"]
    assert len(he["hf_total"]) == 160 and max(he["hf_total"]) > 0.5
    assert len(he["subshells"]["1s"]["hf_curve"]) == 160


@pytest.mark.parametrize("Z,q,ref,tol", [(3, 1, -7.2364, 0.01), (11, 1, -161.677, 0.03), (9, -1, -99.4594, 0.03)])
def test_ion_hf_energy(Z, q, ref, tol):
    _, mf, _ = hf_atom(Z, q)
    assert mf.e_tot == pytest.approx(ref, abs=tol)


@pytest.mark.parametrize("Z,q", [(11, 1), (17, -1), (26, 2), (26, 3), (29, 2)])
def test_ion_density_integrates_to_electrons(Z, q):
    mol, _, dm = hf_atom(Z, q)
    r = 40.0 * np.linspace(0, 1, 1500) ** 3 + 1e-8
    D = radial_density(mol, dm, r)
    N = np.sum((D[1:] + D[:-1]) * (r[1:] - r[:-1])) / 2.0
    assert N == pytest.approx(Z - q, rel=2e-3)


def test_cations_shrink_and_anions_swell_in_hf_and_slater():
    na, na_p = compare_atom(11), compare_atom(11, 1)
    cl, cl_m = compare_atom(17), compare_atom(17, -1)
    assert na_p["mean_r_hf_a0"] < na["mean_r_hf_a0"] and na_p["mean_r_sto_a0"] < na["mean_r_sto_a0"]
    assert cl_m["mean_r_hf_a0"] > cl["mean_r_hf_a0"] and cl_m["mean_r_sto_a0"] > cl["mean_r_sto_a0"]


def test_main_includes_ions(tmp_path):
    main(["--out", str(tmp_path), "--atoms", "11"])
    idx = json.loads((tmp_path / "hf_compare.json").read_text(encoding="utf-8"))
    assert len(idx["atoms"]) == 15 and "11+1" in idx["atoms"] and "17-1" in idx["atoms"]


@pytest.mark.parametrize("Z,q", [(25, 0), (26, 0), (28, 0), (6, 0), (36, 0), (26, 2), (17, -1)])
def test_hf_state_has_the_models_configuration(Z, q):
    rec = compare_atom(Z, q)
    assert rec["config_matches"] and rec["min_l_purity"] > 0.99


def test_hydrogen_subshell_radius():
    s = compare_atom(1)["subshells"]["1s"]
    assert s["mean_r_hf_a0"] == pytest.approx(1.5, abs=0.01)
    assert s["mean_r_sto_a0"] == pytest.approx(1.5)


def test_forcing_the_configuration_lowers_the_fe_error():
    # HF libre converge a 3d7 4s1 en el Fe y daría +25 %; con 3d6 4s2 forzado el error es menor.
    assert 15.0 < compare_atom(26)["mean_r_error_pct"] < 22.0


def test_valence_subshell_errors_have_known_signs():
    assert compare_atom(11)["subshells"]["3s"]["error_pct"] > 5.0
    assert compare_atom(6)["subshells"]["2p"]["error_pct"] < -5.0
    fe = compare_atom(26)["subshells"]
    assert fe["4s"]["error_pct"] > 15.0 and fe["3d"]["error_pct"] > 15.0


def test_subshell_curves_add_up_to_the_total_density():
    for Z in (6, 26):
        rec = compare_atom(Z)
        total = np.sum([v["hf_curve"] for v in rec["subshells"].values()], axis=0)
        assert np.allclose(total, rec["hf_total"], atol=3e-3)
