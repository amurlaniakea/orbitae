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
    main(["--out", str(tmp_path), "--atoms", "1", "2"])
    idx = json.loads((tmp_path / "hf_compare.json").read_text(encoding="utf-8"))
    assert set(idx["atoms"]) == {"1", "2"}
    he = idx["atoms"]["2"]
    assert len(he["hf_total"]) == 160 and max(he["hf_total"]) > 0.5
