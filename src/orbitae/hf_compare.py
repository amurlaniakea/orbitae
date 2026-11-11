# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Compara el modelo de Slater (atoms.py) con Hartree-Fock no restringido (PySCF).

Dependencia opcional: pip install -e ".[hf]". Método de referencia: UHF con la
multiplicidad de Hund de la configuración, base def2-TZVP y densidad promediada
esféricamente por cuadratura angular. Es la referencia HF de este repositorio, no
un valor experimental ni exacto: HF no incluye correlación electrónica.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .atoms import IONS, SYMBOLS, _pdf, atom_record, ion_configuration

BASIS = "def2-tzvp"
BASIS_ION = "def2-tzvpd"  # con funciones difusas, necesarias para aniones
_X, _W = np.polynomial.legendre.leggauss(12)
_CT, _PH = np.meshgrid(_X, np.arange(16) * 2.0 * np.pi / 16, indexing="ij")
_ST = np.sqrt(1.0 - _CT**2)
_DIRS = np.c_[(_ST * np.cos(_PH)).ravel(), (_ST * np.sin(_PH)).ravel(), _CT.ravel()]
_WD = np.repeat(_W, 16) / 32.0  # pesos del promedio angular (suman 1)


def unpaired_electrons(occ: dict) -> int:
    """Electrones desapareados por la regla de Hund para una ocupación {(n, l): e}."""
    u = 0
    for (_, l), o in occ.items():
        c = 2 * l + 1
        u += o if o <= c else 2 * c - o
    return u


def hf_atom(Z: int, charge: int = 0, basis: str | None = None):
    """Devuelve (mol, mf, matriz densidad total) de un cálculo UHF convergido."""
    from pyscf import gto, scf  # importación diferida: dependencia opcional

    occ = ion_configuration(Z, charge)
    basis = basis or (BASIS if charge == 0 else BASIS_ION)
    mol = gto.M(atom=f"{SYMBOLS[Z - 1]} 0 0 0", basis=basis, charge=charge, spin=unpaired_electrons(occ), verbose=0)
    mf = scf.UHF(mol)
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        mf = mf.newton()
        mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"UHF no convergió para Z = {Z}")
    dm = mf.make_rdm1()
    return mol, mf, dm[0] + dm[1]


def radial_density(mol, dm, r, chunk: int = 150):
    """D(r) = 4 pi r^2 <rho>, con <rho> el promedio esférico de la densidad."""
    from pyscf.dft import numint

    r = np.asarray(r, dtype=float)
    out = []
    for i in range(0, len(r), chunk):
        rr = r[i:i + chunk]
        pts = (rr[:, None, None] * _DIRS[None]).reshape(-1, 3)
        rho = numint.eval_rho(mol, numint.eval_ao(mol, pts), dm).reshape(len(rr), -1)
        out.append(rho @ _WD)
    return 4.0 * np.pi * r**2 * np.concatenate(out)


def _trapz(y, x) -> float:
    return float(np.sum((y[1:] + y[:-1]) * (x[1:] - x[:-1])) / 2.0)


def compare_atom(Z: int, charge: int = 0, npts: int = 1500) -> dict:
    mol, mf, dm = hf_atom(Z, charge)
    r = 40.0 * np.linspace(0.0, 1.0, npts) ** 3 + 1e-8
    D = radial_density(mol, dm, r)
    N = _trapz(D, r)
    rec = atom_record(Z, charge)
    ne = Z - charge
    subs = rec["subshells"]
    Dsto = sum(s["occ"] * _pdf(r, s["nstar"], s["zeta"]) for s in subs)

    def cum(y):
        return np.concatenate([[0.0], np.cumsum((y[1:] + y[:-1]) * (r[1:] - r[:-1]) / 2.0)])

    ch, cs = cum(D) / N, cum(Dsto) / ne
    mean_hf = _trapz(r * D, r) / N
    mean_sto = sum(s["occ"] * (2 * s["nstar"] + 1) / (2 * s["zeta"]) for s in subs) / ne
    rg = np.array(rec["radial"]["r"])  # misma rejilla y normalización que atoms.json
    top = sum(s["occ"] * _pdf(rg, s["nstar"], s["zeta"]) for s in subs).max()
    Dg = radial_density(mol, dm, np.maximum(rg, 1e-8))
    return {
        "Z": Z, "key": rec["key"], "charge": charge, "symbol": rec["symbol"],
        "basis": BASIS if charge == 0 else BASIS_ION, "energy_hf": float(mf.e_tot),
        "n_electrons_numeric": N,
        "mean_r_hf_a0": mean_hf, "mean_r_sto_a0": mean_sto,
        "mean_r_error_pct": 100.0 * (mean_sto - mean_hf) / mean_hf,
        "r50_hf_a0": float(np.interp(0.5, ch, r)), "r50_sto_a0": float(np.interp(0.5, cs, r)),
        "r95_hf_a0": float(np.interp(0.95, ch, r)), "r95_sto_a0": float(np.interp(0.95, cs, r)),
        "l1_distance": 0.5 * _trapz(np.abs(D - Dsto), r) / ne,
        "hf_total": np.round(Dg / top, 4).tolist(),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description="Compara Slater con Hartree-Fock (PySCF)")
    ap.add_argument("--out", default="web/data")
    ap.add_argument("--atoms", type=int, nargs="*", default=list(range(1, 37)))
    ap.add_argument("--sin-iones", action="store_true", help="no calcula los iones")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    res = {}
    for Z in a.atoms:
        rec = compare_atom(Z)
        res[str(Z)] = rec
        print(f"{Z:2d} {rec['symbol']:2s} <r> HF={rec['mean_r_hf_a0']:.3f} Slater={rec['mean_r_sto_a0']:.3f} "
              f"({rec['mean_r_error_pct']:+.1f} %) L1={rec['l1_distance']:.3f}")
    if not a.sin_iones:
        for Z, q in IONS:
            rec = compare_atom(Z, q)
            res[rec["key"]] = rec
            print(f"{rec['symbol']:5s} <r> HF={rec['mean_r_hf_a0']:.3f} Slater={rec['mean_r_sto_a0']:.3f} "
                  f"({rec['mean_r_error_pct']:+.1f} %) L1={rec['l1_distance']:.3f}")
    (out / "hf_compare.json").write_text(
        json.dumps({"basis": BASIS, "method": "UHF (spin de Hund), promedio esférico", "atoms": res},
                   ensure_ascii=False),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
