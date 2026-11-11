# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Compara el modelo de Slater (atoms.py) con Hartree-Fock no restringido (PySCF).

Dependencia opcional: pip install -e ".[hf]". Referencia: UHF con la multiplicidad de
Hund de la configuración, base def2-TZVP (def2-TZVPD en iones), densidad promediada
esféricamente por cuadratura angular. HF no incluye correlación electrónica.

Importante: se FUERZA en HF la misma ocupación por subcapa que usa el modelo (por ejemplo
3d6 4s2 en el Fe). Sin esa restricción, HF converge a 3d(n+1) 4s1 en Mn, Fe, Co y Ni, y la
comparación mezclaría el error del modelo con una configuración distinta.
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


def _ao_l(mol):
    return np.concatenate([[mol.bas_angular(i)] * ((2 * mol.bas_angular(i) + 1) * mol.bas_nctr(i))
                           for i in range(mol.nbas)])


def _l_weights(C, S, ao_l):
    """Peso de cada orbital en cada l (exacto: en un átomo los AO de distinto l son ortogonales)."""
    return np.stack([np.einsum("ij,ik,kj->j", C[ao_l == l], S[np.ix_(ao_l == l, ao_l == l)], C[ao_l == l])
                     for l in range(ao_l.max() + 1)], axis=1)


def _target_counts(occ: dict) -> dict:
    """{l: [orbitales alfa ocupados, orbitales beta ocupados]} con espín alto (Hund) por subcapa."""
    t = {}
    for (_, l), e in occ.items():
        a = min(e, 2 * l + 1)
        t.setdefault(l, [0, 0])
        t[l][0] += a
        t[l][1] += e - a
    return t


def hf_atom(Z: int, charge: int = 0, basis: str | None = None):
    """Devuelve (mol, mf, matriz densidad total) de un UHF convergido con la ocupación del modelo."""
    from pyscf import gto, scf  # importación diferida: dependencia opcional

    occ = ion_configuration(Z, charge)
    basis = basis or (BASIS if charge == 0 else BASIS_ION)
    mol = gto.M(atom=f"{SYMBOLS[Z - 1]} 0 0 0", basis=basis, charge=charge,
                spin=unpaired_electrons(occ), verbose=0)
    ao_l, S, target = _ao_l(mol), mol.intor("int1e_ovlp"), _target_counts(occ)
    mf = scf.UHF(mol)
    mf.max_cycle = 300

    def get_occ(mo_energy=None, mo_coeff=None):
        mo_energy = mf.mo_energy if mo_energy is None else mo_energy
        mo_coeff = mf.mo_coeff if mo_coeff is None else mo_coeff
        mo_occ = np.zeros(np.shape(mo_energy))
        for s in (0, 1):
            lab = _l_weights(mo_coeff[s], S, ao_l).argmax(1)
            order = np.argsort(mo_energy[s])
            for l, cnt in target.items():
                members = [j for j in order if lab[j] == l]
                mo_occ[s][members[:cnt[s]]] = 1.0
        return mo_occ

    mf.get_occ = get_occ
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"UHF no convergió para Z = {Z}, carga {charge}")
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


def subshell_densities(mol, mf, r, chunk: int = 150):
    """Densidad radial por subcapa (n, l) de HF, clasificando los orbitales por l y por energía.

    Devuelve ({(n, l): D(r)}, {(n, l): nº de electrones}, pureza mínima en l de los orbitales).
    """
    from pyscf.dft import numint

    S, ao_l = mol.intor("int1e_ovlp"), _ao_l(mol)
    r = np.asarray(r, dtype=float)
    groups, purity = {}, 1.0
    for s in (0, 1):
        idx = np.where(mf.mo_occ[s] > 0.5)[0]
        if len(idx) == 0:
            continue
        Co = mf.mo_coeff[s][:, idx]
        w = _l_weights(Co, S, ao_l)
        lab = w.argmax(1)
        purity = min(purity, float(w.max(1).min()))
        D = []
        for i in range(0, len(r), chunk):
            rr = r[i:i + chunk]
            pts = (rr[:, None, None] * _DIRS[None]).reshape(-1, 3)
            mo = numint.eval_ao(mol, pts) @ Co
            D.append(4.0 * np.pi * rr[:, None] ** 2 *
                     ((mo**2).reshape(len(rr), len(_DIRS), -1) * _WD[None, :, None]).sum(1))
        D = np.concatenate(D)
        order = np.argsort(mf.mo_energy[s][idx])
        for l in np.unique(lab):
            for k, j in enumerate([j for j in order if lab[j] == l]):
                groups.setdefault((int(l) + 1 + k // (2 * int(l) + 1), int(l)), []).append(D[:, j])
    return ({k: np.sum(v, axis=0) for k, v in groups.items()}, {k: len(v) for k, v in groups.items()}, purity)


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
    rg = np.maximum(np.array(rec["radial"]["r"]), 1e-8)  # misma rejilla y normalización que atoms.json
    top = sum(s["occ"] * _pdf(rg, s["nstar"], s["zeta"]) for s in subs).max()
    Dg = radial_density(mol, dm, rg)
    Ds, cnt, purity = subshell_densities(mol, mf, np.concatenate([r, rg]))
    sub_out = {}
    for s in subs:
        key = (s["n"], s["l"])
        d_fine, d_grid = Ds[key][:npts], Ds[key][npts:]
        mr = _trapz(r * d_fine, r) / _trapz(d_fine, r)
        sto = (2 * s["nstar"] + 1) / (2 * s["zeta"])
        sub_out[s["label"]] = {
            "electrons_hf": cnt[key], "mean_r_hf_a0": mr, "mean_r_sto_a0": sto,
            "error_pct": 100.0 * (sto - mr) / mr, "hf_curve": np.round(d_grid / top, 4).tolist(),
        }
    return {
        "Z": Z, "key": rec["key"], "charge": charge, "symbol": rec["symbol"],
        "basis": BASIS if charge == 0 else BASIS_ION, "energy_hf": float(mf.e_tot),
        "n_electrons_numeric": N, "min_l_purity": purity,
        "config_matches": all(cnt.get((s["n"], s["l"])) == s["occ"] for s in subs) and len(cnt) == len(subs),
        "mean_r_hf_a0": mean_hf, "mean_r_sto_a0": mean_sto,
        "mean_r_error_pct": 100.0 * (mean_sto - mean_hf) / mean_hf,
        "r50_hf_a0": float(np.interp(0.5, ch, r)), "r50_sto_a0": float(np.interp(0.5, cs, r)),
        "r95_hf_a0": float(np.interp(0.95, ch, r)), "r95_sto_a0": float(np.interp(0.95, cs, r)),
        "l1_distance": 0.5 * _trapz(np.abs(D - Dsto), r) / ne,
        "hf_total": np.round(Dg / top, 4).tolist(), "subshells": sub_out,
    }


def _line(rec) -> str:
    outer = max(rec["subshells"], key=lambda k: (int(k[0]), "spdf".index(k[1])))
    so = rec["subshells"][outer]
    return (f"{rec['symbol']:5s} <r> HF={rec['mean_r_hf_a0']:.3f} Slater={rec['mean_r_sto_a0']:.3f} "
            f"({rec['mean_r_error_pct']:+.1f} %) | capa {outer}: {so['error_pct']:+.1f} %"
            f"{'' if rec['config_matches'] else '  ¡CONFIGURACIÓN HF DISTINTA!'}")


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
        print(f"{Z:2d} {_line(rec)}")
    if not a.sin_iones:
        for Z, q in IONS:
            rec = compare_atom(Z, q)
            res[rec["key"]] = rec
            print(_line(rec))
    (out / "hf_compare.json").write_text(
        json.dumps({"basis": BASIS, "method": "UHF (spin de Hund, ocupación forzada), promedio esférico",
                    "atoms": res}, ensure_ascii=False),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
