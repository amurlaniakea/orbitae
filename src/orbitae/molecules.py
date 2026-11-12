# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Moléculas pequeñas de capa cerrada con Hartree-Fock restringido (PySCF).

Dependencia opcional: pip install -e ".[hf]". Base def2-TZVP y geometrías experimentales
aproximadas FIJAS (sin optimizar). Se muestran los orbitales de valencia ocupados y el LUMO;
los orbitales internos no se muestrean porque la rejilla de 0,2 bohr no los resuelve.
Los orbitales virtuales de HF no describen estados excitados: son solo una referencia.

Salida para el visor: molecules.json y un .bin por orbital (N posiciones int16 normalizadas
por la extensión de la molécula, seguidas de N bytes de signo; 1 si psi < 0).
También exporta orbitales en formato .cube (VMD, Avogadro, PyMOL).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

BASIS = "def2-tzvp"
A0_ANGSTROM = 0.529177210903
HARTREE_EV = 27.211386245988
GRID_STEP = 0.2   # bohr
GRID_MARGIN = 6.0  # bohr
ENERGY_REF = {"H2O": -76.0675, "NH3": -56.2252, "CO": -112.7909, "N2": -108.9964, "CO2": -187.7248}
DIPOLE_EXP = {"H2O": 1.855, "NH3": 1.471, "CO": 0.122, "N2": 0.0, "CO2": 0.0}  # debye
KEYS = ["H2O", "NH3", "CO", "N2", "CO2"]
INFO = {
    "H2O": ("H₂O", "Agua", "Disolvente de la vida: la geometría doblada y el par solitario del oxígeno explican los enlaces de hidrógeno del agua consigo misma y con las biomoléculas."),
    "NH3": ("NH₃", "Amoníaco", "Se produce en el metabolismo de las proteínas y se elimina en el ciclo de la urea; el par solitario del nitrógeno le da carácter básico."),
    "CO": ("CO", "Monóxido de carbono", "Gas tóxico: su par solitario σ sobre el carbono se une al hierro del grupo hemo con más afinidad que el O₂ y reduce el transporte de oxígeno de la hemoglobina."),
    "N2": ("N₂", "Nitrógeno molecular", "Principal componente del aire (78 %); casi inerte químicamente por su triple enlace."),
    "CO2": ("CO₂", "Dióxido de carbono", "Producto final del metabolismo celular; en sangre se transporta sobre todo como bicarbonato y se elimina por los pulmones."),
}
NOTES = {
    ("H2O", "HOMO"): "Par solitario del oxígeno perpendicular al plano de la molécula; la densidad se anula en ese plano.",
    ("NH3", "HOMO"): "Par solitario del nitrógeno: es el que acepta un protón (NH₃ + H⁺ → NH₄⁺).",
    ("CO", "HOMO"): "Par solitario σ concentrado en el carbono: es el que se une al hierro del grupo hemo.",
    ("CO", "LUMO"): "Par degenerado de orbitales π* antienlazantes; reciben densidad del metal (retrodonación) cuando el CO se une al hierro.",
    ("N2", "HOMO"): "Par degenerado de orbitales π enlazantes: dos de los tres enlaces del triple enlace N≡N.",
    ("N2", "LUMO"): "Par degenerado de orbitales π* antienlazantes.",
    ("CO2", "HOMO"): "Par degenerado de orbitales π no enlazantes localizados en los oxígenos.",
}


def geometries() -> dict:
    """Geometrías experimentales aproximadas en angstrom: {clave: [(símbolo, (x, y, z)), ...]}."""
    g = {}
    r, t = 0.9572, np.radians(104.52 / 2)
    g["H2O"] = [("O", (0, 0, 0)), ("H", (r * np.sin(t), 0, r * np.cos(t))), ("H", (-r * np.sin(t), 0, r * np.cos(t)))]
    r = 1.012
    h = r * np.sqrt((np.cos(np.radians(106.7)) + 0.5) / 1.5)
    rho = np.sqrt(r * r - h * h)
    g["NH3"] = [("N", (0, 0, 0))] + [("H", (rho * np.cos(a), rho * np.sin(a), -h)) for a in (0, 2 * np.pi / 3, 4 * np.pi / 3)]
    g["CO"] = [("C", (0, 0, 0)), ("O", (0, 0, 1.128))]
    g["N2"] = [("N", (0, 0, -0.549)), ("N", (0, 0, 0.549))]
    g["CO2"] = [("O", (0, 0, -1.162)), ("C", (0, 0, 0)), ("O", (0, 0, 1.162))]
    return g


def solve(key: str):
    from pyscf import gto, scf  # importación diferida: dependencia opcional

    mol = gto.M(atom=geometries()[key], basis=BASIS, unit="Angstrom", symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    if not mf.converged:
        raise RuntimeError(f"RHF no convergió para {key}")
    return mol, mf


def select_orbitals(mol, mf) -> list:
    """Orbitales de valencia ocupados y el nivel LUMO, agrupados por degeneración."""
    from pyscf import symm

    e, nocc = mf.mo_energy, mol.nelectron // 2
    ncore = sum(1 for z in mol.atom_charges() if z > 1)
    groups = []
    for i in range(ncore, len(e)):
        if groups and abs(e[i] - e[groups[-1][-1]]) < 1e-4:
            groups[-1].append(i)
        else:
            groups.append([i])
    occ = [g for g in groups if g[-1] < nocc][::-1]
    virt = [g for g in groups if g[0] >= nocc][:1]
    irreps = symm.addons.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mf.mo_coeff)
    out = []
    for k, g in enumerate(occ):
        out += _entries(g, "HOMO" if k == 0 else f"HOMO−{k}", True, e, irreps)
    for g in virt:
        out += _entries(g, "LUMO", False, e, irreps)
    return out


def _entries(g, base, occupied, e, irreps):
    res = []
    for j, i in enumerate(g):
        res.append({"index": int(i), "base_label": base, "occupied": occupied,
                    "label": base if len(g) == 1 else f"{base} ({j + 1}/{len(g)})",
                    "energy_ev": float(e[i] * HARTREE_EV), "irrep": str(irreps[i])})
    return res


def make_grid(mol, step: float = GRID_STEP, margin: float = GRID_MARGIN):
    c = mol.atom_coords()
    lo, hi = c.min(0) - margin, c.max(0) + margin
    ax = [np.arange(a, b + step / 2, step) for a, b in zip(lo, hi)]
    X = np.stack(np.meshgrid(*ax, indexing="ij"), -1).reshape(-1, 3)
    return X, step**3


def orbital_values(mol, C, X, chunk: int = 40000):
    """psi(X) para las columnas de C (orbitales moleculares), en trozos para ahorrar memoria."""
    from pyscf.dft import numint

    C = np.asarray(C).reshape(C.shape[0], -1)
    return np.concatenate([numint.eval_ao(mol, X[i:i + chunk]) @ C for i in range(0, len(X), chunk)])


def sample_orbital_points(mol, C, X, psi, step, size, rng):
    """Puntos distribuidos según psi^2: celdas de la rejilla por probabilidad, con ruido dentro de la celda."""
    w = psi**2
    idx = rng.choice(len(X), size=size, p=w / w.sum())
    pts = X[idx] + (rng.random((size, 3)) - 0.5) * step
    sign = (orbital_values(mol, C, pts)[:, 0] < 0).astype(np.uint8)
    return pts, sign


def export_molecule(key: str, out: Path, points: int = 20000, seed: int = 0) -> dict:
    mol, mf = solve(key)
    sel = select_orbitals(mol, mf)
    X, dV = make_grid(mol)
    cols = [s["index"] for s in sel]
    P = orbital_values(mol, mf.mo_coeff[:, cols], X)
    rng = np.random.default_rng(seed)
    center = mol.atom_coords().mean(0)
    clouds = []
    for j, s in enumerate(sel):
        pts, sign = sample_orbital_points(mol, mf.mo_coeff[:, [s["index"]]], X, P[:, j], dV ** (1 / 3), points, rng)
        clouds.append((pts - center, sign))
        s["norm_grid"] = float(np.sum(P[:, j] ** 2) * dV)
    extent = float(max(np.abs(p).max() for p, _ in clouds)) * 1.001
    for j, (s, (p, sign)) in enumerate(zip(sel, clouds)):
        s["file"] = f"mol_{key}_{j}.bin"
        s["count"] = points
        (out / s["file"]).write_bytes(np.round(p / extent * 32767).astype("<i2").tobytes() + sign.tobytes())
        s["note"] = NOTES.get((key, s["base_label"]))
        del s["index"]
    formula, name, bio = INFO[key]
    return {
        "key": key, "formula": formula, "name": name, "bio": bio, "basis": BASIS,
        "energy_hf": float(mf.e_tot), "energy_ref": ENERGY_REF[key],
        "dipole_hf_debye": float(np.linalg.norm(mf.dip_moment(verbose=0))), "dipole_exp_debye": DIPOLE_EXP[key],
        "extent_a0": extent,
        "atoms": [{"symbol": mol.atom_symbol(i), "xyz_a0": (mol.atom_coords()[i] - center).tolist()} for i in range(mol.natm)],
        "orbitals": sel,
    }


def export_cube(key: str, label: str, path, n: int = 80):
    """Escribe un orbital (por etiqueta, p. ej. HOMO o LUMO (1/2)) en formato .cube."""
    from pyscf.tools import cubegen

    mol, mf = solve(key)
    sel = select_orbitals(mol, mf)
    s = next((x for x in sel if x["label"] == label or x["base_label"] == label), None)
    if s is None:
        raise ValueError(f"Orbital no disponible: {label}. Opciones: {[x['label'] for x in sel]}")
    cubegen.orbital(mol, str(path), mf.mo_coeff[:, s["index"]], nx=n, ny=n, nz=n, margin=GRID_MARGIN)
    return s


def read_cube(path):
    """Lee un .cube: devuelve (nº de átomos, volumen de voxel en bohr^3, valores)."""
    lines = Path(path).read_text().splitlines()
    natm = abs(int(lines[2].split()[0]))
    nx, ny, nz = (int(lines[i].split()[0]) for i in (3, 4, 5))
    vec = np.array([[float(v) for v in lines[i].split()[1:4]] for i in (3, 4, 5)])
    vals = np.array(" ".join(lines[6 + natm:]).split(), dtype=float).reshape(nx, ny, nz)
    return natm, abs(float(np.linalg.det(vec))), vals


def main(argv=None):
    ap = argparse.ArgumentParser(description="Moléculas con HF: datos para el visor y exportación .cube")
    ap.add_argument("--out", default="web/data")
    ap.add_argument("--points", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cube", choices=KEYS, help="exporta un orbital a .cube en lugar de generar datos")
    ap.add_argument("--orbital", default="HOMO")
    ap.add_argument("--cube-out", default=None)
    a = ap.parse_args(argv)
    if a.cube:
        path = a.cube_out or f"{a.cube}_{a.orbital}.cube".replace("−", "-").replace(" ", "")
        s = export_cube(a.cube, a.orbital, path)
        print(f"{path}: {a.cube} {s['label']} ({s['energy_ev']:.2f} eV)")
        return
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    mols = []
    for i, key in enumerate(KEYS):
        rec = export_molecule(key, out, a.points, a.seed + i)
        mols.append(rec)
        print(f"{rec['formula']:4s} E={rec['energy_hf']:.4f} dip={rec['dipole_hf_debye']:.2f} D "
              f"({len(rec['orbitals'])} orbitales: {', '.join(o['label'] for o in rec['orbitals'])})")
    (out / "molecules.json").write_text(
        json.dumps({"method": "RHF, geometrías experimentales aproximadas fijas", "basis": BASIS, "molecules": mols},
                   ensure_ascii=False),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
