# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Moléculas pequeñas con Hartree-Fock (PySCF): capa cerrada (RHF) y con espín (UHF: O2 y NO).

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
DIPOLE_EXP = {"H2O": 1.855, "NH3": 1.471, "CO": 0.122, "N2": 0.0, "CO2": 0.0, "O2": 0.0, "NO": 0.159}  # debye
KEYS = ["H2O", "NH3", "CO", "N2", "CO2"]  # capa cerrada, RHF
OPEN_KEYS = ["O2", "NO"]  # con espín, UHF con la multiplicidad de Hund
ALL_KEYS = KEYS + OPEN_KEYS
SPIN = {"O2": 2, "NO": 1}  # 2S = nº de electrones desapareados
INFO = {
    "H2O": ("H₂O", "Agua", "Disolvente de la vida: la geometría doblada y el par solitario del oxígeno explican los enlaces de hidrógeno del agua consigo misma y con las biomoléculas."),
    "NH3": ("NH₃", "Amoníaco", "Se produce en el metabolismo de las proteínas y se elimina en el ciclo de la urea; el par solitario del nitrógeno le da carácter básico."),
    "CO": ("CO", "Monóxido de carbono", "Gas tóxico: su par solitario σ sobre el carbono se une al hierro del grupo hemo con más afinidad que el O₂ y reduce el transporte de oxígeno de la hemoglobina."),
    "N2": ("N₂", "Nitrógeno molecular", "Principal componente del aire (78 %); casi inerte químicamente por su triple enlace."),
    "O2": ("O₂", "Oxígeno molecular", "Su estado fundamental es un triplete: tiene dos electrones desapareados, por eso es paramagnético y reacciona despacio con la materia orgánica salvo que haya catalizadores como el hierro de las enzimas. En la hemoglobina se une al Fe²⁺ del hemo."),
    "NO": ("NO", "Óxido nítrico", "Radical con un electrón desapareado que actúa como mensajero: relaja el músculo liso de los vasos sanguíneos (vasodilatación) activando la guanilato ciclasa, enzima con grupo hemo al que se une."),
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
    ("O2", "α HOMO"): "Par degenerado de orbitales π* antienlazantes (1πg), cada uno con un electrón α: son los dos electrones desapareados del triplete.",
    ("O2", "β LUMO"): "Par degenerado de orbitales π* antienlazantes vacío para el espín β.",
    ("O2", "SPIN"): "Densidad de espín (ρα − ρβ): los dos electrones desapareados se reparten por igual entre los dos oxígenos; las zonas de color opuesto son polarización de espín.",
    ("NO", "α HOMO"): "Orbital π* antienlazante que contiene el electrón desapareado, con más peso sobre el nitrógeno; su orientación alrededor del eje del enlace es arbitraria.",
    ("NO", "SPIN"): "Densidad de espín (ρα − ρβ): el electrón desapareado se localiza sobre todo en el nitrógeno; las zonas de color opuesto son polarización de espín.",
}

INFO_EN = {
    "H2O": ("Water", "Solvent of life: the bent geometry and the lone pair of oxygen explain the hydrogen bonds of water with itself and with biomolecules."),
    "NH3": ("Ammonia", "Produced in protein metabolism and eliminated through the urea cycle; the lone pair on nitrogen makes it basic."),
    "CO": ("Carbon monoxide", "Toxic gas: its σ lone pair on carbon binds to the iron of the heme group with higher affinity than O₂ and reduces oxygen transport by hemoglobin."),
    "N2": ("Molecular nitrogen", "Main component of air (78 %); almost chemically inert because of its triple bond."),
    "O2": ("Molecular oxygen", "Its ground state is a triplet: it has two unpaired electrons, so it is paramagnetic and reacts slowly with organic matter unless catalysts such as the iron of enzymes are present. In hemoglobin it binds to the Fe²⁺ of the heme."),
    "NO": ("Nitric oxide", "A radical with one unpaired electron that acts as a messenger: it relaxes the smooth muscle of blood vessels (vasodilation) by activating guanylate cyclase, a heme-containing enzyme to which it binds."),
    "CO2": ("Carbon dioxide", "End product of cellular metabolism; in blood it is carried mostly as bicarbonate and eliminated by the lungs."),
}
NOTES_EN = {
    ("H2O", "HOMO"): "Lone pair of oxygen perpendicular to the plane of the molecule; the density vanishes in that plane.",
    ("NH3", "HOMO"): "Lone pair of nitrogen: the one that accepts a proton (NH₃ + H⁺ → NH₄⁺).",
    ("CO", "HOMO"): "σ lone pair concentrated on carbon: the one that binds to the iron of the heme group.",
    ("CO", "LUMO"): "Degenerate pair of antibonding π* orbitals; they receive density from the metal (back-donation) when CO binds to iron.",
    ("N2", "HOMO"): "Degenerate pair of bonding π orbitals: two of the three bonds of the N≡N triple bond.",
    ("N2", "LUMO"): "Degenerate pair of antibonding π* orbitals.",
    ("CO2", "HOMO"): "Degenerate pair of nonbonding π orbitals localized on the oxygens.",
    ("O2", "α HOMO"): "Degenerate pair of antibonding π* orbitals (1πg), each with one α electron: they are the two unpaired electrons of the triplet.",
    ("O2", "β LUMO"): "Degenerate pair of antibonding π* orbitals, empty for β spin.",
    ("O2", "SPIN"): "Spin density (ρα − ρβ): the two unpaired electrons are shared equally between the two oxygens; regions of opposite color are spin polarization.",
    ("NO", "α HOMO"): "Antibonding π* orbital holding the unpaired electron, with more weight on nitrogen; its orientation around the bond axis is arbitrary.",
    ("NO", "SPIN"): "Spin density (ρα − ρβ): the unpaired electron is localized mostly on nitrogen; regions of opposite color are spin polarization.",
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
    g["O2"] = [("O", (0, 0, -0.60375)), ("O", (0, 0, 0.60375))]
    g["NO"] = [("N", (0, 0, 0)), ("O", (0, 0, 1.1508))]
    return g


def solve(key: str):
    from pyscf import gto, scf  # importación diferida: dependencia opcional

    # Con simetría; si el SCF no converge (el NO, con un electrón en un par π* degenerado), sin ella.
    for sym in (True, False):
        mol = gto.M(atom=geometries()[key], basis=BASIS, unit="Angstrom", symmetry=sym, spin=SPIN.get(key, 0), verbose=0)
        if key in SPIN:
            mf = scf.UHF(mol)
            mf.max_cycle = 100
            mf.kernel()
        else:
            mf = scf.RHF(mol).run()
        if mf.converged:
            return mol, mf
    raise RuntimeError(f"El SCF no convergió para {key}")


def _groups(idx, e):
    groups = []
    for i in idx:
        if groups and abs(e[i] - e[groups[-1][-1]]) < 1e-4:
            groups[-1].append(i)
        else:
            groups.append([i])
    return groups


def _entries(g, base, occupied, e, irreps, group, spin=None):
    res = []
    for j, i in enumerate(g):
        res.append({"index": int(i), "base_label": base, "occupied": occupied, "spin": spin, "group": group,
                    "label": base if len(g) == 1 else f"{base} ({j + 1}/{len(g)})",
                    "energy_ev": float(e[i] * HARTREE_EV), "irrep": str(irreps[i])})
    return res


def _irreps(mol, C):
    from pyscf import symm

    if not mol.symmetry:
        return [""] * C.shape[1]
    return symm.addons.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, C)


def select_orbitals(mol, mf) -> list:
    """Orbitales de valencia ocupados y el primer nivel vacío, agrupados por degeneración.
    En UHF (moléculas con espín) se eligen por separado para los espines α y β."""
    ncore = sum(1 for z in mol.atom_charges() if z > 1)
    out = []
    if np.ndim(mf.mo_energy) == 1:  # RHF
        e, nocc = mf.mo_energy, mol.nelectron // 2
        irreps = _irreps(mol, mf.mo_coeff)
        for k, g in enumerate(_groups(range(ncore, nocc), e)[::-1]):
            out += _entries(g, "HOMO" if k == 0 else f"HOMO−{k}", True, e, irreps, "occ")
        for g in _groups(range(nocc, len(e)), e)[:1]:
            out += _entries(g, "LUMO", False, e, irreps, "virt")
        return out
    for s, (tag, sym) in enumerate([("alpha", "α"), ("beta", "β")]):
        e, nocc = mf.mo_energy[s], mol.nelec[s]
        irreps = _irreps(mol, mf.mo_coeff[s])
        suf = "a" if s == 0 else "b"
        for k, g in enumerate(_groups(range(ncore, nocc), e)[::-1]):
            out += _entries(g, f"{sym} HOMO" if k == 0 else f"{sym} HOMO−{k}", True, e, irreps, f"occ_{suf}", tag)
        for g in _groups(range(nocc, len(e)), e)[:1]:
            out += _entries(g, f"{sym} LUMO", False, e, irreps, f"virt_{suf}", tag)
    return out


def orbital_coeff(mf, s):
    """Coeficientes (columna) del orbital descrito por la entrada s, en RHF o UHF."""
    C = mf.mo_coeff
    if s.get("spin") == "alpha":
        return C[0][:, [s["index"]]]
    if s.get("spin") == "beta":
        return C[1][:, [s["index"]]]
    return C[:, [s["index"]]]


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


def _spin_density(mol, mf, X):
    na, nb = mol.nelec
    return ((orbital_values(mol, mf.mo_coeff[0][:, :na], X) ** 2).sum(1)
            - (orbital_values(mol, mf.mo_coeff[1][:, :nb], X) ** 2).sum(1))


def export_molecule(key: str, out: Path, points: int = 20000, seed: int = 0) -> dict:
    mol, mf = solve(key)
    sel = select_orbitals(mol, mf)
    X, dV = make_grid(mol)
    step = dV ** (1 / 3)
    P = orbital_values(mol, np.hstack([orbital_coeff(mf, s) for s in sel]), X)
    rng = np.random.default_rng(seed)
    center = mol.atom_coords().mean(0)
    clouds, extra = [], {}
    for j, s in enumerate(sel):
        pts, sign = sample_orbital_points(mol, orbital_coeff(mf, s), X, P[:, j], step, points, rng)
        clouds.append((pts - center, sign))
        s["norm_grid"] = float(np.sum(P[:, j] ** 2) * dV)
    if key in SPIN:  # densidad de espín: se muestrea |rho_alfa - rho_beta| y se colorea por su signo
        rs = _spin_density(mol, mf, X)
        idx = rng.choice(len(X), size=points, p=np.abs(rs) / np.abs(rs).sum())
        pts = X[idx] + (rng.random((points, 3)) - 0.5) * step
        clouds.append((pts - center, (_spin_density(mol, mf, pts) < 0).astype(np.uint8)))
        sel.append({"index": -1, "base_label": "SPIN", "label": "ρα − ρβ", "occupied": True, "spin": None,
                    "group": "spin", "energy_ev": None, "irrep": "", "norm_grid": None,
                    "spin_integral": float(rs.sum() * dV), "abs_integral": float(np.abs(rs).sum() * dV)})
        s2, _ = mf.spin_square()
        extra = {"s2": float(s2), "s2_ideal": SPIN[key] / 2 * (SPIN[key] / 2 + 1),
                 "spin_pop": [float(v) for v in mf.mulliken_spin_pop(verbose=0)[1]]}
    extent = float(max(np.abs(p).max() for p, _ in clouds)) * 1.001
    for j, (s, (p, sign)) in enumerate(zip(sel, clouds)):
        s["file"] = f"mol_{key}_{j}.bin"
        s["count"] = points
        (out / s["file"]).write_bytes(np.round(p / extent * 32767).astype("<i2").tobytes() + sign.tobytes())
        s["note"] = NOTES.get((key, s["base_label"]))
        s["note_en"] = NOTES_EN.get((key, s["base_label"]))
        del s["index"]
    formula, name, bio = INFO[key]
    return {
        "key": key, "formula": formula, "name": name, "bio": bio, "basis": BASIS,
        "name_en": INFO_EN[key][0], "bio_en": INFO_EN[key][1], "open_shell": key in SPIN,
        "energy_hf": float(mf.e_tot), "energy_ref": ENERGY_REF.get(key),
        "dipole_hf_debye": float(np.linalg.norm(mf.dip_moment(verbose=0))), "dipole_exp_debye": DIPOLE_EXP[key],
        "extent_a0": extent,
        "atoms": [{"symbol": mol.atom_symbol(i), "xyz_a0": (mol.atom_coords()[i] - center).tolist()} for i in range(mol.natm)],
        "orbitals": sel, **extra,
    }


def export_cube(key: str, label: str, path, n: int = 80):
    """Escribe un orbital (por etiqueta, p. ej. HOMO o LUMO (1/2)) en formato .cube."""
    from pyscf.tools import cubegen

    mol, mf = solve(key)
    if label in ("SPIN", "ρα − ρβ") and key in SPIN:
        dma, dmb = mf.make_rdm1()
        cubegen.density(mol, str(path), dma - dmb, nx=n, ny=n, nz=n, margin=GRID_MARGIN)
        return {"label": "ρα − ρβ", "energy_ev": None}
    sel = select_orbitals(mol, mf)
    s = next((x for x in sel if x["label"] == label or x["base_label"] == label), None)
    if s is None:
        raise ValueError(f"Orbital no disponible: {label}. Opciones: {[x['label'] for x in sel]}")
    cubegen.orbital(mol, str(path), orbital_coeff(mf, s)[:, 0], nx=n, ny=n, nz=n, margin=GRID_MARGIN)
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
    ap.add_argument("--cube", choices=ALL_KEYS, help="exporta un orbital a .cube en lugar de generar datos")
    ap.add_argument("--orbital", default="HOMO")
    ap.add_argument("--cube-out", default=None)
    a = ap.parse_args(argv)
    if a.cube:
        path = a.cube_out or f"{a.cube}_{a.orbital}.cube".replace("−", "-").replace(" ", "")
        s = export_cube(a.cube, a.orbital, path)
        en = "" if s["energy_ev"] is None else f" ({s['energy_ev']:.2f} eV)"
        print(f"{path}: {a.cube} {s['label']}{en}")
        return
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    mols = []
    for i, key in enumerate(ALL_KEYS):
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
