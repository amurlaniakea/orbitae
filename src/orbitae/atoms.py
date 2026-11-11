# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Átomos de Z = 1 a 36: reglas de Slater y orbitales tipo Slater (STO).

Aproximación de campo central y didáctica: R(r) = N r^(n*-1) exp(-zeta r), con
zeta = Z_eff / n*. Sin nodos radiales ni correlación electrónica. La densidad
radial r^2 R^2 es una Gamma(2n*+1, escala 1/(2 zeta)), lo que permite muestrear
y validar de forma exacta. La nube total es el promedio esférico de las subcapas.

Salida para el visor: atoms.json, reservas angulares (ang*.bin: direcciones
unitarias int16 + signo) y reservas radiales (rad*.bin: Gamma(2n*+1, 1) en float32).
El visor combina ambas: r = g / (2 zeta).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import gammaln
from scipy.stats import gamma as _gamma

from .export import orbital_label
from .hydrogenic import real_harmonic, sample_orbital

NSTAR = {1: 1.0, 2: 2.0, 3: 3.0, 4: 3.7}
_ORDER = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 0), (3, 2), (4, 1)]
_CAP = {0: 2, 1: 6, 2: 10}
SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
           "Ga Ge As Se Br Kr").split()
NAMES = ("Hidrógeno Helio Litio Berilio Boro Carbono Nitrógeno Oxígeno Flúor Neón Sodio Magnesio "
         "Aluminio Silicio Fósforo Azufre Cloro Argón Potasio Calcio Escandio Titanio Vanadio Cromo "
         "Manganeso Hierro Cobalto Níquel Cobre Zinc Galio Germanio Arsénico Selenio Bromo Criptón").split()
BIO = {
    1: "Con el oxígeno forma el agua y casi toda la química orgánica; el gradiente de H⁺ impulsa la síntesis de ATP.",
    6: "Esqueleto de toda la química orgánica: forma cuatro enlaces estables.",
    7: "Presente en los aminoácidos y en las bases del ADN y el ARN.",
    8: "Aceptor final de electrones en la respiración celular; el O₂ se une al hierro de la hemoglobina.",
    11: "Principal catión extracelular; mantiene el potencial de membrana con la bomba Na⁺/K⁺.",
    12: "Se une al ATP y es cofactor de cientos de enzimas.",
    15: "Forma parte del ATP, el ADN, los fosfolípidos y el hueso.",
    16: "Presente en cisteína y metionina; forma puentes disulfuro.",
    17: "Principal anión extracelular; componente del ácido clorhídrico gástrico.",
    19: "Principal catión intracelular; fija el potencial de reposo y participa en el ritmo cardíaco.",
    20: "Contracción muscular, coagulación, señalización celular y hueso.",
    25: "Cofactor de la superóxido dismutasa mitocondrial.",
    26: "Centro de la hemoglobina y la mioglobina: une y transporta O₂.",
    27: "Centro de la vitamina B12.",
    29: "Cofactor de la citocromo c oxidasa y de la ceruloplasmina.",
    30: "Cofactor de cientos de enzimas, como la anhidrasa carbónica, y de las proteínas de dedos de zinc.",
    34: "Forma parte de las selenoproteínas, como la glutatión peroxidasa.",
}
_SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
_REMOVAL = [(4, 1), (4, 0), (3, 2), (3, 1), (3, 0), (2, 1), (2, 0), (1, 0)]  # orden de pérdida
IONS = [(3, 1), (11, 1), (19, 1), (12, 2), (20, 2), (9, -1), (17, -1), (25, 2), (26, 2),
        (26, 3), (27, 3), (29, 1), (29, 2), (30, 2)]
ION_BIO = {
    (3, 1): "El carbonato de litio se usa como estabilizador del ánimo; el Li⁺ es un ion pequeño que se comporta de forma parecida al Na⁺ y al Mg²⁺.",
    (11, 1): "Principal catión extracelular; su gradiente con el K⁺ genera el potencial de membrana y arrastra agua.",
    (19, 1): "Principal catión intracelular; fija el potencial de reposo y su concentración en sangre afecta al ritmo cardíaco.",
    (12, 2): "Se une al ATP (el sustrato de muchas enzimas es el complejo Mg-ATP) y es cofactor de cientos de enzimas.",
    (20, 2): "Señal intracelular que dispara la contracción muscular y la liberación de neurotransmisores; también coagulación y hueso.",
    (9, -1): "El fluoruro se incorpora al esmalte dental como fluorapatita y lo hace más resistente a los ácidos.",
    (17, -1): "Principal anión extracelular; acompaña al Na⁺ y participa en el equilibrio ácido-base y en el ácido clorhídrico gástrico.",
    (25, 2): "Cofactor de enzimas como la arginasa; la superóxido dismutasa mitocondrial alterna Mn³⁺ y Mn²⁺.",
    (26, 2): "Estado del hierro en el grupo hemo de la hemoglobina que une y transporta O₂.",
    (26, 3): "El hierro oxidado no une O₂: en la metahemoglobina el hemo tiene Fe³⁺. También es la forma que transporta la transferrina.",
    (27, 3): "Forma habitual del cobalto en la vitamina B12 (cobalamina), en el centro del anillo de corrina.",
    (29, 1): "Forma reducida del cobre; alterna con Cu²⁺ en proteínas redox como la citocromo c oxidasa.",
    (29, 2): "Forma oxidada del cobre; la superóxido dismutasa Cu/Zn alterna Cu²⁺ y Cu⁺.",
    (30, 2): "Con la subcapa 3d completa no hace redox: actúa como ácido de Lewis en la anhidrasa carbónica y estructura las proteínas de dedos de zinc.",
}


def configuration(Z: int) -> dict:
    """Ocupación del estado fundamental {(n, l): electrones}, con las excepciones Cr y Cu."""
    if not 1 <= Z <= 36:
        raise ValueError("Solo se admite Z = 1..36")
    occ, left = {}, Z
    for n, l in _ORDER:
        k = min(_CAP[l], left)
        if k:
            occ[(n, l)] = k
            left -= k
    if Z == 24:
        occ[(4, 0)], occ[(3, 2)] = 1, 5
    if Z == 29:
        occ[(4, 0)], occ[(3, 2)] = 1, 10
    return dict(sorted(occ.items()))


def ion_configuration(Z: int, charge: int) -> dict:
    """Ocupación de un ion. Los cationes pierden primero 4p y 4s y luego 3d, 3p...
    (Fe2+ es 3d6, no 3d4 4s2); los aniones ocupan el siguiente hueco del orden aufbau."""
    occ = dict(configuration(Z))
    for _ in range(max(charge, 0)):
        k = next((k for k in _REMOVAL if occ.get(k, 0) > 0), None)
        if k is None:
            raise ValueError("Carga mayor que el número de electrones")
        occ[k] -= 1
        if occ[k] == 0:
            del occ[k]
    for _ in range(max(-charge, 0)):
        k = next((k for k in _ORDER if occ.get(k, 0) < _CAP[k[1]]), None)
        if k is None:
            raise ValueError("No hay hueco para más electrones")
        occ[k] = occ.get(k, 0) + 1
    return dict(sorted(occ.items()))


def ion_symbol(symbol: str, charge: int) -> str:
    if charge == 0:
        return symbol
    mag = "" if abs(charge) == 1 else str(abs(charge)).translate(_SUP)
    return f"{symbol}{mag}{'⁺' if charge > 0 else '⁻'}"


def config_string(occ: dict) -> str:
    return " ".join(f"{n}{'spdf'[l]}{str(o).translate(_SUP)}" for (n, l), o in occ.items())


def _group(n: int, l: int):
    return (n, 0 if l < 2 else l - 1)  # (ns, np) = 0; nd = 1; nf = 2


def slater_zeff(Z: int, occ: dict, n: int, l: int) -> float:
    """Carga nuclear efectiva Z - S con las reglas de Slater (1930)."""
    g = _group(n, l)
    same = sum(o for (n2, l2), o in occ.items() if _group(n2, l2) == g) - 1
    S = (0.30 if n == 1 else 0.35) * same
    for (n2, l2), o in occ.items():
        g2 = _group(n2, l2)
        if g2 == g:
            continue
        if l < 2:
            S += 0.85 * o if n2 == n - 1 else (1.0 * o if n2 <= n - 2 else 0.0)
        elif g2 < g:
            S += 1.0 * o
    return Z - S


def sto_radial(r, nstar: float, zeta: float):
    """R(r) normalizada: integral de r^2 R^2 dr = 1."""
    r = np.asarray(r, dtype=float)
    log_norm = (nstar + 0.5) * np.log(2.0 * zeta) - 0.5 * gammaln(2.0 * nstar + 1.0)
    return np.exp(log_norm) * r ** (nstar - 1.0) * np.exp(-zeta * r)


def sample_sto_radial(nstar: float, zeta: float, size: int, seed: int | None = None):
    return np.random.default_rng(seed).gamma(2.0 * nstar + 1.0, 1.0 / (2.0 * zeta), size)


def _pdf(r, nstar, zeta):
    return _gamma.pdf(r, a=2.0 * nstar + 1.0, scale=1.0 / (2.0 * zeta))


def atom_record(Z: int, charge: int = 0, npts: int = 160) -> dict:
    occ = ion_configuration(Z, charge)
    subs = []
    for (n, l), o in occ.items():
        zeff = slater_zeff(Z, occ, n, l)
        ns = NSTAR[n]
        zeta = zeff / ns
        subs.append({
            "label": f"{n}{'spdf'[l]}", "n": n, "l": l, "occ": o, "nstar": ns,
            "zeff": round(zeff, 4), "zeta": round(zeta, 6),
            "extent_a0": float(_gamma.ppf(0.995, a=2.0 * ns + 1.0, scale=1.0 / (2.0 * zeta))),
            "orbitals": [{"m": m, "label": orbital_label(n, l, m)} for m in range(-l, l + 1)],
        })
    ext = max(s["extent_a0"] for s in subs)
    r = ext * np.linspace(0.0, 1.0, npts) ** 2  # rejilla densa cerca del núcleo
    shells = {s["label"]: s["occ"] * _pdf(r, s["nstar"], s["zeta"]) for s in subs}
    total = sum(shells.values())
    top = total.max()
    ref = ext if charge == 0 else max(s["extent_a0"] for s in atom_record(Z, 0, npts)["subshells"])
    return {
        "Z": Z, "key": str(Z) if charge == 0 else f"{Z}{charge:+d}", "charge": charge,
        "electrons": Z - charge, "symbol": ion_symbol(SYMBOLS[Z - 1], charge), "name": NAMES[Z - 1],
        "config": config_string(occ), "subshells": subs, "extent_a0": ext, "ref_extent_a0": ref,
        "bio": ION_BIO.get((Z, charge)) if charge else BIO.get(Z),
        "radial": {
            "r": np.round(r, 5).tolist(), "total": np.round(total / top, 4).tolist(),
            "shells": {k: np.round(v / top, 4).tolist() for k, v in shells.items()},
        },
    }


def angular_file(l: int, m: int) -> str:
    return f"ang{l}{'m' if m < 0 else 'p'}{abs(m)}.bin"


def export_angular_pool(l: int, m: int, size: int, seed: int) -> bytes:
    """Direcciones unitarias (int16 x32767) y signo del armónico real (1 si < 0)."""
    pts = sample_orbital(l + 1, l, m, size=size, seed=seed)
    u = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    theta = np.arccos(np.clip(u[:, 2], -1.0, 1.0))
    phi = np.arctan2(u[:, 1], u[:, 0])
    sign = (real_harmonic(l, m, theta, phi) < 0).astype(np.uint8)
    return np.round(u * 32767.0).astype("<i2").tobytes() + sign.tobytes()


def main(argv=None):
    ap = argparse.ArgumentParser(description="Exporta átomos 1-36 para el visor web")
    ap.add_argument("--out", default="web/data")
    ap.add_argument("--points", type=int, default=30000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ang, k = [], 0
    for l in range(3):
        for m in range(-l, l + 1):
            f = angular_file(l, m)
            (out / f).write_bytes(export_angular_pool(l, m, a.points, a.seed + 100 + k))
            ang.append({"l": l, "m": m, "file": f})
            k += 1
    rad = []
    for n in sorted(NSTAR):
        g = np.random.default_rng(a.seed + 200 + n).gamma(2.0 * NSTAR[n] + 1.0, 1.0, a.points)
        f = f"rad{n}.bin"
        (out / f).write_bytes(g.astype("<f4").tobytes())
        rad.append({"n": n, "file": f})
    atoms = [atom_record(Z) for Z in range(1, 37)]
    ions = [atom_record(Z, q) for Z, q in IONS]
    (out / "atoms.json").write_text(
        json.dumps({"pools": {"n": a.points, "ang": ang, "rad": rad}, "atoms": atoms, "ions": ions}, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"{len(atoms)} átomos y {len(ions)} iones exportados en {out}")


if __name__ == "__main__":
    main()
