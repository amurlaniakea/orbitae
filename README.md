# Orbitae

**Orbitales atómicos con rigor verificable.** Del hidrógeno exacto a la visualización
educativa, con una suite de validación numérica que comprueba cada resultado.

> Herramienta educativa y científica. No es un dispositivo médico ni una guía clínica.

## Qué hace (v0.1)

- Orbitales hidrogenoides exactos (n, l, m) para cualquier carga nuclear Z.
- Armónicos esféricos **reales** (p_x, p_y, p_z, d_xy, ...) con la convención de los libros,
  y complejos (m = ±1 da densidad toroidal).
- Distribución radial y radios de corte (90 % / 95 %).
- Muestreo de nubes electrónicas (x, y, z) según |ψ|².
- 20 pruebas de validación: normalización, ortogonalidad, ecuación de Schrödinger radial,
  energías E_n = -Z²/(2n²) hartree, forma de los orbitales y momentos del muestreo.

## Instalación

    python3 -m venv .venv && source .venv/bin/activate
    pip install -e ".[dev]"
    pytest -v

Requiere Python 3.10 o superior y SciPy 1.15 o superior (usa `scipy.special.sph_harm_y`).

## Ejemplo

    from orbitae import sample_orbital, radial_cutoff

    puntos = sample_orbital(n=2, l=1, m=1, size=50_000, seed=0)   # 2p_x, forma (N, 3)
    r95 = radial_cutoff(2, 1, p=0.95)                              # radio en bohr (a0)

## Hoja de ruta

Ver `docs/SDD.md`: átomos 1-36, visor web, iones y módulo de biología y medicina,
moléculas pequeñas, exportación `.cube`.

## Licencia

Copyright (C) Pedro Sordo Martínez. AGPL-3.0-or-later (ver `LICENSE`).

## Visor web

    python -m orbitae.export --out web/data
    python3 -m http.server -d web 8000

Abre http://localhost:8000. El visor es WebGL puro, sin librerías de terceros, y solo
se conecta a su propio servidor (política CSP).
