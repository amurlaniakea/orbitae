# Orbitae

**Orbitales atómicos con rigor verificable.** Visualiza nubes electrónicas en 3D, del hidrógeno
exacto a los átomos 1-36, con una suite de pruebas numéricas que comprueba cada resultado.

[![tests](https://github.com/amurlaniakea/orbitae/actions/workflows/tests.yml/badge.svg)](https://github.com/amurlaniakea/orbitae/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24135885.svg)](https://doi.org/10.5281/zenodo.24135885)
**Demo en línea:** https://amurlaniakea.github.io/orbitae/

![Visor de Orbitae](docs/img/orbitae-visor.png)

> Herramienta educativa y científica. No es un dispositivo médico ni una guía clínica.

## Qué hace

- **Hidrógeno exacto** (modo "Hidrógeno exacto"): orbitales hasta n = 4 (s, p, d, f) con armónicos
  esféricos **reales** (p_x, d_xy, ...), fase en dos colores, energía, radio del 95 % y
  distribución radial con sus nodos.
- **Átomos 1-36** (modo "Átomos"): tabla periódica clicable, nube total y orbitales individuales de
  cada subcapa, con carga efectiva de las reglas de Slater y orbitales tipo Slater.
- **Contexto biomédico**: los elementos esenciales en biología humana van marcados y llevan una
  frase de relevancia (hierro en la hemoglobina, sodio y potasio en el potencial de membrana, etc.).
- **Iones de interés biomédico** (Na⁺, K⁺, Mg²⁺, Ca²⁺, Cl⁻, Fe²⁺, Fe³⁺, Cu⁺, Cu²⁺, Zn²⁺, entre otros),
  con su comparación frente a Hartree-Fock.
- **Visor WebGL puro**, sin librerías de terceros y con política CSP que solo permite conectar con
  su propio servidor.

## Límites (a propósito, explícitos)

- Los átomos usan una **aproximación didáctica**: reglas de Slater y orbitales sin nodos radiales,
  sin correlación electrónica. La nube total es el promedio esférico de las subcapas. El error
  frente a Hartree-Fock está medido (ver más abajo).
- Todavía no hay moléculas. Ver la hoja de ruta.
- Las pruebas automáticas cubren el núcleo en Python. El visor se comprueba a mano en el navegador.

## Instalación y uso

    python3 -m venv .venv && source .venv/bin/activate
    pip install -e ".[dev,hf]"
    pytest -q

Requiere Python 3.10 o superior (probado en 3.10, 3.11 y 3.12) y SciPy 1.15 o superior. PySCF es opcional (`.[hf]`) y solo hace falta para la comparación con
Hartree-Fock.

Para ejecutar el visor en local:

    python -m orbitae.export --out web/data
    python -m orbitae.atoms --out web/data
    python -m orbitae.hf_compare --out web/data   # opcional, requiere PySCF
    python3 -m http.server -d web 8000

Abre http://localhost:8000.

Ejemplo en Python:

    from orbitae import sample_orbital, radial_cutoff

    puntos = sample_orbital(n=2, l=1, m=1, size=50_000, seed=0)   # 2p_x, forma (N, 3)
    r95 = radial_cutoff(2, 1, p=0.95)                              # radio en bohr (a0)

## Validación

Cada resultado numérico tiene una prueba: normalización radial y angular, ortogonalidad, ecuación
de Schrödinger radial, energías E_n = -Z²/(2n²), forma de los armónicos s, p, d y f contra sus
fórmulas cartesianas, nodos radiales contra valores analíticos (2s y 3s), muestreo contra
momentos exactos, cargas efectivas de Slater contra valores de referencia y configuraciones
electrónicas de los 36 elementos. La integración continua las ejecuta en cada cambio.

## Error del modelo de Slater frente a Hartree-Fock

Referencia: Hartree-Fock no restringido con PySCF (base def2-TZVP, espín de Hund, densidad
promediada esféricamente). La métrica es el **radio medio de la densidad electrónica total**
(todos los electrones; no mide el tamaño de valencia). HF no incluye correlación electrónica,
así que tampoco es el valor exacto.

| Elemento | Radio medio HF (a₀) | Radio medio Slater (a₀) | Diferencia |
|---|---|---|---|
| H | 1,499 | 1,500 | +0,1 % |
| C | 1,189 | 1,113 | −6,4 % |
| Ne | 0,788 | 0,715 | −9,4 % |
| Na | 0,987 | 0,990 | +0,4 % |
| Ar | 0,892 | 0,861 | −3,5 % |
| Fe | 0,841 | 1,000 | +18,9 % |
| Cu | 0,726 | 0,828 | +14,1 % |
| Kr | 0,729 | 0,746 | +2,4 % |

El modelo funciona bien en los primeros periodos y falla en los metales de transición 3d
(hasta +24 % en el Cr). La tabla de los 36 elementos se genera con
`python -m orbitae.hf_compare` y el visor la muestra sobre la gráfica radial.

### Subcapas de valencia

Mismo método, pero comparando el radio medio de cada subcapa (los orbitales de HF se clasifican por l y por
energía; las densidades de las subcapas suman la densidad total, y eso está comprobado en las pruebas).

| Subcapa | Radio medio HF (a₀) | Radio medio Slater (a₀) | Diferencia |
|---|---|---|---|
| C 2p | 1,701 | 1,538 | −9,5 % |
| Na 3s | 4,215 | 4,773 | +13,2 % |
| Cl 3p | 1,839 | 1,721 | −6,4 % |
| Fe 3d | 1,075 | 1,680 | +56,2 % |
| Fe 4s | 3,255 | 4,144 | +27,3 % |
| Cu 3d | 0,993 | 1,338 | +34,7 % |
| Cu 4s | 3,313 | 4,200 | +26,8 % |
| Kr 4p | 1,950 | 1,884 | −3,4 % |

El modelo de Slater coloca las subcapas 3d y 4s demasiado lejos del núcleo (entre un 27 % y un 56 % en el hierro) y las 2p y 3p
demasiado cerca (entre un 6 % y un 12 %).

> **Corrección (v0.7).** Las versiones 0.5 a 0.6.1 dieron para Mn, Fe, Co y Ni un error del radio medio mayor del real (Fe +25,1 %
> en lugar de +18,9 %). HF sin restricciones converge en esos cuatro átomos a la configuración 3d⁽ⁿ⁺¹⁾ 4s¹ y no a 3d⁽ⁿ⁾ 4s², que es la
> que usa el modelo, y la comparación mezclaba ambos efectos. Desde la v0.7 se fuerza en HF la misma ocupación por subcapa y hay
> una prueba que lo verifica en cada cambio. Los registros anteriores de Zenodo conservan los valores antiguos.

### Iones

Misma métrica y mismo método de referencia (UHF; base def2-TZVPD, con funciones difusas, para los iones).
Los cationes se encogen y los aniones se hinchan respecto al átomo neutro tanto en HF como en el modelo de Slater.

| Ion | Radio medio HF (a₀) | Radio medio Slater (a₀) | Diferencia |
|---|---|---|---|
| Li⁺ | 0,572 | 0,556 | −2,9 % |
| Na⁺ | 0,663 | 0,612 | −7,6 % |
| K⁺ | 0,788 | 0,761 | −3,4 % |
| Mg²⁺ | 0,573 | 0,535 | −6,6 % |
| Ca²⁺ | 0,708 | 0,682 | −3,7 % |
| F⁻ | 0,996 | 0,859 | −13,8 % |
| Cl⁻ | 1,049 | 0,995 | −5,2 % |
| Mn²⁺ | 0,647 | 0,762 | +17,8 % |
| Fe²⁺ | 0,633 | 0,738 | +16,7 % |
| Fe³⁺ | 0,592 | 0,678 | +14,6 % |
| Co³⁺ | 0,580 | 0,662 | +14,1 % |
| Cu⁺ | 0,629 | 0,708 | +12,5 % |
| Cu²⁺ | 0,588 | 0,666 | +13,1 % |
| Zn²⁺ | 0,574 | 0,643 | +12,0 % |

## Hoja de ruta

Ver `docs/SDD.md`. Siguiente: moléculas pequeñas y exportación `.cube`.

## Cómo citar

Usa el botón **Cite this repository** de GitHub (a partir de `CITATION.cff`) o el DOI de Zenodo: https://doi.org/10.5281/zenodo.24135885

## Licencia y contacto

Copyright (C) Pedro Sordo Martínez <amurlaniakea@gmail.com>. AGPL-3.0-or-later (ver `LICENSE`).
