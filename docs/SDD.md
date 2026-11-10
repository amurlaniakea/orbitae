# SDD Orbitae (borrador v0.1)

Autor: Pedro Sordo Martínez. Licencia: AGPL-3.0-or-later.

## 1. Propósito
Plataforma educativa y científica de código abierto para visualizar orbitales y densidad
electrónica de átomos, iones y moléculas sencillas, con validación numérica reproducible
y contenido por nivel de lectura (estudiante, medicina y enfermería, investigador).

## 2. Posicionamiento
Ya existen herramientas que generan orbitales (por ejemplo atorvi, que exporta a formato
XSF, y atomictools). Orbitae no compite en volumen de funciones sino en tres ejes:
1. Validación numérica pública y reproducible de cada resultado.
2. Contenido didáctico por nivel, en español e inglés.
3. Puente hacia biología y medicina (iones, metales, isótopos, moléculas pequeñas).

## 3. Alcance y límites
Herramienta educativa. No ofrece uso clínico, diagnóstico ni dosificación.
Los orbitales atómicos tienen poca aplicación clínica directa; el valor médico viene de
iones fisiológicos (Na+, K+, Ca2+, Mg2+, Cl-), metales en biología (Fe, Zn, Cu),
isótopos de medicina nuclear y moléculas pequeñas.

## 4. Requisitos funcionales
- RF1 Orbitales hidrogenoides reales y complejos (n, l, m).
- RF2 Átomos 1-36 con orbitales tipo Slater tabulados y densidad total.
- RF3 Alternar orbital individual, valencia y densidad total.
- RF4 Fases en dos colores, isosuperficies 90/95 %, curva de transferencia ajustable.
- RF5 Iones cambiando la ocupación.
- RF6 Moléculas pequeñas con PySCF.
- RF7 Fichas de contenido por elemento en tres niveles de lectura (ES/EN).
- RF8 Exportar imagen, datos y formato .cube.

## 5. Requisitos no funcionales
Reproducibilidad (versiones fijadas), accesibilidad (daltonismo, contraste), ES/EN,
AGPL-3.0-or-later, atribución "Pedro Sordo Martínez", cada tabla de datos con fuente
y licencia documentadas.

## 6. Arquitectura
core (Python) -> validation -> export (volúmenes precalculados) -> web estática (Three.js)

Decisión: precalcular en Python y servir como sitio estático; sin backend.
Módulos: core/, validation/, data/, export/, web/, content/, docs/.

## 7. Criterios de aceptación (validación)
- Normalización radial y angular dentro de tolerancia.
- Ortogonalidad radial y angular.
- Ecuación de Schrödinger radial y energías E_n = -Z^2/(2 n^2) hartree.
- Radio de corte contra forma cerrada (1s).
- Muestreo contra momentos analíticos (<r> = 1.5 en 1s; <x^2>/<r^2> = 3/5 en p_x).
- Más adelante: orbitales STO frente a PySCF con error documentado.

## 8. Hoja de ruta
| Versión | Entrega |
|---|---|
| v0.1 | Hidrógeno exacto, orbitales reales, muestreo, validación (hecho: 20 pruebas) |
| v0.2 | Visor web básico y átomos 1-36 con STO |
| v0.3 | Iones y módulo Biología y medicina |
| v0.4 | Moléculas pequeñas con PySCF |
| v0.5 | Exportación .cube, modo investigador, DOI en Zenodo |

## 9. Riesgos
Licencias de tablas de datos, rendimiento en móviles, afirmaciones médicas excesivas,
saltar a DFT demasiado pronto.
