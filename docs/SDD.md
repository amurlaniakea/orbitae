# SDD Orbitae (v0.4)

Autor: Pedro Sordo Martínez. Licencia: AGPL-3.0-or-later.

## 1. Propósito
Plataforma educativa y científica de código abierto para visualizar orbitales y densidad
electrónica de átomos, iones y moléculas sencillas, con validación numérica reproducible
y contenido por nivel de lectura (estudiante, medicina y enfermería, investigador).

## 2. Posicionamiento
Ya existen herramientas que generan orbitales (por ejemplo atorvi, que exporta a formato XSF, y
atomictools). Orbitae no compite en volumen de funciones sino en tres ejes:
1. Validación numérica pública y reproducible de cada resultado, incluido el error de sus aproximaciones.
2. Contenido didáctico por nivel, en español e inglés.
3. Puente hacia biología y medicina (iones, metales, isótopos, moléculas pequeñas).

## 3. Alcance y límites
Herramienta educativa. No ofrece uso clínico, diagnóstico ni dosificación.
Los orbitales atómicos tienen poca aplicación clínica directa; el valor médico viene de iones
fisiológicos (Na+, K+, Ca2+, Mg2+, Cl-), metales en biología (Fe, Zn, Cu), isótopos de medicina
nuclear y moléculas pequeñas. Las aproximaciones se declaran en la propia interfaz.

## 4. Requisitos funcionales
- RF1 Orbitales hidrogenoides reales y complejos (n, l, m), hasta n = 4. HECHO (v0.1, v0.3).
- RF2 Átomos 1-36 con reglas de Slater y orbitales tipo Slater; nube total y orbitales. HECHO (v0.4).
- RF3 Alternar nube total y orbital individual. HECHO (v0.4).
- RF4 Fases en dos colores y escala común; isosuperficies y curva de transferencia: pendientes.
- RF5 Iones cambiando la ocupación: pendiente.
- RF6 Comparación con Hartree-Fock (PySCF) y moléculas pequeñas: pendiente.
- RF7 Fichas de contenido por elemento y nivel (ES/EN): parcial (frase biomédica por elemento).
- RF8 Exportar imagen, datos y formato .cube: pendiente.

## 5. Requisitos no funcionales
Reproducibilidad (semillas fijas, versiones probadas en CI), accesibilidad (daltonismo, contraste),
ES/EN, AGPL-3.0-or-later, atribución "Pedro Sordo Martínez", y cada fuente de datos con su licencia.
Visor sin dependencias de terceros y con política CSP restrictiva.

## 6. Arquitectura
core (Python) -> export (reservas y perfiles precalculados) -> web estática (WebGL puro)
Átomos: el visor combina una reserva angular (direcciones y signo) con una reserva radial
(Gamma) mediante r = g / (2 zeta); no se guarda una nube por orbital.
Despliegue: GitHub Pages desde web/ con integración continua.

## 7. Criterios de aceptación (validación)
- Normalización y ortogonalidad radial y angular.
- Ecuación de Schrödinger radial y energías E_n = -Z^2/(2 n^2) hartree.
- Armónicos reales s, p, d, f contra sus formas cartesianas.
- Nodos radiales contra valores analíticos; radio de corte contra forma cerrada (1s).
- Muestreo contra momentos analíticos (hidrógeno y orbitales de Slater).
- Reglas de Slater contra valores de referencia y configuraciones de Z = 1 a 36.
- Pendiente: orbitales de Slater frente a PySCF con error documentado.

## 8. Hoja de ruta
| Versión | Entrega | Estado |
|---|---|---|
| v0.1 | Hidrógeno exacto, armónicos reales, muestreo y validación | hecho |
| v0.2 | Exportador y visor web WebGL, Pages, CI | hecho |
| v0.3 | Orbitales f (n = 4) y distribución radial con nodos | hecho |
| v0.4 | Átomos 1-36 con Slater, tabla periódica, nube total | hecho |
| v0.5 | Comparación con PySCF: error de la aproximación de Slater | siguiente |
| v0.6 | Iones y módulo Biología y medicina ampliado | |
| v0.7 | Moléculas pequeñas, exportación .cube, DOI en Zenodo | |

## 9. Riesgos
Licencias de tablas de datos, rendimiento en móviles, afirmaciones médicas excesivas,
presentar una aproximación como si fuera exacta, saltar a DFT demasiado pronto.
