# SDD Orbitae (v0.9)

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
- RF5 Iones de interés biomédico con comparación frente a HF: HECHO (v0.6).
- RF6 Comparación con Hartree-Fock (PySCF), total y por subcapa: HECHO (v0.5, v0.7). Moléculas pequeñas de capa cerrada: HECHO (v0.8).
- RF7 Fichas de contenido por elemento y nivel (ES/EN): parcial (frase biomédica por elemento, ya en español e inglés).
- RF8 Exportar orbitales en formato .cube: HECHO para moléculas (v0.8); imagen y datos de átomos: pendiente.

## 5. Requisitos no funcionales
Reproducibilidad (semillas fijas, versiones probadas en CI), accesibilidad (daltonismo, contraste),
ES/EN, AGPL-3.0-or-later, atribución "Pedro Sordo Martínez", y cada fuente de datos con su licencia.
Visor sin dependencias de terceros y con política CSP restrictiva.

## 6. Arquitectura
core (Python) -> export (reservas y perfiles precalculados) -> web estática (WebGL puro)
Átomos: el visor combina una reserva angular (direcciones y signo) con una reserva radial
(Gamma) mediante r = g / (2 zeta); no se guarda una nube por orbital.
Despliegue: GitHub Pages desde web/. Los datos del visor se generan en la CI (no se versionan) y una prueba funcional
del visor en jsdom debe pasar antes de publicar.

## 7. Criterios de aceptación (validación)
- Normalización y ortogonalidad radial y angular.
- Ecuación de Schrödinger radial y energías E_n = -Z^2/(2 n^2) hartree.
- Armónicos reales s, p, d, f contra sus formas cartesianas.
- Nodos radiales contra valores analíticos; radio de corte contra forma cerrada (1s).
- Muestreo contra momentos analíticos (hidrógeno y orbitales de Slater).
- Reglas de Slater contra valores de referencia y configuraciones de Z = 1 a 36.
- Radio medio y densidad radial del modelo de Slater frente a UHF/def2-TZVP (PySCF), con el error
  documentado en el README y mostrado en el visor.
- Subcapas de HF con la misma ocupación que el modelo (comprobado en pruebas) y error por subcapa.
- Visor: prueba funcional en jsdom (modos, tabla periódica, iones, moléculas, español e inglés) antes de cada despliegue.
- Textos: todo texto biomédico de elementos, iones y moléculas tiene su versión en inglés (comprobado en pruebas).
- Moléculas: energías HF frente a referencia, normalización y ortogonalidad de los orbitales de valencia, degeneración de los niveles,
  plano nodal del HOMO del agua, HOMO del CO en el carbono y lectura de un .cube exportado.

## 8. Hoja de ruta
| Versión | Entrega | Estado |
|---|---|---|
| v0.1 | Hidrógeno exacto, armónicos reales, muestreo y validación | hecho |
| v0.2 | Exportador y visor web WebGL, Pages, CI | hecho |
| v0.3 | Orbitales f (n = 4) y distribución radial con nodos | hecho |
| v0.4 | Átomos 1-36 con Slater, tabla periódica, nube total | hecho |
| v0.5 | Comparación con PySCF: error de la aproximación de Slater | hecho |
| v0.6 | Iones (Na⁺, K⁺, Ca²⁺, Fe²⁺/Fe³⁺…) con comparación frente a HF | hecho |
| v0.7 | Subcapas frente a HF y corrección de la referencia (ocupación forzada) | hecho |
| v0.8 | Moléculas pequeñas de capa cerrada y exportación .cube | hecho |
| v0.8.1 | Datos generados en la CI y prueba funcional del visor | hecho |
| v0.9 | Interfaz bilingüe ES/EN del visor y textos biomédicos en inglés | hecho |
| v0.10 | Moléculas con espín (O₂, NO) y más biomoléculas pequeñas | siguiente |

## 9. Riesgos
Licencias de tablas de datos, rendimiento en móviles, afirmaciones médicas excesivas,
presentar una aproximación como si fuera exacta, saltar a DFT demasiado pronto, comparar contra un estado de referencia distinto del modelado (corregido en v0.7), presentar los orbitales virtuales de HF como estados excitados.
