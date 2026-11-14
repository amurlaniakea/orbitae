# Contributing / Contribuir

**English.** Thanks for your interest. The rule of this repository: **every numerical claim has a test**, and a change
that alters published numbers must include a correction note (see the README, "Correction (v0.7)").

    python3 -m venv .venv && source .venv/bin/activate
    pip install -e ".[dev,hf]"
    pytest -q

The viewer data (`web/data`) is generated, not committed; see the README to regenerate it, and
`tests/viewer/viewer.test.js` for the functional test of the viewer (Node 22 or later). Issues and pull requests are
welcome in English or Spanish. By contributing you agree that your contribution is licensed under AGPL-3.0-or-later.

**Español.** Gracias por tu interés. La regla de este repositorio: **toda cifra numérica tiene una prueba**, y un cambio
que modifique cifras ya publicadas debe incluir una nota de corrección (véase el README, "Corrección (v0.7)").
Los comandos son los de arriba. Los datos del visor (`web/data`) se generan, no se suben a git; el README explica cómo
regenerarlos. Se aceptan issues y pull requests en español o en inglés. Al contribuir aceptas que tu aportación se
licencie bajo AGPL-3.0-or-later.
