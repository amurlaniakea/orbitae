---
title: Pathogen Forecast
emoji: 📈
colorFrom: red
colorTo: indigo
sdk: streamlit
pinned: false
license: agpl-3.0
short_description: Visor del desarrollo de pandemia epidemiológica (SIR/SEIR)
sdk_version: 1.66.0
---

Check out the configuration reference at https://huggingface.co/docs/hub/spaces-config-reference
## Visor web

    python -m orbitae.export --out web/data
    python3 -m http.server -d web 8000

Abre http://localhost:8000. El visor es WebGL puro, sin librerías de terceros, y solo
se conecta a su propio servidor (política CSP).
