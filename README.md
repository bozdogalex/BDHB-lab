# BDHB — Big Data in Health & Bioinformatics Labs

[![CI](https://github.com/bozdogalex/BDHB-lab/actions/workflows/ci.yml/badge.svg)](https://github.com/bozdogalex/BDHB-lab/actions/workflows/ci.yml)
[![Open in GitHub Codespaces](https://img.shields.io/badge/Open%20in-Codespaces-000?logo=github)](https://github.com/codespaces/new?hide_repo_select=true&ref=main&repo=bozdogalex%2FBDHB-lab)

> Master’s-level laboratories blending classical bioinformatics with modern ML, networks, and GenAI. The environment is CPU‑only and identical across Codespaces and Docker via the prebuilt image `ghcr.io/bozdogalex/bdhb:base`.


## Labs (index)

- 01 — Databases & GitHub: [labs/01_intro&databases](labs/01_intro&databases)
- 02 — Sequence Alignment: [labs/02_alignment](labs/02_alignment)
- 03 — NGS: [labs/03_formats&NGS](labs/03_formats&NGS)
- 04 — Phylogenetics: [labs/04_phylogenetics](labs/04_phylogenetics) · [Slides (PowerPoint)](docs/presentations/04_phylogenetics.pptx)
- 05 — Clustering: [labs/05_clustering](labs/05_clustering)
- 07— WGCNA + Diseasome: [labs/06_wgcna](labs/06_wgcna)
- 08 — Network Viz & GNN: [labs/07_network_viz](labs/07_network_viz)
- 09 — Federated Learning: [labs/08_ML_flower](labs/08_ML_flower)
- 10 — Drug Repurposing: [labs/09_repurposing](labs/09_repurposing)
- 11 — Integrative + Digital Twin: [labs/10_integrative](labs/10_integrative)
- 12 — Multi‑omics + Quantum : [labs/11_multiomics](labs/11_multiomics)
- 13 — Generative AI : [labs/12_genAI](labs/12_genAI)
- Presentations

---

> Full onboarding (screenshots, tips): **[Onboarding](docs/onboarding.md)**

---

## Repo map

- `labs/` — all weekly lab content
- `docs/` — onboarding, ANIS pack (before/after, one‑pagers, screenshots)
- `mlops/` — MLflow helpers
- `.devcontainer/` — Codespaces/Devcontainer (pulls GHCR image)
- `.github/workflows/` — CI + image publish
- `Dockerfile`, `requirements.txt` — env definition
- `dev.ps1`, `Makefile` — local helpers

---
# Docs

Supporting material and submission pack live under `docs/`:

- [Onboarding](docs/onboarding.md) — Codespaces & Docker setup, smoke test, troubleshooting
- [Before/After](docs/before_after.md) — ANIS submission improvement pack
- [One-pagers](docs/lab_onepagers/) — PDF summaries of labs
- [Screenshots](docs/screens/) — environment/UI captures (MLflow, Codespaces, Argo)
- [Changelog](docs/changelog.md) — changes across versions
- [Policies](docs/policies.md) — third-party license references, repository policies
- [Resources](docs/resources.md) — recommended readings/tutorials
---

## Contributing / Policies / Citation

- `CONTRIBUTING.md` — contribution rules & PR tips  
- `CITATION.cff`  — how to cite this work  
- `docs/changelog.md` — changelog (linked from releases)

