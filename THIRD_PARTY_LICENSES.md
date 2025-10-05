# Third-Party Licenses and Notices

This document lists third-party software used by **NGVU-Multiscale-3D1D0D** and provides links to the original projects and their license files. Please review upstream licenses carefully for your specific versions.

> **Important:** Licensing terms can vary by version. Always check the LICENSE/COPYING files in the exact revision you build.

---

## Contents

* [DUNE ecosystem (Microcirculation)](#dune-ecosystem-microcirculation)
* [Macrocirculation external libraries](#macrocirculation-external-libraries)
* [Upstream research code & references](#upstream-research-code--references)
* [MATLAB](#matlab)
* [Attribution format](#attribution-format)

---

## DUNE ecosystem (Microcirculation)

The microcirculation module builds on the DUNE framework and modules, typically licensed under **GPL-2.0 (or later) with the DUNE Runtime Exception** (often called the “DUNE License”). Verify the exact terms for each module revision you use.

These are fetched/built by `Microcirculation/Code_Upload/DUNE/build.sh`:

* **dune-common** — License: see upstream `COPYING`/`LICENSE`
* **dune-geometry** — License: see upstream `COPYING`/`LICENSE`
* **dune-grid** — License: see upstream `COPYING`/`LICENSE`
* **dune-istl** — License: see upstream `COPYING`/`LICENSE`
* **dune-typetree** — License: see upstream `COPYING`/`LICENSE`
* **dune-localfunctions** — License: see upstream `COPYING`/`LICENSE`
* **dune-functions** — License: see upstream `COPYING`/`LICENSE`
* **dune-foamgrid** — License: see upstream `COPYING`/`LICENSE`

Additional module:

* **dune-angiogenesis** (this repo’s micro application depends on it) — License: see upstream `LICENSE`/`COPYING` in that project.

> **Action for users:** After your first successful DUNE build, record the commit IDs for each module in this file.

---

## Macrocirculation external libraries

Inside `Macrocirculation/external/`:

* **pybind11**

  * Location: `Macrocirculation/external/pybind11`
  * License: **BSD-3-Clause** (see `LICENSE` file included in the directory, or upstream repository)

* **GMM++ (gmm)**

  * Location: `Macrocirculation/external/gmm`
  * License: typically **LGPL**; confirm the exact license/version in the bundled `COPYING`/`LICENSE` file or upstream repository.

---

## Upstream research code & references

Parts of this work build on ideas, formulations, and/or helper code from these research projects. Please consult their repositories for exact licenses and cite as appropriate.

* **Flows1D0D3D (Köppl, Vidotto, Wohlmuth; Fritz et al.)**

  * Repo (release referenced in the paper): [https://github.com/CancerModeling/Flows1D0D3D/releases/tag/v1.0](https://github.com/CancerModeling/Flows1D0D3D/releases/tag/v1.0)
  * License: see upstream `LICENSE` in that release.

* **envy-you (Dormanns et al., Neurovascular coupling)**

  * Repo: [https://github.com/brainstrust/envy-you](https://github.com/brainstrust/envy-you)
  * License: see upstream `LICENSE`.

---

## MATLAB

This repository includes MATLAB scripts (NGVU model and the master orchestrator). **MATLAB is proprietary software**. Users must have a valid MATLAB license to run the MATLAB components. MATLAB itself is **not** redistributed here.

---

## Attribution format

When you pin versions, please update the sections above like this:

```
dune-grid
  Version/commit: v2.9.1 (abcd123)
  License: GPL-2.0-or-later with DUNE Runtime Exception
  Source: https://gitlab.dune-project.org/core/dune-grid
  Local LICENSE: Microcirculation/Code_Upload/DUNE/dune-grid/COPYING
```

---

### Notes

* If you distribute binaries that statically link against GPL/LGPL code, ensure your distribution complies with those licenses (source offer, relinking rights, notices, etc.).
* This project’s top-level license is **GPL-3.0**; ensure your integration is compatible with upstream terms.

---
