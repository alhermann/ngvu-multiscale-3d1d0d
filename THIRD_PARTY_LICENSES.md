# Third-Party Licenses

This project vendors and/or links against several third-party packages. Below is a summary of each component’s license with pointers to the full texts bundled in this repository or upstream. **Always consult the original license text for authoritative terms.**

> **Note on DUNE modules:** Many DUNE core modules are under **GPL-2.0 with a special “runtime exception.”** This exception is intended for template/header libraries and permits using/instantiating the headers without imposing the GPL on your entire executable *solely due to that usage*. See each module’s license for the exact wording.

---

## DUNE modules (vendored under `Microcirculation/Code_Upload/DUNE/…`)

* **dune-common** — GPL-2.0 **with runtime exception**
  Full text mirrored at: `dune-common/LICENSE.md` (also present in the module directory).

* **dune-grid** — GPL-2.0 **with runtime exception**
  Full text mirrored at: `dune-grid/LICENSE.md`.

* **dune-geometry** — GPL-2.0 **with runtime exception**
  Full text mirrored at: `dune-geometry/LICENSE_geometry.md`.

* **dune-istl** — GPL-2.0 **with runtime exception**
  Full text mirrored at: `dune-istl/LICENSE_istl.md`.

* **dune-localfunctions** — GPL-2.0 **with runtime exception**
  Full text mirrored at: `dune-localfunctions/LICENSE_localfunctions.md`.

* **dune-typetree** — **Dual license**: **LGPL-3.0-or-later** *or* **GPL-2.0 with runtime exception**
  Full text mirrored at: `dune-typetree/LICENSE_typetree.md`.

* **dune-foamgrid** — **Dual license**: **LGPL-3.0-or-later** *or* **GPL-2.0 with runtime exception**
  Full text mirrored at: `dune-foamgrid/LICENSE_foamgrid.md`.

* **dune-angiogenesis** (module used by Microcirculation) — **Dual license**: **LGPL-3.0-or-later** *or* **GPL-2.0 with runtime exception**
  Full text mirrored at: `dune-angiogenesis/LICENSE_angiogenesis.md`.

---

## Referenced/related upstream research code (not necessarily vendored here)

* **envy-you (Dormanns et al., Neurovascular coupling)**
  Repo: https://github.com/brainstrust/envy-you

  **Status in this repository:** No source code from *envy-you* is included or redistributed here. Our implementation was developed independently and only draws conceptual inspiration from the published work (Dormanns et al., 2015). At the time of writing, the *envy-you* repository does not contain a clear license file. Because we do **not** ship their code, the licensing of *envy-you* does not impose terms on this repository. If you choose to use *envy-you* directly, consult the authors and/or the upstream repository for licensing details.

* **Flows1D0D3D (Köppl, Vidotto, Wohlmuth; Fritz et al.)**
  Release referenced in the paper: [https://github.com/CancerModeling/Flows1D0D3D/releases/tag/v1.0](https://github.com/CancerModeling/Flows1D0D3D/releases/tag/v1.0)
  **License:** See the **upstream release’s LICENSE**. We base parts of our modeling ideas and comparisons on this work; the code is **not** redistributed here.

---

## Citations (as noted in the README)

* Hermann, Alexander, *et al.* (2025) — *A 3D-1D-0D Multiscale Model of the Neuro-Glial-Vascular Unit for Synaptic and Vascular Dynamics in the Dorsal Vagal Complex.* arXiv:2504.02540.
* Köppl, T., Vidotto, E., Wohlmuth, B. (2020) — *A 3D-1D coupled blood flow and oxygen transport model to generate microvascular networks.* IJNMBE 36(10): e3386.
* Fritz, M., *et al.* (2022) — *A 1D–0D–3D coupled model for simulating blood flow and transport processes in breast tissue.* IJNMBE 38(7): e3612.
* Dormanns, K., *et al.* (2015) — *Neurovascular coupling and the influence of luminal agonists via the endothelium.* Journal of Theoretical Biology 364: 49–70.

---

## The “Runtime Exception” (summary)

Several DUNE modules add a special exception to GPL-2.0 to address template/header usage. In short, **including/instantiating** templates or **linking** object files produced from those headers **does not by itself** make the entire executable a derivative work under GPL-2.0. Other obligations (for non-template source, copied code, or direct modifications) still apply. Always review the exact exception text in the module’s license.

---

## Practical notes

* **We mirror license texts** for DUNE modules at the repo root as `LICENSE_*.md` for convenience; the originals are also present within the corresponding subfolders.
* If you **add or update** third-party code, append a section here and include its license file.

---

**This document is a convenience summary, not legal advice.**
