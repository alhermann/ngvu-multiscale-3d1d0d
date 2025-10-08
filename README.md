# NGVU-Multiscale-3D1D0D

> Research code for our preprint
> **A 3D-1D-0D Multiscale Model of the Neuro-Glial-Vascular Unit for Synaptic and Vascular Dynamics in the Dorsal Vagal Complex**
> Hermann, Alexander, *et al.*, arXiv:2504.02540 (2025)

This repository couples (i) a 1D macrocirculation solver, (ii) a DUNE-based 3D/1D microcirculation and transport module, and (iii) a MATLAB neuro-glial-vascular unit (NGVU) model. It builds the flow inputs, launches the micro model, and time-marches the NGVU while exchanging vessel radii via CSV.

---

## Table of contents

* [Repository layout](#repository-layout)
* [Requirements](#requirements)
* [Quick start](#quick-start)
* [Build & run — details](#build--run--details)
* [Data & inputs](#data--inputs)
* [Data handoff](#data-handoff)
* [Reproducing the paper](#reproducing-the-paper)
* [Citations](#citations)
* [License](#license)
* [Acknowledgments](#acknowledgments)
* [Troubleshooting](#troubleshooting)

---

## Repository layout

```
.
├── Macrocirculation/                # 1D nonlinear solver (C++/CMake)
├── Microcirculation/Code_Upload/    # DUNE modules + dune-angiogenesis (C++)
├── NGVU/                            # MATLAB NGVU model + coupling
│   ├── helpers/                     # functions (all_constants.m, NVC.m, ...)
│   └── quadripartite_NVC_main.m     # merged quadripartite + NVC driver
└── master.m                         # top-level orchestrator (MATLAB)
```

---

## Requirements

**System**

* Linux (tested on Ubuntu 22.04+). macOS may work with GNU toolchain, but not tested.

**Compilers / toolchain**

* GCC/G++/GFortran **9.x** (required by DUNE stack)
* CMake ≥ 3.16
* Make / Ninja

**Libraries**

* DUNE core modules (fetched/built by `Microcirculation/Code_Upload/DUNE/build.sh`)
* Python 3.8+ with standard library (for macrocirculation post-processing)

**MATLAB**

* MATLAB R2022b or newer (tested), toolboxes: base MATLAB only.

---

## Quick start

1. **Clone**

```bash
git clone https://github.com/alhermann/ngvu-multiscale-3d1d0d.git
cd ngvu-multiscale-3d1d0d
```

2. **Open MATLAB in repo root** and run:

```matlab
run('master.m')
```

`master.m`:

* builds and runs **Macrocirculation** (CMake + solver),
* runs a Python helper to write `pica_flow_heart_period.csv`,
* builds **Microcirculation** (DUNE; checks GCC=9),
* copies the CSV into the DUNE build tree,
* launches the **micro** executable (`dune-angiogenesis ArtificialNetwork`) **in the background**,
* runs **NGVU** (`NGVU/quadripartite_NVC_main.m`), streaming vessel radii to CSV files consumed by the micro module,
* finishes when NGVU is done.

> Default configuration is set inside `master.m`. Adjust solver flags there if needed.

---

## Build & run — details

### Macrocirculation (C++)

The master script can do this for you. To do it manually:

```bash
cd Macrocirculation
mkdir -p build && cd build
cmake ..
make -j
# Run (example)
cd bin/macrocirculation
./MacrocirculationNonlinear1DSolver \
  --mesh-file "../../data/1d-meshes/37-vessels.json" \
  --t-end 20 \
  --output-directory "./output"
# Post-process (writes pica_flow_heart_period.csv)
python3 write_flow.py \
  --vessels 36 \
  --filepath "./output/abstract_vessels.json" \
  --output "./output/pica_flow_heart_period.csv" \
  --t-end 20
```

### Microcirculation (DUNE + dune-angiogenesis)

The master script runs the DUNE `build.sh` and then launches:

```
Microcirculation/Code_Upload/DUNE/dune-angiogenesis/build-cmake/src/dune-angiogenesis ArtificialNetwork
```

This program:

* reads parameters from `parameters.ini`,
* waits for NGVU radii streams in `NVU_Vessels/NVU_Vessel_<i>.csv`,
* interpolates radii in time during the 3D/1D simulation.

> GCC/G++/GFortran **must** be version 9. The master script checks and prints a helpful message if not.

### NGVU (MATLAB)

`NGVU/quadripartite_NVC_main.m`:

* runs a merged quadripartite + NVC driver,
* appends per-vessel radii (µm) as `NVU_Vessels/NVU_Vessel_<i>.csv` with rows:

  ```
  time_in_seconds, radius_in_micrometers
  ```
* no headers, one row per exchange step.

---

## Data & inputs

* **Macrocirculation → Microcirculation flow CSV**
  `pica_flow_heart_period.csv` is produced by the macrocirculation post-process and copied by `master.m` to the micro build tree.

* **Reference microvascular network (brain99)**
  The microcirculation module can be configured to use the **brain 1999** 3D capillary network dataset (“brain99”):

  * Data page: [https://sites.arizona.edu/secomb/3d-network-data-brain-1999/](https://sites.arizona.edu/secomb/3d-network-data-brain-1999/)

  Please follow the usage terms provided on the data page and cite appropriately (see **Citations** below). If you use another network, adjust `parameters.ini` and any network loading paths accordingly.

---

## Data handoff

* **Macrocirculation → Microcirculation**: `pica_flow_heart_period.csv` (written by Python postprocess) is copied by `master.m` into
  `Microcirculation/Code_Upload/DUNE/dune-angiogenesis/build-cmake/src/`.

* **NGVU → Microcirculation**: the NGVU writes:

  ```
  Microcirculation/Code_Upload/DUNE/dune-angiogenesis/build-cmake/src/NVU_Vessels/NVU_Vessel_<1..50>.csv
  ```

  Each file is monotonically increasing in time. The micro code blocks until each file has at least the current time or later, then **interpolates** the radius for its current time step.

---

## Reproducing the paper

1. Ensure compilers are GCC=9 toolchain; install MATLAB and Python.
2. Run `master.m` from the repo root.
3. Adjust simulation times and mesh inside `master.m` and the DUNE `parameters.ini` to match the scenarios in the preprint.
4. If using the **brain99** network, make sure the corresponding network/parameter configuration is selected in the microcirculation module.
5. Results (pressures, concentrations, radii) are written into module-specific output folders and CSV files.

---

## Citations

If you use this code, please cite:

* **This work**
  Hermann, Alexander, *et al.* “A 3D-1D-0D Multiscale Model of the Neuro-Glial-Vascular Unit for Synaptic and Vascular Dynamics in the Dorsal Vagal Complex.” *arXiv:2504.02540* (2025).

* **Upstream components & models**

  * Köppl, Vidotto, Wohlmuth (2020). *Int. J. Numer. Meth. Biomed. Eng.*, 36(10):e3386.
  * Fritz, *et al.* (2022). *IJNMBE*, 38(7):e3612. Source: [https://github.com/CancerModeling/Flows1D0D3D/releases/tag/v1.0](https://github.com/CancerModeling/Flows1D0D3D/releases/tag/v1.0)
  * Dormanns, *et al.* (2015). *J. Theor. Biol.* 364:49–70. Repo: [https://github.com/brainstrust/envy-you](https://github.com/brainstrust/envy-you)

* **Reference microvascular network (brain99)**
  Secomb Lab, “3D Network Data – Brain 1999 (brain99).”
  Data and description: [https://sites.arizona.edu/secomb/3d-network-data-brain-1999/](https://sites.arizona.edu/secomb/3d-network-data-brain-1999/)
  *(Please cite the dataset/webpage and any linked publications per the site’s guidance.)*

---

## License

This project is released under **GNU General Public License v3.0** (GPL-3.0).
See [`LICENSE`](LICENSE).

* Third-party components (DUNE modules, upstream repositories) carry their own licenses.
  See `THIRD_PARTY_LICENSES.md`.

---

## Acknowledgments

This repository builds on excellent open-source efforts by the DUNE community and the authors of the cited works and datasets. We are grateful for their contributions and for the availability of the **brain99** network data.

---

## Troubleshooting

* **DUNE build fails**: verify `gcc --version`, `g++ --version`, `gfortran --version` all show **9.x**.
* **CSV not found**: ensure `RUN_MACRO = true` in `master.m`, or provide an existing `pica_flow_heart_period.csv`.
* **Micro waits forever**: confirm NGVU writes `NVU_Vessels/NVU_Vessel_<i>.csv` with increasing time; check permissions and paths.
* **MATLAB errors**: add `NGVU/helpers` to the MATLAB path or run from repo root so relative paths resolve.

---

### How to contribute

Issues and pull requests are welcome. Please include:

* OS and compiler versions,
* minimal reproduction steps,
* logs from `master.m` and the module that failed.

Happy modeling!

