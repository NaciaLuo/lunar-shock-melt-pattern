# lunar-shock-melt-pattern
This repository contains input files and processed output of iSALE-3D simulations. The simulations are performed to investigate impact melt production, melt ejection, and melt deposit patterns under lunar gravity conditions. 

## About iSALE
iSALE-3D is a shock physics code designed for simulating impact processes in planetary science. The code is distributed on a case-by-case basis within the impact community and is restricted to non-commercial use.

For more information about iSALE, including documentation and resources, please visit the  [iSALE Wiki](https://github.com/isale-code/isale-wiki).

## Repository Contents
- Input files used to run iSALE-3D simulations
- Data and Python scripts for producing graphs similar to Figures 4, 6, 9, and 10 in the article (_publication in preparation_)

Simulations were performed with varying impact parameters. To distinguish them, the folders containing the modeling results follow the naming convention: `I(impactor diameter)_deg(impact angle)_V(impact velocity)`, where the units are km, degree, and km/s, respectively. For example, **I06_deg20_V15** stands for **I = 6 km**, **deg = 20 degrees**, **V = 15 km/s**.

## iSALE input files
- **15 CPPR models**: Impact settings of low-resolution ejecta models
- **60 CPPR models**: Impact settings of high-resolution melt models
- **material.inp**: Material parameters used for all simulations

## Figure 4 Melt zones
- `Plot_mzone_xz_yz.py`: plots the xz-plane and yz-plane cross-sections of the shock pressure distribution and the melt zone boundary

- For each model, the projectile (Proj) and target (Tar) grids are stored separately:
  - `(Proj/Tar)_Trp.npy`: peak shock pressure
  - `(Proj/Tar)_mfrac.npy`: melt fraction
  - `(Proj/Tar)_(x/y/z)grid.npy`: 1D grid in x/y/z direction

## Figure 6 Excavation zones
- **Note**: The excavation zone data is generated from the 15-CPPR ejecta model, and the melt zone data is generated from the 60-CPPR melt model.

- `Plot_exzone_mzone.py`: plots the xz-plane cross-section of the excavation zone and the melt zone

- Excavation zone information for ejecta with launch velocities < lunar escape velocity:
  - `ej_(x/y/z)00.npy`: pre-impact location of ejecta
  - `ej_v.npy`: ejecta launch velocity

- **ejecta_escaped**: escaped ejecta (ej_v > lunar escape velocity)

- Melt zone information where the projectile (Proj) and target (Tar) grids are stored separately:
  - `(Proj/Tar)_Trp.npy`: peak shock pressure
  - `(Proj/Tar)_mfrac.npy`: melt fraction
  - `(Proj/Tar)_(x/y/z)grid.npy`: 1D grid in x/y/z direction

- Four additional files that store the (x, y, z) coordinates and the melt fraction of each individual melted tracer as four 1-D arrays:
  - `melt_(x/y/z)0.npy`: pre-impact location of melted tracer
  - `melt_frac.npy`: melt fraction

## Figures 9,10 Melt and ejecta patterns
- `Plot_blankets.py`: plots top-down views of the ejecta thickness pattern, melt thickness pattern, and percentage of melt in ejecta
- Ejecta deposit information:
  - `ej_land(x/y).npy`: landing positions of ejecta
  - `ej_mf.npy`: melt fraction of ejecta
- Transient crater information:
  - `Integrated_surface.npy`: topography of the transient crater
  - `crater_metric.txt`: crater volume, depth, crossrange diameter, and downrange diameter with timestep
 

