
  This repository contains experimental data, simulation data, and scripts to
  reproduce all figures from the paper on stress transmission, force anisotropy,
  and hysteresis in quasi-statically loaded granular materials (glass beads).

 
  ## Requirements

  ### MATLAB (all figure scripts)
  - MATLAB R2019b or later
  - No additional toolboxes required

  ### Python (simulation and data extraction only)
  - Python 3.7+
  - `numpy`, `vtk`
  - [ESyS-Particle](https://launchpad.net/esys-particle) (only to re-run the simulation)

  ```bash
  pip install numpy vtk

  ---
  Reproducing the Figures

  Set the MATLAB working directory to Experiments/ before running experimental
  scripts, and to Simulation/ before running simulation scripts.

  ---
Experimental figures

  Figure 1a — Strain vs. confining pressure

  Script: Experiments/strain_stress_figs1a_c.m
  Reads: Experiments/Glass beadsGlob.mat

  Computes relative strain from LVDT displacement:

  strain (%) = displacement (mm) × 1e⁻³ / 0.5 m × 100

  Plots compression (solid teal) and decompression (dashed gray) branches for each
  cycle, plus an inset with the tangent elastic modulus M = dP/dε per cycle.

  ---
  Figure 1c — Internal stress vs. confining pressure

  Script: Experiments/strain_stress_figs1a_c.m
  Reads: Experiments/Glass beadsGlob.mat

  Plots lateral stress σ_xx (average of 3 sensors) and vertical stress σ_yy (average
  of 3 sensors) against confining pressure, with compression and decompression
  branches separated. Sensor pressures are converted to kPa using the contact area
  of each transducer tip (r = 0.0185 m) relative to the cell area (0.5 × 0.5 m²).

  ---
  Figure 3 — Cumulative strain and dissipated energy per cycle

  Scripts:
  - Experiments/cummulative_strain_fig3.m → cumulative strain panel
  - Experiments/energy_fig3.m → dissipated energy panel

  Reads: All 6 ReadyQuasiStaticData_*.mat files + corresponding dt_*.mat files

  Both scripts combine six experimental datasets spanning ~400 days of loading.
  cummulative_strain_fig3.m computes per-cycle strain increment ΔΕ with error bars.
  energy_fig3.m integrates the pressure–strain hysteresis loop area (W = ∫p dε) for
  each cycle as a measure of dissipated energy. Both figures share a dual x-axis
  (cycle number / elapsed time in days).

  ---
Simulation figures

  All simulation figures are generated from pre-computed data already in
  Simulation/Data/. Set the MATLAB working directory to Simulation/ and run:

  Macroscopic stress–strain curve

  Script: Simulation/Data-Processing/WallForcesAnalysis.m
  Reads: Data/Walls/floorPosition.dat, floorForce.dat, roofForce.dat,
  x±WallForce.dat, z±WallForce.dat
  Output: Figures/FloorForcevsStrain.png

  Computes strain ε = (L₀ − floor displacement) / L₀ with L₀ = 0.4 m and wall
  stress σ = Force / 0.16 m². Plots cycle-averaged compression and decompression
  branches with directional arrows.

  ---
  Mean and maximum contact forces vs. confining pressure

  Script: Simulation/Data-Processing/ForceComponentsAnalysis.m
  Reads: Data/Forces/xyz_Forces_Center.csv, _Medium.csv, _Corner.csv
  Output: Figures/Center/, Figures/Medium/, Figures/Corner/ — each with:

  ┌────────────┬───────────────────────────────────────────────────────────────┐
  │    File    │                            Content                            │
  ├────────────┼───────────────────────────────────────────────────────────────┤
  │ meanfx.png │ Mean radial force ⟨Fx⟩ vs. pressure — top / middle / bottom   │
  ├────────────┼───────────────────────────────────────────────────────────────┤
  │ meanfy.png │ Mean vertical force ⟨Fy⟩ vs. pressure — top / middle / bottom │
  ├────────────┼───────────────────────────────────────────────────────────────┤
  │ Fxmax.png  │ Maximum radial force vs. pressure                             │
  ├────────────┼───────────────────────────────────────────────────────────────┤
  │ Fymax.png  │ Maximum vertical force vs. pressure                           │
  ├────────────┼───────────────────────────────────────────────────────────────┤
  │ Fzmax.png  │ Maximum tangential force vs. pressure                         │
  └────────────┴───────────────────────────────────────────────────────────────┘

  ---
  Internal stress tensor vs. confining pressure

  Script: Simulation/Data-Processing/StressAnalysis.m
  Reads: Data/Stress/data.*.txt, Data/Walls/floorForce.dat
  Output: Screen figures (add saveas calls to export)

  Places 6 virtual sensors in the stress field (3 bottom → σ_yy, 3 lateral → σ_xx),
  applies cycle averaging, and plots internal stress vs. confining pressure with
  hysteresis. Also includes a consistency check comparing the section-integrated
  stress with the measured wall force.

  ---
  Data File Formats

  Glass beadsGlob.mat

  ┌──────────────────┬─────────────────────────────────────────────────────────────────────────────────┐
  │     Variable     │                                   Description                                   │
  ├──────────────────┼─────────────────────────────────────────────────────────────────────────────────┤
  │ DataTotGlob      │ Internal pressures (N × 6): columns 1–3 lateral σ_xx, columns 4–6 vertical σ_yy │
  ├──────────────────┼─────────────────────────────────────────────────────────────────────────────────┤
  │ InputPresureGlob │ Applied confining pressure (kPa), shape (N,)                                    │
  ├──────────────────┼─────────────────────────────────────────────────────────────────────────────────┤
  │ DisplacementGlob │ LVDT displacement (mm), shape (N,)                                              │
  ├──────────────────┼─────────────────────────────────────────────────────────────────────────────────┤
  │ VeldltGlob       │ Loading velocity metric, shape (N,)                                             │
  ├──────────────────┼─────────────────────────────────────────────────────────────────────────────────┤
  │ Indices          │ Timestep indices marking each pressure step boundary                            │
  └──────────────────┴─────────────────────────────────────────────────────────────────────────────────┘

  ReadyQuasiStaticData_*.mat

  ┌──────────────────┬──────────────────────────────────────────────────────────────────────────────┐
  │     Variable     │                                 Description                                  │
  ├──────────────────┼──────────────────────────────────────────────────────────────────────────────┤
  │ ExternalPressure │ Confining pressure time series (kPa)                                         │
  ├──────────────────┼──────────────────────────────────────────────────────────────────────────────┤
  │ Strain           │ Sample strain (%)                                                            │
  ├──────────────────┼──────────────────────────────────────────────────────────────────────────────┤
  │ stops            │ N × 2 matrix with start/end indices of each compression and relaxation phase │
  └──────────────────┴──────────────────────────────────────────────────────────────────────────────┘

  dt_*.mat

  Single scalar — sampling timestep in seconds. Converts sample indices to physical time.

  filtered_partForce.*.dat

  ESyS-Particle RAW_WITH_POS_ID format. One contact per row:
  particle_id_1  particle_id_2  x  y  z  Fx  Fy  Fz
  Indexed by simulation timestep (e.g., filtered_partForce.100000.dat).

  data.*.txt

  Plain ASCII, one grid cell per row:
  x  y  z  sigma_xx  sigma_xy  sigma_xz  sigma_yy  sigma_yz  sigma_zz

  xyz_Forces_*.csv

  One row per (timestep, zone):
  timestep, zone, N, Fy_mean, Fy_std, Fy_max, frac_strong_y,
  Fx_mean, Fx_std, Fx_max, frac_strong_x,
  Fz_mean, Fz_std, Fz_max, frac_strong_z, anis_xz
  - frac_strong: fraction of contacts with F > ⟨F⟩
  - anis_xz: (⟨Fx²⟩ − ⟨Fz²⟩) / (⟨Fx²⟩ + ⟨Fz²⟩)

  Wall files (floorPosition.dat, etc.)

  Two-column ASCII: simulation time and scalar value (force in N or position in m).

  ---
  Re-running the Simulation

  ▎ The pre-computed data in Simulation/Data/ is sufficient to reproduce all
  ▎ simulation figures. Re-running requires a cluster with ESyS-Particle and MPI.

  cd Simulation/Execution
  mpirun -np 18 python blockCompression.py <geometry_file>

  Then convert outputs to analysis-ready formats:

  # stress tensor grid from contact forces
  python stress3vti.py partForce.N.dat output.vti Xmin Xmax Ymin Ymax Zmin Zmax Nx Ny Nz

  # VTK → plain text
  python extract_data.py

  # spatial force statistics → CSV
  python PDFforce_components.py

  Parameters are in Parameters_general.py (material properties) and
  Parameters_blockCompression.py (loading protocol, MPI decomposition).

  ---

  Contact

  For questions about the data or scripts, please open an issue in this repository.
